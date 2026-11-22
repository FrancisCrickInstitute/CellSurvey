"""Neighbourhood niches: groups of cells whose surroundings have a similar cell-type mix.

For each cell, the fraction of each k-means cluster among the cells within a
radius (µm) is its *neighbourhood composition*; k-means on those compositions
gives **niche types** (the CODEX "cellular neighbourhoods" approach, Schürch et
al. 2020). Unlike the Louvain communities (compact patches of the Delaunay
graph), niches follow regions with a characteristic composition, e.g. white
matter vs cortical layers, and one niche type can occur in several places.

Niches are built on the k-means labels, so they change with ``--n-clusters``.
Suitable radius / niche counts depend heavily on the tissue, panel, cell
density and structures of interest: the defaults are starting points to check
with the niche sweep (``python -m cellsurvey.stability --mode niches``).

numpy / scipy / scikit-learn only (no Sopa).
"""

import itertools

import numpy as np
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans

DEFAULT_NICHE_RADIUS_UM = 100.0
DEFAULT_N_NICHES = 7
# Silhouette compares every pair of cells, so it's computed on a fixed random sample.
SILHOUETTE_SAMPLE = 20_000


def table_centroids_um(sdata, pixel_size_um, table_key='table', cell_boundaries=None):
    """``(n_cells, 2)`` cell centroids in µm, in the table's ``obs`` order.

    From ``table.obsm['spatial']`` (Sopa's centroids, in pixels) when present;
    otherwise from the cell shapes' centroids, matched by cell id (``NaN`` for
    cells without a shape).
    """
    table = sdata.tables[table_key]
    if 'spatial' in table.obsm:
        return np.asarray(table.obsm['spatial'], dtype=float)[:, :2] * pixel_size_um
    from cellsurvey.utils import resolve_shape_keys
    shapes = sdata.shapes[cell_boundaries or resolve_shape_keys(sdata)[1]]
    cent = shapes.geometry.centroid
    lookup = {str(i): (p.x, p.y) for i, p in zip(shapes.index, cent)}
    xy = np.array([lookup.get(str(i), (np.nan, np.nan)) for i in table.obs.index], dtype=float)
    return xy * pixel_size_um


def neighbourhood_composition(xy_um, labels, radius_um, n_labels=None):
    """Fraction of each label among the cells within ``radius_um`` of each cell (itself included).

    Args:
        xy_um: ``(n, 2)`` cell centroids in µm.
        labels: ``(n,)`` integer cluster labels; cells with a negative label (or
            non-finite coordinates) are neither counted nor given a composition.
        radius_um: neighbourhood radius in µm.
        n_labels: number of label columns (default: max label + 1).

    Returns:
        ``(n, n_labels)`` float array whose rows sum to 1; ``NaN`` rows for
        excluded cells. Built with one ``cKDTree`` per label and
        ``query_ball_point(return_length=True)``, so memory stays flat at any
        radius (no pair list).
    """
    xy = np.asarray(xy_um, dtype=float)[:, :2]
    labels = np.asarray(labels, dtype=int)
    valid = (labels >= 0) & np.isfinite(xy).all(axis=1)
    if n_labels is None:
        n_labels = int(labels[valid].max()) + 1 if valid.any() else 0
    counts = np.zeros((len(labels), n_labels))
    query_idx = np.flatnonzero(valid)
    for c in range(n_labels):
        pts = xy[valid & (labels == c)]
        if len(pts) == 0:
            continue
        tree = cKDTree(pts)
        counts[query_idx, c] = tree.query_ball_point(xy[query_idx], r=radius_um, return_length=True,
                                                     workers=-1)
    total = counts.sum(axis=1, keepdims=True)
    comp = np.divide(counts, total, out=np.zeros_like(counts), where=total > 0)
    comp[~valid] = np.nan
    return comp


def renumber_by_size(labels):
    """Relabel 0 = largest group, 1 = next, …; negative labels stay -1. Ties keep the original order."""
    labels = np.asarray(labels, dtype=int)
    out = np.full(labels.shape, -1, dtype=int)
    valid = labels >= 0
    if not valid.any():
        return out
    uniq, counts = np.unique(labels[valid], return_counts=True)
    order = np.lexsort((uniq, -counts))
    mapping = np.full(uniq.max() + 1, -1, dtype=int)
    mapping[uniq[order]] = np.arange(len(uniq))
    out[valid] = mapping[labels[valid]]
    return out


def _fit(composition, n_niches, seed, n_init):
    """k-means on the valid rows; returns (labels with -1 for invalid rows, inertia)."""
    valid = np.isfinite(composition).all(axis=1)
    labels = np.full(len(composition), -1, dtype=int)
    if valid.sum() < n_niches:
        raise ValueError(f"Only {int(valid.sum())} cells have a neighbourhood; can't make {n_niches} niches")
    km = KMeans(n_clusters=n_niches, random_state=seed, n_init=n_init)
    labels[valid] = km.fit_predict(composition[valid])
    return labels, float(km.inertia_)


def assign_niches(composition, n_niches, seed=42, n_init=4):
    """Niche label per cell: k-means on the compositions, renumbered by size (0 = largest).

    Rows of ``composition`` with ``NaN`` (excluded cells) get -1.
    """
    labels, _ = _fit(composition, n_niches, seed, n_init)
    return renumber_by_size(labels)


def niche_profiles(composition, niche_labels):
    """``(n_niches, n_labels)`` mean neighbourhood composition of each niche (rows: niche 0, 1, …)."""
    niche_labels = np.asarray(niche_labels)
    ids = np.unique(niche_labels[niche_labels >= 0])
    return np.vstack([np.nanmean(composition[niche_labels == n], axis=0) for n in ids]) if len(ids) else \
        np.zeros((0, composition.shape[1]))


def sweep_niches(xy_um, cluster_labels, radii_um, n_niches_list, seeds=(42, 1, 2), n_init=4,
                 silhouette_sample=SILHOUETTE_SAMPLE, verbose=True):
    """Radius × niche-count sweep with silhouette and seed stability.

    For each radius, compositions are computed once; for each niche count,
    k-means runs once per seed (each with ``n_init`` starts). The lowest-inertia
    fit gives the labels; **seed stability** is the mean pairwise adjusted Rand
    index between the seeds' labels (1 = the same niches every time).
    Silhouette uses one fixed random sample of cells for every setting.

    Note: silhouette rises with radius (bigger windows give smoother, more
    separable compositions), so it can compare niche counts at one radius but
    can't choose the radius; that is a resolution choice.

    Returns a dict: ``params`` (list of ``(radius_um, n_niches)``), ``labels``
    (``(n, n_settings)`` int array, renumbered by size), ``silhouette``,
    ``seed_ari``, ``inertia``, ``min_niche_size`` and ``median_neighbours``
    (per radius, per setting).
    """
    from sklearn.metrics import adjusted_rand_score, silhouette_score

    cluster_labels = np.asarray(cluster_labels, dtype=int)
    n_labels = int(cluster_labels[cluster_labels >= 0].max()) + 1
    rng = np.random.default_rng(seeds[0])
    params, labels_out = [], []
    sil, ari, inertia, min_size, med_nb = [], [], [], [], []
    for radius in radii_um:
        comp = neighbourhood_composition(xy_um, cluster_labels, radius, n_labels)
        valid = np.isfinite(comp).all(axis=1)
        valid_idx = np.flatnonzero(valid)
        sample = (rng.choice(valid_idx, min(silhouette_sample, len(valid_idx)), replace=False)
                  if len(valid_idx) else valid_idx)
        # Median neighbour count (cells within the radius, itself included).
        xy = np.asarray(xy_um, dtype=float)[:, :2]
        tree = cKDTree(xy[valid])
        sub = sample if len(sample) else valid_idx
        neighbours = tree.query_ball_point(xy[sub], r=radius, return_length=True, workers=-1)
        median_nb = float(np.median(neighbours)) if len(neighbours) else float('nan')
        for k in n_niches_list:
            fits = [_fit(comp, k, s, n_init) for s in seeds]
            best_labels, best_inertia = min(fits, key=lambda f: f[1])
            pair_ari = [adjusted_rand_score(a[0][valid], b[0][valid]) for a, b in itertools.combinations(fits, 2)]
            s_val = (float(silhouette_score(comp[sample], best_labels[sample]))
                     if len(np.unique(best_labels[sample])) > 1 else float('nan'))
            renum = renumber_by_size(best_labels)
            params.append((float(radius), int(k)))
            labels_out.append(renum)
            sil.append(s_val)
            ari.append(float(np.mean(pair_ari)) if pair_ari else float('nan'))
            inertia.append(best_inertia)
            min_size.append(int(np.bincount(renum[renum >= 0], minlength=k).min()))
            med_nb.append(median_nb)
            if verbose:
                print(f"radius={radius:g} µm (median {median_nb:.0f} neighbours), niches={k}: "
                      f"silhouette={s_val:.3f}, seed ARI={ari[-1]:.2f}, smallest niche {min_size[-1]:,} cells")
    return {
        'params': params,
        'labels': np.column_stack(labels_out) if labels_out else np.zeros((len(cluster_labels), 0), dtype=int),
        'silhouette': sil, 'seed_ari': ari, 'inertia': inertia, 'min_niche_size': min_size,
        'median_neighbours': med_nb,
    }
