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

**Niche regions** split niches into contiguous places: connected pieces of
same-niche cells on the Delaunay graph, with small pieces merged into a
neighbour (:func:`niche_regions`).

numpy / scipy / scikit-learn only (no Sopa).
"""

import itertools

import numpy as np
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans

DEFAULT_NICHE_RADIUS_UM = 100.0
DEFAULT_N_NICHES = 7
# Niche regions smaller than this (cells) are merged into a neighbouring region.
DEFAULT_MIN_REGION_SIZE = 1000
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


def _sym_adjacency(edges, n):
    """Symmetric CSR adjacency (values 1) from an ``(m, 2)`` edge array."""
    import scipy.sparse as sp
    a, b = edges[:, 0], edges[:, 1]
    A = sp.coo_matrix((np.ones(len(a)), (a, b)), shape=(n, n)).tocsr()
    return (A + A.T).tocsr()


def smooth_labels(labels, edges, n_passes=1):
    """Majority vote of each cell's label over itself and its graph neighbours, ``n_passes`` times.

    Ties keep the cell's own label. ``labels`` must be non-negative integers.
    """
    import scipy.sparse as sp
    labels = np.asarray(labels, dtype=int).copy()
    n = len(labels)
    if n == 0 or n_passes < 1:
        return labels
    k = int(labels.max()) + 1
    M = _sym_adjacency(edges, n) + sp.identity(n, format='csr')
    rows = np.arange(n)
    for _ in range(n_passes):
        onehot = sp.csr_matrix((np.ones(n), (rows, labels)), shape=(n, k))
        counts = (M @ onehot).toarray()
        counts[rows, labels] += 0.5          # tie-break in favour of the cell's own label
        new = counts.argmax(axis=1)
        if np.array_equal(new, labels):
            break
        labels = new
    return labels


def merge_small_regions(regions, edges, min_size, max_rounds=100):
    """Merge regions under ``min_size`` cells into the neighbouring region they share most edges with.

    Each round, every small region moves into the larger neighbouring region
    (ties broken by id) it shares most edges with; rounds repeat until no small
    region can move. Small regions
    with no neighbour (isolated islands) become -1. Returns labels renumbered by
    size (0 = largest).
    """
    import scipy.sparse as sp
    regions = np.asarray(regions, dtype=int).copy()
    a, b = edges[:, 0], edges[:, 1]
    for _ in range(max_rounds):
        valid = regions >= 0
        if not valid.any():
            break
        sizes = np.bincount(regions[valid])
        small = (sizes < min_size) & (sizes > 0)
        if not small.any():
            break
        ra, rb = regions[a], regions[b]
        cross = (ra != rb) & (ra >= 0) & (rb >= 0)
        src = np.concatenate([ra[cross], rb[cross]])
        dst = np.concatenate([rb[cross], ra[cross]])
        # Only into a larger neighbour (ties by id): no merge cycles, and a small region
        # never pulls a bigger one into itself.
        keep = small[src] & ((sizes[dst] > sizes[src]) | ((sizes[dst] == sizes[src]) & (dst > src)))
        if not keep.any():
            break
        shared = sp.coo_matrix((np.ones(int(keep.sum())), (src[keep], dst[keep])),
                               shape=(len(sizes), len(sizes))).tocsr()
        target = np.asarray(shared.argmax(axis=1)).ravel()
        has_neighbour = np.diff(shared.indptr) > 0
        mapping = np.arange(len(sizes))
        movers = np.flatnonzero(small & has_neighbour)
        if not len(movers):
            break
        mapping[movers] = target[movers]
        for _ in range(64):                  # follow chains (a → b → c) to their end
            nxt = mapping[mapping]
            if np.array_equal(nxt, mapping):
                break
            mapping = nxt
        regions = np.where(valid, mapping[np.maximum(regions, 0)], -1)
    valid = regions >= 0
    if valid.any():
        sizes = np.bincount(regions[valid])
        regions[valid & (sizes[np.maximum(regions, 0)] < min_size)] = -1
    return renumber_by_size(regions)


def niche_regions(xy_um, niche_labels, max_edge_distance_um, min_size=DEFAULT_MIN_REGION_SIZE,
                  smooth_passes=1):
    """Split niches into contiguous regions ("this stretch of tract", "this patch of deep cortex").

    1. A Delaunay graph of the cells, edges up to ``max_edge_distance_um`` (as the network).
    2. ``smooth_passes`` majority votes of the niche labels over each cell's
       neighbours, which removes speckle at niche borders.
    3. Connected components of same-niche cells.
    4. Regions under ``min_size`` cells merged into the neighbouring region
       they share most edges with (see :func:`merge_small_regions`).

    Regions are contiguous by construction and numbered by size (0 = largest).
    Cells with niche -1 or non-finite coordinates get -1, as do small islands
    with no neighbour. A merged-in fragment keeps its region's label, so a
    region can hold a few cells of other niches; its niche is the majority one.

    Returns ``(regions, info)``; ``info`` has ``n_components`` (before merging).
    """
    from scipy.sparse.csgraph import connected_components
    from cellsurvey.network_analysis import delaunay_edges

    xy = np.asarray(xy_um, dtype=float)[:, :2]
    niche_labels = np.asarray(niche_labels, dtype=int)
    valid = (niche_labels >= 0) & np.isfinite(xy).all(axis=1)
    out = np.full(len(niche_labels), -1, dtype=int)
    idx = np.flatnonzero(valid)
    if len(idx) == 0:
        return out, {'n_components': 0}
    edges = delaunay_edges(xy[idx], max_edge_distance_um).astype(int)
    lab = smooth_labels(niche_labels[idx], edges, smooth_passes)
    same = lab[edges[:, 0]] == lab[edges[:, 1]]
    n_components, comp = connected_components(_sym_adjacency(edges[same], len(idx)), directed=False)
    out[idx] = merge_small_regions(comp, edges, min_size)
    return out, {'n_components': int(n_components)}


def region_niches(regions, niche_labels):
    """Majority niche of each region: a dict ``region → niche`` (regions ≥ 0 only)."""
    regions = np.asarray(regions, dtype=int)
    niche_labels = np.asarray(niche_labels, dtype=int)
    keep = (regions >= 0) & (niche_labels >= 0)
    if not keep.any():
        return {}
    import scipy.sparse as sp
    ct = sp.coo_matrix((np.ones(int(keep.sum())), (regions[keep], niche_labels[keep]))).tocsr()
    major = np.asarray(ct.argmax(axis=1)).ravel()
    present = np.diff(ct.indptr) > 0
    return {int(r): int(major[r]) for r in np.flatnonzero(present)}


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
