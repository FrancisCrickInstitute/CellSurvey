"""Parameter sweeps over an existing segmented Zarr (read-only).

Two modes (selected with ``--mode``):

- **Community sweep** (default): re-runs Delaunay + Louvain community detection
  across a grid of (``community_resolution``, ``max_edge_distance``) values.
- **Cluster sweep**: re-runs k-means clustering across a range of ``n_clusters``
  values.

Both write per-cell labels to CSV (by default into the Zarr's directory) and
never write back to the Zarr, so the original result is left untouched.

Note: k-means clustering is independent of community detection (Louvain weights
edges from the ``--cluster-on`` features, not the k-means labels).

Run:  python -m cellsurvey.stability --zarr <output>_seg.zarr [--resolutions ...]
      python -m cellsurvey.stability --zarr <output>_seg.zarr --mode clusters
"""
import argparse
import os

import numpy as np
import pandas as pd
import spatialdata

from sklearn.cluster import KMeans

from cellsurvey.network_analysis import compute_louvain_communities, filter_small_communities, n_labels
from cellsurvey.utils import (resolve_shape_keys, resolve_pixel_size, resolve_cluster_on,
                              resolve_min_community_size, extract_feature_tables,
                              build_clustering_features, CLUSTER_ON_CHOICES)


# Silhouette scores compare every pair of cells, so they're computed on a random
# sample (20k cells → 4e8 distances, a few seconds per k).
SILHOUETTE_SAMPLE = 20_000


def _fmt(v):
    """Render a number as a tidy column-name token (drop a trailing ``.0``)."""
    return str(int(v)) if float(v) == int(v) else str(v)


def extract_intensity_matrix(sdata, table_key="table"):
    """Per-cell channel intensities (from ``obsm['intensities']`` when the run had spots)."""
    return extract_feature_tables(sdata.tables[table_key])[0]


def extract_features(sdata, cluster_on=None, pixel_size_um=None, table_key="table"):
    """The scaled per-cell features the pipeline clusters on for ``cluster_on``
    (see :func:`cellsurvey.utils.build_clustering_features`). ``None`` uses the
    mode the run recorded in the Zarr (see :func:`cellsurvey.utils.resolve_cluster_on`).
    Spot modes need the pixel size (``pixel_size_um`` or the value stored in the
    Zarr) to turn cell areas into µm²."""
    cluster_on, _ = resolve_cluster_on(sdata, cluster_on, table_key)
    table = sdata.tables[table_key]
    intensity_df, counts_df = extract_feature_tables(table)
    cell_area_um2 = None
    if cluster_on != "intensities" and counts_df is not None:
        pixel_size_um = resolve_pixel_size(sdata, pixel_size_um)
        cell_area_um2 = table.obs["area"].astype(float) * pixel_size_um ** 2
    return build_clustering_features(intensity_df, counts_df, cluster_on, cell_area_um2)


def extract_coords(sdata, cell_boundaries=None):
    boundaries = sdata.shapes[cell_boundaries or resolve_shape_keys(sdata)[1]]
    centroids = boundaries.geometry.centroid
    return boundaries.index.values, np.column_stack([centroids.x.values, centroids.y.values])


def run_stability_sweep(sdata, resolutions, max_edge_distances,
                        cell_boundaries=None, seed=42, pixel_size_um=None, cluster_on=None,
                        min_community_size=None):
    """Recompute Louvain communities across a parameter grid.

    ``max_edge_distances`` are in µm, converted to pixels with ``pixel_size_um``
    (µm/px) or, if that is None, the pixel size stored in the table's ``.uns``.
    Edges are weighted by the ``cluster_on`` features, as in the pipeline
    (``None``: the mode the run recorded in the Zarr). As in the pipeline,
    communities smaller than ``min_community_size`` cells become -1 and the
    rest are numbered by size (``None``: the value recorded in the Zarr).

    Returns a dict with keys ``cell_ids``, ``coords``, ``labels`` (an
    ``(n_cells, n_sweeps)`` int array, one column per parameter combo),
    ``params`` (the list of ``(resolution, max_edge_distance_um)`` tuples
    matching the columns), ``n_communities_raw`` (Louvain's count per sweep,
    before the minimum size) and ``min_community_size``.
    """
    pixel_size_um = resolve_pixel_size(sdata, pixel_size_um)
    min_community_size, _ = resolve_min_community_size(sdata, min_community_size)
    cell_ids, coords = extract_coords(sdata, cell_boundaries)
    features = extract_features(sdata, cluster_on, pixel_size_um)

    imat = np.zeros((len(coords), features.shape[1]))
    for pos_idx, cell_id in enumerate(cell_ids):
        if cell_id in features.index:
            imat[pos_idx] = features.loc[cell_id].values

    labels, params, n_raw = [], [], []
    for max_edge_distance in max_edge_distances:
        for resolution in resolutions:
            raw_labels, n_edges = compute_louvain_communities(
                coords, imat, comm_detect_res=resolution,
                max_edge_distance=max_edge_distance / pixel_size_um, seed=seed)
            community_labels = filter_small_communities(raw_labels, min_community_size)
            n_comm = n_labels(community_labels)
            n_unassigned = int((community_labels < 0).sum())
            print(f"resolution={resolution}, max_edge_distance={max_edge_distance} µm: "
                  f"{n_edges} edges, {n_comm} communities of >= {min_community_size} cells "
                  f"({n_labels(raw_labels)} before the minimum; {n_unassigned} cells unassigned)")
            labels.append(community_labels)
            params.append((resolution, max_edge_distance))
            n_raw.append(n_labels(raw_labels))

    return {
        "cell_ids": cell_ids,
        "coords": coords,
        "labels": np.column_stack(labels) if labels else np.zeros((len(coords), 0), dtype=int),
        "params": params,
        "n_communities_raw": n_raw,
        "min_community_size": min_community_size,
    }


def sweep_to_csv(result, output_path):
    """Write one row per cell, with a ``community_r<res>_d<dist_um>`` column per sweep."""
    df = pd.DataFrame({
        "cell_id": result["cell_ids"],
        "x": result["coords"][:, 0],
        "y": result["coords"][:, 1],
    })
    for (resolution, max_edge_distance), col in zip(result["params"], result["labels"].T):
        df[f"community_r{_fmt(resolution)}_d{_fmt(max_edge_distance)}"] = col
    df.to_csv(output_path, index=False)
    return df


def sweep_summary(result, output_path):
    """Write one row per sweep with the community counts and unassigned cells."""
    rows = []
    n_raw = result.get("n_communities_raw") or [None] * len(result["params"])
    for (resolution, max_edge_distance), col, raw in zip(result["params"], result["labels"].T, n_raw):
        rows.append({
            "resolution": resolution,
            "max_edge_distance_um": max_edge_distance,
            "n_communities": n_labels(col),
            "n_communities_raw": raw,
            "min_community_size": result.get("min_community_size"),
            "n_unassigned_cells": int((col < 0).sum()),
        })
    pd.DataFrame(rows).to_csv(output_path, index=False)


def _coords_for_cell_ids(sdata, cell_ids, cell_boundaries=None, table_key="table"):
    """Return an (n, 2) array of centroids aligned to ``cell_ids``.

    Uses the table's ``obsm['spatial']`` (Sopa's centroids, one row per table
    cell) when no shape layer is requested, so the polygons needn't be loaded.
    """
    table = sdata.tables.get(table_key) if hasattr(sdata.tables, 'get') else sdata.tables[table_key]
    if cell_boundaries is None and table is not None and "spatial" in table.obsm:
        xy = np.asarray(table.obsm["spatial"], dtype=float)[:, :2]
        pos = pd.Series(np.arange(table.n_obs), index=table.obs.index.astype(str))
        idx = pos.reindex([str(c) for c in cell_ids]).to_numpy()
        out = np.full((len(cell_ids), 2), np.nan)
        ok = ~np.isnan(idx)
        out[ok] = xy[idx[ok].astype(int)]
        return out
    boundaries = sdata.shapes[cell_boundaries or resolve_shape_keys(sdata)[1]]
    lookup = {cid: (g.centroid.x, g.centroid.y)
              for cid, g in zip(boundaries.index.values, boundaries.geometry)}
    return np.array([lookup.get(cid, (np.nan, np.nan)) for cid in cell_ids], dtype=float)


def run_cluster_sweep(sdata, n_clusters_list, cell_boundaries=None, seed=42, cluster_on=None,
                      pixel_size_um=None):
    """Recompute k-means clusters across a range of ``n_clusters``.

    Clusters on the ``cluster_on`` features, scaled as in the pipeline
    (``None``: the mode the run recorded in the Zarr).
    Returns a dict with keys ``cell_ids``, ``coords``, ``labels`` (an
    ``(n_cells, n_sweeps)`` int array, one column per ``k``), ``params`` (the
    list of ``k`` values), ``inertias``, ``silhouettes`` (mean silhouette on a
    random sample of up to ``SILHOUETTE_SAMPLE`` cells; higher = better
    separated clusters) and ``min_cluster_sizes`` (cells in the smallest cluster).
    """
    from sklearn.metrics import silhouette_score

    features = extract_features(sdata, cluster_on, pixel_size_um)
    scaled = features.to_numpy()
    cell_ids = features.index.values
    coords = _coords_for_cell_ids(sdata, cell_ids, cell_boundaries)
    sample = min(SILHOUETTE_SAMPLE, len(scaled))

    labels, inertias, k_values, silhouettes, min_sizes = [], [], [], [], []
    for k in n_clusters_list:
        km = KMeans(n_clusters=k, random_state=seed)
        lab = km.fit_predict(scaled)
        sil = (float(silhouette_score(scaled, lab, sample_size=sample, random_state=seed))
               if 1 < k < len(scaled) else float('nan'))
        labels.append(lab)
        inertias.append(float(km.inertia_))
        silhouettes.append(sil)
        min_sizes.append(int(np.bincount(lab, minlength=k).min()))
        k_values.append(k)
        print(f"n_clusters={k}: inertia={km.inertia_:.1f}, silhouette={sil:.3f} (sample of {sample:,}), "
              f"smallest cluster {min_sizes[-1]:,} cells")

    return {
        "cell_ids": cell_ids,
        "coords": coords,
        "labels": np.column_stack(labels) if labels else np.zeros((len(cell_ids), 0), dtype=int),
        "params": k_values,
        "inertias": inertias,
        "silhouettes": silhouettes,
        "min_cluster_sizes": min_sizes,
    }


def cluster_sweep_to_csv(result, output_path):
    """Write one row per cell, with a ``cluster_k<k>`` column per sweep."""
    df = pd.DataFrame({
        "cell_id": result["cell_ids"],
        "x": result["coords"][:, 0],
        "y": result["coords"][:, 1],
    })
    for k, col in zip(result["params"], result["labels"].T):
        df[f"cluster_k{k}"] = col
    df.to_csv(output_path, index=False)
    return df


def cluster_sweep_summary(result, output_path):
    """Write one row per sweep: ``n_clusters``, k-means inertia, silhouette and the smallest cluster."""
    n_cells = len(result["cell_ids"])
    rows = []
    for i, k in enumerate(result["params"]):
        row = {"n_clusters": k, "inertia": result["inertias"][i]}
        if "silhouettes" in result:
            row["silhouette"] = result["silhouettes"][i]
        if "min_cluster_sizes" in result:
            row["min_cluster_size"] = result["min_cluster_sizes"][i]
            row["min_cluster_pct"] = 100 * result["min_cluster_sizes"][i] / max(n_cells, 1)
        rows.append(row)
    pd.DataFrame(rows).to_csv(output_path, index=False)


def sweep_prefix(zarr_path):
    """Run name for sweep outputs: ``out/run_seg.zarr`` → ``run`` (so runs sharing a folder don't clash)."""
    name = os.path.basename(os.path.normpath(zarr_path))
    for suffix in ("_seg.zarr", ".zarr"):
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return name


def niche_sweep_to_csv(result, cell_ids, xy_um, output_path):
    """Write one row per cell, with a ``niche_r<radius>_n<k>`` column per setting (-1 = no niche)."""
    df = pd.DataFrame({"cell_id": cell_ids, "x_um": xy_um[:, 0], "y_um": xy_um[:, 1]})
    for (radius, k), col in zip(result["params"], result["labels"].T):
        df[f"niche_r{_fmt(radius)}_n{k}"] = col
    df.to_csv(output_path, index=False)
    return df


def niche_sweep_summary(result, n_cells, output_path):
    """One row per (radius, niche count): silhouette, seed stability (ARI), smallest niche, neighbours."""
    rows = []
    for i, (radius, k) in enumerate(result["params"]):
        rows.append({
            "niche_radius_um": radius,
            "n_niches": k,
            "silhouette": result["silhouette"][i],
            "seed_ari": result["seed_ari"][i],
            "min_niche_size": result["min_niche_size"][i],
            "min_niche_pct": 100 * result["min_niche_size"][i] / max(n_cells, 1),
            "median_neighbours": result["median_neighbours"][i],
            "inertia": result["inertia"][i],
        })
    pd.DataFrame(rows).to_csv(output_path, index=False)


def sweep_paths(zarr_path, mode, folder=None):
    """Default ``(per-cell CSV, summary CSV)`` for a sweep of ``zarr_path`` (``mode``: communities|clusters|niches)."""
    folder = folder or os.path.dirname(os.path.abspath(zarr_path))
    stem = {"communities": "community", "clusters": "cluster", "niches": "niche"}[mode]
    base = os.path.join(folder, f"{sweep_prefix(zarr_path)}_{stem}_sweep")
    return f"{base}.csv", f"{base}_summary.csv"


def main():
    parser = argparse.ArgumentParser(
        description="Parameter sweeps over an existing segmented Zarr (read-only).")
    parser.add_argument("--zarr", required=True, help="Path to the segmented Zarr (*_seg.zarr).")
    parser.add_argument("--mode", choices=["communities", "clusters", "niches"], default="communities",
                        help="Which sweep to run (default: communities).")
    parser.add_argument("--niche-radii", default="75,100,150",
                        help="Neighbourhood radii in µm to sweep (niche mode).")
    parser.add_argument("--n-niches", default="4,5,6,7,8",
                        help="Niche counts to sweep (niche mode).")
    parser.add_argument("--resolutions", default="0.1,0.05,0.02,0.01",
                        help="Louvain resolutions to sweep (community mode).")
    parser.add_argument("--max-edge-distances", default="50",
                        help="Max edge distances in µm to sweep (community mode).")
    parser.add_argument("--min-community-size", type=int, default=None,
                        help="Communities with fewer cells become -1 (unassigned), as in the pipeline "
                             "(community mode; default: the value the run used, recorded in the Zarr; "
                             "10 for Zarrs that predate it).")
    parser.add_argument("--pixel-size", type=float, default=None,
                        help="Pixel size in µm/px (default: the value stored in the Zarr by the pipeline; "
                             "required for Zarrs written before pixel sizes were stored).")
    parser.add_argument("--n-clusters", default="5,8,10,12,15,20",
                        help="k-means k values to sweep (cluster mode).")
    parser.add_argument("--cluster-on", choices=CLUSTER_ON_CHOICES, default=None,
                        help="Features for k-means (cluster mode) and the network edge weights (community "
                             "mode), as the pipeline's --cluster-on (default: the mode the run used, "
                             "recorded in the Zarr; 'intensities' for Zarrs that predate it).")
    parser.add_argument("--cell-boundaries", default=None,
                        help="Shape layer to sweep over (default: cell_boundaries if present, else the nuclei layer).")
    parser.add_argument("--output", default=None,
                        help="Output CSV for per-cell labels (default: <zarr_dir>/<run>_<community|cluster>"
                             "_sweep.csv, where <run> is the Zarr name without _seg.zarr).")
    parser.add_argument("--summary", default=None,
                        help="Output CSV for per-sweep summary (default: <zarr_dir>/<run>_<community|cluster>"
                             "_sweep_summary.csv).")
    args = parser.parse_args()

    # The cluster sweep needs only the table (coordinates from obsm['spatial']); skip
    # loading the cell polygons, which takes minutes on whole slides.
    sdata = None
    if args.mode in ("clusters", "niches") and args.cell_boundaries is None:
        try:
            sdata = spatialdata.read_zarr(args.zarr, selection=("tables",))
        except TypeError:
            sdata = None
    if sdata is None:
        sdata = spatialdata.read_zarr(args.zarr)
    default_output, default_summary = sweep_paths(args.zarr, args.mode)
    cluster_on, origin = resolve_cluster_on(sdata, args.cluster_on)
    print(f"Features: --cluster-on {cluster_on} ({origin})")

    if args.mode == "clusters":
        n_clusters_list = [int(x) for x in args.n_clusters.split(",")]
        try:
            result = run_cluster_sweep(sdata, n_clusters_list, cell_boundaries=args.cell_boundaries,
                                       cluster_on=cluster_on, pixel_size_um=args.pixel_size)
        except ValueError as e:
            parser.error(str(e))
        output = args.output or default_output
        summary = args.summary or default_summary
        cluster_sweep_to_csv(result, output)
        cluster_sweep_summary(result, summary)
    elif args.mode == "niches":
        from cellsurvey.niches import sweep_niches, table_centroids_um
        try:
            radii = [float(x) for x in args.niche_radii.split(",")]
            counts = [int(x) for x in args.n_niches.split(",")]
        except ValueError:
            parser.error("--niche-radii must be comma-separated numbers and --n-niches comma-separated integers")
        if not all(r > 0 for r in radii) or not all(k >= 2 for k in counts):
            parser.error("--niche-radii must be positive and --n-niches at least 2")
        table = sdata.tables["table"]
        if "kmeans_cluster" not in table.obs:
            parser.error(f"{args.zarr} has no k-means clusters (obs['kmeans_cluster'])")
        try:
            pixel_size_um = resolve_pixel_size(sdata, args.pixel_size)
        except ValueError as e:
            parser.error(str(e))
        xy_um = table_centroids_um(sdata, pixel_size_um, cell_boundaries=args.cell_boundaries)
        cluster_labels = table.obs["kmeans_cluster"].astype(int).to_numpy()
        print(f"Niches from each cell's mix of {len(np.unique(cluster_labels))} k-means clusters; "
              f"radii {radii} µm × niche counts {counts}, 3 seeds each")
        result = sweep_niches(xy_um, cluster_labels, radii, counts)
        output = args.output or default_output
        summary = args.summary or default_summary
        niche_sweep_to_csv(result, table.obs.index.astype(str).to_numpy(), xy_um, output)
        niche_sweep_summary(result, table.n_obs, summary)
    else:
        resolutions = [float(x) for x in args.resolutions.split(",")]
        max_edge_distances = [float(x) for x in args.max_edge_distances.split(",")]
        try:
            pixel_size_um = resolve_pixel_size(sdata, args.pixel_size)
        except ValueError as e:
            parser.error(str(e))
        if args.min_community_size is not None and args.min_community_size < 1:
            parser.error("--min-community-size must be at least 1")
        min_size, min_origin = resolve_min_community_size(sdata, args.min_community_size)
        print(f"Minimum community size: {min_size} cells ({min_origin})")
        try:
            result = run_stability_sweep(sdata, resolutions, max_edge_distances,
                                         cell_boundaries=args.cell_boundaries,
                                         pixel_size_um=pixel_size_um, cluster_on=cluster_on,
                                         min_community_size=min_size)
        except ValueError as e:
            parser.error(str(e))
        output = args.output or default_output
        summary = args.summary or default_summary
        sweep_to_csv(result, output)
        sweep_summary(result, summary)

    print(f"Wrote {output} and {summary}")


if __name__ == "__main__":
    main()
