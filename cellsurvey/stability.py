"""Parameter sweeps over an existing segmented Zarr (read-only).

Two modes (selected with ``--mode``):

- **Community sweep** (default): re-runs Delaunay + Louvain community detection
  across a grid of (``community_resolution``, ``max_edge_distance``) values.
- **Cluster sweep**: re-runs k-means clustering across a range of ``n_clusters``
  values.

Both write per-cell labels to CSV (by default into the Zarr's directory) and
never write back to the Zarr, so the original result is left untouched.

Note: k-means clustering is independent of community detection (Louvain weights
edges from the raw intensity matrix, not the k-means labels).

Run:  python -m cellsurvey.stability --zarr <output>_seg.zarr [--resolutions ...]
      python -m cellsurvey.stability --zarr <output>_seg.zarr --mode clusters
"""
import argparse
import os

import numpy as np
import pandas as pd
import spatialdata

from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

from cellsurvey.network_analysis import compute_louvain_communities
from cellsurvey.utils import remove_channel_suffix


def _fmt(v):
    """Render a number as a tidy column-name token (drop a trailing ``.0``)."""
    return str(int(v)) if float(v) == int(v) else str(v)


def extract_intensity_matrix(sdata, table_key="table"):
    measurements = sdata.tables[table_key]
    X = measurements.X.toarray() if hasattr(measurements.X, "toarray") else measurements.X
    df = pd.DataFrame(X, index=measurements.obs.index, columns=measurements.var.index)
    df.columns = [remove_channel_suffix(c) for c in df.columns]
    return df.loc[:, ~df.columns.duplicated(keep="first")]


def extract_coords(sdata, cell_boundaries="stardist_boundaries"):
    boundaries = sdata.shapes[cell_boundaries]
    centroids = boundaries.geometry.centroid
    return boundaries.index.values, np.column_stack([centroids.x.values, centroids.y.values])


def run_stability_sweep(sdata, resolutions, max_edge_distances,
                        cell_boundaries="stardist_boundaries", seed=42):
    """Recompute Louvain communities across a parameter grid.

    Returns a dict with keys ``cell_ids``, ``coords``, ``labels`` (an
    ``(n_cells, n_sweeps)`` int array, one column per parameter combo), and
    ``params`` (the list of ``(resolution, max_edge_distance)`` tuples matching
    the columns).
    """
    cell_ids, coords = extract_coords(sdata, cell_boundaries)
    intensity_df = extract_intensity_matrix(sdata)

    imat = np.zeros((len(coords), intensity_df.shape[1]))
    for pos_idx, cell_id in enumerate(cell_ids):
        if cell_id in intensity_df.index:
            imat[pos_idx] = intensity_df.loc[cell_id].values

    labels, params = [], []
    for max_edge_distance in max_edge_distances:
        for resolution in resolutions:
            community_labels, n_edges = compute_louvain_communities(
                coords, imat, comm_detect_res=resolution,
                max_edge_distance=max_edge_distance, seed=seed)
            n_comm = len(set(community_labels)) - (1 if -1 in community_labels else 0)
            print(f"resolution={resolution}, max_edge_distance={max_edge_distance}: "
                  f"{n_edges} edges, {n_comm} communities")
            labels.append(community_labels)
            params.append((resolution, max_edge_distance))

    return {
        "cell_ids": cell_ids,
        "coords": coords,
        "labels": np.column_stack(labels) if labels else np.zeros((len(coords), 0), dtype=int),
        "params": params,
    }


def sweep_to_csv(result, output_path):
    """Write one row per cell, with a ``community_r<res>_d<dist>`` column per sweep."""
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
    """Write one row per sweep with the resulting community count."""
    rows = []
    for (resolution, max_edge_distance), col in zip(result["params"], result["labels"].T):
        n_comm = len(set(col)) - (1 if -1 in col else 0)
        rows.append({
            "resolution": resolution,
            "max_edge_distance": max_edge_distance,
            "n_communities": n_comm,
        })
    pd.DataFrame(rows).to_csv(output_path, index=False)


def _coords_for_cell_ids(sdata, cell_ids, cell_boundaries="stardist_boundaries"):
    """Return an (n, 2) array of centroids aligned to ``cell_ids``."""
    boundaries = sdata.shapes[cell_boundaries]
    lookup = {cid: (g.centroid.x, g.centroid.y)
              for cid, g in zip(boundaries.index.values, boundaries.geometry)}
    return np.array([lookup.get(cid, (np.nan, np.nan)) for cid in cell_ids], dtype=float)


def run_cluster_sweep(sdata, n_clusters_list, cell_boundaries="stardist_boundaries", seed=42):
    """Recompute k-means clusters across a range of ``n_clusters``.

    Returns a dict with keys ``cell_ids``, ``coords``, ``labels`` (an
    ``(n_cells, n_sweeps)`` int array, one column per ``k``), ``params`` (the
    list of ``k`` values), and ``inertias``.
    """
    intensity_df = extract_intensity_matrix(sdata)
    scaled = StandardScaler().fit_transform(intensity_df)
    cell_ids = intensity_df.index.values
    coords = _coords_for_cell_ids(sdata, cell_ids, cell_boundaries)

    labels, inertias, k_values = [], [], []
    for k in n_clusters_list:
        km = KMeans(n_clusters=k, random_state=seed)
        lab = km.fit_predict(scaled)
        labels.append(lab)
        inertias.append(float(km.inertia_))
        k_values.append(k)
        print(f"n_clusters={k}: inertia={km.inertia_:.1f}")

    return {
        "cell_ids": cell_ids,
        "coords": coords,
        "labels": np.column_stack(labels) if labels else np.zeros((len(cell_ids), 0), dtype=int),
        "params": k_values,
        "inertias": inertias,
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
    """Write one row per sweep with ``n_clusters`` and the k-means inertia."""
    rows = [{"n_clusters": k, "inertia": inertia}
            for k, inertia in zip(result["params"], result["inertias"])]
    pd.DataFrame(rows).to_csv(output_path, index=False)


def main():
    parser = argparse.ArgumentParser(
        description="Parameter sweeps over an existing segmented Zarr (read-only).")
    parser.add_argument("--zarr", required=True, help="Path to the segmented Zarr (*_seg.zarr).")
    parser.add_argument("--mode", choices=["communities", "clusters"], default="communities",
                        help="Which sweep to run (default: communities).")
    parser.add_argument("--resolutions", default="0.1,0.05,0.02,0.01",
                        help="Louvain resolutions to sweep (community mode).")
    parser.add_argument("--max-edge-distances", default="1000",
                        help="Max edge distances to sweep (community mode).")
    parser.add_argument("--n-clusters", default="5,8,10,12,15,20",
                        help="k-means k values to sweep (cluster mode).")
    parser.add_argument("--output", default=None,
                        help="Output CSV for per-cell labels (default: <zarr_dir>/<mode>_sweep.csv).")
    parser.add_argument("--summary", default=None,
                        help="Output CSV for per-sweep summary (default: <zarr_dir>/<mode>_sweep_summary.csv).")
    args = parser.parse_args()

    sdata = spatialdata.read_zarr(args.zarr)
    zarr_dir = os.path.dirname(args.zarr) or "."

    if args.mode == "clusters":
        n_clusters_list = [int(x) for x in args.n_clusters.split(",")]
        result = run_cluster_sweep(sdata, n_clusters_list)
        output = args.output or os.path.join(zarr_dir, "cluster_sweep.csv")
        summary = args.summary or os.path.join(zarr_dir, "cluster_sweep_summary.csv")
        cluster_sweep_to_csv(result, output)
        cluster_sweep_summary(result, summary)
    else:
        resolutions = [float(x) for x in args.resolutions.split(",")]
        max_edge_distances = [float(x) for x in args.max_edge_distances.split(",")]
        result = run_stability_sweep(sdata, resolutions, max_edge_distances)
        output = args.output or os.path.join(zarr_dir, "community_sweep.csv")
        summary = args.summary or os.path.join(zarr_dir, "community_sweep_summary.csv")
        sweep_to_csv(result, output)
        sweep_summary(result, summary)

    print(f"Wrote {output} and {summary}")


if __name__ == "__main__":
    main()
