import os
import json
import numpy as np
from scipy.spatial import Delaunay
import networkx as nx

from cellsurvey.utils import resolve_shape_keys, resolve_pixel_size


def delaunay_edges(coords, max_edge_distance):
    """Return the unique Delaunay edges of ``coords`` no longer than ``max_edge_distance``.

    ``max_edge_distance`` is in the same units as ``coords``. Each edge appears
    once as ``(a, b)`` with ``a < b``; an edge is shared by two triangles, so
    without sorting it would be counted twice.
    """
    simplices = Delaunay(coords).simplices
    pairs = np.vstack([simplices[:, [0, 1]], simplices[:, [1, 2]], simplices[:, [0, 2]]])
    pairs = np.unique(np.sort(pairs, axis=1), axis=0)
    lengths = np.linalg.norm(coords[pairs[:, 0]] - coords[pairs[:, 1]], axis=1)
    return pairs[lengths <= max_edge_distance]


def compute_louvain_communities(coords, imat, comm_detect_res=0.1, max_edge_distance=None, seed=42):
    """Build a Delaunay graph of ``coords``, weighted by expression similarity in
    ``imat``, and return per-node Louvain community labels.

    ``max_edge_distance`` is in the same units as ``coords`` (pixels for
    CellSurvey shapes); ``None`` keeps every edge. ``imat`` is an
    ``(n_cells, n_markers)`` array aligned to ``coords``; if None, the graph is
    left unweighted. Returns ``(community_labels, n_edges)``.
    """
    n_cells = len(coords)
    edges = delaunay_edges(coords, np.inf if max_edge_distance is None else max_edge_distance)
    edges = [(int(a), int(b)) for a, b in edges]

    G = nx.Graph()
    G.add_nodes_from(range(n_cells))
    if imat is not None:
        for a, b in edges:
            corr = np.corrcoef(imat[a], imat[b])[0, 1]
            if not np.isfinite(corr):
                corr = 0.0
            G.add_edge(a, b, weight=1 + corr)
    else:
        G.add_edges_from(edges)

    communities = nx.community.louvain_communities(
        G, weight='weight' if imat is not None else None,
        resolution=comm_detect_res, seed=seed)

    community_labels = np.full(n_cells, -1, dtype=int)
    for comm_idx, comm in enumerate(communities):
        for node in comm:
            community_labels[node] = comm_idx

    return community_labels, len(edges)


def run_network_analysis(sdata, intensity_matrix=None, cell_boundaries=None, index_name='cell_id',
               output_dir='.', comm_detect_res=0.1, max_edge_distance=50, pixel_size_um=None):
    """Delaunay network + Louvain communities on the cell layer.

    ``max_edge_distance`` is in µm and converted to pixels with
    ``pixel_size_um`` (µm/px); if that is None, the value stored in the table's
    ``.uns`` is used (see :func:`cellsurvey.utils.store_pixel_size`).
    """
    if cell_boundaries is None:
        cell_boundaries = resolve_shape_keys(sdata)[1]
    pixel_size_um = resolve_pixel_size(sdata, pixel_size_um)
    max_edge_distance_px = max_edge_distance / pixel_size_um
    boundaries = sdata.shapes[cell_boundaries]
    boundaries.index.name = index_name

    centroids = boundaries.geometry.centroid
    coords = np.column_stack([centroids.x.values, centroids.y.values])

    table = sdata.tables['table']
    n_cells = len(boundaries)

    cluster_labels = np.full(n_cells, -1, dtype=int)
    for idx, cell_id in enumerate(boundaries.index.values):
        if cell_id in table.obs.index:
            cluster_labels[idx] = int(table.obs.loc[cell_id, 'kmeans_cluster'])

    print("\nBuilding Delaunay triangulation...")
    imat = None
    if intensity_matrix is not None:
        cell_ids = boundaries.index.values
        imat = np.zeros((n_cells, intensity_matrix.shape[1]))
        for pos_idx, cell_id in enumerate(cell_ids):
            if cell_id in intensity_matrix.index:
                imat[pos_idx] = intensity_matrix.loc[cell_id].values

    community_labels, n_edges = compute_louvain_communities(
        coords, imat, comm_detect_res=comm_detect_res,
        max_edge_distance=max_edge_distance_px, seed=42)
    print(f"Delaunay network: {n_edges} edges (filtered to max distance "
          f"{max_edge_distance} µm = {max_edge_distance_px:.1f} px)")

    n_communities = len(set(community_labels)) - (1 if -1 in community_labels else 0)
    print(f"Found {n_communities} communities")

    # Embed community and cluster labels into the shapes GeoDataFrame
    boundaries['kmeans_cluster'] = cluster_labels
    boundaries['community'] = community_labels

    # Write summary statistics
    stats = {
        'n_cells': int(n_cells),
        'n_clusters': int(len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)),
        'n_communities': int(n_communities),
        'n_edges': int(n_edges),
        'max_edge_distance_um': max_edge_distance,
        'max_edge_distance_px': float(max_edge_distance_px),
        'pixel_size_um': pixel_size_um,
        'community_resolution': comm_detect_res,
        'cluster_sizes': {int(k): int(v) for k, v in zip(*np.unique(cluster_labels[cluster_labels >= 0], return_counts=True))},
        'community_sizes': {int(k): int(v) for k, v in zip(*np.unique(community_labels[community_labels >= 0], return_counts=True))},
    }
    with open(os.path.join(output_dir, 'summary.json'), 'w') as f:
        json.dump(stats, f, indent=2)
    print(f"\nSummary statistics written to {os.path.join(output_dir, 'summary.json')}")

    return {
        'cell_ids': boundaries.index.values,
        'community_labels': community_labels,
        'cluster_labels': cluster_labels,
    }
