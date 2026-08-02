import os
import json
import numpy as np
from scipy.spatial import Delaunay
import networkx as nx


def run_muspan(sdata, intensity_matrix=None, cell_boundaries='stardist_boundaries', index_name='cell_id',
               output_dir='.', cell_colour='table: kmeans_cluster', comm_detect_res=0.1,
               max_edge_distance=1000, fig_size=20):
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
    tri = Delaunay(coords)
    edges = set()
    for simplex in tri.simplices:
        for i in range(3):
            for j in range(i + 1, 3):
                a, b = simplex[i], simplex[j]
                dist = np.linalg.norm(coords[a] - coords[b])
                if dist <= max_edge_distance:
                    edges.add((a, b))

    print(f"Delaunay network: {len(edges)} edges (filtered to max distance {max_edge_distance})")

    print("\nDetecting Louvain communities...")
    G = nx.Graph()
    G.add_nodes_from(range(n_cells))

    if intensity_matrix is not None:
        cell_ids = boundaries.index.values
        imat = np.zeros((n_cells, intensity_matrix.shape[1]))
        for pos_idx, cell_id in enumerate(cell_ids):
            if cell_id in intensity_matrix.index:
                imat[pos_idx] = intensity_matrix.loc[cell_id].values
        print("Computing edge weights from expression similarity...")
        for a, b in edges:
            corr = np.corrcoef(imat[a], imat[b])[0, 1]
            weight = 1 + corr
            G.add_edge(a, b, weight=weight)
    else:
        G.add_edges_from(edges)

    communities = nx.community.louvain_communities(G, weight='weight' if intensity_matrix is not None else None,
                                                  resolution=comm_detect_res, seed=42)

    community_labels = np.full(n_cells, -1, dtype=int)
    for comm_idx, comm in enumerate(communities):
        for node in comm:
            community_labels[node] = comm_idx

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
        'n_edges': len(edges),
        'max_edge_distance': max_edge_distance,
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
