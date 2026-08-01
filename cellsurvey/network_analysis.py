import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
from scipy.spatial.distance import pdist, squareform
import networkx as nx

fig_kwargs = dict(figsize=(20, 20))


def run_muspan(sdata, intensity_matrix=None, cell_boundaries='stardist_boundaries', index_name='cell_id',
               output_dir='.', cell_colour='table: kmeans_cluster', comm_detect_res=0.1,
               max_edge_distance=1000):
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

    print("\nGenerating cells plot...")
    fig, ax = plt.subplots(**fig_kwargs)
    scatter = ax.scatter(coords[:, 0], coords[:, 1], c=cluster_labels, s=3.0, cmap='tab10')
    ax.set_title("Cells by k-means cluster")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cells.png'))
    plt.close()

    print("\nGenerating Delaunay network plot...")
    fig, ax = plt.subplots(**fig_kwargs)
    ax.scatter(coords[:, 0], coords[:, 1], c=cluster_labels, s=1.0, cmap='tab10')
    edge_coords = np.array([[coords[a], coords[b]] for a, b in edges])
    for a, b in edge_coords:
        ax.plot([a[0], b[0]], [a[1], b[1]], 'k-', linewidth=0.25, alpha=0.3)
    ax.set_title(f"Delaunay network ({len(edges)} edges)")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'delaunay_network.png'))
    plt.close()

    print("\nGenerating communities plot...")
    fig, ax = plt.subplots(**fig_kwargs)
    ax.scatter(coords[:, 0], coords[:, 1], c=community_labels, s=5.0, cmap='tab20')
    ax.set_title(f"Louvain communities ({n_communities} communities)")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'communities_network.png'))
    plt.close()

    return {
        'cell_ids': boundaries.index.values,
        'community_labels': community_labels,
        'cluster_labels': cluster_labels,
    }
