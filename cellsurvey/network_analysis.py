import os
import numpy as np
import matplotlib.pyplot as plt

fig_kwargs = dict(figsize=(20, 20))


def run_muspan(sdata, cell_boundaries='stardist_boundaries', index_name='cell_id', output_dir='.',
               cell_colour='table: kmeans_cluster', comm_detect_res=0.1, max_edge_distance=1000):
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

    community_labels = np.zeros(n_cells, dtype=int)

    print("\nGenerating placeholder cells plot...")
    fig, ax = plt.subplots(**fig_kwargs)
    scatter = ax.scatter(coords[:, 0], coords[:, 1], c=cluster_labels, s=3.0, cmap='tab10')
    ax.set_title("Cells by k-means cluster (placeholder — no network analysis)")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cells.png'))
    plt.close()

    print("\nGenerating placeholder Delauney network plot...")
    fig, ax = plt.subplots(**fig_kwargs)
    ax.scatter(coords[:, 0], coords[:, 1], c=cluster_labels, s=1.0, cmap='tab10')
    ax.set_title("Delaunay network (placeholder — no network analysis)")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'delaunay_network.png'))
    plt.close()

    print("\nGenerating placeholder communities plot...")
    fig, ax = plt.subplots(**fig_kwargs)
    ax.scatter(coords[:, 0], coords[:, 1], c=community_labels, s=5.0, cmap='tab10')
    ax.set_title("Communities (placeholder — all cells in community 0)")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'communities_network.png'))
    plt.close()

    return {
        'cell_ids': boundaries.index.values,
        'community_labels': community_labels,
        'cluster_labels': cluster_labels,
    }
