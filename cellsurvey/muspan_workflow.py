import os
import numpy as np
import muspan as ms
import matplotlib.pyplot as plt

fig_kwargs = dict(figsize=(20, 20))


def get_colors_for_communities(n_communities):
    if n_communities <= 10:
        cmap = plt.cm.get_cmap('tab10')
    elif n_communities <= 20:
        cmap = plt.cm.get_cmap('tab20')
    else:
        cmap = plt.cm.get_cmap('hsv')

    colors = []
    for i in range(n_communities):
        rgba = cmap(i / max(n_communities - 1, 1))
        rgb = [int(rgba[0] * 255), int(rgba[1] * 255), int(rgba[2] * 255)]
        colors.append(rgb)

    return colors


def run_muspan(sdata, cell_boundaries='stardist_boundaries', index_name='cell_id', output_dir='.',
               cell_colour='table: kmeans_cluster', comm_detect_res=0.1, max_edge_distance=1000):
    boundaries = sdata.shapes[cell_boundaries]
    boundaries.index.name = index_name

    centroids = boundaries.geometry.centroid
    coords = np.column_stack([centroids.x.values, centroids.y.values])

    # Build MuSpAn domain directly with just points
    muspan_domain = ms.domain('cell_domain')
    point_ids = muspan_domain.add_points(points=coords, return_IDs=True)

    # Map cell IDs to MuSpAn point IDs (for later use by export)
    cell_id_to_point = dict(zip(boundaries.index.values, map(str, point_ids)))

    # Attach cluster labels for visualisation using a simple dict,
    # keyed by point index to avoid MuSpAn's add_labels internals
    table = sdata.tables['table']
    cluster_map = {}
    for idx, cell_id in enumerate(boundaries.index.values):
        if cell_id in table.obs.index:
            cluster_val = int(table.obs.loc[cell_id, 'kmeans_cluster'])
        else:
            cluster_val = -1
        cluster_map[idx] = float(cluster_val)

    muspan_domain.add_labels(label_name='kmeans_cluster', labels=cluster_map,
                              add_labels_to=None, label_type='categorical')

    print("\nVisualising cells...")
    ms.visualise.visualise(muspan_domain, color_by='kmeans_cluster', marker_size=3.0, figure_kwargs=fig_kwargs)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cells.png'))

    print("\nGenerating Delauney network...")
    ms.networks.generate_network(
        muspan_domain,
        network_name='Centroid Delaunay',
        network_type='Delaunay',
        max_edge_distance=max_edge_distance
    )

    print("\nVisualising Delauney network...")
    ms.visualise.visualise_network(
        muspan_domain,
        network_name='Centroid Delaunay',
        edge_width=0.25,
        visualise_kwargs=dict(color_by='kmeans_cluster', marker_size=1.0),
        figure_kwargs=fig_kwargs
    )
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'delaunay_network.png'))

    print(f'\nDetecting communities at res {comm_detect_res}...')
    communities_res_1 = ms.networks.community_detection(
        muspan_domain,
        network_name='Centroid Delaunay',
        edge_weight_name=None,
        community_method='louvain',
        community_method_parameters=dict(resolution=comm_detect_res),
        community_label_name='Communities'
    )

    print(f'\nVisualising communities at res {comm_detect_res}...')
    ms.visualise.visualise(muspan_domain, color_by='Communities', marker_size=5.0, figure_kwargs=fig_kwargs)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'communities_network.png'))

    # Store mappings on domain for export use
    muspan_domain._cell_id_map = cell_id_to_point
    muspan_domain._boundaries = boundaries

    return muspan_domain
