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

    # Build MuSpAn domain manually to avoid spatialdata_to_domain indexing bugs
    muspan_domain = ms.domain
    muspan_domain.add_points(muspan_domain, points=coords, labels=boundaries.index.values, label_name=index_name)

    # Attach table columns as labels
    table = sdata.tables['table']
    for col in table.obs.columns:
        values = dict(zip(table.obs['cell_id'].astype(str), table.obs[col].values))
        muspan_domain.add_labels(muspan_domain, labels=values, label_name=col, label_key=index_name)

    print("\nVisualising cells...")
    ms.visualise.visualise(muspan_domain, color_by=cell_colour, marker_size=3.0, figure_kwargs=fig_kwargs)
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
        visualise_kwargs=dict(color_by=cell_colour, marker_size=1.0),
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
        community_label_name=f'Communities'
    )

    print(f'\nVisualising communities at res {comm_detect_res}...')
    ms.visualise.visualise(muspan_domain, color_by='Communities', marker_size=5.0, figure_kwargs=fig_kwargs)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'communities_network.png'))

    return muspan_domain
