import ctypes.util

ctypes.CDLL("/nemo/stp/lm/working/barryd/hpc/pixi/sopaspan/.pixi/envs/sopaspan/lib/libstdc++.so.6")

import sopa
import argparse
import spatialdata
import re
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import muspan as ms
import json
import os
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from shapely.geometry import mapping
import numpy as np
from bioio import BioImage
from spotiflow.model.spotiflow import Spotiflow as sp
from spatialdata.models import PointsModel
from spatialdata.transformations import Identity
from dask_image.imread import imread
import geopandas as gpd
from skimage.feature import blob_log

plt.rcParams['font.size'] = 20
plt.rcParams['axes.linewidth'] = 3
fig_kwargs = dict(figsize=(20, 20))


# Strip the _ch_N suffix from column names
def remove_channel_suffix(name):
    return re.sub(r'_ch_\d+$', '', name)


def cluster_data(data, n_clusters=10, random_seed=42):
    # Perform k-means clustering
    print("\nPerforming k-means clustering...")

    # Standardize the data (important for clustering)
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(data)

    kmeans = KMeans(n_clusters=n_clusters, random_state=random_seed)
    cluster_labels = kmeans.fit_predict(scaled_data)

    print(f"Clustering complete. Found {n_clusters} clusters.")
    print(f"Cluster distribution:")
    for i in range(n_clusters):
        count = np.sum(cluster_labels == i)
        print(f"  Cluster {i}: {count} cells ({100 * count / len(cluster_labels):.1f}%)")

    return cluster_labels


def get_colors_for_communities(n_communities):
    """Generate distinct colors for communities."""
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


def run_muspan(spatial_data, cell_boundaries='stardist_boundaries', index_name='cell_id', output_dir='.',
               cell_colour='table: kmeans_cluster', comm_detect_res=0.1):
    # Set the index name (as we learned earlier)
    spatial_data.shapes[cell_boundaries].index.name = index_name

    # Create clean version without image_patches
    sdata_clean = spatialdata.SpatialData(
        shapes={cell_boundaries: spatial_data.shapes[cell_boundaries]},
        tables=spatial_data.tables
    )

    # Convert to muspan domain
    muspan_domain = ms.io.spatialdata_to_domain(sdata_clean, import_shapes_as_points=True)

    print("\nVisualising cells...")
    ms.visualise.visualise(muspan_domain, color_by=cell_colour, marker_size=3.0, figure_kwargs=fig_kwargs)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'cells.png'))

    print("\nGenerating Delauney network...")
    ms.networks.generate_network(
        muspan_domain,
        network_name='Centroid Delaunay',
        network_type='Delaunay',
        max_edge_distance=1000
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


def export_to_qupath(domain, communities, clusters, output_path, cell_id='cell_id', spots_with_cells=None):
    print("Exporting to QuPath GeoJSON format...")

    # Get cell IDs and community labels from the muspan domain
    cell_ids = domain.labels[cell_id]['labels']
    communities_labels = domain.labels[communities]['labels']
    cluster_labels_from_domain = domain.labels[clusters]['labels']

    # Create a mapping from cell_id to community labels
    cell_to_community = dict(zip(cell_ids, communities_labels))
    cell_to_cluster = dict(zip(cell_ids, cluster_labels_from_domain))

    # Get unique communities for coloring
    unique_communities = np.unique(communities_labels)
    community_colors = get_colors_for_communities(len(unique_communities))

    features = []

    print("Exporting cell boundaries with community labels...")
    boundaries = sdata.shapes['stardist_boundaries']
    for idx, row in boundaries.iterrows():
        # Get community assignments
        community = cell_to_community.get(idx, -1)
        cluster_id = cell_to_cluster.get(idx, -1)

        # Get color based on community
        if 0 <= community < len(community_colors):
            color = community_colors[int(community)]
            classification_name = f"Community_{int(community)}"
        else:
            color = [128, 128, 128]  # Gray for unassigned
            classification_name = "Unassigned"

        # Build measurements dictionary
        measurements_dict = {
            "Community ID": float(community),
            "Cluster ID": float(cluster_id)
        }

        # Add intensity measurements if available
        if idx in intensity_df.index:
            for col in intensity_df.columns:
                if col != 'Cluster':
                    try:
                        val = intensity_df.loc[idx, col]
                        if pd.notna(val):
                            measurements_dict[f"Cell: {col} mean"] = float(val)
                    except (ValueError, TypeError, KeyError):
                        pass

        feature = {
            "type": "Feature",
            "id": str(idx),
            "geometry": mapping(row.geometry),
            "properties": {
                "objectType": "cell",
                "classification": {
                    "name": classification_name,
                    "color": color
                },
                "measurements": measurements_dict
            }
        }
        features.append(feature)

    print(f"Exported {len(boundaries)} cell boundaries")

    # Export spots as point detections with parent cell assignment
    if spots_with_cells is not None:
        print("Exporting spot detections with cell assignments...")
        for idx, row in spots_with_cells.iterrows():
            cell_id_val = row['cell_id']  # already a string hash or None

            feature = {
                "type": "Feature",
                "geometry": {
                    "type": "Point",
                    "coordinates": [float(row['x']), float(row['y'])]
                },
                "properties": {
                    "objectType": "detection",
                    "parent_id": str(cell_id_val) if cell_id_val is not None else None,
                    "classification": {
                        "name": row['gene'],
                        "color": [0, 255, 0]  # Green for spots
                    },
                    "measurements": {
                        "prob": float(row['prob']),
                    }
                }
            }
            features.append(feature)

        assigned = spots_with_cells['cell_id'].notna().sum()
        print(f"Exported {len(spots_with_cells)} spots ({assigned} assigned to cells, "
              f"{len(spots_with_cells) - assigned} unassigned)")

    # Create final GeoJSON
    geojson = {
        "type": "FeatureCollection",
        "features": features
    }

    # Write to file
    with open(output_path, 'w') as f:
        json.dump(geojson, f, indent=2)

    print(f"\n✓ Exported {len(features)} total features")
    print(f"Output file: {output_path}")


def detect_blobs_tiled(image, tile_size=1024, overlap=50, min_sigma=1, max_sigma=5, num_sigma=5, threshold=0.1):
    print(f"Detecting blobs in tiles of {tile_size}x{tile_size} with {overlap}px overlap...")

    height, width = image.shape
    all_blobs = []

    # Calculate tile starts
    y_starts = list(range(0, height, tile_size - overlap))
    x_starts = list(range(0, width, tile_size - overlap))

    total_tiles = len(y_starts) * len(x_starts)
    tile_count = 0

    # Compute once outside the loop for efficiency
    image_min = image.min()
    image_max = image.max()

    for y in y_starts:
        for x in x_starts:
            tile_count += 1

            # Calculate tile boundaries (clamp to image size)
            y_end = min(y + tile_size, height)
            x_end = min(x + tile_size, width)

            # Extract and normalise tile
            tile = image[y:y_end, x:x_end].astype(float)
            tile = (tile - image_min) / (image_max - image_min + 1e-8)
            print(f'Finding blobs in tile {tile.shape}')
            # Detect blobs in tile
            blobs = blob_log(
                tile,
                min_sigma=min_sigma,
                max_sigma=max_sigma,
                num_sigma=num_sigma,
                threshold=threshold
            )

            if len(blobs) > 0:
                # Convert tile coordinates back to global image coordinates
                blobs[:, 0] += y
                blobs[:, 1] += x

                # Remove blobs detected in overlap region (keep only those in the
                # non-overlapping part of the tile, except at image edges)
                half_overlap = overlap // 2
                y_min_keep = y + half_overlap if y > 0 else y
                x_min_keep = x + half_overlap if x > 0 else x
                y_max_keep = y_end - half_overlap if y_end < height else y_end
                x_max_keep = x_end - half_overlap if x_end < width else x_end

                mask = (
                        (blobs[:, 0] >= y_min_keep) & (blobs[:, 0] < y_max_keep) &
                        (blobs[:, 1] >= x_min_keep) & (blobs[:, 1] < x_max_keep)
                )
                blobs = blobs[mask]
                all_blobs.append(blobs)

            if tile_count % 10 == 0:
                print(f"  Processed {tile_count}/{total_tiles} tiles...")

    if all_blobs:
        all_blobs = np.concatenate(all_blobs, axis=0)
    else:
        all_blobs = np.zeros((0, 3))

    print(f"Found {len(all_blobs)} blobs total")
    return all_blobs


def assign_spots_to_cells(spatial_data, spots_key='spots', cell_boundaries='stardist_boundaries'):
    print("Computing spot-to-cell assignments...")

    spots_df = spatial_data.points[spots_key].compute()
    spots_gdf = gpd.GeoDataFrame(
        spots_df,
        geometry=gpd.points_from_xy(spots_df['x'], spots_df['y']),
        crs=spatial_data.shapes[cell_boundaries].crs
    )

    spots_with_cells = gpd.sjoin(
        spots_gdf,
        spatial_data.shapes[cell_boundaries][['geometry']],
        how='left',
        predicate='within'
    )

    # cell_id is a string hash - keep as string, use None for unassigned
    spots_with_cells['cell_id'] = spots_with_cells['cell_id'].where(
        spots_with_cells['cell_id'].notna(), other=None
    )

    assigned = spots_with_cells['cell_id'].notna().sum()
    total = len(spots_with_cells)
    print(f"Assigned {assigned} of {total} spots to cells ({100 * assigned / total:.1f}%)")
    print(f"Unassigned spots: {total - assigned} ({100 * (total - assigned) / total:.1f}%)")

    return spots_with_cells


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', help='Path to input image',
                        default='/nemo/stp/lm/working/barryd/hpc/projects/stps/ehp/2026.01/comet_lunaphore/data/test_crop_2.ome.tiff')
    parser.add_argument('-o', '--output_file', help='Path to output Zarr',
                        default='/nemo/stp/lm/working/barryd/hpc/projects/stps/ehp/2026.01/comet_lunaphore/data/test_crop_2')
    parser.add_argument('-p', '--plot_dir', help='Output directory for data plots', default='.')
    args = parser.parse_args()

    imagepath = args.input_file
    zarr_path = args.output_file

    channel_index = 9
    # Get a BioImage object
    img = BioImage(imagepath)  # selects the first scene found
    channel_names = img.channel_names

    dataset = sopa.io.ome_tif(imagepath, as_image=False)

    image = imread(imagepath)

    # result = sp().predict(img=image[channel_index], device="cuda", prob_thresh=0.5, subpix=True, peak_mode='skimage')
    # spots = result[0]  # np.ndarray, shape (N, 2), order (row, col)
    # details = result[1]
    # probs = details.prob  # np.ndarray, shape (N,)
    #
    # # spots is (N, 2) in (row, col) order, so row=y, col=x
    # df = pd.DataFrame({
    #     "x": spots[:, 1],
    #     "y": spots[:, 0],
    #     "prob": probs,
    #     "gene": np.array([channel_names[channel_index]] * len(spots)),
    # })

    blobs = detect_blobs_tiled(image[channel_index], tile_size=1024, overlap=50, threshold=0.01)

    df = pd.DataFrame({
        "x": blobs[:, 1],
        "y": blobs[:, 0],
        "sigma": blobs[:, 2],  # sigma as proxy for confidence
        "gene": np.array([channel_names[channel_index]] * len(blobs)),
    })

    points = PointsModel.parse(
        df,
        coordinates={"x": "x", "y": "y"},
        transformations={"global": Identity()},
    )

    dataset["spots"] = points

    orig_zarr_path = f'{zarr_path}.zarr'

    print(f'Saving Zarr to {orig_zarr_path}')

    dataset.write(f'{orig_zarr_path}')

    print("Done")

    print("Loading Zarr...")

    dataset = spatialdata.read_zarr(orig_zarr_path)  # we can read the data back

    image_name = list(dataset.images.keys())[0]

    print("Make image patches...")

    sopa.make_image_patches(dataset)

    print("Set backend to None (will use GPU)...")

    sopa.settings.parallelization_backend = None

    print("Get channel names...")

    channels = sopa.utils.get_channel_names(dataset)

    print(channels)

    unique_channels = [f"{ch}_ch_{i}" for i, ch in enumerate(channels)]
    print("Fixed channels:", unique_channels)

    for scale_name in dataset.images[image_name].children:
        scale_node = dataset.images[image_name][scale_name]
        # Update the dataset with new coordinates
        scale_node.ds = scale_node.ds.assign_coords(c=unique_channels)

    print("Run stardist...")

    sopa.segmentation.stardist(dataset, model_type='2D_versatile_fluo', channels=unique_channels[0])

    print("Aggregating...")

    sopa.aggregate(dataset, aggregate_genes=True, points_key='spots', gene_column='gene')

    seg_zarr_path = f'{zarr_path}_seg.zarr'

    print(f'Saving Zarr to {seg_zarr_path}')

    dataset.write(f'{seg_zarr_path}', overwrite=True)

    print("Done")

    # Set random seed for reproducibility
    np.random.seed(42)

    # Load the data
    path_to_spatialData_file = seg_zarr_path
    sdata = spatialdata.read_zarr(path_to_spatialData_file)

    # Get the boundaries and measurements
    measurements = sdata.tables['table']

    # Convert AnnData to DataFrame
    intensity_df = pd.DataFrame(
        measurements.X.toarray(),
        index=measurements.obs.index,
        columns=measurements.var.index
    )

    intensity_df.columns = [remove_channel_suffix(col) for col in intensity_df.columns]

    # Keep only first occurrence of duplicate columns
    intensity_df = intensity_df.loc[:, ~intensity_df.columns.duplicated(keep='first')]
    print(f"Columns after removing duplicates: {list(intensity_df.columns)}")
    print(f"Number of cells: {len(intensity_df)}")

    # ===== ADD CLUSTERING RESULTS TO SPATIALDATA =====
    cluster_labels = cluster_data(intensity_df)

    # Add cluster labels to the AnnData obs (as categorical for efficiency)
    sdata.tables['table'].obs['kmeans_cluster'] = pd.Categorical(cluster_labels)

    # Optionally, add the cluster labels as a string for better visualization
    sdata.tables['table'].obs['kmeans_cluster_label'] = pd.Categorical(
        [f'Cluster_{i}' for i in cluster_labels]
    )

    print("\n✓ Cluster labels added to sdata.tables['table'].obs['kmeans_cluster']")
    print("✓ Cluster labels added to sdata.tables['table'].obs['kmeans_cluster_label']")

    # Verify it was added
    print(f"\nUpdated obs columns: {list(sdata.tables['table'].obs.columns)}")

    example_domain = run_muspan(sdata)

    spots_with_cells = assign_spots_to_cells(sdata)

    export_to_qupath(example_domain, 'Communities', 'table: kmeans_cluster', output_path='./qupath_export.geojson',
                     spots_with_cells=spots_with_cells)

    adata = sdata.tables['table']
    #    idx = np.random.choice(adata.n_obs, size=30000, replace=False)
    adata_subset = adata
    sopa.spatial.spatial_neighbors(adata, radius=(0, 1000))
    cell_type_to_cell_type = sopa.spatial.mean_distance(adata, "kmeans_cluster", "kmeans_cluster")
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.linewidth'] = 2
    plt.figure(figsize=(10, 10))
    heatmap_kwargs = {"cmap": sns.cm.rocket_r, "cbar_kws": {"label": "Mean hop distance"}}
    sns.heatmap(cell_type_to_cell_type, **heatmap_kwargs)
    plt.tight_layout()
    plt.savefig(os.path.join(args.plot_dir, 'cell_type_to_cell_type.png'))

    sc.pp.normalize_total(adata_subset)
    sc.pp.log1p(adata_subset)
    sc.pp.neighbors(adata_subset)
    sc.tl.umap(adata_subset)
    sc.tl.leiden(adata_subset)
    plt.figure(figsize=(10, 10))
    sc.pl.umap(adata_subset, color="kmeans_cluster")
    sc.pl.umap(adata_subset, color="leiden")
