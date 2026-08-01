import ctypes.util
import pathlib as _pl

# Preload the pixi environment's libstdc++ to avoid ABI conflicts with the
# system library. Locate it relative to this script's pixi env directory.
_pixi_env_lib = _pl.Path(__file__).resolve().parent.parent / ".pixi" / "envs" / "default" / "lib"
_libstdcpp = _pixi_env_lib / "libstdc++.so.6"
if _libstdcpp.exists():
    ctypes.CDLL(str(_libstdcpp))

import argparse
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import spatialdata
import sopa
import tensorflow as tf
from bioio import BioImage
from dask_image.imread import imread
from spatialdata.models import PointsModel
from spatialdata.transformations import Identity

from cellsurvey.utils import remove_channel_suffix, cluster_data, assign_spots_to_cells
from cellsurvey.blob_detection import detect_blobs_tiled
from cellsurvey.network_analysis import run_muspan
from cellsurvey.export import export_to_qupath


def main():
    plt.rcParams['font.size'] = 20
    plt.rcParams['axes.linewidth'] = 3

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', help='Path to input image', required=True)
    parser.add_argument('-o', '--output_file', help='Path to output Zarr', required=True)
    parser.add_argument('-p', '--plot_dir', help='Output directory for data plots', default='.')
    parser.add_argument('--detect-blobs', action='store_true',
                        help='Enable RNA spot blob detection on the specified channels')
    parser.add_argument('--use-gpu', action='store_true',
                        help='Force GPU usage for Stardist segmentation (auto-detected if not specified)')
    parser.add_argument('--channels', help='Comma-separated channel indices for blob detection',
                        default='9,10,11,12')
    parser.add_argument('--thresholds', help='Comma-separated blob detection thresholds (one per channel)',
                        default='0.01,0.1,0.1,0.1')
    parser.add_argument('--tile-size', type=int, default=2048,
                        help='Tile size for blob detection (default: 2048)')
    parser.add_argument('--overlap', type=int, default=50,
                        help='Tile overlap for blob detection (default: 50)')
    parser.add_argument('--workers', type=int, default=14,
                        help='Number of worker threads for blob detection (default: 14)')
    parser.add_argument('--n-clusters', type=int, default=10,
                        help='Number of clusters for k-means (default: 10)')
    parser.add_argument('--community-resolution', type=float, default=0.1,
                        help='Resolution for Louvain community detection (default: 0.1)')
    parser.add_argument('--max-edge-distance', type=float, default=1000,
                        help='Max edge distance for Delaunay network (default: 1000)')
    parser.add_argument('--radius-min', type=float, default=0,
                        help='Min radius for spatial neighbors graph (default: 0)')
    parser.add_argument('--radius-max', type=float, default=1000,
                        help='Max radius for spatial neighbors graph (default: 1000)')
    parser.add_argument('--resume-from',
                        help='Path to existing Zarr to resume from (skips image loading and spot detection)')
    args = parser.parse_args()

    os.makedirs(args.plot_dir, exist_ok=True)

    imagepath = args.input_file
    zarr_path = args.output_file
    if not zarr_path.endswith('.zarr'):
        zarr_path += '.zarr'

    if args.resume_from:
        orig_zarr_path = args.resume_from
        if not orig_zarr_path.endswith('.zarr'):
            orig_zarr_path += '.zarr'
        if not os.path.exists(orig_zarr_path):
            parser.error(f"Resume Zarr not found: {orig_zarr_path}")
        print(f"Resuming from existing Zarr: {orig_zarr_path}")

        channel_names = BioImage(imagepath).channel_names
    elif args.detect_blobs:
        channel_index = list(map(int, args.channels.split(',')))
        thresholds = list(map(float, args.thresholds.split(',')))

        if len(channel_index) != len(thresholds):
            parser.error(
                f"Number of channels ({len(channel_index)}) must match number of thresholds ({len(thresholds)})")

        channel_names = BioImage(imagepath).channel_names
        image = imread(imagepath)[channel_index]

        dataset = sopa.io.ome_tif(imagepath, as_image=False)

        frames = []

        for i, (ci, t) in enumerate(zip(channel_index, thresholds)):
            ch_name = channel_names[ci]
            print(f'Finding blobs in channel {ch_name} (index {ci})...')
            blobs = detect_blobs_tiled(image[i], tile_size=args.tile_size, overlap=args.overlap, threshold=t,
                                       n_workers=args.workers)

            frames.append(pd.DataFrame({
                "x": blobs[:, 1],
                "y": blobs[:, 0],
                "sigma": blobs[:, 2],
                "gene": np.array([ch_name] * len(blobs)),
            }))

        df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["x", "y", "sigma", "gene"])

        points = PointsModel.parse(
            df,
            coordinates={"x": "x", "y": "y"},
            transformations={"global": Identity()},
        )

        dataset["spots"] = points

        orig_zarr_path = zarr_path

        print(f'Saving Zarr to {orig_zarr_path}')

        try:
            dataset.write(orig_zarr_path, overwrite=True)
        except Exception as e:
            print(f"ERROR: Failed to write Zarr to {orig_zarr_path}: {e}", file=sys.stderr)
            sys.exit(1)

        print("Done")
    else:
        channel_names = BioImage(imagepath).channel_names
        dataset = sopa.io.ome_tif(imagepath, as_image=False)
        orig_zarr_path = zarr_path

        print(f'Saving Zarr to {orig_zarr_path}')

        try:
            dataset.write(orig_zarr_path, overwrite=True)
        except Exception as e:
            print(f"ERROR: Failed to write Zarr to {orig_zarr_path}: {e}", file=sys.stderr)
            sys.exit(1)

        print("Done")

    seg_zarr_path = zarr_path.replace('.zarr', '_seg.zarr')
    needs_stardist = True
    needs_aggregation = True

    if os.path.exists(seg_zarr_path):
        try:
            dataset = spatialdata.read_zarr(seg_zarr_path)
            if 'table' in dataset.tables and dataset.tables['table'] is not None:
                print(f"Segmented Zarr with table found, skipping to clustering: {seg_zarr_path}")
                needs_stardist = False
                needs_aggregation = False
                sdata = dataset
            elif 'stardist_boundaries' in dataset.shapes:
                print(f"Segmented Zarr has boundaries but no table, re-running aggregation only")
                needs_stardist = False
                # dataset is loaded, proceed to aggregation
        except Exception:
            import shutil
            shutil.rmtree(seg_zarr_path, ignore_errors=True)

    if needs_stardist:
        print("Loading Zarr...")

        dataset = spatialdata.read_zarr(orig_zarr_path)

        image_name = list(dataset.images.keys())[0]

        print("Make image patches...")

        sopa.make_image_patches(dataset)

        print("Set backend to None (will use GPU)...")

        gpus = tf.config.list_physical_devices('GPU')
        if args.use_gpu or gpus:
            if gpus:
                print(f"Found {len(gpus)} GPU(s), using GPU backend for Stardist")
            else:
                print("No GPU detected but --use-gpu specified, attempting GPU backend")
            sopa.settings.parallelization_backend = None
        else:
            print("WARNING: No GPU detected. Stardist segmentation will run on CPU and may be very slow.")
            sopa.settings.parallelization_backend = None

        print("Get channel names...")

        channels = sopa.utils.get_channel_names(dataset)

        print(channels)

        unique_channels = [f"{ch}_ch_{i}" for i, ch in enumerate(channels)]
        print("Fixed channels:", unique_channels)

        for scale_name in dataset.images[image_name].children:
            scale_node = dataset.images[image_name][scale_name]
            scale_node.ds = scale_node.ds.assign_coords(c=unique_channels)

        print("Run stardist...")

        sopa.segmentation.stardist(dataset, model_type='2D_versatile_fluo', channels=unique_channels[0])

    if needs_aggregation:
        print("Aggregating...")

        # Force pandas to use plain object dtype for strings, not ArrowStringArray,
        # which can't be written to Zarr backing stores by anndata
        with pd.option_context('future.infer_string', False):
            if "spots" in dataset.points:
                sopa.aggregate(dataset, aggregate_genes=True, points_key='spots', gene_column='gene')
            else:
                sopa.aggregate(dataset)

        seg_zarr_path = zarr_path.replace('.zarr', '_seg.zarr')

        print(f'Saving Zarr to {seg_zarr_path}')

        try:
            dataset.write(seg_zarr_path, overwrite=True)
        except Exception as e:
            print(f"ERROR: Failed to write segmented Zarr to {seg_zarr_path}: {e}", file=sys.stderr)
            sys.exit(1)

        print("Done")

        np.random.seed(42)

        path_to_spatialData_file = seg_zarr_path
        sdata = spatialdata.read_zarr(path_to_spatialData_file)

    measurements = sdata.tables['table']

    X_data = measurements.X.toarray() if hasattr(measurements.X, 'toarray') else measurements.X

    intensity_df = pd.DataFrame(
        X_data,
        index=measurements.obs.index,
        columns=measurements.var.index
    )

    intensity_df.columns = [remove_channel_suffix(col) for col in intensity_df.columns]

    intensity_df = intensity_df.loc[:, ~intensity_df.columns.duplicated(keep='first')]
    print(f"Columns after removing duplicates: {list(intensity_df.columns)}")
    print(f"Number of cells: {len(intensity_df)}")

    cluster_labels = cluster_data(intensity_df, n_clusters=args.n_clusters)

    sdata.tables['table'].obs['kmeans_cluster'] = pd.Categorical(cluster_labels)

    sdata.tables['table'].obs['kmeans_cluster_label'] = pd.Categorical(
        [f'Cluster_{i}' for i in cluster_labels]
    )

    # Force obs index to plain string dtype to avoid ArrowStringArray
    # write errors with Zarr backing stores
    sdata.tables['table'].obs.index = sdata.tables['table'].obs.index.astype(str)

    print("\n\u2713 Cluster labels added to sdata.tables['table'].obs['kmeans_cluster']")
    print("\u2713 Cluster labels added to sdata.tables['table'].obs['kmeans_cluster_label']")

    print(f"\nUpdated obs columns: {list(sdata.tables['table'].obs.columns)}")

    result = run_muspan(sdata, comm_detect_res=args.community_resolution,
                                max_edge_distance=args.max_edge_distance)

    spots_with_cells = assign_spots_to_cells(sdata)

    export_to_qupath(result['cell_ids'], result['community_labels'], result['cluster_labels'],
                     output_path='./qupath_export.geojson',
                     sdata=sdata, intensity_df=intensity_df,
                     spots_with_cells=spots_with_cells)

    adata = sdata.tables['table']
    sopa.spatial.spatial_neighbors(adata, radius=(args.radius_min, args.radius_max))
    cell_type_to_cell_type = sopa.spatial.mean_distance(adata, "kmeans_cluster", "kmeans_cluster")
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.linewidth'] = 2
    plt.figure(figsize=(10, 10))
    heatmap_kwargs = {"cmap": sns.cm.rocket_r, "cbar_kws": {"label": "Mean hop distance"}}
    sns.heatmap(cell_type_to_cell_type, **heatmap_kwargs)
    plt.tight_layout()
    plt.savefig(os.path.join(args.plot_dir, 'cell_type_to_cell_type.png'))

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, flavor='igraph', n_iterations=2, directed=False)
    sc.pl.umap(adata, color="kmeans_cluster", show=False)
    plt.savefig(os.path.join(args.plot_dir, 'umap_kmeans_cluster.png'))
    plt.close()
    sc.pl.umap(adata, color="leiden", show=False)
    plt.savefig(os.path.join(args.plot_dir, 'umap_leiden.png'))
    plt.close()


if __name__ == '__main__':
    main()
