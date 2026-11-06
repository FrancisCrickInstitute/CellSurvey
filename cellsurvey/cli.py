import ctypes.util
import pathlib as _pl
import sys

# Preload the active environment's libstdc++ (e.g. .pixi/envs/<env>/lib) to
# avoid ABI conflicts with the system library. sys.prefix is the running env,
# so this works for both the default (Stardist) and instanseg pixi envs.
_libstdcpp = _pl.Path(sys.prefix) / "lib" / "libstdc++.so.6"
if _libstdcpp.exists():
    ctypes.CDLL(str(_libstdcpp))

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import spatialdata
import sopa
from bioio import BioImage
from dask_image.imread import imread
from spatialdata.models import PointsModel
from spatialdata.transformations import Identity

from cellsurvey.utils import (remove_channel_suffix, cluster_data, assign_spots_to_cells,
                              resolve_shape_keys, STARDIST_NUCLEI_KEY, NUCLEI_KEY, CELLS_KEY)
from cellsurvey.blob_detection import detect_blobs_tiled
from cellsurvey.network_analysis import run_network_analysis
from cellsurvey.segmentation import get_pixel_size, expand_nuclei, run_stardist, run_instanseg
from cellsurvey.export import export_to_qupath


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', help='Path to input image', required=True)
    parser.add_argument('-o', '--output_file', help='Path to output Zarr', required=True)
    parser.add_argument('-p', '--plot_dir', help='Output directory for data plots', default='.')
    parser.add_argument('--detect-blobs', action='store_true',
                        help='Enable RNA spot blob detection on the specified channels')
    parser.add_argument('--segmenter', choices=['stardist', 'instanseg'], default='stardist',
                        help='Segmentation method (default: stardist). instanseg needs the instanseg '
                             'pixi environment: pixi run -e instanseg python run.py ...')
    parser.add_argument('--seg-channels', default=None,
                        help='Comma-separated channel indices to segment on '
                             '(default: 0 for stardist, all channels for instanseg)')
    parser.add_argument('--instanseg-model', default='fluorescence_nuclei_and_cells',
                        help='InstanSeg pretrained model name or path (default: fluorescence_nuclei_and_cells)')
    parser.add_argument('--pixel-size', type=float, default=None,
                        help='Pixel size in µm/px, overriding the image metadata')
    parser.add_argument('--use-gpu', action='store_true',
                        help='Request GPU segmentation (auto-detected if not specified)')
    parser.add_argument('--cell-expansion', type=float, default=5.0,
                        help='Stardist only: expand nuclei into cells by this radius in µm; 0 = nuclei only '
                             '(default: 5.0). Ignored for instanseg, which segments cells directly')
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
    parser.add_argument('--min-sigma', type=float, default=2,
                        help='Minimum blob radius for spot detection (default: 2)')
    parser.add_argument('--max-sigma', type=float, default=5,
                        help='Maximum blob radius for spot detection (default: 5)')
    parser.add_argument('--num-sigma', type=int, default=5,
                        help='Number of sigma steps for blob detection (default: 5)')
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
    parser.add_argument('--geojson-path', default='./qupath_export.geojson',
                        help='Output path for QuPath GeoJSON (default: ./qupath_export.geojson)')
    parser.add_argument('--font-size', type=int, default=20,
                        help='Font size for plots (default: 20)')
    parser.add_argument('--axes-linewidth', type=int, default=3,
                        help='Axes line width for plots (default: 3)')
    args = parser.parse_args()

    plt.rcParams['font.size'] = args.font_size
    plt.rcParams['axes.linewidth'] = args.axes_linewidth

    os.makedirs(args.plot_dir, exist_ok=True)

    imagepath = args.input_file

    if args.pixel_size is not None:
        if args.pixel_size <= 0:
            parser.error("--pixel-size must be positive")
        pixel_size_um = args.pixel_size
        print(f"Physical pixel size (from --pixel-size): {pixel_size_um:.4f} µm/px")
    else:
        pixel_size_um = get_pixel_size(imagepath)
        if pixel_size_um is None:
            print("WARNING: no physical pixel size in image metadata; pass --pixel-size "
                  "(needed for --cell-expansion and InstanSeg rescaling)")
        else:
            print(f"Physical pixel size: {pixel_size_um:.4f} µm/px")

    # Shape layers the chosen segmenter writes; used to decide whether a resumed
    # _seg.zarr already holds this segmenter's output. (cell_boundaries alone is
    # not enough for InstanSeg: Stardist + expansion also writes it.)
    seg_layers = [STARDIST_NUCLEI_KEY] if args.segmenter == 'stardist' else [NUCLEI_KEY, CELLS_KEY]

    zarr_path = args.output_file
    if not zarr_path.endswith('.zarr'):
        zarr_path += '.zarr'

    if args.resume_from:
        zarr_path = args.resume_from
        if not zarr_path.endswith('.zarr'):
            zarr_path += '.zarr'
        if not os.path.exists(zarr_path):
            parser.error(f"Resume Zarr not found: {zarr_path}")
        print(f"Resuming from existing Zarr: {zarr_path}")

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
                                       min_sigma=args.min_sigma, max_sigma=args.max_sigma,
                                       num_sigma=args.num_sigma, n_workers=args.workers)

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

        print(f'Saving Zarr to {zarr_path}')

        try:
            dataset.write(zarr_path, overwrite=True)
        except Exception as e:
            print(f"ERROR: Failed to write Zarr to {zarr_path}: {e}", file=sys.stderr)
            sys.exit(1)

        print("Done")
    else:
        dataset = sopa.io.ome_tif(imagepath, as_image=False)
        print(f'Saving Zarr to {zarr_path}')

        try:
            dataset.write(zarr_path, overwrite=True)
        except Exception as e:
            print(f"ERROR: Failed to write Zarr to {zarr_path}: {e}", file=sys.stderr)
            sys.exit(1)

        print("Done")

    seg_zarr_path = zarr_path.replace('.zarr', '_seg.zarr')
    needs_segmentation = True
    needs_aggregation = True

    if args.resume_from and os.path.exists(seg_zarr_path):
        try:
            dataset = spatialdata.read_zarr(seg_zarr_path)
            if 'table' in dataset.tables and dataset.tables['table'] is not None:
                print(f"Segmented Zarr with table found, skipping to clustering: {seg_zarr_path}")
                needs_segmentation = False
                needs_aggregation = False
                sdata = dataset
            elif all(k in dataset.shapes for k in seg_layers):
                print(f"Segmented Zarr has {', '.join(seg_layers)} but no table, re-running aggregation only")
                needs_segmentation = False
        except Exception:
            import shutil
            shutil.rmtree(seg_zarr_path, ignore_errors=True)

    if needs_segmentation:
        print("Loading Zarr...")

        dataset = spatialdata.read_zarr(zarr_path)

        image_name = list(dataset.images.keys())[0]

        print("Make image patches...")

        sopa.make_image_patches(dataset)

        print("Get channel names...")

        channels = sopa.utils.get_channel_names(dataset)

        print(channels)

        unique_channels = [f"{ch}_ch_{i}" for i, ch in enumerate(channels)]
        print("Fixed channels:", unique_channels)

        for scale_name in dataset.images[image_name].children:
            scale_node = dataset.images[image_name][scale_name]
            scale_node.ds = scale_node.ds.assign_coords(c=unique_channels)

        if args.seg_channels is not None:
            try:
                seg_channels = [unique_channels[int(i)] for i in args.seg_channels.split(',')]
            except (ValueError, IndexError):
                parser.error(f"--seg-channels must be comma-separated indices in 0..{len(unique_channels) - 1}")
        elif args.segmenter == 'stardist':
            seg_channels = [unique_channels[0]]
        else:
            seg_channels = unique_channels

        if args.segmenter == 'stardist':
            if len(seg_channels) != 1:
                parser.error("Stardist segments a single channel; pass one index to --seg-channels")
            run_stardist(dataset, seg_channels[0], use_gpu=args.use_gpu)
        else:
            run_instanseg(dataset, seg_channels, pixel_size_um,
                          model_type=args.instanseg_model, use_gpu=args.use_gpu)

    if needs_aggregation:
        # Stardist gives nuclei only; approximate whole cells by watershed expansion.
        # InstanSeg already produced cell_boundaries.
        if args.segmenter == 'stardist' and args.cell_expansion > 0:
            if STARDIST_NUCLEI_KEY not in dataset.shapes:
                print(f"WARNING: no {STARDIST_NUCLEI_KEY} to expand; skipping cell expansion")
            elif pixel_size_um is None:
                print("WARNING: no physical pixel size available; cannot expand cells in µm; skipping")
            else:
                nuclei = dataset.shapes[STARDIST_NUCLEI_KEY]
                dataset[CELLS_KEY] = expand_nuclei(nuclei, args.cell_expansion, pixel_size_um)
                print(f"Expanded {len(nuclei)} nuclei into cell boundaries (--cell-expansion={args.cell_expansion} µm)")
        elif args.segmenter == 'instanseg' and args.cell_expansion != 5.0:
            print("NOTE: --cell-expansion is ignored for instanseg (cells are segmented directly)")

        # Aggregate over whole cells when available, else nuclei.
        agg_shapes_key = resolve_shape_keys(dataset)[1]
        print(f"Aggregating over {agg_shapes_key}...")

        # Force pandas to use plain object dtype for strings, not ArrowStringArray,
        # which can't be written to Zarr backing stores by anndata
        with pd.option_context('future.infer_string', False):
            if "spots" in dataset.points:
                sopa.aggregate(dataset, aggregate_genes=True, points_key='spots', gene_column='gene',
                               shapes_key=agg_shapes_key)
            else:
                sopa.aggregate(dataset, shapes_key=agg_shapes_key)

        np.random.seed(42)

        sdata = dataset

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

    cluster_labels, cluster_distances = cluster_data(intensity_df, n_clusters=args.n_clusters)

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

    # Build the network on the cell layer (not the nuclei), so cell ids match the
    # table produced by aggregation over that layer.
    nuclei_key, network_shapes_key = resolve_shape_keys(sdata)
    result = run_network_analysis(sdata, intensity_matrix=intensity_df,
                                cell_boundaries=network_shapes_key,
                                comm_detect_res=args.community_resolution,
                                max_edge_distance=args.max_edge_distance,
                                output_dir=args.plot_dir)

    spots_with_cells = assign_spots_to_cells(sdata, cell_boundaries=network_shapes_key)

    export_to_qupath(result['cell_ids'], result['community_labels'], result['cluster_labels'],
                     output_path=args.geojson_path,
                     sdata=sdata, intensity_df=intensity_df,
                     spots_with_cells=spots_with_cells,
                     cell_boundaries=network_shapes_key)

    # Persist everything to Zarr.
    # Write to a temporary path first, then replace the target path atomically
    # to avoid "path in use" errors when the segmented Zarr was loaded for resumption.
    import tempfile
    import shutil as _shutil
    # Attach cluster/community to the network shapes layer and the table (they
    # share the same index/cell ids).
    sdata.shapes[network_shapes_key]['kmeans_cluster'] = result['cluster_labels']
    sdata.shapes[network_shapes_key]['community'] = result['community_labels']
    sdata.tables['table'].obs['community'] = pd.Categorical(result['community_labels'])

    # Also attach cluster/community to the nuclei layer via the nucleus_id
    # recorded on each cell (survives sopa.aggregate's re-indexing). Nuclei not
    # linked to a cell get -1.
    if nuclei_key is not None and nuclei_key != network_shapes_key \
            and 'nucleus_id' in sdata.shapes[network_shapes_key].columns:
        cell_gdf = sdata.shapes[network_shapes_key]
        linked = cell_gdf['nucleus_id'].notna().to_numpy()
        nucleus_ids = cell_gdf['nucleus_id'].to_numpy()[linked]
        nucleus_to_cluster = dict(zip(nucleus_ids, np.asarray(result['cluster_labels'])[linked]))
        nucleus_to_community = dict(zip(nucleus_ids, np.asarray(result['community_labels'])[linked]))
        nuclei = sdata.shapes[nuclei_key]
        nuclei['kmeans_cluster'] = [nucleus_to_cluster.get(cid, -1) for cid in nuclei.index]
        nuclei['community'] = [nucleus_to_community.get(cid, -1) for cid in nuclei.index]
    print(f'Saving Zarr to {seg_zarr_path}')
    try:
        tmpdir = tempfile.mkdtemp(dir=os.path.dirname(seg_zarr_path) or '.')
        tmp_zarr = os.path.join(tmpdir, os.path.basename(seg_zarr_path))
        sdata.write(tmp_zarr)
        if os.path.exists(seg_zarr_path):
            _shutil.rmtree(seg_zarr_path)
        os.rename(tmp_zarr, seg_zarr_path)
        os.rmdir(tmpdir)
    except Exception as e:
        print(f"ERROR: Failed to write segmented Zarr to {seg_zarr_path}: {e}", file=sys.stderr)
        sys.exit(1)
    print("Done")

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

    # Clean up Dask workers that may keep the process alive
    import dask.distributed as dd
    try:
        client = dd.get_client()
        client.close()
    except ValueError:
        pass

    # Per-cluster cell density maps with distance-based soft weighting.
    # Closer to cluster center → higher weight in that cluster's density.
    obs = sdata.tables['table'].obs
    obs_ids = [str(cid) for cid in obs.index]

    # Cell centroids aligned to the table (a filtered subset of the cell layer).
    cell_gdf = sdata.shapes[network_shapes_key]
    centroid_map = {str(cid): (geom.centroid.x, geom.centroid.y)
                    for cid, geom in zip(cell_gdf.index, cell_gdf.geometry)}
    x = np.array([centroid_map.get(cid, (np.nan, np.nan))[0] for cid in obs_ids])
    y = np.array([centroid_map.get(cid, (np.nan, np.nan))[1] for cid in obs_ids])
    bins = min(200, int(np.sqrt(len(x))))

    cluster_labels_arr = obs['kmeans_cluster'].astype(int).to_numpy()
    community_labels_arr = obs['community'].astype(int).to_numpy()

    n_clusters = args.n_clusters
    n_communities = len(set(community_labels_arr)) - (1 if -1 in community_labels_arr else 0)

    exp_weights = np.exp(-cluster_distances)
    exp_weights /= exp_weights.sum(axis=1, keepdims=True)

    n_rows = n_clusters + n_communities
    n_cols = 1
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(10, 6 * n_rows), squeeze=False)

    for c in range(n_clusters):
        weights = exp_weights[:, c]
        density, xedges, yedges = np.histogram2d(x, y, bins=bins, weights=weights)
        ax = axes[c, 0]
        im = ax.imshow(np.log1p(density.T * 1000), origin='lower', aspect='auto', cmap='inferno',
                       extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
        plt.colorbar(im, ax=ax, label='log(weighted density + 1)')
        count = np.sum(cluster_labels_arr == c)
        ax.set_title(f'Cluster {c} density ({count} cells)')

    for comm in range(n_communities):
        mask = community_labels_arr == comm
        if not mask.any():
            continue
        density, xedges, yedges = np.histogram2d(x[mask], y[mask], bins=bins)
        ax = axes[n_clusters + comm, 0]
        im = ax.imshow(np.log1p(density.T), origin='lower', aspect='auto', cmap='inferno',
                       extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
        plt.colorbar(im, ax=ax, label='log(cell count + 1)')
        count = mask.sum()
        ax.set_title(f'Community {comm} density ({count} cells)')

    plt.tight_layout()
    plt.savefig(os.path.join(args.plot_dir, 'cell_density.png'))
    plt.close()

    # Channel intensity heatmap per k-means cluster
    cluster_means = intensity_df.groupby(sdata.tables['table'].obs['kmeans_cluster']).mean()
    plt.figure(figsize=(max(8, cluster_means.shape[1] * 0.4), max(6, cluster_means.shape[0] * 0.5)))
    sns.heatmap(cluster_means, cmap='RdBu_r', center=0, xticklabels=True, yticklabels=True,
                cbar_kws={'label': 'Mean intensity (z-score)'})
    plt.title('Mean channel intensity per k-means cluster')
    plt.xlabel('Channel')
    plt.ylabel('Cluster')
    plt.tight_layout()
    plt.savefig(os.path.join(args.plot_dir, 'cluster_intensity_heatmap.png'))
    plt.close()

    # Morphological features per cluster — cell area
    obs_df = sdata.tables['table'].obs
    morph_cols = [c for c in obs_df.columns if c in ('area',)]
    if morph_cols:
        import matplotlib.gridspec as gridspec
        n_morph = len(morph_cols)
        fig = plt.figure(figsize=(max(8, n_morph * 4), 6))
        for i, col in enumerate(morph_cols):
            ax = fig.add_subplot(1, n_morph, i + 1)
            sns.swarmplot(data=obs_df, x='kmeans_cluster', y=col, ax=ax, size=1)
            ax.set_title(f'{col} by cluster')
            ax.set_xlabel('Cluster')
        plt.tight_layout()
        plt.savefig(os.path.join(args.plot_dir, 'morphology_by_cluster.png'))
        plt.close()


if __name__ == '__main__':
    main()
