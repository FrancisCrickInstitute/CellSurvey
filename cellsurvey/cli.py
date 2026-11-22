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

from cellsurvey.utils import (cluster_data, assign_spots_to_cells, extract_feature_tables,
                              build_clustering_features, resolve_shape_keys, store_pixel_size,
                              store_table_meta, resolve_pixel_size, CLUSTER_ON_CHOICES,
                              DEFAULT_MIN_COMMUNITY_SIZE,
                              STARDIST_NUCLEI_KEY, NUCLEI_KEY, CELLS_KEY)
from cellsurvey.blob_detection import detect_spots
from cellsurvey.image_source import open_image_source, stored_pixel_size, is_omero
from cellsurvey.network_analysis import run_network_analysis
from cellsurvey.segmentation import expand_nuclei, run_stardist, run_instanseg
from cellsurvey.export import export_to_qupath
from cellsurvey.report import write_report, code_version
from cellsurvey.niches import (neighbourhood_composition, assign_niches, table_centroids_um,
                               DEFAULT_N_NICHES, DEFAULT_NICHE_RADIUS_UM)

# cell_density.png shows only this many of the largest communities.
MAX_DENSITY_COMMUNITIES = 20


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', required=True,
                        help='Input image: a multichannel (OME-)TIFF path, or an OMERO image as '
                             'omero://<host>[:<port>]/Image:<id> (log in first with `omero login`)')
    parser.add_argument('--omero-readers', type=int, default=8,
                        help='OMERO input only: parallel tile readers (server connections) used when '
                             'copying the image into <output>.zarr (default: 8)')
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
    parser.add_argument('--keep-unnucleated', action='store_true',
                        help='InstanSeg only: keep cells that contain no nucleus (dropped by default)')
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
    parser.add_argument('--cluster-on', choices=CLUSTER_ON_CHOICES, default='intensities',
                        help='Per-cell features for k-means, the network edge weights and UMAP/Leiden: '
                             'channel intensities, spot densities (needs --detect-blobs), or both, '
                             'weighted equally (default: intensities)')
    parser.add_argument('--community-resolution', type=float, default=0.1,
                        help='Resolution for Louvain community detection (default: 0.1)')
    parser.add_argument('--min-community-size', type=int, default=DEFAULT_MIN_COMMUNITY_SIZE,
                        help='Communities with fewer cells than this are labelled -1 (unassigned); the rest '
                             f'are numbered by size, 0 = largest (default: {DEFAULT_MIN_COMMUNITY_SIZE}; '
                             '1 keeps every community)')
    parser.add_argument('--n-niches', type=int, default=DEFAULT_N_NICHES,
                        help='Number of neighbourhood niches: groups of cells whose surroundings (within '
                             '--niche-radius) have a similar mix of k-means clusters. 0 skips niches '
                             f'(default: {DEFAULT_N_NICHES}; a starting point - check with '
                             'python -m cellsurvey.stability --mode niches)')
    parser.add_argument('--niche-radius', type=float, default=DEFAULT_NICHE_RADIUS_UM,
                        help='Neighbourhood radius in µm for the niches (default: '
                             f'{DEFAULT_NICHE_RADIUS_UM:g}; suitable values depend on the tissue - check with the '
                             'niche sweep)')
    parser.add_argument('--max-edge-distance', type=float, default=50,
                        help='Max edge length in µm for the Delaunay network (default: 50)')
    parser.add_argument('--radius-min', type=float, default=0,
                        help='Min radius in µm for the spatial neighbors graph (default: 0)')
    parser.add_argument('--radius-max', type=float, default=280,
                        help='Max radius in µm for the spatial neighbors graph (default: 280)')
    parser.add_argument('--resume-from',
                        help='Path to existing Zarr to resume from (never re-reads the input image; '
                             'reuses spots and a valid _seg.zarr)')
    parser.add_argument('--geojson-path', default=None,
                        help='Output path for QuPath GeoJSON (default: <output_file>.geojson, '
                             'next to the Zarr outputs)')
    parser.add_argument('--font-size', type=int, default=20,
                        help='Font size for plots (default: 20)')
    parser.add_argument('--axes-linewidth', type=int, default=3,
                        help='Axes line width for plots (default: 3)')
    parser.add_argument('--no-report', action='store_true',
                        help='Skip the HTML overview report (<plot_dir>/report.html), written by default at '
                             'the end of the run. It can be made later with: '
                             'python -m cellsurvey.report --zarr <output>_seg.zarr -p <plot_dir>')
    args = parser.parse_args()

    plt.rcParams['font.size'] = args.font_size
    plt.rcParams['axes.linewidth'] = args.axes_linewidth

    os.makedirs(args.plot_dir, exist_ok=True)

    imagepath = args.input_file

    if args.omero_readers < 1:
        parser.error("--omero-readers must be at least 1")
    if args.min_community_size < 1:
        parser.error("--min-community-size must be at least 1")
    if args.n_niches < 0:
        parser.error("--n-niches must be 0 (no niches) or more")
    if args.n_niches > 0 and not args.niche_radius > 0:
        parser.error("--niche-radius must be positive")

    # Spot-based clustering needs spots; a new run only gets them with --detect-blobs.
    # (A resumed run is checked once its table is loaded.)
    if args.cluster_on != 'intensities' and not args.detect_blobs and not args.resume_from:
        parser.error(f"--cluster-on {args.cluster_on} needs spots: add --detect-blobs")

    # Check blob-detection arguments before reading anything.
    if args.detect_blobs:
        try:
            channel_index = [int(c) for c in args.channels.split(',')]
            thresholds = [float(t) for t in args.thresholds.split(',')]
        except ValueError:
            parser.error("--channels must be comma-separated integers and --thresholds comma-separated numbers")
        if len(channel_index) != len(thresholds):
            parser.error(
                f"Number of channels ({len(channel_index)}) must match number of thresholds ({len(thresholds)})")

    resume_zarr_path = None
    if args.resume_from:
        resume_zarr_path = args.resume_from
        if not resume_zarr_path.endswith('.zarr'):
            resume_zarr_path += '.zarr'
        if not os.path.exists(resume_zarr_path):
            parser.error(f"Resume Zarr not found: {resume_zarr_path}")

    def open_source():
        if not is_omero(imagepath) and not os.path.exists(imagepath):
            hint = (" (an OMERO URL must start with omero://; check the host part isn't an empty "
                    "shell variable)" if imagepath.lower().lstrip('/').startswith('image:') else "")
            parser.error(f"Input image not found: {imagepath}{hint}")
        try:
            return open_image_source(imagepath, readers=args.omero_readers)
        except (ValueError, RuntimeError, FileNotFoundError, ImportError) as e:
            parser.error(str(e))
        except Exception as e:
            if is_omero(imagepath):
                raise   # unexpected OMERO/Ice error: keep the traceback
            parser.error(f"Can't read input image {imagepath}: {type(e).__name__}: {e}. "
                         f"Expected a multichannel (OME-)TIFF, or an omero://<host>/Image:<id> URL.")

    # New runs open the input image once, here (lazily: nothing is read yet);
    # stage 1 writes it to <prefix>.zarr and every later stage reads only that.
    # Resumed runs never open the source unless an older Zarr lacks a stored
    # pixel size.
    source = None if args.resume_from else open_source()

    def fail(message):
        """parser.error, closing the source first (e.g. the OMERO connection)."""
        if source is not None:
            source.close()
        parser.error(message)

    if args.pixel_size is not None:
        if not args.pixel_size > 0:  # also rejects NaN
            fail("--pixel-size must be positive")
        pixel_size_um = args.pixel_size
        print(f"Physical pixel size (from --pixel-size): {pixel_size_um:.4f} µm/px")
    else:
        pixel_size_origin = "image metadata"
        if source is not None:
            pixel_size_um = source.pixel_size_um
        else:
            pixel_size_um = stored_pixel_size(spatialdata.read_zarr(resume_zarr_path))
            if pixel_size_um is not None:
                pixel_size_origin = f"stored in {resume_zarr_path}"
            else:  # Zarr written before pixel sizes were stored: read the source's metadata only
                old_source = open_source()
                pixel_size_um = old_source.pixel_size_um
                old_source.close()
        if pixel_size_um is None:
            # Fail before any loading/segmentation: every distance parameter is in µm.
            fail("No physical pixel size in the image metadata; pass --pixel-size (µm/px). "
                 "It is needed because --max-edge-distance, --radius-min/--radius-max and "
                 "--cell-expansion are in µm.")
        print(f"Physical pixel size ({pixel_size_origin}): {pixel_size_um:.4f} µm/px")

    # Shape layers the chosen segmenter writes; used to decide whether a resumed
    # _seg.zarr already holds this segmenter's output. (cell_boundaries alone is
    # not enough for InstanSeg: Stardist + expansion also writes it.)
    seg_layers = [STARDIST_NUCLEI_KEY] if args.segmenter == 'stardist' else [NUCLEI_KEY, CELLS_KEY]

    zarr_path = args.output_file
    if not zarr_path.endswith('.zarr'):
        zarr_path += '.zarr'
    # Default the GeoJSON to sit beside the Zarrs (from -o, even when resuming),
    # so runs with different -o don't overwrite each other's export.
    if args.geojson_path is None:
        args.geojson_path = zarr_path[:-len('.zarr')] + '.geojson'

    if args.resume_from:
        zarr_path = resume_zarr_path
        print(f"Resuming from existing Zarr: {zarr_path}")
    else:
        # Stage 1: write the image (only) to <prefix>.zarr. Blob detection runs
        # afterwards on this local copy, so the source is read once.
        if args.detect_blobs:
            bad = [ci for ci in channel_index if not 0 <= ci < source.data.shape[0]]
            if bad:
                fail(f"--channels {bad} out of range: the image has {source.data.shape[0]} channels")
        dataset = source.to_spatialdata(pixel_size_um)
        print(f'Saving Zarr to {zarr_path}')

        try:
            source.write(dataset, zarr_path)
        except Exception as e:
            print(f"ERROR: Failed to write Zarr to {zarr_path}: {e}", file=sys.stderr)
            source.close()
            sys.exit(1)
        # The source is never read again (e.g. this closes the OMERO connection).
        source.close()
        del dataset

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
                # Segmentation/expansion used the stored pixel size; keep using it
                # so the µm conversions stay consistent with the existing shapes.
                try:
                    stored_px = resolve_pixel_size(sdata)
                except ValueError:
                    stored_px = None  # older Zarr, or invalid stored value
                if stored_px is not None and not np.isclose(stored_px, pixel_size_um, rtol=1e-3):
                    print(f"WARNING: pixel size {pixel_size_um:.4f} µm/px differs from the "
                          f"{stored_px:.4f} µm/px stored in {seg_zarr_path}, which segmentation "
                          f"used; continuing with the stored value. Re-segment (without "
                          f"--resume-from) to use the new one.")
                    pixel_size_um = stored_px
            elif all(k in dataset.shapes for k in seg_layers):
                print(f"Segmented Zarr has {', '.join(seg_layers)} but no table, re-running aggregation only")
                needs_segmentation = False
        except Exception:
            import shutil
            shutil.rmtree(seg_zarr_path, ignore_errors=True)

    # Blob detection, on the local <prefix>.zarr (never the source). Only when
    # segmentation will re-read that Zarr; a reused _seg.zarr already has its spots.
    if args.detect_blobs:
        if needs_segmentation:
            image_sdata = spatialdata.read_zarr(zarr_path)
            if 'spots' in image_sdata.points:
                print(f"{zarr_path} already has spots; skipping blob detection")
            else:
                try:
                    points = detect_spots(sopa.utils.get_spatial_image(image_sdata), channel_index, thresholds,
                                          tile_size=args.tile_size, overlap=args.overlap,
                                          min_sigma=args.min_sigma, max_sigma=args.max_sigma,
                                          num_sigma=args.num_sigma, n_workers=args.workers)
                except ValueError as e:
                    parser.error(str(e))
                if points is None:
                    # Nothing to count: carry on as a run without spots.
                    if args.cluster_on != 'intensities':
                        parser.error(f"No spots were detected in channels {channel_index} "
                                     f"(--thresholds {args.thresholds}), so --cluster-on "
                                     f"{args.cluster_on} has nothing to use. Lower --thresholds, or "
                                     f"use --cluster-on intensities.")
                    print(f"WARNING: no spots detected in channels {channel_index} (--thresholds "
                          f"{args.thresholds}); continuing without spots. Lower --thresholds to find some.")
                else:
                    image_sdata['spots'] = points
                    print(f'Adding spots to {zarr_path}')
                    try:
                        image_sdata.write_element('spots')
                    except Exception as e:
                        print(f"ERROR: Failed to write spots to {zarr_path}: {e}", file=sys.stderr)
                        sys.exit(1)
                    print("Done")
        elif 'spots' not in dataset.points:
            print(f"WARNING: --detect-blobs ignored: the reused {seg_zarr_path} has no spots. "
                  f"Re-run without --resume-from to detect them.")

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
                          model_type=args.instanseg_model, use_gpu=args.use_gpu,
                          keep_unnucleated=args.keep_unnucleated)

    if needs_aggregation:
        # Stardist gives nuclei only; approximate whole cells by watershed expansion.
        # InstanSeg already produced cell_boundaries.
        if args.segmenter == 'stardist' and args.cell_expansion > 0:
            if STARDIST_NUCLEI_KEY not in dataset.shapes:
                print(f"WARNING: no {STARDIST_NUCLEI_KEY} to expand; skipping cell expansion")
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

    # Force obs index to plain string dtype to avoid ArrowStringArray
    # write errors with Zarr backing stores (done first, so every per-cell
    # table below shares this index).
    sdata.tables['table'].obs.index = sdata.tables['table'].obs.index.astype(str)
    measurements = sdata.tables['table']

    # Per-cell channel intensities and, with spots, spot counts. (With spots,
    # sopa.aggregate puts the counts in .X and the intensities in obsm.)
    intensity_df, counts_df = extract_feature_tables(measurements)
    print(f"Intensity channels: {list(intensity_df.columns)}")
    if counts_df is not None:
        print(f"Spot channels: {list(counts_df.columns)}")
    print(f"Number of cells: {len(intensity_df)}")

    # The features k-means, the network edge weights and UMAP/Leiden use.
    if args.cluster_on != 'intensities' and counts_df is None:
        parser.error(f"--cluster-on {args.cluster_on} needs spots, but {zarr_path} has none: "
                     f"add --detect-blobs (and re-run without --resume-from if a _seg.zarr is reused)")
    cell_area_um2 = None
    if counts_df is not None and 'area' in measurements.obs:
        cell_area_um2 = measurements.obs['area'].astype(float) * pixel_size_um ** 2   # px² → µm²
    features = build_clustering_features(intensity_df, counts_df, args.cluster_on, cell_area_um2)
    print(f"Clustering on {args.cluster_on}: {features.shape[1]} features")

    cluster_labels, cluster_distances = cluster_data(features, n_clusters=args.n_clusters, scale=False)

    sdata.tables['table'].obs['kmeans_cluster'] = pd.Categorical(cluster_labels)

    sdata.tables['table'].obs['kmeans_cluster_label'] = pd.Categorical(
        [f'Cluster_{i}' for i in cluster_labels]
    )

    print("\n\u2713 Cluster labels added to sdata.tables['table'].obs['kmeans_cluster']")
    print("\u2713 Cluster labels added to sdata.tables['table'].obs['kmeans_cluster_label']")

    print(f"\nUpdated obs columns: {list(sdata.tables['table'].obs.columns)}")

    # Neighbourhood niches: k-means on each cell's local mix of k-means clusters
    # (independent of the Louvain network below). Labels follow the table's obs order.
    niche_summary = {}
    obs_niche = None
    if args.n_niches > 0:
        xy_um = table_centroids_um(sdata, pixel_size_um)
        composition = neighbourhood_composition(xy_um, cluster_labels, args.niche_radius,
                                                n_labels=int(args.n_clusters))
        niche_labels = assign_niches(composition, args.n_niches)
        sdata.tables['table'].obs['niche'] = pd.Categorical(niche_labels)
        obs_niche = dict(zip(sdata.tables['table'].obs.index.astype(str), niche_labels))
        sizes = np.bincount(niche_labels[niche_labels >= 0], minlength=args.n_niches)
        niche_summary = {'n_niches': int(args.n_niches), 'niche_radius_um': float(args.niche_radius),
                         'niche_n_clusters': int(args.n_clusters),
                         'niche_sizes': {int(i): int(s) for i, s in enumerate(sizes)}}
        print(f"Niches: {args.n_niches} niches from each cell's mix of the {args.n_clusters} k-means clusters "
              f"within {args.niche_radius:g} µm (numbered by size): "
              + ', '.join(f'{i}: {s:,}' for i, s in enumerate(sizes)))
    else:
        print("Niches skipped (--n-niches 0)")
        # Don't leave a resumed run's old niche labels looking current.
        if 'niche' in sdata.tables['table'].obs:
            del sdata.tables['table'].obs['niche']
        for key in resolve_shape_keys(sdata):
            if key is not None and 'niche' in sdata.shapes[key].columns:
                del sdata.shapes[key]['niche']

    # Build the network on the cell layer (not the nuclei), so cell ids match the
    # table produced by aggregation over that layer.
    nuclei_key, network_shapes_key = resolve_shape_keys(sdata)
    result = run_network_analysis(sdata, intensity_matrix=features,
                                cell_boundaries=network_shapes_key,
                                comm_detect_res=args.community_resolution,
                                max_edge_distance=args.max_edge_distance,
                                pixel_size_um=pixel_size_um,
                                min_community_size=args.min_community_size,
                                output_dir=args.plot_dir,
                                summary_extra={'cluster_on': args.cluster_on,
                                               'n_features': int(features.shape[1]),
                                               **niche_summary})
    # Niche per cell shape (aligned with result['cell_ids']); -1 for shapes not in the table.
    shape_niches = (np.array([obs_niche.get(str(cid), -1) for cid in result['cell_ids']], dtype=int)
                    if obs_niche is not None else None)

    spots_with_cells = assign_spots_to_cells(sdata, cell_boundaries=network_shapes_key)

    export_to_qupath(result['cell_ids'], result['community_labels'], result['cluster_labels'],
                     output_path=args.geojson_path,
                     sdata=sdata, intensity_df=intensity_df, counts_df=counts_df,
                     spots_with_cells=spots_with_cells,
                     cell_boundaries=network_shapes_key,
                     niche_labels=obs_niche)

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
    if shape_niches is not None:
        sdata.shapes[network_shapes_key]['niche'] = shape_niches
    # Persist the pixel size so read-only consumers (stability sweep) can convert µm params.
    store_pixel_size(sdata, pixel_size_um)
    # Run parameters, so the report (and anyone reading the Zarr) needn't rely on summary.json.
    store_table_meta(sdata, cluster_on=args.cluster_on, min_community_size=args.min_community_size,
                     n_clusters=int(args.n_clusters), community_resolution=float(args.community_resolution),
                     max_edge_distance_um=float(args.max_edge_distance), segmenter=args.segmenter,
                     **({'instanseg_model': args.instanseg_model} if args.segmenter == 'instanseg'
                        else {'cell_expansion_um': float(args.cell_expansion)}),
                     code_version=code_version() or 'unknown',
                     n_niches=int(args.n_niches),
                     **({'niche_radius_um': float(args.niche_radius), 'niche_n_clusters': int(args.n_clusters)}
                        if args.n_niches > 0 else {}))

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
        if shape_niches is not None:
            nucleus_to_niche = dict(zip(nucleus_ids, shape_niches[linked]))
            nuclei['niche'] = [nucleus_to_niche.get(cid, -1) for cid in nuclei.index]
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
    # Radius flags are in µm; cell coordinates (obsm['spatial']) are in pixels.
    sopa.spatial.spatial_neighbors(adata, radius=(args.radius_min / pixel_size_um,
                                                  args.radius_max / pixel_size_um))
    cell_type_to_cell_type = sopa.spatial.mean_distance(adata, "kmeans_cluster", "kmeans_cluster")
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.linewidth'] = 2
    plt.figure(figsize=(10, 10))
    heatmap_kwargs = {"cmap": sns.cm.rocket_r, "cbar_kws": {"label": "Mean hop distance"}}
    sns.heatmap(cell_type_to_cell_type, **heatmap_kwargs)
    plt.tight_layout()
    plt.savefig(os.path.join(args.plot_dir, 'cell_type_to_cell_type.png'))

    # UMAP/Leiden on the same scaled features as k-means and the network (not
    # Scanpy's count normalisation of .X, which holds intensities or spot counts).
    adata.obsm['X_cellsurvey'] = features.reindex(adata.obs.index).fillna(0).to_numpy()
    sc.pp.neighbors(adata, use_rep='X_cellsurvey')
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

    exp_weights = np.exp(-cluster_distances)
    exp_weights /= exp_weights.sum(axis=1, keepdims=True)

    plot_cell_density(x, y, cluster_labels_arr, community_labels_arr, exp_weights, args.n_clusters,
                      bins, os.path.join(args.plot_dir, 'cell_density.png'))

    # Heatmap of the clustering features per k-means cluster: each feature
    # z-scored across all cells (as k-means sees it, but without the --cluster-on
    # both group weighting, so values are comparable), then averaged per cluster;
    # 0 = the all-cell mean. Spot features are log spot densities, '(spots)'.
    heatmap_features = build_clustering_features(intensity_df, counts_df, args.cluster_on, cell_area_um2,
                                                 weight_groups=False)
    cluster_means = heatmap_features.groupby(sdata.tables['table'].obs['kmeans_cluster']).mean()
    plt.figure(figsize=(max(8, cluster_means.shape[1] * 0.4), max(6, cluster_means.shape[0] * 0.5)))
    # Cap the colour scale at ±3 SD so one outlier cluster (e.g. a single debris
    # cell at z≈20) can't wash out the differences between the others.
    sns.heatmap(cluster_means, cmap='RdBu_r', center=0, vmin=-3, vmax=3,
                xticklabels=True, yticklabels=True,
                cbar_kws={'label': 'Mean z-score across cells (capped at ±3)', 'extend': 'both'})
    plt.title(f'Clustering features per k-means cluster (--cluster-on {args.cluster_on})')
    plt.xlabel('Feature')
    plt.ylabel('Cluster')
    plt.tight_layout()
    plt.savefig(os.path.join(args.plot_dir, 'cluster_intensity_heatmap.png'))
    plt.close()

    # Morphology per cluster: cell area. Sopa's obs['area'] is in px² (shape
    # coordinates are pixels), so convert to µm² for the plot.
    obs_df = sdata.tables['table'].obs
    if 'area' in obs_df.columns:
        morph_df = pd.DataFrame({
            'kmeans_cluster': obs_df['kmeans_cluster'],
            'area_um2': obs_df['area'].astype(float) * pixel_size_um ** 2,
        })
        fig, ax = plt.subplots(figsize=(8, 6))
        # Violin, not swarm: swarmplot places every point individually and
        # takes many hours at whole-slide scale (100k+ cells per cluster).
        sns.violinplot(data=morph_df, x='kmeans_cluster', y='area_um2', ax=ax, cut=0, inner='quartile')
        ax.set_title('Cell area by cluster')
        ax.set_xlabel('Cluster')
        ax.set_ylabel('Cell area (µm²)')
        plt.tight_layout()
        plt.savefig(os.path.join(args.plot_dir, 'morphology_by_cluster.png'))
        plt.close(fig)

    # Overview report, from the _seg.zarr just written. Last, and never fatal: a
    # failure here mustn't lose a finished run (it can be re-made standalone).
    if not args.no_report:
        try:
            write_report(seg_zarr_path, args.plot_dir)
        except Exception as e:
            print(f"WARNING: overview report failed ({type(e).__name__}: {e}). Re-make it with: "
                  f"python -m cellsurvey.report --zarr {seg_zarr_path} -p {args.plot_dir}", file=sys.stderr)


def plot_cell_density(x, y, cluster_labels, community_labels, exp_weights, n_clusters, bins, path,
                      max_communities=MAX_DENSITY_COMMUNITIES):
    """Write ``cell_density.png``: one soft-weighted density panel per k-means
    cluster, plus plain cell-count panels for the ``max_communities`` largest
    communities (-1, unassigned, is never plotted). The title says how many
    communities and unassigned cells were left out."""
    # Pick communities by actual size (not label order).
    comm_ids, comm_counts = np.unique(community_labels[community_labels >= 0], return_counts=True)
    plotted_comms = comm_ids[np.lexsort((comm_ids, -comm_counts))][:max_communities]
    n_comm_omitted = len(comm_ids) - len(plotted_comms)

    # Grid layout: clusters first, then communities starting on a new row.
    n_cols = 4
    cluster_rows = -(-n_clusters // n_cols)
    comm_rows = -(-len(plotted_comms) // n_cols)
    n_rows = max(cluster_rows + comm_rows, 1)
    # Constrained layout (not tight_layout) leaves room for the suptitle.
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows), squeeze=False,
                             layout='constrained')
    for ax in axes.ravel():
        ax.axis('off')

    def _density_panel(ax, density, xedges, yedges, label, title):
        ax.axis('on')
        # Image coordinates: y increases downwards (as in QuPath/napari), so
        # the smallest y goes at the top.
        im = ax.imshow(np.log1p(density.T), origin='upper', aspect='auto', cmap='inferno',
                       extent=[xedges[0], xedges[-1], yedges[-1], yedges[0]])
        plt.colorbar(im, ax=ax, label=label)
        ax.set_title(title)

    # Every panel spans the whole tissue, so communities are shown in place
    # rather than each zoomed to its own extent.
    finite = np.isfinite(x) & np.isfinite(y)
    extent = [[x[finite].min(), x[finite].max()], [y[finite].min(), y[finite].max()]] if finite.any() else None

    for c in range(n_clusters):
        density, xedges, yedges = np.histogram2d(x, y, bins=bins, range=extent, weights=exp_weights[:, c])
        _density_panel(axes[c // n_cols, c % n_cols], density * 1000, xedges, yedges,
                       'log(weighted density + 1)',
                       f'Cluster {c} density ({np.sum(cluster_labels == c)} cells)')

    for i, comm in enumerate(plotted_comms):
        mask = community_labels == comm
        density, xedges, yedges = np.histogram2d(x[mask], y[mask], bins=bins, range=extent)
        _density_panel(axes[cluster_rows + i // n_cols, i % n_cols], density, xedges, yedges,
                       'log(cell count + 1)', f'Community {comm} density ({mask.sum()} cells)')

    title = f'Cell density: {n_clusters} k-means clusters'
    if len(plotted_comms):
        title += f', the {len(plotted_comms)} largest communities'
    if n_comm_omitted:
        title += f' ({n_comm_omitted} smaller communities not shown)'
    n_unassigned = int((community_labels < 0).sum())
    if n_unassigned:
        title += f'; {n_unassigned} unassigned cells not shown'
    fig.suptitle(title)
    fig.savefig(path)
    plt.close(fig)


if __name__ == '__main__':
    main()
