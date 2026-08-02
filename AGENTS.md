# AGENTS.md

## Project overview

CellSurvey is a spatial biology/omics analysis pipeline built on [Sopa](https://gustaveroussy.github.io/sopa/) (segmentation, aggregation). It processes multichannel microscopy images (OME-TIFF) into spatial data objects, segments nuclei with Stardist, detects RNA spots via blob detection, clusters cells with k-means, builds Delaunay networks, detects Louvain communities, and exports GeoJSON for QuPath visualization.

**MuSpAn removed (complete):** MuSpAn has been fully replaced with open-source libraries. Delaunay triangulation uses `scipy.spatial.Delaunay`, Louvain community detection uses `networkx`, and visualisation uses matplotlib. No private dependencies remain.

## Environment and package management

This project uses **pixi** (via `pixi.toml`) for environment management targeting `linux-64` only. The lockfile is `pixi.lock`.

Key dependency constraints:
- **Python**: `>=3.12, <3.13` on `linux-64`
- **TensorFlow**: `>=2.18` on `linux-64` (with `and-cuda` extras for GPU)
- **CUDA/cuDNN**: TF >=2.18 bundles its own CUDA 12/cuDNN 9 libraries; no separate conda CUDA/cuDNN packages are needed
- **`tf_keras`** is a required pypi dependency — TF >=2.16 defaults to Keras 3, but Stardist needs legacy Keras 2 API to avoid cuDNN autotuner failures on CUDA 12/cuDNN 9
- **`scipy` (`>=1.14, <2`)** and **`networkx` (`>=3.4, <4`)** for Delaunay triangulation and Louvain community detection (replacing MuSpAn)
- **`python-igraph`** for fast Leiden clustering in Scanpy's spatial neighborhood analysis
- **Windows and macOS are not supported** via pixi — only `linux-64` is in the platforms list.

A **Dockerfile** is provided: Ubuntu 24.04 base, installs pixi, copies `pixi.toml`, sets `TF_USE_LEGACY_KERAS=1`, entrypoint is `pixi run python run.py`.

There is no Makefile, no CI/CD, and no tests. The source code is split across 6 files under the `cellsurvey/` package, with `run.py` as the entry-point shim.

**TODO**: Set up linting and formatting (Ruff, mypy) with a `pyproject.toml` config and pre-commit hooks.

## Commands

**Development environment setup (pixi):**
```bash
pixi install
```

**Run the pipeline (pixi):**
```bash
pixi run python run.py -i <input_tiff> -o <output_zarr_prefix> -p <plot_output_dir>
```

**Run the pipeline (standalone):**
```bash
python run.py -i <input_tiff> -o <output_zarr_prefix> -p <plot_output_dir>
```

Three required arguments:
- `-i`: Path to input OME-TIFF image
- `-o`: Path prefix for output Zarr file (`.zarr` suffix is appended automatically)
- `-p`: Directory for output plots

**Docker build and run:**
```bash
docker build -t cellsurvey .
docker run --gpus all -v /path/to/data:/data cellsurvey -i /data/input.tiff -o /data/output -p /data/plots
```

There is no build step, no test command, and no linting configured.

## Architecture and data flow

The pipeline is split into modules under the `cellsurvey/` package. `run.py` is a shim that sets environment variables (TF_USE_LEGACY_KERAS), preloads libstdc++, and delegates to `cellsurvey.cli.main()`. The pipeline runs these stages sequentially:

1. **Image loading** (`cli.py`): Reads the OME-TIFF via `BioImage` (from `bioio`) to get channel names, and via `sopa.io.ome_tif()` as a SpatialData dataset. Three code paths at this stage:
   - `--resume-from`: Skips image loading entirely — reads channel names from the existing image via `BioImage` but loads the Zarr directly.
   - `--detect-blobs`: Loads only the subset of channels needed for blob detection via `dask_image.imread` to limit memory usage.
   - Neither flag: Loads the full dataset via `sopa.io.ome_tif()` and writes Zarr immediately (no blob detection).

2. **Spot/blob detection** (`cli.py` → `blob_detection.py`): Only runs when `--detect-blobs` is passed. For each configured channel, runs tiled Laplacian-of-Gaussian blob detection (`skimage.feature.blob_log`) parallelized with `dask.delayed` and a thread pool (`dask.compute(..., scheduler='threads')`). Overlapping tiles with overlap region filtering prevent duplicate detections. Results are assembled into a `PointsModel` and stored in `dataset["spots"]`.

3. **Initial Zarr write** (`cli.py`): The dataset (with or without spots) is written to disk with try/except error handling.

4. **Stardist segmentation** (`cli.py`): Reads back the Zarr (materialized checkpoint), creates image patches via `sopa.make_image_patches()`, detects GPU availability with `tf.config.list_physical_devices('GPU')` and warns if absent, renames channel coordinates to include `_ch_` suffixes (e.g., `DAPI_ch_0`), and runs `sopa.segmentation.stardist()` with the `2D_versatile_fluo` model. Only the first unique channel is passed to Stardist.

   **Segmented Zarr reuse**: Before running Stardist, checks if `_seg.zarr` already exists:
   - Has `tables['table']` → skips both Stardist and aggregation, jumps to clustering
   - Has `stardist_boundaries` but no table → skips Stardist, re-runs aggregation only
   - Missing or corrupt → full Stardist + aggregation

5. **Channel aggregation** (`cli.py`): Runs `sopa.aggregate()` to compute per-cell mean intensities for each channel (genes). If `"spots"` exists in the dataset points, passes `aggregate_genes=True`, `points_key='spots'`, and `gene_column='gene'` to assign spots to cells. Otherwise runs plain aggregation. Wraps the aggregation in `pd.option_context('future.infer_string', False)` to prevent ArrowStringArray errors on Zarr write. Writes segmented Zarr with try/except error handling. The segmented Zarr replaces `.zarr` with `_seg.zarr`.

6. **K-means clustering** (`cli.py` → `utils.py`): Extracts the intensity matrix from the AnnData table, standardizes with `StandardScaler`, runs k-means, and attaches cluster labels to `sdata.tables['table'].obs`.

7. **Network analysis** (`cli.py` → `network_analysis.py`): Extracts centroids from the cell boundaries GeoDataFrame. Builds a `scipy.spatial.Delaunay` triangulation, filters edges by `max_edge_distance`. Constructs a `networkx.Graph` from the filtered edges and runs `nx.community.louvain_communities()` with the `community_resolution` parameter and fixed seed 42. Returns a dict with `cell_ids`, `community_labels`, and `cluster_labels` arrays. Embeds `kmeans_cluster` and `community` labels into both the `stardist_boundaries` GeoDataFrame and the AnnData table obs. The segmented Zarr is written at this stage (single write after all labels are computed).

8. **Spot-to-cell assignment** (`cli.py` → `utils.py`): Spatial join of spots to cell boundaries using GeoPandas (`gpd.sjoin` with `predicate='within'`). Returns `None` if no spots are present in the dataset (no guard needed in `cli.py` — `export_to_qupath` handles `None`).

9. **QuPath GeoJSON export** (`cli.py` → `export.py`): Exports cell boundaries and spot detections as GeoJSON features with community/cluster assignments and intensity measurements for QuPath visualization. Takes `cell_ids`, `community_labels`, and `cluster_labels` as direct arrays. Output is always `./qupath_export.geojson` (hardcoded).

10. **Spatial neighborhood analysis** (`cli.py`): Computes spatial neighbors radius graph, mean hop distance heatmap between clusters (`cell_type_to_cell_type.png`), UMAP embedding with k-means coloring (`umap_kmeans_cluster.png`), and Leiden clustering with `igraph` backend (`umap_leiden.png`). UMAP plots are saved to `--plot_dir`.

### CLI arguments

All analysis parameters are exposed as command-line flags with sensible defaults:

| Flag | Default | Description |
|---|---|---|
| `-i`, `--input_file` | *(required)* | Path to input OME-TIFF image |
| `-o`, `--output_file` | *(required)* | Path to output Zarr (`.zarr` appended if missing) |
| `-p`, `--plot_dir` | `.` | Output directory for plots |
| `--detect-blobs` | — | Enable RNA spot blob detection on specified channels (default: off) |
| `--use-gpu` | — | Force GPU usage for Stardist (auto-detected by default) |
| `--channels` | `9,10,11,12` | Comma-separated channel indices for blob detection |
| `--thresholds` | `0.01,0.1,0.1,0.1` | Comma-separated blob detection thresholds (one per channel) |
| `--tile-size` | `2048` | Tile size for blob detection |
| `--overlap` | `50` | Tile overlap for blob detection |
| `--workers` | `14` | Worker threads for blob detection |
| `--min-sigma` | `2` | Minimum blob radius for spot detection |
| `--max-sigma` | `5` | Maximum blob radius for spot detection |
| `--num-sigma` | `5` | Number of sigma steps for blob detection |
| `--n-clusters` | `10` | Number of k-means clusters |
| `--community-resolution` | `0.1` | Louvain community detection resolution |
| `--max-edge-distance` | `1000` | Max edge distance for Delaunay network |
| `--radius-min` | `0` | Min radius for spatial neighbors graph |
| `--radius-max` | `1000` | Max radius for spatial neighbors graph |
| `--resume-from` | — | Path to existing Zarr to resume from (skips image loading and spot detection) |
| `--geojson-path` | `./qupath_export.geojson` | Output path for QuPath GeoJSON |
| `--fig-size` | `20` | Figure size for plots |
| `--font-size` | `20` | Font size for plots |
| `--axes-linewidth` | `3` | Axes line width for plots |

## Key gotchas

- **`TF_USE_LEGACY_KERAS='1'`**: Set in `run.py` before imports. TF >=2.16 defaults to Keras 3, which compiles Stardist's model with XLA JIT, triggering a cuDNN autotuner failure on 1x1 convolutions with CUDA 12/cuDNN 9. Legacy Keras 2 uses the non-XLA cuDNN path and retains GPU acceleration. Requires the `tf_keras` pip package.

- **System-specific shared library**: `run.py` preloads the pixi environment's `libstdc++.so.6` (resolved relative to the script's `.pixi/` directory) to avoid ABI conflicts with the system library. If the file doesn't exist (e.g., on a non-pixi setup), it silently skips. Note: this preload code is duplicated in `cli.py` (lines 3-9) so that `cli.py` can also be run standalone via `python -m cellsurvey.cli`.

- **GPU detection and `--use-gpu` flag**: GPU is detected at runtime via `tf.config.list_physical_devices('GPU')`. If no GPU is found and `--use-gpu` is not set, a warning is printed but execution continues — Stardist will run on CPU and be very slow. The `--use-gpu` flag forces GPU backend even without auto-detection. Both paths set `sopa.settings.parallelization_backend = None`.

- **Zarr write-then-read pattern**: The pipeline writes the initial Zarr to disk then immediately reads it back before continuing to segmentation. This is intentional — it materializes the spots-including dataset as a clean checkpoint.

- **Segmented Zarr incremental reuse**: The segmented Zarr (`_seg.zarr`) is checked before Stardist. If it already has a valid table, both Stardist and aggregation are skipped (jump to clustering). If it has boundaries but no table, only aggregation is re-run. If the file exists but can't be read, it's deleted via `shutil.rmtree` and a full re-run proceeds. This enables crash recovery at finer granularity than `--resume-from`.

- **Pandas 2.x + anndata ArrowStringArray compatibility**: Two workarounds in the aggregation stage:
  1. `pd.option_context('future.infer_string', False)` wraps the `sopa.aggregate()` call to force plain object dtype for strings, avoiding ArrowStringArray which can't be written to Zarr backing stores.
  2. `obs.index.astype(str)` is called on the AnnData table's obs index after clustering to force plain string dtype (same ArrowStringArray issue).

- **Channel naming**: Channel names with suffixes like `_ch_0`, `_ch_1` are constructed in the Stardist stage. The `remove_channel_suffix()` function strips these for clean column headers in the intensity DataFrame.

- **Duplicate column handling**: After stripping channel suffixes, duplicate column names are dropped keeping the first occurrence. This happens because channels with the same name at different indices collapse to the same column name.

- **`.zarr` suffix is auto-appended**: If `-o` doesn't end in `.zarr`, it's appended automatically. The segmented output replaces `.zarr` with `_seg.zarr`. This applies to both `--output_file` and `--resume-from` values.

- **Blob detection must be explicitly enabled**: The `--detect-blobs` flag is required to perform RNA spot/blob detection. Without it, the pipeline writes the raw OME-TIFF as a Zarr and proceeds directly to segmentation/aggregation without any spot data. The `sopa.aggregate()` call checks for `"spots" in dataset.points` to decide whether to aggregate genes.

- **Louvain community detection uses `networkx`**: `nx.community.louvain_communities()` with fixed seed 42 (for reproducibility). Requires `networkx>=3.4` in `pixi.toml`. The resolution parameter from `--community-resolution` is passed directly.

- **Zarr writes have basic error handling**: Both Zarr write points are wrapped in try/except — if a write fails, the error and path are printed to stderr and the script exits with code 1.

- **`--resume-from` enables crash recovery**: If a Zarr already exists at the expected path, you can skip image loading and spot detection and resume from Stardist segmentation. When used with the segmented Zarr reuse logic, this provides two levels of checkpoint restart.

- **Channel subset loading**: The full multichannel image is not loaded into memory for blob detection. Only the channels specified by `--channels` are loaded via `dask_image.imread`, reducing memory footprint. The full image is available on disk via the Zarr for Stardist segmentation.

- **AnnData `.X` can be sparse or dense**: `sopa.aggregate()` may produce either a scipy sparse matrix or a dense numpy array depending on the input data size and sopa version. The intensity extraction at `cli.py:238` handles both with `hasattr(measurements.X, 'toarray')`. Never assume `.X` is sparse.

- **Hardcoded output paths**: The GeoJSON export always writes to `./qupath_export.geojson` regardless of the `-o` or `-p` flags.

- **Matplotlib rcParams are set twice**: Global `font.size=20` and `axes.linewidth=3` at the start of `main()`, then overridden to `font.size=10` and `axes.linewidth=2` before the heatmap/UMAP plots. The `Agg` non-interactive backend is set at import time (`matplotlib.use('Agg')` before `import matplotlib.pyplot as plt`) to prevent plot windows from appearing on headless systems.

- **Plot directory auto-created**: `os.makedirs(args.plot_dir, exist_ok=True)` is called at the start of `main()` to ensure the output directory exists before any plots are saved.

- **`assign_spots_to_cells` returns `None` when no spots exist**: If `spots_key` is not in `spatial_data.points` (no `--detect-blobs` used, or no spots detected), the function returns `None` and `export_to_qupath` skips spot export.

- **Leiden clustering uses `igraph` backend**: `sc.tl.leiden(adata, flavor='igraph', n_iterations=2, directed=False)` — orders of magnitude faster than `leidenalg` for large datasets. Requires `python-igraph` in `pixi.toml` pypi dependencies. `show=False` passed to `sc.pl.umap()` to prevent `plt.show()` calls on the `Agg` backend.

## Code patterns and conventions

- Modular package structure under `cellsurvey/` — functions grouped by concern (blob detection, network analysis, export, utilities, CLI orchestration)
- Matplotlib global rcParams are set inside `cli.main()` at startup; the `Agg` non-interactive backend is forced at import time to prevent display on headless systems
- Random seeds (42) are used at multiple points for reproducibility (k-means clustering, Louvain communities)
- Print-based logging with no logging framework
- Dask is used for parallel blob detection but the scheduler is explicitly set to `'threads'` (not the default multiprocessing)
- NumPy, pandas, GeoPandas, and AnnData/Scanpy are the primary data structures
- `run.py` is the entry-point shim — all logic lives in the `cellsurvey/` package modules. However, `cli.py` can also be run standalone (`python -m cellsurvey.cli`) since it duplicates the libstdc++ preload.
- `export_to_qupath` takes `cell_ids`, `community_labels`, `cluster_labels`, `sdata`, `intensity_df`, and `spots_with_cells` as explicit parameters (no implicit closure on module globals)
- The pipeline writes Zarr files at two points: the initial Zarr after image loading (and optionally blob detection), and the segmented Zarr after Stardist + aggregation. Both writes are wrapped in try/except and exit with code 1 on failure.
- The Zarr write-then-immediate-read pattern between stages 3 and 4 materializes a clean checkpoint. If segmented Zarr reuse kicks in (stage 4), the read of the initial Zarr is skipped entirely.
- The `sopa.segmentation.stardist()` call passes only `unique_channels[0]` as the channels argument, not all channel names.
- Dynamic imports inside except blocks: `import shutil` is imported inside the segmented Zarr corruption handler to avoid pulling it in unnecessarily.
- `run_network_analysis()` in `network_analysis.py` returns a plain dict with keys `cell_ids`, `community_labels`, `cluster_labels`.
- `run_network_analysis()` accepts `output_dir` parameter for plot output (set to `args.plot_dir` from `cli.py`).

## File structure

```
.
├── run.py                            # Entry-point shim: TF_USE_LEGACY_KERAS + libstdc++ guard + delegates to cli.main()
├── cellsurvey/
│   ├── __init__.py                   # Re-exports all public symbols
│   ├── cli.py                        # main() with argparse and pipeline orchestration
│   ├── blob_detection.py             # detect_blobs_in_tile, detect_blobs_tiled
│   ├── network_analysis.py           # run_network_analysis (scipy Delaunay + networkx Louvain)
│   ├── export.py                     # export_to_qupath
│   └── utils.py                      # remove_channel_suffix, cluster_data, assign_spots_to_cells, get_colors_for_communities
├── Dockerfile                        # Ubuntu 24.04 + pixi + GPU-ready container
├── pixi.toml                         # Pixi environment config (linux-64 only)
├── pixi.lock                         # Pixi lockfile (generated)
├── requirements.txt                  # Minimal pip requirements (sopa)
├── README.md                         # User-facing installation and usage docs
└── .pixi/                            # Pixi environment directory (gitignored)
```
