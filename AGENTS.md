# AGENTS.md

## Project overview

SopaSpan is a single-script spatial biology/omics analysis pipeline combining [Sopa](https://gustaveroussy.github.io/sopa/) (segmentation, aggregation) and [MuSpAn](https://www.muspan.co.uk/) (network analysis, community detection). It processes multichannel microscopy images (OME-TIFF) into spatial data objects, segments nuclei with Stardist, detects RNA spots via blob detection, clusters cells with k-means, builds Delaunay networks, detects Louvain communities, and exports GeoJSON for QuPath visualization.

## Environment and package management

This project uses **pixi** (via `pixi.toml`) for environment management targeting `win-64`, `linux-64`, and `osx-64`. The lockfile is `pixi.lock` (marked as binary/generated in `.gitattributes`).

Key dependency constraints:
- **Python < 3.11** (pixi), though the README documents Python 3.12 with conda — both work; the pixi config is the authoritative development environment
- **TensorFlow < 2.11** (for Stardist compatibility)
- **CUDA 11.2 / cuDNN 8.1.0** for GPU acceleration
- **MuSpAn** is distributed as a password-protected zip from `https://docs.muspan.co.uk/code/latest.zip` and installed from `./latest.zip` via `pixi.toml`

There is no Makefile, no CI/CD, and no tests. The source code is split across 6 files under the `sopaspan/` package, with `sopaspan.py` as the entry-point shim.

**TODO**: Set up linting and formatting (Ruff, mypy) with a `pyproject.toml` config and pre-commit hooks.

## Commands

**Development environment setup (pixi):**
```bash
pixi install
```

**Run the pipeline:**
```bash
python sopaspan.py -i <input_tiff> -o <output_zarr_prefix> -p <plot_output_dir>
```

Three arguments:
- `-i`: Path to input OME-TIFF image
- `-o`: Path prefix for output Zarr file (`.zarr` suffix is appended automatically)
- `-p`: Directory for output plots

There is no build step, no test command, and no linting configured.

## Architecture and data flow

The pipeline is split into modules under the `sopaspan/` package. `sopaspan.py` is a shim that guards the libstdc++ preload and delegates to `sopaspan.cli.main()`. The pipeline runs these stages sequentially:

1. **Image loading** (`cli.py`): Reads the OME-TIFF via `BioImage` (from `bioio`) to get channel names, and via `sopa.io.ome_tif()` as a SpatialData dataset. Only loads the subset of channels needed for blob detection via `dask_image.imread` to limit memory usage.

2. **Spot/blob detection** (`cli.py` → `blob_detection.py`): For each configured channel, runs tiled Laplacian-of-Gaussian blob detection (`skimage.feature.blob_log`) parallelized with `dask.delayed` and a thread pool (`dask.compute(..., scheduler='threads')`). Overlapping tiles with overlap region filtering prevent duplicate detections. Results are assembled into a `PointsModel` and stored in `dataset["spots"]`.

3. **Initial Zarr write** (`cli.py`): The dataset (containing spots) is written to disk with try/except error handling.

4. **Stardist segmentation** (`cli.py`): Reads back the Zarr, creates image patches via `sopa.make_image_patches()`, detects GPU availability with `tf.config.list_physical_devices('GPU')` and warns if absent, fixes channel names, and runs `sopa.segmentation.stardist()` with the `2D_versatile_fluo` model.

5. **Channel aggregation** (`cli.py`): Runs `sopa.aggregate()` to compute per-cell mean intensities for each channel (genes), assigning spots to cells. Writes segmented Zarr with try/except error handling.

6. **K-means clustering** (`cli.py` → `utils.py`): Extracts the intensity matrix from the AnnData table, standardizes with `StandardScaler`, runs k-means, and attaches cluster labels to `sdata.tables['table'].obs`.

7. **MuSpAn network analysis** (`cli.py` → `muspan_workflow.py`): Converts SpatialData to a MuSpAn domain (shapes as points), builds a Delaunay triangulation network, detects Louvain communities, and generates visualization plots.

8. **Spot-to-cell assignment** (`cli.py` → `utils.py`): Spatial join of spots to cell boundaries using GeoPandas.

9. **QuPath GeoJSON export** (`cli.py` → `export.py`): Exports cell boundaries and spot detections as GeoJSON features with community/cluster assignments and intensity measurements for QuPath visualization.

10. **Spatial neighborhood analysis** (`cli.py`): Computes spatial neighbors radius graph, mean hop distance heatmap between clusters, UMAP embedding, and Leiden clustering via Scanpy. UMAP plots are saved to the plot output directory.

### CLI arguments

All analysis parameters are exposed as command-line flags with sensible defaults:

| Flag | Default | Description |
|---|---|---|
| `-i`, `--input_file` | *(required)* | Path to input OME-TIFF image |
| `-o`, `--output_file` | *(required)* | Path to output Zarr (`.zarr` appended if missing) |
| `-p`, `--plot_dir` | `.` | Output directory for plots |
| `--channels` | `9,10,11,12` | Comma-separated channel indices for blob detection |
| `--thresholds` | `0.01,0.1,0.1,0.1` | Comma-separated blob detection thresholds |
| `--tile-size` | `2048` | Tile size for blob detection |
| `--overlap` | `50` | Tile overlap for blob detection |
| `--workers` | `14` | Worker threads for blob detection |
| `--n-clusters` | `10` | Number of k-means clusters |
| `--community-resolution` | `0.1` | Louvain community detection resolution |
| `--max-edge-distance` | `1000` | Max edge distance for Delaunay network |
| `--radius-min` | `0` | Min radius for spatial neighbors graph |
| `--radius-max` | `1000` | Max radius for spatial neighbors graph |
| `--resume-from` | — | Path to existing Zarr to resume from (skips image loading and spot detection) |

## Key gotchas

- **System-specific shared library**: `sopaspan.py` preloads the pixi environment's `libstdc++.so.6` (resolved relative to the script's `.pixi/` directory) to avoid ABI conflicts with the system library. If the file doesn't exist (e.g., on Windows or a non-pixi setup), it silently skips.

- **GPU vs CPU**: GPU is detected at runtime via `tf.config.list_physical_devices('GPU')`. If no GPU is found, a warning is printed but execution continues — Stardist will run on CPU and be very slow.

- **Zarr write-then-read pattern**: The pipeline writes the initial Zarr to disk then immediately reads it back before continuing to segmentation. This is intentional — it materializes the spots-including dataset as a clean checkpoint.

- **Channel naming**: Channel names with suffixes like `_ch_0`, `_ch_1` are constructed in the Stardist stage. The `remove_channel_suffix()` function strips these for clean column headers in the intensity DataFrame.

- **Duplicate column handling**: After stripping channel suffixes, duplicate column names are dropped keeping the first occurrence. This happens because channels with the same name at different indices collapse to the same column name.

- **`.zarr` suffix is auto-appended**: If `-o` doesn't end in `.zarr`, it's appended automatically. The segmented output replaces `.zarr` with `_seg.zarr`.

- **MuSpAn dependency is not public**: `muspan` requires a username/password obtained from https://www.muspan.co.uk/get-the-code. The `latest.zip` in the repo root is the MuSpAn distribution file referenced in `pixi.toml` as a path dependency.

**TODO**: Evaluate replacing MuSpAn with an alternative network analysis / community detection library (e.g., Scanpy's native spatial tools, Squidpy, or NetworkX + GeoPandas) to eliminate the private-dependency friction. The pipeline currently uses MuSpAn only for Delaunay triangulation + Louvain community detection + visualization — all of which have equivalents in open, pip-installable libraries.

- **Zarr writes have basic error handling**: Both Zarr write points are wrapped in try/except — if a write fails, the error and path are printed to stderr and the script exits with code 1.

- **`--resume-from` enables crash recovery**: If a Zarr already exists at the expected path, you can skip image loading and spot detection and resume from Stardist segmentation. This saves re-running expensive steps after a crash in later stages.

- **Channel subset loading**: The full multichannel image is not loaded into memory. Only the channels specified by `--channels` are loaded for blob detection, reducing memory footprint for large images. The full image is still available on disk via the Zarr for Stardist segmentation.

## Code patterns and conventions

- Modular package structure under `sopaspan/` — functions grouped by concern (blob detection, MuSpAn workflow, export, utilities, CLI orchestration)
- Matplotlib global rcParams are set inside `cli.main()` at startup
- Random seeds (42) are used at multiple points for reproducibility
- Print-based logging with no logging framework
- Dask is used for parallel blob detection but the scheduler is explicitly set to `'threads'` (not the default multiprocessing)
- NumPy, pandas, GeoPandas, and AnnData/Scanpy are the primary data structures
- `sopaspan.py` is the entry-point shim — all logic lives in the `sopaspan/` package modules
- `export_to_qupath` takes `sdata` and `intensity_df` as explicit parameters (no implicit closure on module globals)

## File structure

```
.
├── sopaspan.py                       # Entry-point shim: libstdc++ guard + delegates to cli.main()
├── sopaspan/
│   ├── __init__.py                   # Re-exports all public symbols
│   ├── cli.py                        # main() with argparse and pipeline orchestration
│   ├── blob_detection.py             # detect_blobs_in_tile, detect_blobs_tiled
│   ├── muspan_workflow.py            # run_muspan, get_colors_for_communities, fig_kwargs
│   ├── export.py                     # export_to_qupath
│   └── utils.py                      # remove_channel_suffix, cluster_data, assign_spots_to_cells
├── pixi.toml                         # Pixi environment config (primary)
├── pixi.lock                         # Pixi lockfile (generated)
├── requirements.txt                  # Minimal pip requirements (sopa, muspan)
├── latest.zip                        # MuSpAn distribution (password-protected)
├── README.md                         # User-facing installation and usage docs
├── .gitignore                        # Ignores .pixi/* except config.toml
├── .gitattributes                    # Marks pixi.lock as binary/YAML-generated
└── .pixi/                            # Pixi environment directory (gitignored)
```
