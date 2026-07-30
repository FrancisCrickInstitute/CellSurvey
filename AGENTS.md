# AGENTS.md

## Project overview

CellSurvey is a spatial biology/omics analysis pipeline combining [Sopa](https://gustaveroussy.github.io/sopa/) (segmentation, aggregation) and [MuSpAn](https://www.muspan.co.uk/) (network analysis, community detection). It processes multichannel microscopy images (OME-TIFF) into spatial data objects, segments nuclei with Stardist, detects RNA spots via blob detection, clusters cells with k-means, builds Delaunay networks, detects Louvain communities, and exports GeoJSON for QuPath visualization.

## Environment and package management

This project uses **pixi** (via `pixi.toml`) for environment management targeting `win-64`, `linux-64`, and `osx-64`. The lockfile is `pixi.lock` (marked as binary/generated in `.gitattributes`).

Key dependency constraints:
- **Python < 3.11** on `win-64`, `>=3.11,<3.12` on `linux-64` and `osx-64` — the README documents Python 3.12 with conda; both work; the pixi config is the authoritative development environment
- **TensorFlow**: Constraints are platform-specific in `pixi.toml` — `<2.11` on `win-64` (with `and-cuda` extras), `>=2.18` on `linux-64`, `>=2.11` on `osx-64`
- **CUDA 11.2 / cuDNN 8.1.0** for GPU acceleration (on win-64 and linux-64; osx-64 has no CUDA deps)
- **MuSpAn** is distributed as a password-protected zip from `https://docs.muspan.co.uk/code/latest.zip` and installed from `./latest.zip` via `pixi.toml`

There is no Makefile, no CI/CD, and no tests. The source code is split across 6 files under the `cellsurvey/` package, with `run.py` as the entry-point shim.

**TODO**: Set up linting and formatting (Ruff, mypy) with a `pyproject.toml` config and pre-commit hooks.

## Commands

**Development environment setup (pixi):**
```bash
pixi install
```

**Run the pipeline:**
```bash
python run.py -i <input_tiff> -o <output_zarr_prefix> -p <plot_output_dir>
```

Three required arguments:
- `-i`: Path to input OME-TIFF image
- `-o`: Path prefix for output Zarr file (`.zarr` suffix is appended automatically)
- `-p`: Directory for output plots

With pixi:
```bash
pixi run python run.py -i <input_tiff> -o <output_zarr_prefix> -p <plot_output_dir>
```

There is no build step, no test command, and no linting configured.

## Architecture and data flow

The pipeline is split into modules under the `cellsurvey/` package. `run.py` is a shim that guards the libstdc++ preload and delegates to `cellsurvey.cli.main()`. The pipeline runs these stages sequentially:

1. **Image loading** (`cli.py`): Reads the OME-TIFF via `BioImage` (from `bioio`) to get channel names, and via `sopa.io.ome_tif()` as a SpatialData dataset. Three code paths at this stage:
   - `--resume-from`: Skips image loading entirely — reads channel names from the existing image via `BioImage` but loads the Zarr directly.
   - `--detect-blobs`: Loads only the subset of channels needed for blob detection via `dask_image.imread` to limit memory usage.
   - Neither flag: Loads the full dataset via `sopa.io.ome_tif()` and writes Zarr immediately (no blob detection).

2. **Spot/blob detection** (`cli.py` → `blob_detection.py`): Only runs when `--detect-blobs` is passed. For each configured channel, runs tiled Laplacian-of-Gaussian blob detection (`skimage.feature.blob_log`) parallelized with `dask.delayed` and a thread pool (`dask.compute(..., scheduler='threads')`). Overlapping tiles with overlap region filtering prevent duplicate detections. Results are assembled into a `PointsModel` and stored in `dataset["spots"]`.

3. **Initial Zarr write** (`cli.py`): The dataset (with or without spots) is written to disk with try/except error handling.

4. **Stardist segmentation** (`cli.py`): Reads back the Zarr (materialized checkpoint), creates image patches via `sopa.make_image_patches()`, detects GPU availability with `tf.config.list_physical_devices('GPU')` and warns if absent, renames channel coordinates to include `_ch_` suffixes (e.g., `DAPI_ch_0`), and runs `sopa.segmentation.stardist()` with the `2D_versatile_fluo` model. Only the first unique channel is passed to Stardist.

5. **Channel aggregation** (`cli.py`): Runs `sopa.aggregate()` to compute per-cell mean intensities for each channel (genes). If `"spots"` exists in the dataset points, passes `aggregate_genes=True`, `points_key='spots'`, and `gene_column='gene'` to assign spots to cells. Otherwise runs plain aggregation. Writes segmented Zarr with try/except error handling. The segmented Zarr replaces `.zarr` with `_seg.zarr`.

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
| `--detect-blobs` | — | Enable RNA spot blob detection on specified channels (default: off) |
| `--channels` | `9,10,11,12` | Comma-separated channel indices for blob detection |
| `--thresholds` | `0.01,0.1,0.1,0.1` | Comma-separated blob detection thresholds (one per channel) |
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

- **`POT_BACKEND=numpy` workaround**: `run.py` sets `POT_BACKEND=numpy` before any imports from `cellsurvey`. MuSpAn depends on POT (Python Optimal Transport) which by default tries to import TensorFlow as a backend. Since TF 2.10 is compiled against numpy 1.x, it crashes when numpy 2.x is present. Forcing the numpy backend avoids the TF import entirely — TF is only needed later for Stardist segmentation, where it initializes after numpy is already loaded. This variable MUST be set before importing `cellsurvey.cli`.

- **System-specific shared library**: `run.py` preloads the pixi environment's `libstdc++.so.6` (resolved relative to the script's `.pixi/` directory) to avoid ABI conflicts with the system library. If the file doesn't exist (e.g., on Windows or a non-pixi setup), it silently skips. Note: this preload code is duplicated in `cli.py` (line 3-9) so that `cli.py` can also be run standalone.

- **GPU vs CPU**: GPU is detected at runtime via `tf.config.list_physical_devices('GPU')`. If no GPU is found, a warning is printed but execution continues — Stardist will run on CPU and be very slow.

- **Zarr write-then-read pattern**: The pipeline writes the initial Zarr to disk then immediately reads it back before continuing to segmentation. This is intentional — it materializes the spots-including dataset as a clean checkpoint.

- **Channel naming**: Channel names with suffixes like `_ch_0`, `_ch_1` are constructed in the Stardist stage. The `remove_channel_suffix()` function strips these for clean column headers in the intensity DataFrame.

- **Duplicate column handling**: After stripping channel suffixes, duplicate column names are dropped keeping the first occurrence. This happens because channels with the same name at different indices collapse to the same column name.

- **`.zarr` suffix is auto-appended**: If `-o` doesn't end in `.zarr`, it's appended automatically. The segmented output replaces `.zarr` with `_seg.zarr`. This applies to both `--output_file` and `--resume-from` values.

- **Blob detection must be explicitly enabled**: The `--detect-blobs` flag is required to perform RNA spot/blob detection. Without it, the pipeline writes the raw OME-TIFF as a Zarr and proceeds directly to segmentation/aggregation without any spot data. The `sopa.aggregate()` call checks for `"spots" in dataset.points` to decide whether to aggregate genes.

- **MuSpAn dependency is not public**: `muspan` requires a username/password obtained from https://www.muspan.co.uk/get-the-code. The `latest.zip` in the repo root is the MuSpAn distribution file referenced in `pixi.toml` as a path dependency.

**TODO**: Evaluate replacing MuSpAn with an alternative network analysis / community detection library (e.g., Scanpy's native spatial tools, Squidpy, or NetworkX + GeoPandas) to eliminate the private-dependency friction. The pipeline currently uses MuSpAn only for Delaunay triangulation + Louvain community detection + visualization — all of which have equivalents in open, pip-installable libraries.

- **Zarr writes have basic error handling**: Both Zarr write points are wrapped in try/except — if a write fails, the error and path are printed to stderr and the script exits with code 1.

- **`--resume-from` enables crash recovery**: If a Zarr already exists at the expected path, you can skip image loading and spot detection and resume from Stardist segmentation. This saves re-running expensive steps after a crash in later stages.

- **Channel subset loading**: The full multichannel image is not loaded into memory. Only the channels specified by `--channels` are loaded for blob detection, reducing memory footprint for large images. The full image is still available on disk via the Zarr for Stardist segmentation.

## Code patterns and conventions

- Modular package structure under `cellsurvey/` — functions grouped by concern (blob detection, MuSpAn workflow, export, utilities, CLI orchestration)
- Matplotlib global rcParams are set inside `cli.main()` at startup
- Random seeds (42) are used at multiple points for reproducibility
- Print-based logging with no logging framework
- Dask is used for parallel blob detection but the scheduler is explicitly set to `'threads'` (not the default multiprocessing)
- NumPy, pandas, GeoPandas, and AnnData/Scanpy are the primary data structures
- `run.py` is the entry-point shim — all logic lives in the `cellsurvey/` package modules. However, `cli.py` can also be run standalone (`python -m cellsurvey.cli`) since it duplicates the libstdc++ preload.
- `export_to_qupath` takes `sdata`, `intensity_df`, and `spots_with_cells` as explicit parameters (no implicit closure on module globals)
- The pipeline writes Zarr files at two points: the initial Zarr after image loading (and optionally blob detection), and the segmented Zarr after Stardist + aggregation. Both writes are wrapped in try/except and exit with code 1 on failure. The Zarr write-then-immediate-read pattern between stages 3 and 4 materializes a clean checkpoint.
- The `sopa.segmentation.stardist()` call passes only `unique_channels[0]` as the channels argument, not all channel names.

## File structure

```
.
├── run.py                            # Entry-point shim: libstdc++ guard + delegates to cli.main()
├── cellsurvey/
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
