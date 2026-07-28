# AGENTS.md

## Project overview

SopaSpan is a single-script spatial biology/omics analysis pipeline combining [Sopa](https://gustaveroussy.github.io/sopa/) (segmentation, aggregation) and [MuSpAn](https://www.muspan.co.uk/) (network analysis, community detection). It processes multichannel microscopy images (OME-TIFF) into spatial data objects, segments nuclei with Stardist, detects RNA spots via blob detection, clusters cells with k-means, builds Delaunay networks, detects Louvain communities, and exports GeoJSON for QuPath visualization.

## Environment and package management

This project uses **pixi** (via `pixi.toml`) for environment management targeting `win-64`. The lockfile is `pixi.lock` (marked as binary/generated in `.gitattributes`).

Key dependency constraints:
- **Python < 3.11** (pixi), though the README documents Python 3.12 with conda — both work; the pixi config is the authoritative development environment
- **TensorFlow < 2.11** (for Stardist compatibility)
- **CUDA 11.2 / cuDNN 8.1.0** for GPU acceleration
- **MuSpAn** is distributed as a password-protected zip from `https://docs.muspan.co.uk/code/latest.zip` and installed from `./latest.zip` via `pixi.toml`

There is no Makefile, no CI/CD, and no tests. The only source file is `sopaspan.py`.

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

`sopaspan.py` is a single flat script (~550 lines) with no classes — just functions and a `__main__` block. The pipeline runs these stages sequentially:

1. **Image loading**: Reads the OME-TIFF via `BioImage` (from `bioio`) to get channel names, and via `sopa.io.ome_tif()` as a SpatialData dataset. Also loads the image array via `dask_image.imread`.

2. **Spot/blob detection** (lines 403–415): For each pre-configured channel index in `channel_index` (hardcoded as `[9, 10, 11, 12]` at line 378), runs tiled Laplacian-of-Gaussian blob detection (`skimage.feature.blob_log`) parallelized with `dask.delayed` and a thread pool (`dask.compute(..., scheduler='threads')`). Each channel gets its own threshold from the `thresholds` list (line 379). Overlapping tiles with overlap region filtering prevent duplicate detections. Results are assembled into a `PointsModel` and stored in `dataset["spots"]`.

3. **Initial Zarr write** (lines 425–429): The dataset (containing spots) is written to disk as `{zarr_path}.zarr`.

4. **Stardist segmentation** (lines 441–463): Reads back the Zarr, creates image patches via `sopa.make_image_patches()`, disables backend parallelization (`sopa.settings.parallelization_backend = None` to use GPU), fixes channel names, and runs `sopa.segmentation.stardist()` with the `2D_versatile_fluo` model.

5. **Channel aggregation** (lines 467): Runs `sopa.aggregate()` to compute per-cell mean intensities for each channel (genes), assigning spots to cells. Writes segmented Zarr as `{zarr_path}_seg.zarr`.

6. **K-means clustering** (lines 496–510): Extracts the intensity matrix from the AnnData table in the SpatialData object, standardizes with `StandardScaler`, runs k-means (n=10), and attaches cluster labels to `sdata.tables['table'].obs`.

7. **MuSpAn network analysis** (lines 518): Converts SpatialData to a MuSpAn domain (shapes as points), builds a Delaunay triangulation network, detects Louvain communities, and generates visualization plots.

8. **Spot-to-cell assignment** (lines 520): Spatial join of spots to cell boundaries using GeoPandas.

9. **QuPath GeoJSON export** (lines 522–523): Exports cell boundaries and spot detections as GeoJSON features with community/cluster assignments and intensity measurements for QuPath visualization.

10. **Spatial neighborhood analysis** (lines 528–545): Computes spatial neighbors radius graph, mean hop distance heatmap between clusters, UMAP embedding, and Leiden clustering via Scanpy.

### Hardcoded configuration

Several important values are hardcoded in `__main__` rather than exposed as CLI arguments:
- `channel_index = [9, 10, 11, 12]` — which channels to process for blob detection
- `thresholds = [0.01, 0.1, 0.1, 0.1]` — per-channel blob detection thresholds
- `tile_size = 2048`, `overlap = 50`, `n_workers = 14` — blob detection tiling params
- `n_clusters = 10` — k-means cluster count
- `comm_detect_res = 0.1` — Louvain community detection resolution
- `max_edge_distance = 1000` — Delaunay network max edge distance
- `radius = (0, 1000)` — spatial neighbors radius

## Key gotchas

- **System-specific shared library**: The very first line of `sopaspan.py` loads a C++ standard library from a hardcoded HPC path (`/nemo/stp/lm/working/barryd/hpc/pixi/sopaspan/.pixi/envs/sopaspan/lib/libstdc++.so.6`). This will fail on any non-HPC system. Remove or conditionally skip this line when running elsewhere.

- **GPU vs CPU**: `sopa.settings.parallelization_backend = None` at line 445 disables Sopas parallel backend, relying on GPU via TensorFlow. If GPU is unavailable, Stardist will fall back to CPU and be extremely slow.

- **Zarr write-then-read pattern**: The pipeline writes the initial Zarr to disk then immediately reads it back (line 435) before continuing to segmentation. This is intentional — it materializes the spots-including dataset as a clean checkpoint.

- **Channel naming**: Channel names with suffixes like `_ch_0`, `_ch_1` are constructed at line 453. The `remove_channel_suffix()` function strips these for clean column headers in the intensity DataFrame.

- **Duplicate column handling**: After stripping channel suffixes, duplicate column names are dropped keeping the first occurrence (line 497). This happens because channels with the same name at different indices collapse to the same column name.

- **`.zarr` suffix is auto-appended**: Providing `-o /path/to/output` creates `/path/to/output.zarr`, not `/path/to/output`.

- **MuSpAn dependency is not public**: `muspan` requires a username/password obtained from https://www.muspan.co.uk/get-the-code. The `latest.zip` in the repo root is the MuSpAn distribution file referenced in `pixi.toml` as a path dependency.

- **No error handling**: The script has no try/except blocks. Failures at any stage will crash with a stack trace. Stardist segmentation and blob detection are the most likely failure points (memory, GPU availability, image dimensions).

- **Large memory footprint**: The entire multichannel image is loaded into memory via `dask_image.imread` at line 386. For large whole-slide images, this may cause OOM errors.

## Code patterns and conventions

- Script style with global scope — no classes, no modules, no `if __name__ == '__main__'` isolation beyond the argparse block
- Matplotlib global rcParams are set at module level (lines 30–32)
- Random seeds (42) are used at multiple points for reproducibility
- Print-based logging with no logging framework
- Dask is used for parallel blob detection but the scheduler is explicitly set to `'threads'` (not the default multiprocessing)
- NumPy, pandas, GeoPandas, and AnnData/Scanpy are the primary data structures
- Default argument values in argparse point to HPC paths — these will not exist on other systems

## File structure

```
.
├── sopaspan.py          # The entire pipeline
├── pixi.toml            # Pixi environment config (primary)
├── pixi.lock            # Pixi lockfile (generated)
├── requirements.txt     # Minimal pip requirements (sopa, muspan)
├── latest.zip           # MuSpAn distribution (password-protected)
├── README.md            # User-facing installation and usage docs
├── .gitignore           # Ignores .pixi/* except config.toml
├── .gitattributes       # Marks pixi.lock as binary/YAML-generated
└── .pixi/               # Pixi environment directory (gitignored)
```

## Implementation roadmap

Grouped in dependency order — each chunk can ship independently. P0 blocks everything else.

### Chunk 1: Make pipeline runnable anywhere (P0) ✅ DONE

- **Remove/guard `ctypes.CDLL` call** (line 1–3): wrapped in `os.path.exists()` check so the HPC C++ library is only loaded when present.
- **Remove HPC defaults from argparse** (lines 366–372): `-i` and `-o` are now required with no defaults.
- **Detect GPU and warn on fallback** (line 445): checks `tf.config.list_physical_devices('GPU')` before setting backend. Prints warning if no GPU found.

### Chunk 2: Memory safety (P1) ✅ DONE

- **Lazy-load the image** (line 386): `dask_image.imread(imagepath)[channel_index]` now only loads the 4 channels needed for blob detection, not all channels. `BioImage` is still used for lightweight metadata (channel names).

### Chunk 3: Correctness fixes (P1)

- **Fix `pd.concat` in loop** (lines 401–415): collect DataFrames in a list, call `pd.concat(frames)` once at the end.
- **Remove dead Spotiflow code** (lines 388–399): delete the commented block.
- **Remove useless `adata_subset` alias** (line 527): delete the variable and the commented `np.random.choice` line. Use `adata` directly.
- **Save UMAP/Leiden plots to `plot_dir`** (lines 543–545): wrap `sc.pl.umap` calls with `plt.savefig(os.path.join(args.plot_dir, ...))` and `plt.close()`.
- **Drop the ignored `plt.figure()` before `sc.pl.umap`** (line 543).

### Chunk 4: Remove module-level global dependency (P2)

- **Pass `intensity_df` as a parameter** to `export_to_qupath` instead of relying on implicit closure (line 174). Add it as a required argument. This also means `sdata` accessed at line 153 needs to become a parameter.

### Chunk 5: Expose hardcoded parameters as CLI flags (P2)

- Add argparse arguments for: `--channels`, `--thresholds`, `--tile-size`, `--overlap`, `--workers`, `--n-clusters`, `--community-resolution`, `--max-edge-distance`, `--radius-min`, `--radius-max`.
- Sensible defaults match current hardcoded values.
- `--channels` and `--thresholds` accept comma-separated lists parsed with `list(map(int, ...))` / `list(map(float, ...))`.
- Make `-o` behavior explicit: if the path doesn't end in `.zarr`, append it, and document this.

### Chunk 6: Basic error resilience (P1)

- **Wrap the two Zarr write blocks** (lines 429, 473) in try/except — if the write fails, print the error and the path so the user can debug disk space/permissions.
- **Add a `--resume-from` flag**: if a Zarr already exists at the expected path, skip image loading and spot detection and resume from Stardist segmentation. This avoids re-running the most expensive steps after a crash in later stages.

### Chunk 7: Structural refactor (P3)

- Split `sopaspan.py` into:
  - `sopaspan/blob_detection.py` — `detect_blobs_in_tile`, `detect_blobs_tiled`
  - `sopaspan/muspan_workflow.py` — `run_muspan`, `get_colors_for_communities`
  - `sopaspan/export.py` — `export_to_qupath`
  - `sopaspan/cli.py` — `__main__` and argparse
  - `sopaspan/utils.py` — `remove_channel_suffix`, `cluster_data`, `assign_spots_to_cells`
- Keep `sopaspan.py` as a re-export shim for backward compatibility.
