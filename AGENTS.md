# AGENTS.md

## Project overview

CellSurvey is a spatial biology/omics analysis pipeline built on [Sopa](https://gustaveroussy.github.io/sopa/) (segmentation, aggregation). It processes multichannel microscopy images (OME-TIFF) into spatial data objects, segments nuclei with Stardist, detects RNA spots via blob detection, clusters cells with k-means, builds Delaunay networks, detects Louvain communities, and exports GeoJSON for QuPath visualization.

**MuSpAn removed (complete):** MuSpAn has been fully replaced with open-source libraries. Delaunay triangulation uses `scipy.spatial.Delaunay`, Louvain community detection uses `networkx`, and visualisation uses matplotlib. No private dependencies remain.

## Environment and package management

This project uses **pixi** (via `pixi.toml`) for environment management targeting `linux-64` (GPU) and `win-64` (CPU-only, for local testing). The lockfile is `pixi.lock`.

Key dependency constraints:
- **Python**: `>=3.12, <3.13` on `linux-64`
- **TensorFlow**: `>=2.18` on `linux-64` (with `and-cuda` extras for GPU)
- **CUDA/cuDNN**: TF >=2.18 bundles its own CUDA 12/cuDNN 9 libraries; no separate conda CUDA/cuDNN packages are needed
- **`tf_keras`** is a required pypi dependency — TF >=2.16 defaults to Keras 3, but Stardist needs legacy Keras 2 API to avoid cuDNN autotuner failures on CUDA 12/cuDNN 9
- **`scipy` (`>=1.14, <2`)** and **`networkx` (`>=3.4, <4`)** for Delaunay triangulation and Louvain community detection (replacing MuSpAn)
- **`python-igraph`** for fast Leiden clustering in Scanpy's spatial neighborhood analysis
- **`bioio` (`>=3.4`) and `bioio-ome-tiff`** for reading channel names from OME-TIFF metadata (via `BioImage`); `setuptools` is pinned as a pypi dependency
- **Windows is CPU-only for local testing** via pixi (`win-64`, TensorFlow CPU wheel — no GPU, no `and-cuda`). **macOS is not supported** — only `linux-64` and `win-64` are in the platforms list.
- **Sopa is a single-maintainer dependency risk**: Sopa is maintained almost entirely by one person (Quentin Blampey), who has left academia. CellSurvey currently uses Sopa for image reading, patching, the Stardist wrapper, aggregation, and spatial-neighbour analysis. **Strategic aim**: gradually reduce Sopa coupling — own the segmentation/aggregation/spatial steps directly where practical — while continuing to use it in the short term. The cell-segmentation work (and eventual Cellpose integration) is a natural seam to start decoupling.

A **Dockerfile** is provided: Ubuntu 24.04 base, installs pixi, copies `pixi.toml`, sets `TF_USE_LEGACY_KERAS=1`, entrypoint is `pixi run python run.py`.

There is no Makefile, no CI/CD, and no tests. The source code is split across 6 files under the `cellsurvey/` package, with `run.py` as the entry-point shim.

**TODO**: Set up linting and formatting (Ruff, mypy) with a `pyproject.toml` config and pre-commit hooks.

**TODO (docs)**: Add a *tutorial* page and more visual illustrations to the MkDocs site. All eight nav pages are now written (`index`, `installation`, `getting-started`, `pipeline`, `parameters`, `outputs`, `visualization`, `faq`), but a step-by-step tutorial (with example inputs/outputs) and richer diagrams/example-plot images would help first-time users. Revisit `docs/` once example datasets and output images are available.

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

   **Segmented Zarr reuse**: Only when `--resume-from` is specified, checks if `_seg.zarr` already exists before running Stardist:
   - Has `tables['table']` → skips both Stardist and aggregation, jumps to clustering
   - Has `stardist_boundaries` but no table → skips Stardist, re-runs aggregation only
   - Missing or corrupt → full Stardist + aggregation

5. **Channel aggregation** (`cli.py`): Runs `sopa.aggregate()` to compute per-cell mean intensities for each channel (genes). If `"spots"` exists in the dataset points, passes `aggregate_genes=True`, `points_key='spots'`, and `gene_column='gene'` to assign spots to cells. Otherwise runs plain aggregation. Wraps the aggregation in `pd.option_context('future.infer_string', False)` to prevent ArrowStringArray errors on Zarr write. Writes segmented Zarr with try/except error handling. The segmented Zarr replaces `.zarr` with `_seg.zarr`.

6. **K-means clustering** (`cli.py` → `utils.py`): Extracts the intensity matrix from the AnnData table, standardizes with `StandardScaler`, runs k-means, and attaches cluster labels to `sdata.tables['table'].obs`.

7. **Network analysis** (`cli.py` → `network_analysis.py`): Extracts centroids from the cell boundaries GeoDataFrame. Builds a `scipy.spatial.Delaunay` triangulation, filters edges by `max_edge_distance`. Constructs a `networkx.Graph` from the filtered edges; when an intensity matrix is supplied, edges are weighted by expression similarity (`1 + Pearson corr`). Runs `nx.community.louvain_communities()` with the `community_resolution` parameter and fixed seed 42. Returns a dict with `cell_ids`, `community_labels`, and `cluster_labels` arrays. Embeds `kmeans_cluster` and `community` labels into both the `stardist_boundaries` GeoDataFrame and the AnnData table obs, and writes `summary.json` (cell/cluster/community/edge counts) to `--plot_dir`. The segmented Zarr is written at this stage (single write after all labels are computed).

8. **Spot-to-cell assignment** (`cli.py` → `utils.py`): Spatial join of spots to cell boundaries using GeoPandas (`gpd.sjoin` with `predicate='within'`). Returns `None` if no spots are present in the dataset (no guard needed in `cli.py` — `export_to_qupath` handles `None`).

9. **QuPath GeoJSON export** (`cli.py` → `export.py`): Exports cell boundaries and spot detections as GeoJSON features with community/cluster assignments and intensity measurements for QuPath visualization. Takes `cell_ids`, `community_labels`, and `cluster_labels` as direct arrays. Output path is `--geojson-path` (default `./qupath_export.geojson`).

10. **Spatial neighborhood analysis** (`cli.py`): Computes spatial neighbors radius graph, mean hop distance heatmap between clusters (`cell_type_to_cell_type.png`), UMAP embedding with k-means coloring (`umap_kmeans_cluster.png`), and Leiden clustering with `igraph` backend (`umap_leiden.png`). Also emits per-cluster/per-community density maps (`cell_density.png`), a mean channel-intensity heatmap per cluster (`cluster_intensity_heatmap.png`), and a morphology-by-cluster plot (`morphology_by_cluster.png`, only if an `area` column is present). All plots are saved to `--plot_dir`.

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
| `--font-size` | `20` | Font size for plots |
| `--axes-linewidth` | `3` | Axes line width for plots |

## Key gotchas

- **`TF_USE_LEGACY_KERAS='1'`**: Set in `run.py` before imports. TF >=2.16 defaults to Keras 3, which compiles Stardist's model with XLA JIT, triggering a cuDNN autotuner failure on 1x1 convolutions with CUDA 12/cuDNN 9. Legacy Keras 2 uses the non-XLA cuDNN path and retains GPU acceleration. Requires the `tf_keras` pip package.

- **System-specific shared library**: `run.py` preloads the pixi environment's `libstdc++.so.6` (resolved relative to the script's `.pixi/` directory) to avoid ABI conflicts with the system library. If the file doesn't exist (e.g., on a non-pixi setup), it silently skips. Note: this preload code is duplicated in `cli.py` (lines 3-9) so that `cli.py` can also be run standalone via `python -m cellsurvey.cli`.

- **GPU detection and `--use-gpu` flag**: GPU is detected at runtime via `tf.config.list_physical_devices('GPU')`. If no GPU is found and `--use-gpu` is not set, a warning is printed but execution continues — Stardist will run on CPU and be very slow. The `--use-gpu` flag forces GPU backend even without auto-detection. Both paths set `sopa.settings.parallelization_backend = None`.

- **Zarr write-then-read pattern**: The pipeline writes the initial Zarr to disk then immediately reads it back before continuing to segmentation. This is intentional — it materializes the spots-including dataset as a clean checkpoint.

- **Segmented Zarr incremental reuse (gated on `--resume-from`)**: The segmented Zarr (`_seg.zarr`) is only checked for reuse when `--resume-from` is set. If it then has a valid table, both Stardist and aggregation are skipped (jump to clustering). If it has boundaries but no table, only aggregation is re-run. If the file exists but can't be read, it's deleted via `shutil.rmtree` and a full re-run proceeds. This enables crash recovery at finer granularity than `--resume-from`; without the flag, any existing `_seg.zarr` is simply overwritten.

- **Pandas 2.x + anndata ArrowStringArray compatibility**: Two workarounds in the aggregation stage:
  1. `pd.option_context('future.infer_string', False)` wraps the `sopa.aggregate()` call to force plain object dtype for strings, avoiding ArrowStringArray which can't be written to Zarr backing stores.
  2. `obs.index.astype(str)` is called on the AnnData table's obs index after clustering to force plain string dtype (same ArrowStringArray issue).

- **Channel naming**: Channel names with suffixes like `_ch_0`, `_ch_1` are constructed in the Stardist stage. The `remove_channel_suffix()` function strips these for clean column headers in the intensity DataFrame.

- **Duplicate column handling**: After stripping channel suffixes, duplicate column names are dropped keeping the first occurrence. This happens because channels with the same name at different indices collapse to the same column name.

- **`.zarr` suffix is auto-appended**: If `-o` doesn't end in `.zarr`, it's appended automatically. The segmented output replaces `.zarr` with `_seg.zarr`. This applies to both `--output_file` and `--resume-from` values.

- **Blob detection must be explicitly enabled**: The `--detect-blobs` flag is required to perform RNA spot/blob detection. Without it, the pipeline writes the raw OME-TIFF as a Zarr and proceeds directly to segmentation/aggregation without any spot data. The `sopa.aggregate()` call checks for `"spots" in dataset.points` to decide whether to aggregate genes.

- **Input format is multi-channel TIFF only**: The image is loaded via `sopa.io.ome_tif()`, which is `dask_image.imread` (tifffile/pims) under the hood — it does **not** go through `bioio`. So it accepts 3D `(C, Y, X)` or 4D `(1, C, Y, X)` TIFFs (`.ome.tiff`/`.tiff`/`.tif`) but **not** ND2/CZI/LIF/DV or single-channel 2D images (those raise `ValueError: Number of dimensions not supported`). Channel names are read from the OME-XML when present; otherwise they fall back to numeric names. **Future aim**: broaden input support via `sopa.io.bioio()` (or `bioio` directly), which handles ND2/CZI/LIF/DV — the only blockers are confirming the channel-name/coordinate handling and adding a CLI switch to select the reader.

- **Louvain community detection uses `networkx`**: `nx.community.louvain_communities()` with fixed seed 42 (for reproducibility). Requires `networkx>=3.4` in `pixi.toml`. The resolution parameter from `--community-resolution` is passed directly.

- **Segmented Zarr write uses tempfile + atomic rename**: The segmented Zarr is written to a temporary directory and then atomically renamed into place. This prevents "path in use" errors that occur when `spatialdata.read_zarr()` has the backing store open (during `--resume-from` with a pre-existing `_seg.zarr`). After the rename, the temporary directory is cleaned up.

- **`--resume-from` enables crash recovery**: If a Zarr already exists at the expected path, you can skip image loading and spot detection and resume from Stardist segmentation. When used with the segmented Zarr reuse logic, this provides two levels of checkpoint restart.

- **Channel subset loading**: The full multichannel image is not loaded into memory for blob detection. Only the channels specified by `--channels` are loaded via `dask_image.imread`, reducing memory footprint. The full image is available on disk via the Zarr for Stardist segmentation.

- **AnnData `.X` can be sparse or dense**: `sopa.aggregate()` may produce either a scipy sparse matrix or a dense numpy array depending on the input data size and sopa version. The intensity extraction at `cli.py:238` handles both with `hasattr(measurements.X, 'toarray')`. Never assume `.X` is sparse.

- **GeoJSON output path**: The QuPath GeoJSON path is configurable via `--geojson-path` (default `./qupath_export.geojson`). It is independent of `-o` and `-p`; all plot/PNG artifacts and `summary.json` go to `--plot_dir`.

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

## Planned: Parameter Stability Sweep

**Goal**: Quantify how robust community assignments are under parameter variation, producing per-cell confidence scores and consensus niches.

### Parameters to explore
| Parameter | Range | Rationale |
|---|---|---|
| `community_resolution` (Louvain) | 0.05–1.0 | Directly controls community granularity |
| `max_edge_distance` | 500–2000 | Changes which cells are neighbors in the Delaunay graph |

> **`n_clusters` is NOT part of the community sweep.** k-means clustering and
> Louvain community detection are independent: Louvain weights its Delaunay edges
> from the raw intensity matrix, never from the k-means labels. So `n_clusters`
> only affects the `kmeans_cluster` labels, not the `community` labels.

(Re-running Stardist or blob detection with varied parameters is not in scope — too expensive. The sweep operates on an existing `_seg.zarr`.)

### Phases

**Phase 1 — Sweep**: For each combination (or random sample, e.g. 50–100 draws), re-run Delaunay → Louvain using the existing segmented Zarr. Stack results into an `(n_cells, n_sweeps)` assignment matrix.

**Phase 2 — Stability metrics**:
- **Co-occurrence matrix**: `(n_cells, n_cells)` — fraction of sweeps where cells A and B share a community
- **Per-cell entropy**: how uniformly is a cell assigned across different community labels (low entropy = stable, high entropy = boundary/transitional)
- **Switching probability**: for each pair, how often they switch community together vs. independently

**Phase 3 — Consensus communities**: Hierarchical clustering on the co-occurrence matrix → final high-confidence niches. Export alongside per-cell confidence scores to GeoJSON.

### Implementation (Phase 1 done)
New module `cellsurvey/stability.py` with `run_stability_sweep(sdata, resolutions, max_edge_distances)` returning a dict with `cell_ids`, `coords`, an `(n_cells, n_sweeps)` label matrix, and the `(resolution, max_edge_distance)` params. It reads an existing `_seg.zarr` **read-only** (never writes back), and `sweep_to_csv`/`sweep_summary` write per-cell labels and per-sweep community counts to CSVs. Standalone entry point: `python -m cellsurvey.stability --zarr <output>_seg.zarr --resolutions ... --max-edge-distances ...`.

Phases 2–3 (stability metrics + consensus communities) are still TODO:
- `run_stability_sweep` should gain co-occurrence, per-cell entropy, consensus labels, and confidence scores (the `(n_cells, n_sweeps)` matrix is the input to these).

### Outputs
- `stability_map.png` — spatial heatmap of per-cell entropy (uncertainty)
- `co_occurrence_heatmap.png` — clustered co-occurrence matrix
- `stability_scores` and `consensus_community` columns in GeoJSON export

### Risks
- Full grid search is `O(n_resolutions × n_distances)` — random sampling is more practical
- Co-occurrence matrix is `O(n_cells²)` memory — sparse storage or chunking needed for large datasets

## Planned: Cell segmentation (nucleus expansion → whole cell)

**Status**: in progress — **v1 + v2 (non-overlapping) implemented and wired**: `--cell-expansion` builds non-overlapping `cell_boundaries` (buffer + custom Voronoi clip) and aggregation targets them. **Pending**: the separate `nucleus` table and `--cluster-regions` downstream threading.

**Goal**: The pipeline currently segments **nuclei only** (Stardist → `stardist_boundaries`). Marker intensity and RNA spots are mostly cytoplasmic, so aggregating over nuclei misses the cytoplasm. Add a **whole-cell** boundary as an expanded version of each nucleus (QuPath-style). A distinct **cytoplasm** compartment is deferred to proper whole-cell segmentation later (Cellpose) — it is *not* approximated arithmetically in this phase.

### Why this matters
- `sopa.aggregate()` measures mean intensity inside each shape. Using nuclei alone undercounts cytoplasmic/membrane markers and conflates nuclear vs cytoplasmic signal.
- Many markers are compartment-specific (nuclear transcription factors vs membrane/cytoplasmic proteins); quantifying them in the *correct* compartment is a core goal.
- `assign_spots_to_cells()` (`gpd.sjoin` with `predicate='within'`) should assign spots to the **whole cell**, and additionally classify each spot as **nuclear vs cytoplasmic**.
- Downstream cell `area` (for `morphology_by_cluster.png`) is biologically meaningful only at whole-cell scale.

### The regions

| Region | Shape layer (proposed) | Definition |
|---|---|---|
| Nucleus | `stardist_boundaries` (existing) | Stardist output, unchanged |
| Whole cell | `cell_boundaries` (new) | Nucleus expanded by `--cell-expansion` µm |

Both share the same `cell_id` index, so their measurements are joinable. A third compartment — **cytoplasm** (whole cell minus nucleus) — is intentionally **deferred** to proper whole-cell segmentation (see Future: Cellpose); it is not derived arithmetically.

### Units: microns, not pixels

- `--cell-expansion` is expressed in **µm** (default `5.0`, matching QuPath's `cellExpansion`).
- Physical pixel size (µm/px) is read from OME metadata via `BioImage(imagepath).physical_pixel_sizes` (bioio). For 2D images this is `(Y, X)`; guard for missing metadata.
- Geometry operations run in pixel coordinates (as Stardist/SOPA produce), so the µm radius is converted to px internally: `expansion_px = expansion_um / pixel_size_um`.
- If physical pixel size is unavailable, warn and fall back to an explicit `--cell-expansion-px` (or refuse with a clear error).
- **Related cleanup (separate, deferred)**: the existing pixel-based params (`--max-edge-distance`, `--radius-min/max`) should eventually move to µm too, but that touches Delaunay/neighbour logic and existing outputs — keep it out of this change.

### Approach (QuPath-style expansion)
QuPath detects nuclei, then dilates each nucleus by a `cellExpansion` distance, using a distance-transform/watershed expansion so neighbouring cells stop where they meet (no overlap).

Two implementation tiers:

1. **v1 — simple buffer (start here)**: `shapely.buffer(nucleus_geometry, expansion_px)`. Fast, trivially parallel over rows, but produces **overlapping** polygons where cells are dense; overlap double-counts pixels in aggregation.
2. **v2 — non-overlapping (QuPath-faithful)**: clip each buffered cell to the Voronoi cell of its centroid so adjacent cells stop where they meet. **Implemented with our own `_resolve_overlap()`** (shapely `voronoi_polygons` + STRtree). **Note**: `sopa.shapes.expand.remove_overlap()` — and `sopa.aggregate(no_overlap=True)`, which calls it — were found to empty most cells on dense data (1742/2107 → NaN geometries), so they are **not** used.

v1 + v2 are implemented and wired through `cli.py` (buffer + custom Voronoi clip).

### Aggregation & measurements (two tables)

`cli.py` currently calls `sopa.aggregate()` once (over `stardist_boundaries`). With whole-cell expansion, aggregate **two solid regions** (no holes):

| Compartment | Source | How |
|---|---|---|
| Whole cell | `cell_boundaries` | `sopa.aggregate(..., shapes_key='cell_boundaries')` → `table` |
| Nucleus | `stardist_boundaries` | `sopa.aggregate(..., shapes_key='stardist_boundaries')` → `nucleus` |

Each table's `.obs` index is the same `cell_id`, so `table` and `nucleus` join cleanly. Decide whether the *primary* `table` should be whole-cell (recommended — most biologically meaningful) or nucleus (current behaviour); open question below.

**Cytoplasm** is *not* computed in this phase — it requires a proper cell/cytoplasm segmentation, deferred to the Cellpose path. (This also sidesteps a known Sopa limitation: `sopa.aggregate()` ignores polygon holes, so a donut-shaped cytoplasm couldn't be aggregated correctly anyway.)

### Downstream consumption

- **Region selection for clustering/network** — add `--cluster-regions` (comma-separated region table keys, default `cell,nucleus`). The selected regions' feature matrices are **concatenated** column-wise into one `(n_cells, n_regions × n_markers)` matrix that drives k-means (`cluster_data`) and Louvain edge weights (`run_network_analysis`). Clustering on *both* nuclear and whole-cell signal is fully supported — no need to pick one. (Note: whole-cell already includes the nuclear signal, so the two are partly redundant; the genuinely complementary pair is nucleus + cytoplasm, which arrives with proper segmentation.)
- **Spot assignment** — `assign_spots_to_cells()` assigns to `cell_boundaries` (whole cell), then classifies each spot as nuclear vs cytoplasmic via a second `within` predicate against `stardist_boundaries` (nucleus `within` → nuclear, else cytoplasmic). Works with the two solid geometries — no donut needed.
- **GeoJSON export** — `export_to_qupath()` emits per-region measurements, prefixed by compartment: `Nucleus: <channel> mean`, `Cell: <channel> mean`. Cell boundaries are the primary `objectType: "cell"` features; nuclei optionally exported as `objectType: "nucleus"` with a `parent_id` link to the cell.

### Implementation plan

1. **Read physical pixel size** — **done**: `get_pixel_size(imagepath)` in `cellsurvey/segmentation.py` (returns µm/px or `None`).
2. **Expand nuclei (v1 + v2)** — **done**: `expand_nuclei(nuclei_gdf, expansion_um, pixel_size_um)` buffers then clips to a custom Voronoi tiling (`_resolve_overlap`), producing non-overlapping cells.
3. **`--cell-expansion` flag** — **done** (µm, default `5.0`). `--cluster-regions` deferred to step 6 (only meaningful once the second table exists).
4. **Build & store `cell_boundaries`** — **done**: inserted after Stardist, before aggregation; `stardist_boundaries` kept.
5. **Aggregate over cells** — **partially done**: `sopa.aggregate` now targets `cell_boundaries` → `table` (whole-cell). The separate `stardist_boundaries` → `nucleus` table is **pending**.
6. **Thread region through downstream** — **pending**: replace hardcoded `'stardist_boundaries'` in `export.py`, `assign_spots_to_cells`, `network_analysis.py`, `stability.py`, and add `--cluster-regions`.

### Future: Cellpose (proper cytoplasm segmentation)
Integrate **Cellpose** for true whole-cell segmentation (cell + cytoplasm), replacing the nucleus-expansion approximation. Cellpose is the dominant pretrained whole-cell segmenter for fluorescent microscopy, and would also yield a genuine **cytoplasm** compartment (cell minus nucleus) — the piece this phase deliberately does not approximate. Note: Cellpose is PyTorch-based (see the TF/Torch co-existence note) — recommend a separate pixi environment or a deferred, opt-in `--segmenter cellpose` path. Recorded here for later; not in scope for the expansion work.

### Open questions
- **Primary `table` region**: whole cell (recommended) vs nucleus (current). Changing `table` to whole-cell alters existing results and makes prior `_seg.zarr` outputs non-comparable — acceptable if flagged as a breaking change.
- **Always-on vs opt-in**: default `--cell-expansion 5.0` (on) vs `0` (off, nuclear-only). A non-zero default is the point of the feature, but changes results vs today.
- **Fixed vs per-nucleus radius**: QuPath uses a fixed `cellExpansion`; a per-nucleus radius (scaled to nucleus area) is more accurate but adds a parameter.

### Risks
- ~~v1 overlaps double-count pixels and inflate aggregated intensity~~ — resolved: v2 clips to a Voronoi tiling, so cells no longer overlap.
- Two tables and a region-selection flag widen the downstream surface (clustering, network, export, sweeps) — easy to miss a hardcoded `'stardist_boundaries'`/`'table'` reference; grep thoroughly.
- Missing OME physical pixel size breaks µm conversion; need an explicit, well-documented fallback.

## Planned: Normalisation (DAPI reference)

**Status**: planned — not implemented.

**Goal**: Correct for technical sources of variation — uneven illumination, stitching artefacts, inhomogeneous antibody labelling, tissue depth/exposure — by normalising each cell's marker measurements to that cell's **DAPI** signal. DAPI is a nuclear counterstain that stains nuclei roughly uniformly, so a cell's DAPI intensity is a proxy for the local technical bias at that location.

### Why this matters
- Raw mean intensities are confounded by illumination/stitching/depth gradients; marker differences can be artefact rather than biology.
- Dividing by DAPI yields a *relative* intensity per cell that is far more comparable across the tissue.

### Approach (simple ratio)
For each cell and each marker channel:

```
x'_cell,marker = x_cell,marker / dapi_cell
```

where `dapi_cell` is that cell's DAPI mean intensity. Applied **after aggregation (stage 5), before clustering/network (stage 6)** — the normalised matrix feeds k-means and Louvain edge weights, while the raw matrix is retained for export/visualisation.

### Implementation plan
1. **CLI flags** — `--normalise-dapi` (off by default) and an optional `--dapi-channel <name>` override. Channel names are already read from the image metadata (OME-XML via `BioImage`/`sopa`), so the DAPI channel is **auto-detected by name** (substring `DAPI`, case-insensitive) by default; the flag only overrides this when the counterstain is named differently (Hoechst, Sytox, DRAQ5, etc.).
2. **Identify the DAPI column** — after `remove_channel_suffix` + duplicate-drop, find the column matching the auto-detected / `--dapi-channel` name; error clearly if absent, and if the match is **ambiguous** (multiple channels contain the token).
3. **Normalise** — divide every marker column by the DAPI column (skip the DAPI column itself). Guard against zero/near-zero DAPI (`max(dapi, eps)`, or skip + warn).
4. **Feed downstream** — the normalised matrix becomes `intensity_df` for `cluster_data()` and `run_network_analysis()`; keep the raw matrix for `export_to_qupath()` and `cluster_intensity_heatmap.png`.

### Interaction with cell segmentation
- DAPI is nuclear, so the reference is the **nuclear** DAPI column regardless of which region(s) `--cluster-regions` selects. With the cell plan, normalisation divides every region's measurements by the same nuclear DAPI reference.
- Composes cleanly with region concatenation: concatenate the selected regions first, then divide the whole matrix by the DAPI column.

### Open questions
- **Ratio vs log-ratio**: a plain ratio can be skewed when DAPI is low; a `log1p` of the ratio is a common follow-up. Keep a simple ratio for now, revisit if needed.
- **What to export**: raw, normalised, or both? Recommend both, e.g. `Cell: <marker> mean` and `Cell: <marker> DAPI-normalised`.
- **Zero/missing DAPI fallback**: error vs skip-cell vs clip — decide before implementation.

### Risks
- Normalisation amplifies noise in low-DAPI cells; the zero/epsilon guard is essential.
- Changing the matrix changes clustering/community results, so existing `_seg.zarr` outputs become non-comparable once enabled (same breaking-change consideration as cell segmentation).

## Reference: PANORAMIC (plevritis-lab)

**Note (for future consideration):** [PANORAMIC](https://github.com/plevritis-lab/panoramic) is an R/Bioconductor package for **multi-sample meta-analysis of spatial colocalization**. It is NOT integrated into CellSurvey yet — this section records the concepts worth borrowing or adopting downstream.

### What it does
- Takes pre-segmented single-cell spatial data (`SpatialExperiment` objects with a `cell_type` label).
- Computes within-sample, cell-type-pair spatial statistics: default `local_comp_enrichment` (edge-corrected, bootstrapped percentage-point enrichment within radius `r`), plus L/K-function alternatives (`Lcross`, `Kcross`, etc.).
- Pools sample-level effects with **multilevel random-effects meta-analysis** (`metafor::rma.mv`) to test **group-level differential colocalization** (case vs. control), producing `beta_diff`, `p_diff`, `fdr_diff`.
- `create_spatial_network()` builds an igraph network of cell-type pairs (edge weight `|z_diff|`, FDR-filtered) with **Leiden** community detection and centrality metrics.

### Relevance to CellSurvey
- **Complementary, not overlapping**: CellSurvey is single-sample and ends at per-cell community/cluster labeling. PANORAMIC adds the **cross-sample statistical hypothesis testing** layer that CellSurvey lacks. It would run *after* CellSurvey.
- **Integration path**: CellSurvey's GeoJSON/AnnData output would need conversion to `SpatialExperiment` (cell coordinates + `cell_type` label in `colData`). Modest adapter only.
- **Methodology worth borrowing** (already conceptually aligned with our Planned Stability Sweep):
  - Bootstrap + uncertainty pooling mirrors the sweep's co-occurrence/entropy/consensus goals.
  - Cell-type-granular Leiden network clustering (PANORAMIC) vs. our per-cell Louvain (`network_analysis.py`).
  - K/L-function edge-corrected enrichment as a principled alternative to our Delaunay `max_edge_distance` filtering.
- **Caveats**: R-only (R ≥ 4.6; `spatstat`, `metafor`, `igraph`) — would require an R sidecar/`rpy2`, or a Python port (scipy/numpy for K/L functions + `statsmodels` for meta-analysis). Not a segmentation tool (assumes cells already segmented). Early-stage (v0.99.3, API may shift).

## Reference: Spatial Permutation & Normalization (plevritis-lab)

**Note (for future consideration):** [Spatial_Permutation_and_Normalization](https://github.com/plevritis-lab/Spatial_Permutation_and_Normalization) is an R script for **significance-testing and normalizing cell-cell colocalization** (colocation quotient, CLQ) on multiplexed immunofluorescence data. Not integrated into CellSurvey — concepts retained for potential reuse.

### What it does
- Computes the **colocation quotient (CLQ)** for each cell-type pair over a fixed k-nearest-neighbor set (k=20, via `spdep::knearneigh`): `CLQ_{b→a} = (C_{b→a}/N_a) / (N_b/(N−1))`.
- **Permutation testing**: spatial coordinates stay fixed; cell-type labels are permuted (500 iterations, preserving proportions) to build a null CLQ distribution per pair. Observed CLQs outside the 5th/95th percentile tails are deemed significant positive/negative colocalizations.
- **Normalization**: tail-clipped Z-score (default right 0.05 / left 0) to make CLQs comparable across samples/conditions, especially for rare cell types whose null distributions are naturally wider.
- Batch-oriented: globs all `*_cell_type_assignment.csv` files and processes each sample.

### Relevance to CellSurvey
- **Different spatial statistic family**: CLQ (k-NN co-occurrence quotient) vs. our Delaunay `max_edge_distance` graph + Louvain. CLQ + permutation gives a **p-value per cell-type pair**, which CellSurvey's deterministic Louvain labeling does not provide.
- **Directly complementary to our Planned Stability Sweep**: Panoramic's bootstrap and this tool's permutation null both answer "is this spatial association significant?" — the same uncertainty question the sweep targets per-cell.
- **Portable to Python**: the core logic is small — k-NN via `scipy.spatial.cKDTree`, CLQ matrix via numpy, permutation null via numpy label shuffling (`numpy.random.default_rng`), and normalization via tail-clipped Z-scoring. No heavy dependencies.
- **Inputs**: requires per-cell cell-type assignments + X/Y coordinates (exactly what CellSurvey outputs via GeoJSON/AnnData obs), though it assumes CELESTA's CSV format upstream.
- **Caveats**: R-only (needs `spdep`, `ggplot2`, `dplyr`), single-script architecture with a duplicated function definition quirk, and rare-population handling (cells with ≤5 of a type get CLQ=0) is heuristic. Not a segmentation tool.

## Reference: CELESTA (plevritis-lab)

**Note (for future consideration):** [CELESTA](https://github.com/plevritis-lab/CELESTA) (CELl typE identification with SpaTiAl information; Zhang & Li et al., Nature Methods 2022) is an R package for **unsupervised, spatial-aware cell-type identification** in multiplexed in situ imaging (CODEX, MIBI/IMC). Not integrated into CellSurvey — concepts retained for potential reuse.

### What it does
- Consumes **already-segmented cells** (X/Y coordinates + per-marker expression columns); does NOT segment.
- Assigns cell types with **no training labels**: fits a per-marker Gaussian Mixture Model (`Rmixmod`) → activation probability, then combines expression-based scoring with **spatial neighborhood context** via EM-style mean-field propagation.
- Works in **hierarchical rounds** (coarse lineage → fine subtype), with "anchor" vs. "index" cell assignment, iterative prior-matrix updates, and a distance-decaying `beta` spatial term.
- Optional `FilterCells()` QC removes doublets/artifacts (all markers uniformly high/low).
- Output: per-cell `*_cell_type_assignment.csv` with per-round and final labels; needs a user-defined marker-signature matrix (1/0/NA per marker per type).

### Relationship across the plevritis-lab toolkit (sequential pipeline)
```
Segmented imaging (XY + markers)
   → CELESTA           : cell-type assignment (*_cell_type_assignment.csv)
   → Spatial_Perm...   : per-sample CLQ colocalization + permutation testing
   → PANORAMIC         : cross-sample/group meta-analysis of colocalization
```

### Relevance to CellSurvey
- **Overlaps CellSurvey's cell-typing intent, different method**: CellSurvey types cells implicitly via k-means on aggregated intensities + Louvain communities. CELESTA type-calls *with a spatial prior* and explicit marker signatures — more interpretable, unsupervised, and lineage-aware.
- **Spatial propagation is thematically aligned** with our Delaunay/Louvain network analysis and the Planned Stability Sweep (both use neighbors to refine assignments; CELESTA's `beta` distance-decay is a cleaner alternative to `max_edge_distance`).
- **Potential role**: a post-Stardist cell-type annotation step between aggregation (stage 5) and clustering/network (stages 6-7), replacing the generic k-means label with marker-informed, spatially propagated types.
- **Portable but heavier than the other two tools**: GMM (`sklearn.mixture.GaussianMixture`), k-NN (`scipy.spatial.cKDTree`), and the EM mean-field loop are all reproducible in Python, but the CELESTA R code is a single large file (`CELESTA_functions.R`, ~25-slot S4 object) with non-trivial logic.
- **Caveats**: requires a user-defined marker-signature/lineage matrix (domain input); R-only (`Rmixmod`, `spdep`, `ggplot2`, `zeallot`); heuristic thresholds (`max_iteration`, `cell_change_threshold`, anchor high/low) need tuning.

## Reference: VALIS (MathOnco)

**Note (for future consideration):** [VALIS](https://github.com/MathOnco/VALIS) (Gatenbee & Anderson, Moffitt; *Nature Communications* 14, 4502, 2023) is a **CPU-only whole-slide image registration** pipeline ("Virtual Alignment of pathoLogy Image Series"). Not integrated into CellSurvey — but the highest-traction *upstream* candidate for multi-round/cyclic imaging. Not integrated yet.

### What it does
- Fully automatic **registration (rigid + non-rigid)** of serial sections and repeated-cycling IF/IHC slides into a common coordinate frame, with no reference image required (auto-selected from the stack center; auto slide ordering by feature similarity).
- Reads 322+ formats via Bio-Formats/OpenSlide; writes warped full-resolution slides as **OME-TIFF pyramids**; merges non-RGB channels into a single **highly-multiplexed OME-TIFF** (`warp_and_merge_slides` + `channel_name_dict` — e.g. 32-channel CyCIF from 11 cycles).
- Key hook: can **warp point-coordinate data** (cell centroids, ROI/polygon vertices) into the registered frame. Also cross-modal registration (H&E ↔ IF/DAPI) for transferring annotations.
- Stack: PyTorch 2.7 + kornia/torchvision (feature detectors/matchers: LightGlue, SuperPoint, BRISK, KAZE, DISK, DeDoDe), pyvips/SimpleITK for warping, Java/Bio-Formats via scyjava/jpype. **CPU-only** (no GPU for registration).

### Relevance to CellSurvey
- **Does not segment** — no overlap with Stardist. Its role is *upstream alignment*, exactly the assumption CellSurvey currently takes for granted (input OME-TIFF is already in one frame).
- **Direct seam**: VALIS's merged output is a **multiplexed OME-TIFF** — precisely CellSurvey's `-i` input. For repeated-cycling data (CyCIF/cyclic IHC), VALIS would pre-align rounds/cycles before CellSurvey's segmentation → aggregation → clustering → network stages.
- **Alternative usage**: run Stardist per-frame first, then use VALIS's **point-warping API** to map cell centroids into the common frame for downstream k-means/Delaunay/Louvain — avoiding a full re-segmentation on the merged image.
- **Caveats**: PyTorch stack (see TF/Torch co-existence note — but registration is CPU-only, so a separate light env is feasible); large-image memory/time during optional micro-registration is a known cost; `error_df` per-pair registration error is useful QC to gate downstream spatial analyses.

## Reference: LazySlide (rendeirolab)

**Note (for future consideration):** [LazySlide](https://github.com/rendeirolab/LazySlide) (Zheng, Abila, Rendeiro et al., CeMM; *Nature Methods* 2026; bioRxiv 2025.05.28.656548) is a **PyTorch whole-slide image (WSI) analysis framework** for histopathology, interoperable with scverse via SpatialData. Not integrated into CellSurvey — concepts retained for potential reuse.

### What it does
- Tile-centric WSI pipeline (`zs.pp` → `zs.seg` → `zs.tl` → `zs.pl`) for H&E/histopathology slides: tissue detection, tiling at a chosen MPP, per-tile model inference, slide-level summaries.
- Container is **`WSIData`**, backed by **SpatialData (Zarr)** — results slot into `wsi.shapes` (tiles/cells/annotations), `wsi.tables` (feature embeddings as AnnData), `wsi.images`, `wsi.attrs`.
- Large model zoo: pathology foundation models (UNI, CONCH, GigaPath, Virchow, H-Optimus, CHIEF, phikon), cell segmentation (**InstanSeg** default, Cellpose, NuLite, HistoPLUS — *no Stardist*), tissue/artifact segmentation (GrandQC), `tl.spatial_domain` (unsupervised domain segmentation from tile embeddings), `tl.virtual_stain`.
- PyTorch/timm + wsidata + SpatialData; CLI entry-point `lazyslide`/`zs`; Python 3.11–3.13; pip + conda-forge.

### Relevance to CellSurvey
- **Shared format, different assay**: both use SOPA-adjacent SpatialData (Zarr) storage, but LazySlide targets **H&E whole-slide histopathology** (tile-first, `.svs`-centric) while CellSurvey targets **multiplexed fluorescent OME-TIFF** (StarDist nuclear segmentation → per-cell marker aggregation → k-means → Delaunay/Louvain).
- **Complementary, not overlapping**: LazySlide has no native k-means/Delaunay/Louvain on cells (graph/community analysis is delegated to Squidpy/Scanpy), no StarDist backend, and OME-TIFF is not an advertised first-class input. Conversely CellSurvey has no foundation-model feature extraction or tile-level pathology models.
- **Potential role**: a *front-end* for extracting pathology-foundation-model or tile-prediction features, written back to a SpatialData Zarr that CellSurvey (or a shared downstream) could open — while CellSurvey keeps ownership of OME-TIFF segmentation/aggregation/clustering/network stages.
- **Caveats**: PyTorch + timm + transformers stack (a second DL framework alongside TensorFlow — see the TF/Torch co-existence note); models are gated and live in the separate `lazyslide-models` package (UNI/Virchow etc. require access); v0.12.0 alpha, API in flux.

## Reference: CellVoyager (zou-lab)

**Note (for future consideration):** [CellVoyager](https://github.com/zou-group/CellVoyager) (Salber, Chen, Sun, Isakova, Wilk, Zou; *Nature Methods* 23, 749–759, 2026; bioRxiv 2025.06.03.657517) is an **LLM-agent for autonomous single-cell RNA-seq analysis** from the Zou Lab (Stanford). It is NOT a segmentation or spatial-omics tool — integrated here because it consumes the same AnnData seam CellSurvey emits. Not integrated into CellSurvey.

### What it does
- Ingests a `.h5ad` AnnData file plus a text summary (paper abstract / dataset summary / biological context / "past analyses tried" / focus directions) and API keys.
- An LLM agent (LiteLLM/OpenAI/Anthropic/Claude Agent SDK) **writes and runs its own scanpy analyses incrementally** inside a live Jupyter notebook, self-critiques the results, and iterates to generate/test new biological hypotheses.
- Output: a live notebook in `outputs/` with plots, statistics, findings, and a research narrative. Two execution modes: `claude` (default; interactive notebook) and `legacy` (programmatic idea executor).
- Stack: LLM orchestration via LiteLLM + scanpy/AnnData computing (env files: numpy, scipy, pandas, matplotlib, seaborn, anndata, scanpy, python-igraph/leidenalg, h5py, litellm, streamlit GUI, celltypist, claude-agent-sdk). Documentation injection is limited to `sc.*`/`scanpy.*` namespaces.

### Relevance to CellSurvey
- **Orthogonal, not overlapping**: CellVoyager operates purely in transcriptomic space (UMAP/tsNE/PCA, neighbor graphs, Leiden/Louvain, DE/markers, pseudotime). It has **zero** references to spatial coordinates, images, segmentation, or SOPA/SpatialData — so it does not replace any CellSurvey stage.
- **Natural integration seam = AnnData**: CellSurvey's `sopa.aggregate()` output is an AnnData table (`sdata.tables['table']`) with `kmeans_cluster` and `community` labels already in `.obs`. Exporting that to `.h5ad` and feeding it (plus a written summary of the completed spatial analysis) into CellVoyager lets the agent autonomously explore hypotheses *on top of* CellSurvey's finished segmentation/clustering/network results.
- **Complementary role**: CellSurvey produces deterministic spatial structure (segmentation → k-means → Delaunay/Louvain); CellVoyager would add open-ended, LLM-driven downstream interrogation on the same cells.
- **Caveats**: MIT but research-grade (no tagged releases); requires paid LLM API keys; expensive/agentic (nondeterministic, needs a sensible "past analyses" prelude to be useful); its environment does not install sopa/spatialdata/stardist, so the hand-off must be through an exported `.h5ad`, not a live SpatialData object.

## Reference: WassersteinWormhole (dpeerlab)

**Note (for future consideration):** [WassersteinWormhole](https://github.com/dpeerlab/WassersteinWormhole) (Haviv & Pe'er lab et al., ICML 2024; arXiv:2404.09411) learns a **Transformer autoencoder embedding of point-clouds** such that Euclidean distance in latent space approximates **optimal-transport / Wasserstein distance** between the original point-clouds. Not integrated into CellSurvey — concepts retained for potential reuse.

### What it does
- Python 3 library (JAX/Flax + OTT-JAX); two classes:
  - `Wormhole` — embeds general weighted point-clouds (per-point features), with an encoder (embedding) and decoder (reconstruct point-clouds for barycenter/interpolation).
  - `SpatialWormhole` — operates on `AnnData` with spatial coords in `.obsm['spatial']`; treats each cell's **k-NN spatial "niche"** as a point-cloud of expression profiles and embeds niches so Euclidean distance ≈ OT distance between their expression distributions.
- Supports OT variants: W1/S1 (Sinkhorn), W2/S2, Gromov-Wasserstein (GW/GS), plus Riemannian (`_R`) variants; automatic Sinkhorn iteration count and distance scaling for numerical stability.
- Enables **O(n) Wasserstein-distance approximation** via embedding, plus learned **Wasserstein barycenters / OT interpolation** through the decoder.

### Relevance to CellSurvey
- **Same input convention as CellSurvey**: `SpatialWormhole` natively consumes AnnData + `.obsm['spatial']`, exactly what `sopa.aggregate()` produces (`sdata.tables['table']`). Low-friction integration point.
- **Niches ≈ CellSurvey's communities**: embedding each cell's spatial k-NN neighborhood is conceptually parallel to our Delaunay/Louvain network analysis — Wormhole gives a **continuous, OT-principled niche distance** instead of discrete Louvain labels. Could complement or validate the deterministic community assignments.
- **Potential uses**: (1) niche/domain annotation as an alternative to k-means + Louvain; (2) a principled distance metric for the Planned Stability Sweep's co-occurrence/consensus clustering; (3) cross-sample comparison (OT distance between tissue niches) that parallels PANORAMIC's cross-sample intent.
- **Trade-offs**: brings a heavy JAX/Flax/OTT-JAX stack (GPU-recommended) on top of TensorFlow already present in CellSurvey — a second DL framework in one environment. Requires model training per dataset (not a drop-in analytical step).
- **Caveats**: research code (early API); needs a tuned `k` for niche size; `SpatialWormhole` save/load re-supplies AnnData at load; not a trajectory-inference tool itself (OT distance supports ordering/interpolation but no pseudotime module).

## Reference: PhenoGraph (dpeerlab)

**Note (for future consideration):** [PhenoGraph](https://github.com/dpeerlab/PhenoGraph) (Levine et al., Cell 2015) is a **graph-based clustering method for high-dimensional single-cell data** — a k-NN similarity graph (Jaccard or Gaussian kernel) followed by **Louvain/Leiden community detection**. Not integrated into CellSurvey — concepts retained for potential reuse.

### What it does
- `phenograph.cluster(data)` takes an `(n_cells × d_markers)` array (or a precomputed sparse kNN graph) and returns `(communities, graph, Q)` where `communities` is a per-cell integer label array (`-1` = outlier) and `Q` is the graph modularity.
- Pipeline: kNN search (k=30) → Jaccard/Gaussian affinity graph → symmetrize → **Louvain** (bundled C++ binaries) or optional **Leiden** (`leidenalg`) modularity optimization → small clusters (`min_cluster_size`=10) relabeled as outliers.
- Also ships `classify()` — semi-supervised label propagation (random-walk/Laplacian) for assigning unlabeled cells.
- Lightweight Python stack: `numpy`, `scipy`, `scikit-learn`, `python-igraph`/`leidenalg`, `psutil`.

### Relevance to CellSurvey
- **Direct overlap with our clustering stage**: CellSurvey's k-means (stage 6) and networkx Louvain (stage 7) are two separate steps; PhenoGraph does a unified **graph-based phenotype clustering** that returns both communities *and* a modularity score `Q` we currently don't compute.
- **Extremely low-friction integration**: pure Python, and it already depends on `python-igraph`/`leidenalg` — same family as CellSurvey's existing `python-igraph` (Leiden backend) and `networkx` (Louvain). No new DL framework.
- **Potential role**: a drop-in alternative to k-means for cell-type assignment (marker-intensity-based, no `n_clusters` to guess — communities emerge from the graph), and a way to quantify clustering quality via modularity.
- **Caveats**: operates on marker/expression space only — **ignores spatial coordinates** (unlike our Delaunay spatial graph). For spatial-aware clustering you'd feed coordinates as features or chain it with our neighbor graph. Uses its own bundled C++ Louvain binaries (vs. our `networkx` Louvain) unless the Leiden backend is chosen.
- **Overlap note re: CELESTA**: PhenoGraph (graph clustering, marker-only) and CELESTA (GMM + spatial propagation) are alternative cell-typing approaches — PhenoGraph is simpler and coordinate-agnostic; CELESTA is spatial-aware and lineage-guided.

## Reference: segger (dpeerlab)

**Note (for future consideration):** [segger](https://github.com/dpeerlab/segger) (Heidari et al., bioRxiv 2025.03.14.643160; Pe'er & Gerstung labs) is a **GNN-based cell segmentation tool for imaging-based spatial transcriptomics (IST)** — Xenium/CosMx/MERSCOPE. Not integrated into CellSurvey — concepts retained for potential reuse.

### What it does
- **Transcript-centric, non-image** segmentation: treats each transcript as a graph node and segmentation as **transcript→cell link prediction** on a heterogeneous graph (`tx` transcript nodes, `bd` cell/boundary nodes, GATv2 attention layers). Assigns transcripts to their cell of origin, then aggregates into cells.
- Needs only **transcript coordinates + nucleus masks** — no pixel-level imaging.
- Trains per-dataset (optionally leveraging scRNA-seq gene-correlation references); metric-learning (L2-normalized embeddings = cosine) with triplet + segmentation losses.
- GPU-native (PyTorch Geometric / PyTorch Lightning + RAPIDS cuDF/cuML/cuGraph/cuSpatial/CuPy); atlas-scale speed via tiling.
- **Exports to SOPA / SpatialData conventions** (`export` subcommand → `anndata.h5ad`, `transcripts.parquet` with `segger_cell_id`, `cell_boundaries.parquet`).

### Relevance to CellSurvey
- **High interoperability**: both use SOPA/SpatialData + pixi; segger's output (cell-by-gene AnnData + boundary polygons) is the natural input to CellSurvey's aggregation/clustering/network stages. Could slot in as an alternative segmentation front-end.
- **Different segmentation paradigm**: CellSurvey uses **Stardist** (image-based, star-convex nuclei on a DAPI channel). Segger is for **probe/target-based IST** where transcripts (not just nuclei) define cells — irrelevant to CellSurvey's microscopy/OME-TIFF DAPI workflow but directly relevant if the project ever ingests Xenium/CosMx data.
- **Key conceptual asset — "transcript-to-cell assignment"**: segger explicitly solves the assignment-accuracy problem that CellSurvey handles heuristically via `gpd.sjoin(predicate='within')` spot-to-cell assignment (stage 8 / `assign_spots_to_cells`). Segger's GNN/link-prediction approach is a more principled alternative when spots lie near cell boundaries.
- **Heavy stack trade-off**: requires PyTorch + PyG + full RAPIDS/CuPy GPU toolchain — a *third* ML framework on top of CellSurvey's TensorFlow, and a second segmenter. High integration cost; only justified if IST data becomes a target.
- **Caveats**: per-dataset training (not a pretrained drop-in); very thin README (algorithm lives in the preprint + external docs site); v0.2.0 research code.

## Reference: cellina (PMBio)

**Note (for future consideration):** [cellina](https://github.com/PMBio/cellina) is a **dual-encoder VAE for spatial transcriptomics** built on scvi-tools. It models how a cell's transcription changes when its local neighborhood is altered — "tissue graph counterfactuals." Not integrated into CellSurvey — concepts retained for potential reuse.

### What it does
- Splits each cell into an **intrinsic latent `z`** (cell identity) and a **spatial-context latent `s`** (neighborhood/microenvironment), then reconstructs counts from `[z; s]` under a Negative Binomial likelihood.
- Two variants: `Cellina` (MLP spatial encoder over degree-normalized neighbor pseudobulk) and `CellinaGCN` (GATv2/GCN message-passing over the spatial connectivity graph).
- **Supervised disentanglement**: cell-type classifier anchors `z`; an adversarial discriminator predicts spatial *domain* from `z` to force microenvironment signal into `s`; optional graph-contrastive loss on `s`.
- **Counterfactual inference** (the key feature): `get_counterfactual_expression` (edge perturbation — rewire a cell's neighbors) and `get_perturbed_expression` (node perturbation — modify neighbor gene expression in silico, e.g. ligand knockout/overexpression), to read out downstream effects on the focal cell.
- Input: `AnnData` counts + spatial connectivity (`obsp`) / neighbor features (`obsm`); `spatial_neighbors()` builds squidpy/mistyR-style kNN graphs. Output: latent arrays + counterfactual count matrices.

### Relevance to CellSurvey
- **Complementary, sits after CellSurvey's core**: CellSurvey produces an aggregated AnnData (cell-by-gene + centroids + community labels). Cellina consumes exactly that shape and answers a **different question** — "what would this cell's expression be under a different neighborhood?" (mechanistic signaling/perturbation screen), which CellSurvey doesn't attempt.
- **No overlap with segmentation/blobs**: cellina is transcriptomics-only (no image/stain deconvolution, no segmentation). It is a *downstream* consumer of the same kind of AnnData CellSurvey emits.
- **Potential use**: turning CellSurvey's community/niche labels and spot-to-cell assignments into perturbation experiments — e.g. knock out a ligand in one community and predict transcriptional response in neighboring cells (biomarker/signaling discovery).
- **Trade-offs**: brings scvi-tools + PyTorch Geometric + torch-scatter/sparse — another DL stack alongside TensorFlow. Requires per-dataset training. Spatial context is graph/coordinate based (not image).
- **Caveats**: research code (v0.7.1/v1.1.0 paths); needs cell-type + domain labels for the disentanglement objectives to work well; CPU version available but GPU expected for scale.

## Reference: GBM_analysis (PMBio)

**Note (for future consideration):** [GBM_analysis](https://github.com/PMBio/GBM_analysis) is the **analysis-code companion to the GBM-Space atlas** (single-cell snRNA+snATAC multi-omics of 12 IDH-wildtype glioblastomas). It is the interpretation layer on top of **scDoRI** (bioFAM/scDoRI), which infers enhancer-mediated gene regulatory networks (eGRNs) as "topics." Mostly *not* aligned with CellSurvey — retained mainly for methodological reference.

### What it does
- `python_scripts/topic_regulation.py` computes **Topic Activation Potential (TAP)** and **Topic Repression Score (TRS)** between scDoRI topics — a "regulation potential" of TF→target-topic links, weighted by epigenetic priming (ATAC accessibility) and significance-tested against precomputed permutation nulls (1000 per topic pair).
- `plasticity_analysis.ipynb` measures **epigenetic plasticity** as the entropy of an ATAC state-classifier's predicted probabilities.
- `tf_screen/` is a **55-TF gain-of-function screening pipeline** (Harmony batch correction, LogisticRegression state/topic classifiers, fold-change/percentile consensus differential testing, dose-response metacells, Wilcoxon DE).
- Stack: `numpy`, `pandas`, `scikit-learn`, `scipy`, `statsmodels`, `scanpy`/`harmonypy`; deterministic with seed 42.

### Relevance to CellSurvey
- **Low direct overlap**: this is RNA/ATAC regulatory-network analysis for a cancer atlas — no imaging, no segmentation, no cell-boundary/spatial-community logic matching CellSurvey's pipeline.
- **Borrowable methodology** (most valuable for our Planned Stability Sweep):
  - **Permutation-null significance with precomputed nulls**: cell-by-cell "is this association real?" — the exact pattern CellSurvey's stability sweep could adopt (precompute shuffled nulls once, then threshold cheaply). Same family as PANORAMIC's bootstrap and Spatial_Permutation's label shuffle.
  - **Entropy-of-classifier-probabilities as a "plasticity/uncertainty" score** — conceptually identical to the sweep's per-cell entropy confidence metric.
  - **Epigenetic-priming-weighted regulation (TAP/TRS)** — a principled way to combine a signal with a per-regulator confidence weight, analogous to weighting CLQ/Delaunay edges.
- **Not worth integrating**: domain-specific (GBM topics from scDoRI), requires scDoRI output as input, and no path through CellSurvey's data flow.
- **Caveats**: pandas <3 required (chained-assignment reliance); `scale_topic_regulation_target_topic` mutates in place (over-normalizes if called twice); large precomputed null files.

## Reference: IMAXT (Cancer Grand Challenge)

**Note (for future consideration):** [IMAXT](https://github.com/IMAXT) ("Imaging and Molecular Annotation of Xenografts and Tumours") is the code org for a **Cancer Research UK Cancer Grand Challenge** (Hannon lab, CRUK Cambridge) that built 3D single-cell molecular tumour maps combining imaging mass cytometry (IMC), MERFISH, and serial two-photon tomography. Not integrated into CellSurvey — one repo is worth noting, the rest are off-topic.

### Relevant repos
- **`imc-nuclear-segmentation`** — full **IMC analysis pipeline**: reads IMC images → watershed segmentation → per-cell channel-intensity catalog (positions, shapes, per-antibody intensities). Same goal as CellSurvey's StarDist → aggregation stage, but via **watershed** instead of StarDist. A useful *reference* for watershed-based segmentation and intensity cataloging, not something to adopt wholesale.
- **`mcdlib` / `imdlib`** — C++ parsers for raw IMC **`.mcd` / `.imd`** file formats (Fluidigm Hyperion output). Only relevant if CellSurvey ever ingests raw Hyperion IMC files directly instead of pre-converted OME-TIFF.
- **`stardist`** (fork) — IMAXT's copy of StarDist; already used by CellSurvey, nothing new.

### Not relevant (off-topic for CellSurvey)
- `MERlin` (MERFISH decoding), `stpt-mosaic-pipeline` (serial two-photon tomography), `Bressan_etal_2021_code` (3D VR tumour models), `owl-pipeline-client/server` (Kubernetes job scheduler), `imaxt-image` (generic image utilities).

### Relevance to CellSurvey
- **Low adoption value, some reference value**: CellSurvey already does StarDist + `sopa.aggregate()`. The only genuinely useful concepts are (1) watershed segmentation as an alternative to StarDist (a lighter-weight option for non-nuclear cell structures), and (2) raw `.mcd`/`.imd` ingestion via `mcdlib`/`imdlib` if direct Hyperion IMC input is ever needed.
- **Caveats**: original IMAXT code is largely astronomy-institute-owned and not actively maintained as a general-purpose library; most repos are forks or publication-specific.

## Reference: novae (prism-oncology)

**Note (for future consideration):** [novae](https://github.com/prism-oncology/novae) is a **graph-based foundation model for spatial domain / niche assignment** on spatial transcriptomics data (Nature Methods 2025). It is the **highest-priority integration candidate** reviewed so far — same lab and ecosystem as Sopa, and it overlaps (rather than merely complements) CellSurvey's clustering/network stages. Not integrated yet.

### What it does
- Self-supervised deep clustering on graphs (SwAV / Sinkhorn-Knopp prototyping): a GAT-style `GraphEncoder` learns per-cell representations **within their local spatial neighborhood**, then assigns cells to hierarchical **spatial domains** (niches), not cell types.
- **Zero-shot**: pretrained models on Hugging Face (`novae-human-0`, `-mouse-0`, `-brain-0`); inference on a new slide needs no training (`compute_representations(adata, zero_shot=True)`); optional short `fine_tune`.
- **Native batch-effect correction** across slides/panels/technologies (`batch_effect_correction`).
- Built-in downstream utilities: spatially variable genes, pathway scores, PAGA domain architecture/trajectory, domain proportions, and **LLM-based niche labeling** (`label_domains`).
- Multimodal extension: fuses H&E histology embeddings (CONCH) with transcriptomics (`compute_histo_embeddings`).

### Relevance to CellSurvey
- **Same stack and input convention**: Novae consumes `AnnData` + `.obsm['spatial']`, i.e. exactly what CellSurvey's `sopa.aggregate()` produces (stage 5). It is part of the same scverse/Sopa/SpatialData ecosystem CellSurvey already builds on.
- **Direct overlap with stages 6–7**: Novae's spatial-domain assignment replaces/upgrades CellSurvey's ad-hoc k-means (`cluster_data`) + Delaunay/Louvain community detection (`run_network_analysis`) with a pretrained, hierarchical, biologically meaningful niche labeling that needs no `n_clusters` guess.
- **Enables cross-sample coherence**: Novae's native batch correction + consistent cross-slide labels directly supports a future multi-sample mode, aligning with the PANORAMIC / cross-sample goals already noted.
- **Synergy with the Planned Stability Sweep**: the sweep's co-occurrence/entropy machinery could quantify Novae's domain stability, giving confidence scores on top of a black-box foundation model.
- **Trade-offs**: heavy PyTorch + PyTorch Geometric + Lightning stack (a second DL framework alongside TensorFlow/Stardist), though it is the most natural addition since Sopa shares the scverse ecosystem. Foundation model is less transparent than deterministic k-means/Louvain. Primary target modality is transcriptomics (Xenium/MERSCOPE/CosMx); antibody/OME-TIFF + blob-detected transcripts is not the canonical use case and should be validated.
- **Caveats**: research/foundation model (v1.1.1); may trail dataset-specific methods (GraphST/STAGATE) on tightly-tuned single-sample benchmarks, but wins on generality, cross-slide transfer, and integrated downstream analysis.

## Reference: SACCELERATOR (SpatialHackathon)

**Note (for future consideration):** [SACCELERATOR](https://github.com/SpatialHackathon/SACCELERATOR) ("SA" = spatially-aware, *not* a GPU/rasterization accelerator) is a **Snakemake benchmarking + consensus framework for spatially aware clustering (SAC) methods** (Sun et al., Nature Methods 2026). Not integrated into CellSurvey — most thematically aligned with our Planned Stability Sweep.

### What it does
- Wraps **~24 SAC methods** (BANKSY, STAGATE, GraphST, SpaceFlow, CellCharter, BayesSpace, etc.) over ~28 datasets, scores with ~17 metrics, then produces a **consensus labeling**.
- Signature **consensus module** (3 steps): aggregate per-method labels → **base-clustering (BC) selection** (automatic via cross-method ARI / smoothness-entropy, or **expert-in-the-loop** manual) → combine via three algorithms: **k-modes** (`dicer`), **LCA** (Latent Class Analysis, `poLCA`), and **weighted** (`igraph` + `future.apply`).
- GPU use is **delegated to the individual method modules** (some PyTorch/TF); the orchestration + consensus layers are CPU R/Python. No segmentation, no rasterization kernels.

### Relevance to CellSurvey
- **Directly parallels the Planned Stability Sweep**: the sweep's *consensus communities* phase is a single-method consensus (vary k-means/Louvain params); SACCELERATOR generalizes this to **cross-method** consensus. Its LCA/k-modes/weighted aggregation and **base-clustering selection** logic are directly borrowable.
- **Metric catalog is a goldmine**: includes `cross-method entropy` and `smoothness-entropy` — the exact "per-cell stability/uncertainty" metric the sweep targets, plus spatial metrics (CHAOS, LISI, PAS) we don't currently compute.
- **Not a library to integrate**: it's a benchmarking harness (Snakemake, R+Python, 24 per-method conda envs). Extract *algorithms and ideas*, not the framework.
- **Caveats**: R/Python mix; MIT-0; consensus quality depends on good base-clustering selection (the expert step), which is hard to automate well.

## Reference: TF/Torch co-existence (integration note)

**Note:** Multiple surveyed tools (Novae, segger, cellina) require PyTorch/PyG while CellSurvey currently uses TensorFlow (for StarDist). Co-installing TF + Torch in one environment is *usually fine* (both dlopen their own CUDA/cuDNN pieces at runtime), but the costs are real and worth avoiding unless a stage is genuinely in-pipeline:

- **Footprint/build time**: TF (~2 GB) + Torch (~2–3 GB) + PyG/RAPIDS-class extras → very large, slow-to-solve pixi env. CellSurvey's `pixi.lock` is already TF-heavy.
- **CUDA/cuDNN coupling**: CellSurvey pins TF ≥2.18 (bundled CUDA 12/cuDNN 9); any Torch dep must resolve a matching cu12 build. Pixi makes this *more* tractable than pip/conda, not less.
- **Environment-variable surface**: CellSurvey already fights TF/Keras ABI issues via `TF_USE_LEGACY_KERAS=1` (process-wide, harmless to Torch) — a second DL stack doubles this class of risk.

**Recommendation**: prefer a **separate pixi environment per DL framework** — run CellSurvey (TF) → write `AnnData`/`SpatialData` (`_seg.zarr`) → run the downstream tool (Torch) in its own env. This mirrors the Sopa→Novae modularity the authors themselves chose, and is especially cheap for zero-shot consumers like Novae. Only co-install if the tool becomes a first-class in-pipeline stage.

## Reference: SmartHisto (Vijendran et al.)

**Note (for future consideration):** [SmartHisto: Bayesian active learning for histology images](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013611) (Vijendran, Arruda, Anderson, Eulenstein; *PLoS Comput Biol* 22(9):e1013611, 2026; code [github.com/flu-crew/histology_segmentation](https://github.com/flu-crew/histology_segmentation), pulmonary dataset doi:10.5281/zenodo.18421739, CC0) is a **Bayesian-active-learning framework for training semantic segmentation models with far less expert annotation**. Not integrated into CellSurvey — relevant only if CellSurvey ever adopts a trainable tissue-region segmenter (it currently uses Stardist, which is pretrained drop-in). Retained chiefly for its **uncertainty-decomposition methodology**, which is conceptually aligned with our Planned Stability Sweep.

### What it does
- Trains a **Bayesian U-Net** (smaller UNet, all weights as independent Gaussians via **Bayes by Backprop**, loss = **DiceBCE + scaled KL**) so per-pixel predictive variance measures uncertainty — an *ensemble-like* effect without multiple models.
- **Active-learning sampling** that selects *informative regions* rather than whole images: pixels are grouped by **SLIC superpixels** (default 1000 segs, compactness 28) and ranked by the **average per-pixel divergence** within each superpixel, so experts label only the highlighted uncertain regions (add 5% per "active epoch").
- Explicitly decomposes predictive variance into **epistemic** (reducible, from lack of data) vs **aleatoric** (irreducible, intrinsic noise) uncertainty, fractionally down-weighting aleatoric to prioritize learning what reduces epistemic uncertainty — a distinction point-estimate models cannot make.
- Validated on **GlaS** (colorectal glands), a custom **pulmonary** (pig lung) dataset, and **TIGER ROI** (breast cancer, 6 tissue classes). Mean IoU 0.75 vs. ~0.60 for baselines; on the hardest (TIGER) benchmark, reached peak mIoU with ~44% of the annotation pool, and no baseline matched SmartHisto's peak at *any* annotation level.

### Relevance to CellSurvey
- **Limited direct applicability to the current pipeline**: CellSurvey's segmentation is Stardist (pretrained, no training data needed); SmartHisto addresses the *training-data-scarce* *custom segmenter* problem, which is a different (though plausible future) use case — e.g. if we ever segment tissue regions or non-nuclear structures that Stardist can't handle.
- **The uncertainty methodology is the durable lesson** and directly parallels the Planned Stability Sweep:
  - **Aleatoric vs. epistemic separation** — the sweep's per-cell entropy/confidence score could adopt this framing: distinguish "this cell's community is inherently ambiguous" (aleatoric) from "we simply haven't sampled enough parameters" (epistemic, reducible by more sweep draws).
  - **Bayesian predictive variance as a confidence signal** — the same idea as SACCELERATOR's entropy metrics and PANORAMIC's bootstrapped uncertainty; reinforces the sweep's core design of *quantifying* uncertainty rather than reporting a single deterministic label.
- **Region-aware sampling** (superpixel divergence) is a mining heuristic for *where to spend annotation budget* — relevant if CellSurvey ever builds a supervised cell-type classifier (cf. CellSighter, RIBCA in the shortlist) and needs to label training data efficiently.

### Caveats
- H&E-only, fixed/consistent magnification, mutually-exclusive classes — the authors explicitly flag these as untested generalizations (varying magnification, alternative stains, and overlapping classes need study).
- Uncorrected proof at time of reading.
- Heavy(r) per-dataset training is required; not a drop-in (contrast Stardist).

## Reference: ST "Ten Quick Tips" (Kurogi et al.)

**Note (for future consideration):** [Ten quick tips for spatial transcriptomics analysis](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1014757) (Kurogi, Shimbara, Koreeda, Tsuyuzaki; *PLoS Comput Biol* 22(9):e1014757, 2026) is a practical, platform-neutral review of the entire spatial transcriptomics (ST) workflow. It is not a tool to integrate — it is methodological guidance. Relevant to CellSurvey because it names several tools and practices that overlap or replace CellSurvey's pipeline stages, and because its statistical cautions (pseudoreplication, spatial autocorrelation) bear directly on our Planned Stability Sweep.

### What it does
Ten tips spanning experimental design → platform selection → data structure → analysis → visualization → interpretation → multi-omics integration → open science → limitations. Two platform families (sequencing-based capture vs. imaging-based ISH) are contrasted (resolution vs. genome coverage; FFPE vs. fresh-frozen), with walkthroughs of data layout (e.g. Visium `filtered_feature_bc_matrix.h5` + `Spatial/` subdir) and recommended standard formats (**AnnData/H5AD**, **Seurat objects**).

### Relevance to CellSurvey (stage-by-stage)
- **Visualization/segmentation (Tips 4–5)**: cites **QuPath** (CellSurvey's export target), **napari**, and **Cellpose** for shifting sequencing-based ST toward single-cell resolution. Notably, CellSurvey already uses QuPath for GeoJSON export and Stardist (rather than Cellpose) for segmentation — same ecosystem, different segmenter.
- **Spatial analysis tooling (Tip 6)**: **Scanpy + Squidpy** (spatial graphs, neighborhood enrichment, ligand–receptor) and **BayesSpace** (Bayesian subspot upsampling). CellSurvey already uses Scanpy; Squidpy's neighborhood/enrichment machinery and BayesSpace's subspot upsampling are directly relevant alternatives to our Delaunay/Louvain approach. Workflow managers (**Nextflow/Snakemake**) recommended for reproducibility — CellSurvey currently has none.
- **Interpretation discipline (Tip 7)**: co-localization / ligand–receptor co-expression is *consistent with* but *does not demonstrate* interaction — inferred cell–cell communication must be treated as **hypothesis generation**. Directly applicable to how CellSurvey's community/niche labels should be reported.
- **Multi-omics integration (Tip 8)**: **cell2location** and **Tangram** for scRNA-seq-atlas mapping; spatial proteomics (**CODEX**, **IMC**) as complementary readouts.
- **3D reconstruction (Tip 10)**: **PASTE** (optimal-transport alignment) and image-registration tools (**Fiji/ImageJ, ANTs, Elastix**); a 12-method benchmark found *no single multi-slice integration method dominates*. Relevant if CellSurvey ever goes from 2D sections to 3D tissue.

### Methodology most worth borrowing (for the Planned Stability Sweep)
- **Pseudoreplication** (Tip 4): spots/cells are spatially autocorrelated and *not* independent replicates; treating thousands of spots as sample size inflates significance. This is the statistical foundation for why our sweep must count replication at the biological (sample) level, not the per-cell level.
- **Robustness-to-parameter checks** (Tip 4): "check results for robustness to method/parameter choice" — exactly the sweep's premise, now with a cited endorsement.
- **LLM cell-type annotation caution** (Tip 7): LLMs show only moderate accuracy for cell-type annotation and need expert oversight — a warning against over-relying on the CellVoyager/LLM seam.
- **Segmentation-error propagation** (Tip 5): segmentation errors (merged/split cells, misassigned transcripts) propagate into per-cell estimates and must be validated — dovetails with the Bruhns et al. shortlist entry.

### Caveats
- Not code — a review/roadmap. No library to install; extract guidance, not a dependency.
- Focused on ST proper (Visium/Xenium/MERFISH etc.); CellSurvey is multiplexed-antibody OME-TIFF + blob-detected transcripts, so platform-specific specifics transfer only loosely, but the analysis/statistical principles transfer directly.
- Balanced platform-neutral framing means it favors generality over CellSurvey-specific prescriptiveness.

## Shortlist: additional candidates (not yet deep-dived)

**Note:** flagged as future candidates from a library scan; full deep-dive analyses deferred.

### Spatial statistics / community robustness (plugs into the Planned Stability Sweep)
- **Bruhns et al., "Effects of segmentation errors on downstream analysis in highly-multiplexed tissue imaging"** (*PLoS Comput Biol*, 2025) — perturbs segmentation via affine transforms and measures degradation in k-means/Leiden clustering and GMM phenotyping. Closest empirical validation of our stability-sweep concern: downstream robustness to *upstream* error.
- **SpatialMNN** (Zhou, Hicks; *Bioinformatics*, 2025) — mutual-nearest-neighbor graph + Louvain for cross-sample spatial-domain integration/batch correction.
- **SPF** (Vu, Ghosh; *PLoS Comput Biol*, 2022) — K-function variants + functional Cox regression linking cell-interaction patterns to survival. Principled alternative to `max_edge_distance`.
- **cytoNet** (Mahadevan, Qutub; *PLoS Comput Biol*, 2022) — network-science features of cell communities + cell-cell interaction effects.
- **spicyR** (Canete, Patrick; *Bioinformatics*, 2022) — cross-group colocalization-change inference (R analogue of PANORAMIC's statistical question).

### Cell-type phenotyping (alternatives to k-means)
- **RIBCA — Robust Image-Based Cell Annotator** (Sun, Murphy; *Cell Systems*, 2025) — training-free, reference-based cell-type annotation for multiplexed images (>3M cells, >40 tissues).
- **CellSighter** (Amitay, Keren; *Nat Commun*, 2023) — deep-learning cell classification on multiplexed images with per-cell **prediction confidence**.

### Spot detection (replaces LoG blob detection)
- **Spotiflow** (Mantes, Weigert; *Nat Methods*, 2025) — subpixel-accurate, deep-learning spot detection for spatial transcriptomics; generalizes across chemistries; drop-in upgrade candidate for `blob_log`.

### Statistical rigor for confidence-score reporting
- **Morgan, "Alternative to the statistical mass confusion of testing for 'no effect'"** (*J Cell Biol*, 2025) — replace p-values with effect sizes/confidence intervals.
- **Kitanovski et al., "Uncertainty-aware quantitative analysis"** (*PLoS Comput Biol*, 2026) — Bayesian hierarchical uncertainty quantification, same philosophy as the sweep (quantify uncertainty, avoid NHST pitfalls).

## File structure

```
.
├── run.py                            # Entry-point shim: TF_USE_LEGACY_KERAS + libstdc++ guard + delegates to cli.main()
├── cellsurvey/
│   ├── __init__.py                   # Re-exports all public symbols
│   ├── cli.py                        # main() with argparse and pipeline orchestration
│   ├── blob_detection.py             # detect_blobs_in_tile, detect_blobs_tiled
│   ├── network_analysis.py           # run_network_analysis + compute_louvain_communities (scipy Delaunay + networkx Louvain, expr-similarity weights, summary.json)
│   ├── stability.py                  # run_stability_sweep (community) + run_cluster_sweep (k-means), read-only → CSV
│   ├── export.py                     # export_to_qupath
│   └── utils.py                      # remove_channel_suffix, cluster_data, assign_spots_to_cells, get_colors_for_communities
├── docs/                             # MkDocs Material site (ReadTheDocs-hosted)
│   ├── index.md                      # Landing / overview + pipeline diagram
│   ├── installation.md               # pixi + Docker setup
│   ├── getting-started.md            # First-run walkthrough
│   ├── pipeline.md                   # 10-stage walkthrough
│   ├── parameters.md                 # Full CLI reference + tuning guidance
│   ├── outputs.md                    # Every output file explained
│   ├── visualization.md              # Odon / QuPath / TissUUmaps / napari
│   └── faq.md                        # Question-phrased how-to + troubleshooting
├── Dockerfile                        # Ubuntu 24.04 + pixi + GPU-ready container
├── pixi.toml                         # Pixi environment config (linux-64 only)
├── pixi.lock                         # Pixi lockfile (generated)
├── mkdocs.yml                        # MkDocs + Material config
├── .readthedocs.yaml                 # ReadTheDocs build config
├── requirements.txt                  # Minimal pip requirements (sopa)
├── requirements-docs.txt             # Docs build deps (mkdocs-material)
├── README.md                         # User-facing installation and usage docs
├── LICENSE                           # License
├── .gitignore / .gitattributes       # Git config
└── .pixi/                            # Pixi environment directory (gitignored)
```
