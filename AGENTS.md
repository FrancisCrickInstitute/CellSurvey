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

   **Segmented Zarr reuse**: Only when `--resume-from` is specified, checks if `_seg.zarr` already exists before running Stardist:
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

- **Segmented Zarr write uses tempfile + atomic rename**: The segmented Zarr is written to a temporary directory and then atomically renamed into place. This prevents "path in use" errors that occur when `spatialdata.read_zarr()` has the backing store open (during `--resume-from` with a pre-existing `_seg.zarr`). After the rename, the temporary directory is cleaned up.

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

## Planned: Parameter Stability Sweep

**Goal**: Quantify how robust community assignments are under parameter variation, producing per-cell confidence scores and consensus niches.

### Parameters to explore
| Parameter | Range | Rationale |
|---|---|---|
| `n_clusters` (k-means) | 5–20 | Different cluster resolutions change the expression-feature space fed into Louvain |
| `community_resolution` (Louvain) | 0.05–1.0 | Directly controls community granularity |
| `max_edge_distance` | 500–2000 | Changes which cells are neighbors in the Delaunay graph |

(Re-running Stardist or blob detection with varied parameters is not in scope — too expensive. The sweep operates on an existing `_seg.zarr`.)

### Phases

**Phase 1 — Sweep**: For each combination (or random sample, e.g. 50–100 draws), re-run k-means → Delaunay → Louvain using the existing segmented Zarr. Stack results into an `(n_cells, n_sweeps)` assignment matrix.

**Phase 2 — Stability metrics**:
- **Co-occurrence matrix**: `(n_cells, n_cells)` — fraction of sweeps where cells A and B share a community
- **Per-cell entropy**: how uniformly is a cell assigned across different community labels (low entropy = stable, high entropy = boundary/transitional)
- **Switching probability**: for each pair, how often they switch community together vs. independently

**Phase 3 — Consensus communities**: Hierarchical clustering on the co-occurrence matrix → final high-confidence niches. Export alongside per-cell confidence scores to GeoJSON.

### Implementation sketch
New module `cellsurvey/stability.py` with `run_stability_sweep(sdata, intensity_df, param_grid)` returning a dict of assignment matrix, co-occurrence, entropy, consensus labels, and confidence scores. CLI flag: `--stability-sweep` with optional `--sweep-iterations` (default 50).

### Outputs
- `stability_map.png` — spatial heatmap of per-cell entropy (uncertainty)
- `co_occurrence_heatmap.png` — clustered co-occurrence matrix
- `stability_scores` and `consensus_community` columns in GeoJSON export

### Risks
- Full grid search is `O(n_clusters × n_resolutions × n_distances)` — random sampling is more practical
- Co-occurrence matrix is `O(n_cells²)` memory — sparse storage or chunking needed for large datasets

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
