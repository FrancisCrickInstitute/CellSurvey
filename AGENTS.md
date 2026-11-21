# AGENTS.md

## Project overview

CellSurvey is a spatial biology/omics analysis pipeline built on [Sopa](https://gustaveroussy.github.io/sopa/) (segmentation, aggregation). It processes multichannel microscopy images (OME-TIFF) into spatial data objects, segments nuclei with Stardist, detects RNA spots via blob detection, clusters cells with k-means, builds Delaunay networks, detects Louvain communities, and exports GeoJSON for QuPath visualization.

**MuSpAn removed (complete):** MuSpAn has been fully replaced with open-source libraries. Delaunay triangulation uses `scipy.spatial.Delaunay`, Louvain community detection uses `networkx`, and visualisation uses matplotlib. No private dependencies remain.

## Environment and package management

This project uses **pixi** (via `pixi.toml`) for environment management targeting `linux-64` (GPU) and `win-64` (CPU-only, for local testing). `pixi.lock` is committed (from `feature/omero` on). The omero feature's conda pins are fragile, so the lock keeps every install on the tested versions. Use `pixi install --locked` on the HPC to install exactly what was tested, and commit the updated lock whenever `pixi.toml` changes.

**One pixi environment per segmenter**, so TensorFlow and PyTorch are never co-installed:
- `default` (alias `spatial`) — feature `stardist`: Stardist + TensorFlow + `tf_keras`.
- `instanseg` — feature `instanseg`: `instanseg-torch` (PyTorch; PyPI torch wheels bundle CUDA on linux-64, CPU-only on win-64). No TensorFlow.
- Shared deps (sopa, bioio, scipy, networkx, igraph, …) are in the top-level `[pypi-dependencies]`.
- **Long-term aim**: migrate fully to InstanSeg and drop Stardist/TensorFlow. Keep TF imports lazy and Stardist-specific logic isolated in `segmentation.run_stardist` so it is easy to delete.

Key dependency constraints:
- **Python**: `>=3.12, <3.13`
- **Sopa**: `>=2.3.0` — the version InstanSeg support was developed against. The InstanSeg path only uses `sopa.segmentation.custom_staining_based` (not `sopa.segmentation.instanseg`), which earlier 2.x releases also have, so the pin could be relaxed if needed
- **TensorFlow** (`stardist` feature only): `>=2.18` on `linux-64` (with `and-cuda` extras for GPU)
- **CUDA/cuDNN**: TF >=2.18 bundles its own CUDA 12/cuDNN 9 libraries; no separate conda CUDA/cuDNN packages are needed
- **`tf_keras`** is a required pypi dependency — TF >=2.16 defaults to Keras 3, but Stardist needs legacy Keras 2 API to avoid cuDNN autotuner failures on CUDA 12/cuDNN 9
- **`scipy` (`>=1.14, <2`)** and **`networkx` (`>=3.4, <4`)** for Delaunay triangulation and Louvain community detection (replacing MuSpAn)
- **`python-igraph`** for fast Leiden clustering in Scanpy's spatial neighborhood analysis
- **`bioio` (`>=3.4`) and `bioio-ome-tiff`** for reading channel names from OME-TIFF metadata (via `BioImage`); `setuptools` is pinned as a pypi dependency
- **Windows is CPU-only for local testing** via pixi (`win-64`, TensorFlow CPU wheel — no GPU, no `and-cuda`). **macOS is not supported** — only `linux-64` and `win-64` are in the platforms list.
- **Sopa is a single-maintainer dependency risk**: Sopa is maintained almost entirely by one person (Quentin Blampey), who has left academia. CellSurvey currently uses Sopa for image reading, patching, the Stardist wrapper, aggregation, and spatial-neighbour analysis. **Strategic aim**: gradually reduce Sopa coupling — own the segmentation/aggregation/spatial steps directly where practical — while continuing to use it in the short term. The cell-segmentation work (and eventual Cellpose integration) is a natural seam to start decoupling.

A **Dockerfile** is provided: Ubuntu 24.04 base, installs pixi, installs the pixi env named by build arg `PIXI_ENV` (default `default`; `instanseg` for InstanSeg), copies `run.py` + `cellsurvey/`, sets `TF_USE_LEGACY_KERAS=1`, entrypoint is `pixi run -e $PIXI_ENV python run.py`.

There is no Makefile, no CI/CD, and no tests. The source code is split across 7 modules (plus `__init__.py`) under the `cellsurvey/` package, with `run.py` as the entry-point shim.

**TODO**: Set up linting and formatting (Ruff, mypy) with a `pyproject.toml` config and pre-commit hooks.

**TODO (docs)**: Add a *tutorial* page and more visual illustrations to the MkDocs site. All eight nav pages are now written (`index`, `installation`, `getting-started`, `pipeline`, `parameters`, `outputs`, `visualization`, `faq`), but a step-by-step tutorial (with example inputs/outputs) and richer diagrams/example-plot images would help first-time users. Revisit `docs/` once example datasets and output images are available.

## Commands

**Development environment setup (pixi):**
```bash
pixi install                  # default env (Stardist / TensorFlow)
pixi install -e instanseg     # InstanSeg env (PyTorch)
```

**Run the pipeline (pixi):**
```bash
pixi run python run.py -i <input_tiff> -o <output_zarr_prefix> -p <plot_output_dir>
pixi run -e instanseg python run.py --segmenter instanseg -i <input_tiff> -o <output_zarr_prefix> -p <plot_output_dir>
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

docker build --build-arg PIXI_ENV=instanseg -t cellsurvey-instanseg .
docker run --gpus all -v /path/to/data:/data cellsurvey-instanseg --segmenter instanseg -i /data/input.tiff -o /data/output -p /data/plots
```

There is no build step, no test command, and no linting configured.

## Architecture and data flow

The pipeline is split into modules under the `cellsurvey/` package. `run.py` is a shim that sets environment variables (TF_USE_LEGACY_KERAS), preloads libstdc++, and delegates to `cellsurvey.cli.main()`. The pipeline runs these stages sequentially:

1. **Image loading + initial Zarr** (`cli.py` → `image_source.py`): New runs call `open_image_source(path)` once at startup. This gives a lazy `(c, y, x)` image, channel names (`sopa.io.ome_tif(..., as_image=True)`) and the pixel size (`get_pixel_size`, bioio). `ImageSource.to_spatialdata()` + `write()` then write `<prefix>.zarr` (image only, with `attrs['cellsurvey'] = {pixel_size_um, source}`). **The source is read only here**; every later stage reads the Zarr. `--resume-from` never opens the source: the pixel size comes from `--pixel-size`, else `attrs['cellsurvey']` in the resume Zarr (`image_source.stored_pixel_size`), else (older Zarrs only) `get_pixel_size(imagepath)`. Blob `--channels`/`--thresholds` are parsed and checked (lengths, then index range against the source) before anything is written.

2. **Spot/blob detection** (`cli.py` → `blob_detection.detect_spots`): Only runs when `--detect-blobs` is passed, **after** the initial Zarr write, on the image read back from `<prefix>.zarr` (`sopa.utils.get_spatial_image`); spots are added with `write_element('spots')`. `detect_spots` passes each selected channel as a `(y, x)` array to `detect_blobs_tiled` (tiled Laplacian-of-Gaussian, `skimage.feature.blob_log`, parallelised with `dask.delayed` on threads; overlap-region filtering prevents duplicates), so a different detector (e.g. Spotiflow) can be swapped in there. Gene names are the Zarr's channel names. The `gene` column is built as strings (`np.full(n, name, dtype=object)`); with 0 blobs, `np.array([])` would be float and `sopa.aggregate`'s spot counting crashes (`.str accessor … not floating`; a pre-existing bug, found on the HPC 2026-11-20). If **no** spots are found in any channel, `detect_spots` returns `None` and no `spots` element is added. The run then continues as one without spots (with a warning), or stops straight away if `--cluster-on spots|both` needs them. It runs only when segmentation will re-read `<prefix>.zarr`. On `--resume-from`, existing spots are reused. Missing spots are detected then, unless a `_seg.zarr` is being reused, in which case it warns. Before 2026-11 detection ran on the TIFF before a single write (a March prototype choice for Spotiflow); on a synthetic OME-TIFF both orders give identical spots.

3. **Initial Zarr write** (`cli.py`): The dataset (with or without spots) is written to disk with try/except error handling.

4. **Segmentation** (`cli.py` → `segmentation.py`): Reads back the Zarr (materialized checkpoint), creates image patches via `sopa.make_image_patches()`, renames channel coordinates to include `_ch_` suffixes (e.g., `DAPI_ch_0`), then runs the segmenter chosen by `--segmenter` on the `--seg-channels`:
   - **`stardist`** (default) → `run_stardist`: lazily imports TF, detects GPU via `tf.config.list_physical_devices('GPU')`, runs `sopa.segmentation.stardist()` (`2D_versatile_fluo`) on one channel (default index 0) → `stardist_boundaries` (nuclei).
   - **`instanseg`** → `run_instanseg`: needs the `instanseg` pixi env. Runs InstanSeg (`--instanseg-model`, default `fluorescence_nuclei_and_cells`) through `sopa.segmentation.custom_staining_based` with our own per-patch callable `_instanseg_paired_patch`. That callable runs `eval_small_image(target='all_outputs')`, so InstanSeg pairs nuclei with cells, and returns one channel. It runs **twice** because Sopa stitches one mask per run: once returning nuclei → `nucleus_boundaries`, once returning cells → `cell_boundaries`. Default channels: all (InstanSeg is channel-invariant). Passes `pixel_size` so the image is rescaled to the model's training resolution. Then `link_nuclei_to_cells` adds `nucleus_id` to each cell (largest nucleus whose representative point it contains; NaN if none), and cells with no nucleus are dropped unless `--keep-unnucleated` is set.

   **Segmented Zarr reuse**: Only when `--resume-from` is specified, checks if `_seg.zarr` already exists before segmenting:
   - Has `tables['table']` → skips segmentation and aggregation, jumps to clustering
   - Has the chosen segmenter's layers (stardist: `stardist_boundaries`; instanseg: `nucleus_boundaries` **and** `cell_boundaries`) but no table → re-runs aggregation only
   - Missing or corrupt → full segmentation + aggregation

5. **Cell expansion + channel aggregation** (`cli.py`): For Stardist only, if `--cell-expansion` > 0, `expand_nuclei` grows nuclei into non-overlapping `cell_boundaries` (marker-controlled watershed; records `nucleus_id`). Aggregation targets the cell layer from `resolve_shape_keys` (`cell_boundaries` if present, else nuclei). Runs `sopa.aggregate(shapes_key=...)` to compute per-cell mean intensities for each channel (genes). If `"spots"` exists in the dataset points, passes `aggregate_genes=True`, `points_key='spots'`, and `gene_column='gene'` to assign spots to cells. Otherwise runs plain aggregation. Wraps the aggregation in `pd.option_context('future.infer_string', False)` to prevent ArrowStringArray errors on Zarr write. Writes segmented Zarr with try/except error handling. The segmented Zarr replaces `.zarr` with `_seg.zarr`.

6. **Features + k-means clustering** (`cli.py` → `utils.py`):
   - **Feature tables**: `utils.extract_feature_tables(table)` returns `(intensity_df, counts_df)`. Without spots, `.X` holds the mean intensities. With spots (`sopa.aggregate(aggregate_genes=True)`), `.X` holds the **spot counts** and the intensities are in `obsm['intensities']`. Never read `.X` as intensities directly.
   - **Clustering features**: `utils.build_clustering_features(intensity_df, counts_df, --cluster-on, cell_area_um2)` builds the scaled per-cell features (`intensities`: z-scored; `spots`: spots per 100 µm² of cell area → `log1p` → z-scored, columns suffixed `' (spots)'`; `both`: the two groups, each divided by √(number of its columns)). k-means runs on these with `cluster_data(..., scale=False)`, and cluster labels are attached to `obs`.
   - **Same features everywhere**: the network edge weights (stage 7), UMAP/Leiden (stage 10, via `obsm['X_cellsurvey']` + `use_rep`) and the stability sweeps use them too. `cluster_on` is recorded in `summary.json` and `table.uns['cellsurvey']`.

7. **Network analysis** (`cli.py` → `network_analysis.py`): Extracts centroids from the cell boundaries GeoDataFrame. Builds a `scipy.spatial.Delaunay` triangulation, filters edges by `max_edge_distance`. Constructs a `networkx.Graph` from the filtered edges; edges are weighted by the similarity of the two cells' clustering features (`1 + Pearson corr` across features). Runs `nx.community.louvain_communities()` with the `community_resolution` parameter and fixed seed 42. `filter_small_communities` then relabels communities smaller than `--min-community-size` (default 10) as `-1` ("unassigned") and renumbers the rest by size (0 = largest); the stability sweep does the same, using the value stored in `table.uns['cellsurvey']['min_community_size']` by default (`utils.resolve_min_community_size`). Returns a dict with `cell_ids`, `community_labels`, and `cluster_labels` arrays. Embeds `kmeans_cluster` and `community` labels into the cell-layer GeoDataFrame and the AnnData table obs, and maps them onto the nuclei layer via `nucleus_id` (unlinked nuclei get -1), and writes `summary.json` (cell/cluster/community/edge counts; `n_communities` counts communities at or above the minimum, plus `n_communities_raw`, `min_community_size` and `n_unassigned_cells`) to `--plot_dir`. The segmented Zarr is written at this stage (single write after all labels are computed).

8. **Spot-to-cell assignment** (`cli.py` → `utils.py`): Spatial join of spots to cell boundaries using GeoPandas (`gpd.sjoin` with `predicate='within'`). Returns `None` if no spots are present in the dataset (no guard needed in `cli.py` — `export_to_qupath` handles `None`).

9. **QuPath GeoJSON export** (`cli.py` → `export.py`): Exports cell boundaries and spot detections as GeoJSON features with community/cluster assignments and measurements (`Cell: <channel> mean` intensities and, with spots, `Cell: <channel> spot count`) for QuPath visualization. Takes `cell_ids`, `community_labels`, and `cluster_labels` as direct arrays. Output path is `--geojson-path` (default `<output_file>.geojson`, i.e. beside the Zarrs; derived from `-o` even when resuming).

10. **Spatial neighborhood analysis** (`cli.py`): Computes spatial neighbors radius graph, mean hop distance heatmap between clusters (`cell_type_to_cell_type.png`), UMAP embedding with k-means coloring (`umap_kmeans_cluster.png`), and Leiden clustering with `igraph` backend (`umap_leiden.png`); neighbours/UMAP/Leiden use the clustering features (`use_rep='X_cellsurvey'`), not Scanpy's count normalisation of `.X`. Also emits density maps (`cell_density.png`, `cli.plot_cell_density`: a 4-column grid of one panel per k-means cluster, then the `MAX_DENSITY_COMMUNITIES` = 20 largest communities by size; the suptitle states how many communities and unassigned cells were left out; every panel shares the whole-tissue extent and is drawn with y increasing downwards, like the image), a heatmap of the clustering features per cluster (`cluster_intensity_heatmap.png`, z-scored without the `both` group weighting), and a morphology-by-cluster plot (`morphology_by_cluster.png`, cell area converted from Sopa's px² to µm²; only if an `area` column is present). All plots are saved to `--plot_dir`.

### CLI arguments

All analysis parameters are exposed as command-line flags with sensible defaults:

| Flag | Default | Description |
|---|---|---|
| `-i`, `--input_file` | *(required)* | Input OME-TIFF path, or `omero://host[:port]/Image:<id>` |
| `--omero-readers` | `8` | OMERO input only: parallel tile readers (connections) for the stage-1 copy |
| `-o`, `--output_file` | *(required)* | Path to output Zarr (`.zarr` appended if missing) |
| `-p`, `--plot_dir` | `.` | Output directory for plots |
| `--detect-blobs` | — | Enable RNA spot blob detection on specified channels (default: off) |
| `--segmenter` | `stardist` | `stardist` (default env) or `instanseg` (needs `pixi run -e instanseg`) |
| `--seg-channels` | *(segmenter default)* | Comma-separated channel indices to segment on (stardist: `0`; instanseg: all) |
| `--instanseg-model` | `fluorescence_nuclei_and_cells` | InstanSeg pretrained model name or path |
| `--keep-unnucleated` | — | InstanSeg only: keep cells containing no nucleus (dropped by default) |
| `--pixel-size` | *(from OME metadata)* | µm/px override; **required** if the image has none (run stops at startup) — all distance flags are in µm |
| `--use-gpu` | — | Request GPU segmentation (auto-detected by default) |
| `--cell-expansion` | `5.0` | Stardist only: expand nuclei into cells by this radius (µm); `0` = nuclei only |
| `--channels` | `9,10,11,12` | Comma-separated channel indices for blob detection |
| `--thresholds` | `0.01,0.1,0.1,0.1` | Comma-separated blob detection thresholds (one per channel) |
| `--tile-size` | `2048` | Tile size for blob detection |
| `--overlap` | `50` | Tile overlap for blob detection |
| `--workers` | `14` | Worker threads for blob detection |
| `--min-sigma` | `2` | Minimum blob radius for spot detection |
| `--max-sigma` | `5` | Maximum blob radius for spot detection |
| `--num-sigma` | `5` | Number of sigma steps for blob detection |
| `--n-clusters` | `10` | Number of k-means clusters |
| `--cluster-on` | `intensities` | Features for k-means, network edge weights and UMAP/Leiden: `intensities`, `spots` (spot density; needs spots) or `both` (groups weighted equally). Also in the stability sweep, where it defaults to the mode recorded in the Zarr (`utils.resolve_cluster_on`; `intensities` for older Zarrs) |
| `--community-resolution` | `0.1` | Louvain community detection resolution |
| `--min-community-size` | `10` | Communities with fewer cells are labelled `-1` (unassigned); the rest are numbered by size, 0 = largest. `1` keeps all. Also in the stability sweep, where it defaults to the value stored in the Zarr (10 for older Zarrs) |
| `--max-edge-distance` | `50` | Max edge length (µm) for the Delaunay network |
| `--radius-min` | `0` | Min radius (µm) for spatial neighbors graph |
| `--radius-max` | `280` | Max radius (µm) for spatial neighbors graph |
| `--resume-from` | — | Path to existing Zarr to resume from (never opens the input image; reuses spots, detects them only if missing and `--detect-blobs` is set) |
| `--geojson-path` | `<output_file>.geojson` | Output path for QuPath GeoJSON |
| `--font-size` | `20` | Font size for plots |
| `--axes-linewidth` | `3` | Axes line width for plots |

## Key gotchas

- **`TF_USE_LEGACY_KERAS='1'`**: Set in `run.py` before imports. TF >=2.16 defaults to Keras 3, which compiles Stardist's model with XLA JIT, triggering a cuDNN autotuner failure on 1x1 convolutions with CUDA 12/cuDNN 9. Legacy Keras 2 uses the non-XLA cuDNN path and retains GPU acceleration. Requires the `tf_keras` pip package.

- **TF CUDA libraries are preloaded** (`segmentation._preload_tf_cuda_libs`, called at the top of `run_stardist`): the TF 2.21 wheel's RUNPATH lists only some `site-packages/nvidia/*/lib` folders (no `cusolver`, `nvjitlink`, `curand` or `cuda_nvrtc`). TF then fails to find `libcusolver.so.11`, logs "Cannot dlopen some GPU libraries" and silently runs on CPU. Also, `LD_LIBRARY_PATH` is searched before RUNPATH, so a system CUDA (the Crick GPU nodes put CUDA 13.0's `/usr/local/cuda/lib64` on it) supplies a mismatched `libcusparse.so.12`. Loading the bundled libraries by full path in dependency order fixes both, because later lookups by name reuse the loaded copies. Found on an A100 node (driver 580, 2026-11-18). The InstanSeg env is unaffected: Torch finds its CUDA 13 wheels itself. To diagnose: `TF_CPP_MAX_VLOG_LEVEL=1 … 2>&1 | grep dso_loader`. Note that a bare `python -c "import tensorflow"` skips both this preload and the `libstdc++` one, so it can wrongly suggest the GPU or numpy is broken.

- **HPC: a long-lived GPU allocation can lose CUDA**: on a Slurm interactive job running for more than a day (Crick, 2026-11-19), TF began failing with `failed call to cuInit: CUDA_ERROR_NOT_INITIALIZED`, after hanging about 200 s, and silently ran Stardist on CPU. `nvidia-smi` still showed a healthy, idle GPU. A bare TF check with only the CUDA preload failed the same way, so it was not the pipeline. A fresh allocation fixed it. If "No GPU detected" appears in a job that previously used the GPU, start a new job before debugging the code.

- **System-specific shared library**: `run.py` preloads the active environment's `libstdc++.so.6` (`sys.prefix/lib`, so it works for any pixi env) to avoid ABI conflicts with the system library. If the file doesn't exist (e.g., Windows or a venv), it silently skips. Note: this preload code is duplicated at the top of `cli.py` so that `cli.py` can also be run standalone via `python -m cellsurvey.cli` — change both together.

- **TensorFlow is imported lazily**: only inside `segmentation.run_stardist`. Importing `cellsurvey` (or running the sweeps) must not import TF, otherwise the `instanseg` env breaks. Likewise `torch`/InstanSeg are imported only in `run_instanseg`.

- **GPU detection and `--use-gpu` flag**: Stardist detects GPUs via `tf.config.list_physical_devices('GPU')`; InstanSeg via `torch.cuda.is_available()` (and passes `device='cuda'`/`'cpu'`). If no GPU is found a warning is printed but execution continues on CPU (slow). `--use-gpu` only changes the messages for Stardist; for InstanSeg it warns if CUDA is unavailable. Both set `sopa.settings.parallelization_backend = None`.

- **Shape layers — never hard-code one key**: nuclei are `stardist_boundaries` (Stardist) or `nucleus_boundaries` (InstanSeg); whole cells are `cell_boundaries` (Stardist + watershed expansion, or InstanSeg cells). `utils.resolve_shape_keys(sdata)` returns `(nuclei_key, cells_key)`, with `cells_key` falling back to the nuclei layer for nuclei-only Zarrs (`--cell-expansion 0`, or older outputs written when a missing pixel size skipped expansion). `run_network_analysis`, `assign_spots_to_cells`, `export_to_qupath`, and the sweeps all default to it.

- **Distances are in µm; coordinates are in pixels**: Sopa's OME-TIFF reader uses an identity transform, so all shape/centroid coordinates are pixels. `--max-edge-distance`, `--radius-min/max` and `--cell-expansion` are µm and are converted at the point of use (`value_um / pixel_size_um`). The pixel size comes from `get_pixel_size` (OME metadata) or `--pixel-size`; `cli.main()` calls `parser.error` at startup if neither exists. It is also stored in `table.uns['cellsurvey']['pixel_size_um']` (`utils.store_pixel_size`) so read-only consumers like the stability sweep can convert µm parameters; `utils.resolve_pixel_size` reads it (explicit value wins) and raises for older Zarrs, which need `--pixel-size`. `run_network_analysis(max_edge_distance=µm)` converts internally; `compute_louvain_communities`/`delaunay_edges` take the threshold in coordinate units (pixels).

- **Delaunay edges are deduplicated**: `delaunay_edges` sorts each pair so an edge shared by two triangles is counted once (`n_edges` ≈ 3× cells at most). Before this, `(a, b)`/`(b, a)` duplicates inflated `n_edges` in `summary.json` by ~⅓ and computed each correlation twice (networkx had already merged them, so communities were unaffected).

- **InstanSeg needs a correct pixel size**: it rescales the image to the model's training resolution; a wrong value can make it find few or no cells.

- **InstanSeg must be run with `target='all_outputs'`, not via `sopa.segmentation.instanseg`**: InstanSeg only pairs nuclei with cells (`resolve_cell_and_nucleus_boundaries`) when both outputs are requested in one call. Sopa's wrapper asks for one target per call, which skips the pairing. On a 29-plex COMET image that left about 25% of cells touching no nucleus and 19k nuclei in no cell. Hence our own per-patch callable `_instanseg_paired_patch`. Sopa still stitches each layer separately and the polygon ids do not correspond, so nuclei and cells are matched spatially in `link_nuclei_to_cells`. In a multinucleate cell only the largest nucleus is linked; the others get cluster/community `-1`. Pass `clip_limit=0, gaussian_sigma=0` to `custom_staining_based`; its defaults apply CLAHE and smoothing, which `sopa.segmentation.instanseg` does not.

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

- **`--resume-from` enables crash recovery**: If a Zarr already exists at the expected path, you can skip image loading (and reuse its spots) and resume from segmentation. When used with the segmented Zarr reuse logic, this provides two levels of checkpoint restart.

- **Channel subset loading**: The full multichannel image is not loaded into memory for blob detection. Only the channels specified by `--channels` are read (lazily, from the 1024²-chunked `<prefix>.zarr`); `detect_blobs_tiled` materialises one channel at a time.

- **AnnData `.X` can be sparse or dense**: `sopa.aggregate()` may produce either a scipy sparse matrix or a dense numpy array depending on the input data size and sopa version. The intensity extraction in `cli.main()` (and `stability.extract_intensity_matrix`) handles both with `hasattr(measurements.X, 'toarray')`. Never assume `.X` is sparse.

- **GeoJSON output path**: The QuPath GeoJSON path is configurable via `--geojson-path`. By default it is `<output_file>.geojson` (`.zarr` stripped from `-o`), next to the Zarrs, so runs with different `-o` don't overwrite each other. All plot/PNG artifacts and `summary.json` go to `--plot_dir`.

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
- The pipeline writes Zarr files at two points: the initial Zarr after image loading (spots, if any, are then added to it with `write_element`), and the segmented Zarr after segmentation + aggregation. All writes are wrapped in try/except and exit with code 1 on failure.
- The Zarr write-then-immediate-read pattern between stages 3 and 4 materializes a clean checkpoint. If segmented Zarr reuse kicks in (stage 4), the read of the initial Zarr is skipped entirely.
- Segmentation channels come from `--seg-channels` (indices into `unique_channels`): Stardist takes exactly one (default `0`); InstanSeg defaults to all channels.
- Dynamic imports inside except blocks: `import shutil` is imported inside the segmented Zarr corruption handler to avoid pulling it in unnecessarily.
- `run_network_analysis()` in `network_analysis.py` returns a plain dict with keys `cell_ids`, `community_labels`, `cluster_labels`.
- `run_network_analysis()` accepts `output_dir` parameter for plot output (set to `args.plot_dir` from `cli.py`).

## Known issues (code review, 2026-11-11, untested)

Open:

- **Slow edge weighting** (`network_analysis.compute_louvain_communities`): calls `np.corrcoef` once per edge in a Python loop (~3 edges per cell, so millions of calls on whole slides). The stability sweep repeats this for every parameter combination. Fix: z-score each row of `imat` once, then compute `(z[a] * z[b]).mean(axis=1)` over the edge arrays.
- **Slow per-cell lookups**: `run_network_analysis` (cluster labels and `imat`) and `stability.run_stability_sweep` (`imat`) do one `.loc[...]` per cell. Replace with `intensity_matrix.reindex(boundaries.index).fillna(0).to_numpy()` and the same for `obs['kmeans_cluster']`.
- **k-means makes near-empty clusters and is unstable** (`utils.cluster_data`, k=10 on standardised raw mean intensities): on `test_crop` the GPU and CPU Stardist runs differed by only 2 cells, yet the cluster sizes reshuffled substantially and one cluster held a single cell. The InstanSeg run had two clusters of 2 cells each. Outlier cells are getting their own clusters. On identical input k-means is reproducible (seed 42): resuming the same `_seg.zarr` gave identical cluster sizes. So the reshuffling comes from small segmentation differences amplified by k-means, not from randomness. Relevant to the Planned Normalisation, the `--mode clusters` sweep, and PIPEX-style robust marker scores. The z-scored `cluster_intensity_heatmap.png` (GPU Stardist run on `test_crop`) shows two causes:
  - **Outlier/debris cells**: the 1-cell cluster 8 sits at z≈22 for CD56 and PD-1, and z≈10–12 for FOXP3, CD68 and PD-L1. A cell bright in almost every channel is likely debris or an artefact, not a cell type. Options: clip or `log1p` intensities before scaling, winsorise each channel (e.g. at the 99.5th percentile), or flag and exclude cells that are bright in most channels before clustering (keep them in the table with a label such as `-1`).
  - **Background channels used as markers**: cluster 3 (17 cells) is high in TRITC and Cy5, unstained background channels that appear twice in the COMET channel list, and raised across many markers, so it is likely autofluorescent cells. All aggregated channels currently feed k-means and the Louvain edge weights. Add an option to exclude channels from clustering by name or index (e.g. `--exclude-channels TRITC,Cy5,AF`), still aggregating them so they stay in the table. This would also cover `AF` from the OMERO qptiff findings.
  - Also consider a minimum cluster size, or choosing k from the cluster sweep (inertia plus a subsampled silhouette score).

- **With `--detect-blobs`, clustering uses spot counts instead of marker intensities** (`cli.py`, pre-existing; the old detect-then-write flow made the same `sopa.aggregate` call): when `sopa.aggregate` gets spots (`aggregate_genes=True`), the table's `.X` holds the **spot counts per cell** (one column per spot channel) and the channel intensities move to `obsm['intensities']` (Sopa `compute_table`). `cli.py` builds `intensity_df` from `.X`, so k-means, the Louvain edge weights, the cluster heatmap, UMAP/Leiden and the GeoJSON measurements all run on a few count columns. On `test_crop` with the default `--channels 9,10,11,12`: `Columns after removing duplicates: ['CD20-C6', 'CD4-C5', 'CD56-C6', 'CD8-C5']` instead of 21 channels. Symptoms: `invalid value encountered in divide` (correlation of all-zero count rows), Scanpy `Some cells have zero counts`, UMAP `Spectral initialisation failed`. Fix: always take the intensity matrix from `obsm['intensities']` when present (else `.X`), and keep the counts in the table as extra per-cell data. **Resolved (2026-11-19).** HPC-checked on `test_crop` with `--detect-blobs`:
  - the log lists 21 intensity channels and 4 spot channels; `--cluster-on intensities` gives 21 features and the divide, zero-count and spectral warnings are gone;
  - on resume, `spots` (4 features) and `both` (25 features) work;
  - `--cluster-on spots` without `--detect-blobs` stops at startup;
  - both sweep modes work with `--cluster-on`.

  Synthetic tests (both Sopa table layouts, through a Zarr round trip):
  - intensities and counts are extracted correctly, and the `obsm` DataFrame survives the Zarr round trip;
  - `--cluster-on intensities` gives k-means labels identical to the old path;
  - spot features are standardised, and in `both` each group's total variance is exactly 1;
  - both sweep modes run in all three modes;
  - the GeoJSON gets `mean` and `spot count` measurements.

  Design:
  - **New flag `--cluster-on {intensities,spots,both}`** (default `intensities`), also in the stability sweep, which defaults to the mode recorded in `table.uns['cellsurvey']['cluster_on']` so a sweep uses the same features as the run it tests (`utils.resolve_cluster_on`; falls back to `intensities` for Zarrs that predate it; an explicit flag overrides). `spots`/`both` need spots: a new run without `--detect-blobs` stops at startup.
  - **Spot features**: spot *density* (spots per 100 µm² of cell area, so clusters don't just track cell size), then `log1p`, then z-scored.
  - **Intensity features**: z-scored per channel.
  - **`both`**: each group (intensity channels, spot channels) is weighted by 1/√(number of its columns), so the groups contribute equally.
  - **Same features everywhere**: k-means, the Louvain edge weights and UMAP/Leiden (`use_rep`) all use the same scaled features. This changes default community and UMAP results: edges were previously weighted by raw-intensity correlation, and UMAP used Scanpy's `normalize_total` + `log1p` on `.X`.
  - **GeoJSON**: exports both `Cell: <channel> mean` and `Cell: <channel> spot count` (raw counts).
- **Correlation edge weights are weak with few features** (`network_analysis.compute_louvain_communities`): each edge is weighted by `1 + Pearson r` between the two cells' feature vectors, computed *across features*. With few features this carries little information: with 2 features r is always ±1, so `--cluster-on spots` with two spot channels gave 59 communities on 300 synthetic cells, against 8 for `intensities`. It matters most for `spots` with few spot channels. Options: a distance-based similarity (e.g. `exp(-‖a−b‖² / σ²)`) instead of correlation, or unweighted edges below a minimum feature count. Revisit together with the vectorised edge weighting.
- **Whole-slide run findings** (OMERO Image 602778, 617,527 cells, 2026-11-20; see the OMERO HPC test results):
  - **`cell_density.png` has one panel per community**: `n_rows = n_clusters + n_communities` (`cli.py`), so 1,723 communities gave a **1,000 × 1,039,800 px, 65 MB PNG** that took ~4 min and can't be viewed. Fix: plot only the largest communities (e.g. the top 20 by size) and note how many were left out.
  - **Most communities are isolated cells**: 1,723 communities, but **52 hold 99.1% of cells**; 1,182 are single cells and 1,556 have ≤5 cells. These are cells with no Delaunay neighbour within 50 µm (tissue edges, debris, gaps), so `n_communities` in `summary.json` is misleading. Options: report counts above a minimum size (e.g. `n_communities_ge_10`), and/or relabel communities below a minimum size as `-1` ("unassigned"), as for unlinked nuclei.
  - **`<prefix>.zarr` also gets `shapes` and `tables`**, not only the image (and spots): `nucleus_boundaries`, `cell_boundaries`, `image_patches` and the table are written into it. Sopa writes new elements straight into the Zarr a dataset was read from (`dataset = spatialdata.read_zarr(zarr_path)` before segmentation). This is pre-existing, and the same data is written again to `_seg.zarr`, so it doubles the disk use for shapes and table. To check: whether a later re-run without `--resume-from` cleans these up (stage 1 overwrites the Zarr), and whether reading the image Zarr into memory or a separate store avoids it. Docs that say `<prefix>.zarr` holds only the image + spots are inaccurate until then.
  - **GeoJSON is 1.9 GB** for 617k cells (8 channels): heavy for QuPath. See the PIPEX note on compressed `.pbf` export, or offer gzip / cells-only / no-measurements options.
  - **InstanSeg runs the model twice** (once per layer, each computing both outputs; see `_instanseg_paired_patch`): ~26 + ~48 min here. Caching each patch's cell mask from the first pass (e.g. in the Sopa cache dir) and returning it in the second would roughly halve the GPU time.
  - **27% of InstanSeg cells had no nucleus** and were dropped (226,499), against 9% on `test_crop`. Plausible for brain tissue (fragments of processes), but check visually (e.g. the GeoJSON in QuPath) before relying on the default.
  - ~~**`morphology_by_cluster.png` plots area in px²**~~: now µm² (to-do step 3).

  **To-do plan for these findings (2026-11-20; not started).** Steps 1–4 are small and planned *before* fast-forwarding `feature/omero` into `main`; steps 5–8 are larger and can go on a new branch after the merge. Each step lists how to check it.

  *Before the merge:*
  1. **Minimum community size + size-ordered labels** (`network_analysis.py`; also the stability sweep). **Done (2026-11-21), HPC check pending**: `--min-community-size` (default 10), `filter_small_communities`, the new `summary.json` / sweep-summary fields, the value stored in `uns` (the sweep reads it), and the GeoJSON colour count no longer counting `-1`. Synthetic test (3 tissue blobs + 40 isolated cells, through a Zarr round trip): 50 raw communities → 3, 47 cells unassigned, sizes descending, the GeoJSON's `Unassigned` count matches, and the sweep at the run's parameters reproduces its labels exactly. Decisions: minimum 10, renumber by size. Original plan:
     - add `--min-community-size`; proposed default 10 cells;
     - relabel smaller communities as `-1` ("unassigned", grey in the GeoJSON, which already handles `-1`; nuclei mapping and plots must handle it);
     - renumber the rest by size (0 = largest);
     - `summary.json`: report the communities above the minimum, the number of unassigned cells, and the raw Louvain count;
     - the sweep applies the same minimum, so counts match.

     Check: resume `omero_602778` → ~50 communities instead of 1,723.
  2. **`cell_density.png`** (`cli.py`; after step 1). **Done (2026-11-21), HPC check pending**: moved into `plot_cell_density` (plain arrays, so it can be tested without a run); the 20 largest communities chosen by size (not label order), `-1` never plotted; a 4-column grid with constrained layout instead of one column. Synthetic test at slide scale (617,527 cells, 10 clusters, 52 communities + 5,109 unassigned): 2000 × 3200 px, ~6 s; labels out of size order and an all-unassigned run also render correctly. Original plan: plot the k-means clusters plus only the largest communities (proposed default 20), state the number left out in the title, and loop over actual labels rather than `range(n_communities)`. Check: normal image size, drawn in seconds.
  3. **Area in µm²** (`cli.py`, morphology plot): `area_px × pixel_size_um²`, axis labelled µm². **Done (2026-11-21), HPC check pending.** Also fixed in `cell_density.png` (found while testing step 2):
     - the maps were **upside down** relative to the image (`origin='lower'`, so y increased upwards); now `origin='upper'` with a flipped extent;
     - each community panel **zoomed to its own extent**, hiding where it sits; all panels now share the whole-tissue range, which also skips NaN centroids that would have crashed `histogram2d`'s automatic range.

     Checked with a synthetic layout (a top-left blob and an L along the left and bottom edges): it renders as in QuPath.
  4. **Docs**: this file, `docs/outputs.md` and `docs/parameters.md` (new flag(s), new `summary.json` fields, the "unassigned" label).

  → One resume run on `omero_602778` checks steps 1–3 together (resumes skip segmentation, ~20 min). Then fast-forward into `main`.

  **Decisions (2026-11-21)**: minimum community size 10; renumber by size; `cell_density.png` capped at a fixed 20 communities.

  *After the merge:*
  5. **Stop Sopa writing into `<prefix>.zarr`**:
     - first confirm the mechanism (Sopa writes new elements when the dataset is Zarr-backed) and where Sopa's segmentation cache goes;
     - likely fix: detach after reading (`dataset.path = None`), so new elements stay in memory until the single `_seg.zarr` write while the image is still read lazily from the Zarr;
     - risks to measure on `omero_602778`: memory, and the cache location.

     Check: `<prefix>.zarr` holds only `images` (+ `points`) after a run; `test_crop` results are identical; all resume paths still work, especially "shapes but no table".
  6. **GeoJSON size**:
     - measure what takes up the 1.9 GB first;
     - cheap savings: round coordinates to 2 decimals (currently full float precision) and simplify polygons slightly (~0.5 px), plausibly halving the file;
     - then, if needed: options for no measurements, cells only, or gzip (check whether QuPath imports compressed GeoJSON);
     - longer term: a `.pbf` export (see the PIPEX reference).

     Check: file size, and that it loads and looks right in QuPath.
  7. **Run the InstanSeg model once**:
     - in the nuclei pass, save each patch's cell mask to a temporary cache keyed by a hash of the patch pixels (Sopa's callable gets no patch coordinates);
     - the cells pass returns the cached mask instead of running the model;
     - nuclei and cells then come from identical model output, which should reduce unlinked nuclei;
     - ~16 MB per patch (~10 GB for this slide): node-local temporary disk, deleted afterwards.

     Check: identical results on `test_crop`; ~half the segmentation time on `omero_602778` (~35 min saved).
  8. **27% unnucleated cells**:
     - mainly a visual review: the GeoJSON in QuPath, or a `--keep-unnucleated` run;
     - code support: record the dropped count in `summary.json`, and optionally export the dropped cells as a separate GeoJSON layer for review;
     - then decide whether the default should change.
- **Windows: `UnicodeEncodeError` when output is redirected**: with stdout piped or sent to a file, Windows Python encodes with cp1252, and printing channel names such as `αSMA` crashes the run. Workaround: `PYTHONIOENCODING=utf-8`. Fix: `sys.stdout.reconfigure(encoding='utf-8')` (and stderr) at the top of `run.py`. Interactive consoles and Linux are unaffected.
- **`python -m cellsurvey.stability` warns `'cellsurvey.stability' found in sys.modules…`** (harmless, pre-existing): `cellsurvey/__init__.py` imports `stability` (and `cli`), so running the module with `-m` loads it twice. Fix: import lazily in `__init__`, or give the sweep its own entry point.
- **Windows: Stardist on CPU is extremely slow**: a local run on `test_crop` used more than an hour of CPU time in segmentation (4 patches), against about 14 s on an HPC CPU node. Not investigated; win-64 is only for local testing.
- **TF gets imported late in the default env** (not a bug): after the network step, Scanpy's UMAP imports umap-learn, which optionally imports TensorFlow (for `ParametricUMAP`) when it is installed. It costs some import time and logs TF's CPU-feature message. The instanseg env has no TF.

Checked and resolved / not a bug:

- ~~Unused source-image read on resume~~: the `--resume-from` branch no longer opens the input with `BioImage` just to read channel names it never used. Since the stage-1 reorder it doesn't read the source for the pixel size either, unless the Zarr predates stored pixel sizes.
- ~~Delaunay on too few cells~~: `delaunay_edges` no longer raises. Fewer than 2 points give no edges. With 2 points, or when Qhull can't triangulate because all points are collinear or identical, consecutive points along the main axis are joined, which is the Delaunay graph's limit, and the distance filter still applies. The full pipeline still needs enough cells for k-means (`--n-clusters`), Scanpy neighbours and UMAP, so tiny inputs fail later anyway; this fix protects the sweeps and library callers.
- ~~Pixel size can change on resume~~: when `--resume-from` reuses a `_seg.zarr` with a table, `cli.py` reads the stored pixel size. If it differs from the current one (`rtol=1e-3`), the run warns and continues with the stored value, so the network, radius conversions and re-stored value stay consistent with the existing shapes. Older Zarrs with no stored value keep the current one. The startup check still needs a pixel size from metadata or `--pixel-size`, even when resuming.
- ~~Unreachable InstanSeg warning~~: `run_instanseg` now raises `ValueError` at the start, before importing torch, if `pixel_size_um` is missing or not positive. This protects library callers; the CLI already stops at startup.
- ~~Sweep `--pixel-size` not validated~~: `utils.resolve_pixel_size` now raises `ValueError` for zero, negative or NaN values, whether passed in or stored in the Zarr. The sweep reports it via `parser.error`. `cli.py`'s own check also rejects NaN.
- ~~Misleading cluster intensity heatmap~~: `cluster_intensity_heatmap.png` now z-scores each channel across all cells (as k-means sees it), then averages per cluster, so the "z-score" label and the diverging map centred on 0 are correct.
- ~~GeoJSON default path ignored `-o`~~: the default is now `<output_file>.geojson`, so consecutive runs no longer overwrite `./qupath_export.geojson`.
- ~~Cell ids may not match between shapes and table~~: `sopa.aggregate` → `parse_table` sets the table's `obs_names` and the cell layer's index to the same string ids (`str_cell_id`, sopa 2.3.0), so `cell_id in table.obs.index` matches. This could only break if the shapes were re-indexed after aggregation.
- ~~Fallback layer looked up eagerly (`shapes.get('cell_boundaries', shapes['stardist_boundaries'])`)~~: replaced by `resolve_shape_keys`.
- ~~`n_edges` inflated by double-counted Delaunay edges~~: `delaunay_edges` sorts and de-duplicates pairs.
- ~~Mixed units (`--max-edge-distance`, `--radius-*` in px)~~: all distance flags are now in µm.

## Planned: Parameter Stability Sweep

**Goal**: Quantify how robust community assignments are under parameter variation, producing per-cell confidence scores and consensus niches.

### Parameters to explore
| Parameter | Range | Rationale |
|---|---|---|
| `community_resolution` (Louvain) | 0.05–1.0 | Directly controls community granularity |
| `max_edge_distance` | 25–100 µm | Changes which cells are neighbors in the Delaunay graph |

> **`n_clusters` is NOT part of the community sweep.** k-means clustering and
> Louvain community detection are independent: Louvain weights its Delaunay edges
> from the `--cluster-on` features, never from the k-means labels. So `n_clusters`
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

**Status**: implemented and wired — `--cell-expansion` builds non-overlapping `cell_boundaries` via a marker-controlled watershed (`_watershed_expand`: rasterize → EDT → watershed → vectorize); aggregation, clustering, the Louvain network, spot assignment, export, and the sweeps all target the cell layer, and cluster/community are mapped back onto `stardist_boundaries` (nuclei) via a recorded `nucleus_id`. **Remaining**: the separate `nucleus` table and `--cluster-regions` multi-region clustering.

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
- If physical pixel size is unavailable and `--pixel-size` isn't given, the run stops at startup (all distance flags are µm).
- **Related cleanup — done**: `--max-edge-distance` (default 50 µm, ≈ p90 of Delaunay edge lengths on a dense COMET image: median 9 µm, p90 48 µm, p99 110 µm) and `--radius-min/max` (default 0–280 µm, equivalent to the old 0–1000 px at 0.28 µm/px) are now in µm.

### Approach (QuPath-style expansion)
QuPath detects nuclei, then dilates each nucleus by a `cellExpansion` distance, using a distance-transform/watershed expansion so neighbouring cells stop where they meet (no overlap).

Two implementation tiers:

1. **v1 — simple buffer (start here)**: `shapely.buffer(nucleus_geometry, expansion_px)`. Fast, trivially parallel over rows, but produces **overlapping** polygons where cells are dense; overlap double-counts pixels in aggregation.
2. **v2 — non-overlapping (QuPath-faithful)**: clip each buffered cell to the Voronoi cell of its centroid so adjacent cells stop where they meet. **Implemented with our own `_resolve_overlap()`** (shapely `voronoi_polygons` + STRtree). **Note**: `sopa.shapes.expand.remove_overlap()` — and `sopa.aggregate(no_overlap=True)`, which calls it — were found to empty most cells on dense data (1742/2107 → NaN geometries), so they are **not** used.

v1 + v2 are implemented and wired through `cli.py` (buffer + custom Voronoi clip).

### v3 — marker-controlled watershed (implemented)

**Problem**: v2's Voronoi clip produced straight, "severe" boundaries between
adjacent cells — Voronoi lines at the midpoint between centroids, which ignore
cell morphology. v3 replaces it with a raster-based watershed.

**Implemented** in `_watershed_expand` (`cellsurvey/segmentation.py`): rasterize
nuclei → `distance_transform_edt` → `skimage.segmentation.watershed` (capped at
`cellExpansion`) → vectorize back to polygons. **Gotcha**: `find_contours` on a
region that exactly fills its bounding box returns nothing (only corner cells
survive), so each region is padded with a background border before contouring.

**Reference**: QuPath's `WatershedCellDetection` (ImageJ-based; see
[`qupath.imagej.detect.cells.WatershedCellDetection`](https://qupath.github.io/javadoc/docs/qupath/imagej/detect/cells/WatershedCellDetection.html)).
The algorithm:
1. Label the nuclei (integer markers).
2. Compute the **Euclidean Distance Transform (EDT)** of the nucleus binary mask (distance to nearest background).
3. **Negate** it (nucleus interiors become minima).
4. Run a **marker-controlled watershed**, flood **capped at `cellExpansion`** (pixels whose distance ≥ `cellExpansion` are never flooded).
5. Non-overlap is intrinsic (watershed ridges where two floods meet).

**Python equivalent** (deps already present: `scipy`, `skimage`):
```python
edt = scipy.ndimage.distance_transform_edt(nucleus_mask)
labels = skimage.segmentation.watershed(
    -edt, markers=nucleus_labels,
    mask=(edt <= cell_expansion_px), connectivity=1,
)
# then vectorize labels -> polygons
```

**Key change**: this is **raster-based** (rasterize nuclei → watershed →
vectorize), unlike the current **vector-based** buffer → Voronoi clip. The
watershed boundary follows the distance-transform ridges (curving with nucleus
shape/size), much less severe than straight Voronoi lines — but still not
membrane-aware.

**Implementation plan** (done):
1. Rasterize `stardist_boundaries` to a labeled nucleus image (markers) at image resolution.
2. `distance_transform_edt` + `watershed` with the `cellExpansion` cap.
3. Vectorize the labeled result back to `cell_boundaries` polygons (preserve `cell_id`).
4. Replace `_resolve_overlap` (Voronoi) in `expand_nuclei` with the watershed.

**Considerations / open questions**:
- **Resolution**: watershed on the image raster gives boundaries at pixel resolution (good); µm/px is already known via `get_pixel_size`.
- **Tiling**: large images need a tile-based watershed (with overlap, like Stardist patches) to bound memory.
- **Watershed vs Cellpose**: watershed is deterministic, no training, no new heavy deps, and mirrors QuPath — a cheap immediate improvement. Cellpose learns *actual membranes* (better boundaries) but needs PyTorch + a pretrained model (see the TF/Torch note). Recommendation: do the watershed next, keep Cellpose as the longer-term goal.

### Aggregation & measurements (two tables)

`cli.py` currently calls `sopa.aggregate()` once (over `stardist_boundaries`). With whole-cell expansion, aggregate **two solid regions** (no holes):

| Compartment | Source | How |
|---|---|---|
| Whole cell | `cell_boundaries` | `sopa.aggregate(..., shapes_key='cell_boundaries')` → `table` |
| Nucleus | `stardist_boundaries` | `sopa.aggregate(..., shapes_key='stardist_boundaries')` → `nucleus` |

Each table's `.obs` index is the same `cell_id`, so `table` and `nucleus` join cleanly. Decide whether the *primary* `table` should be whole-cell (recommended — most biologically meaningful) or nucleus (current behaviour); open question below.

**Cytoplasm** is *not* computed in this phase — it requires a proper cell/cytoplasm segmentation, deferred to the Cellpose path. (This also sidesteps a known Sopa limitation: `sopa.aggregate()` ignores polygon holes, so a donut-shaped cytoplasm couldn't be aggregated correctly anyway.)

### Downstream consumption

- **Region selection for clustering/network** — **done** (single region): clustering and the Louvain network run on the whole-cell `table` / `cell_boundaries`. `--cluster-regions` (multi-region concatenation) remains a future option.
- **Spot assignment** — **done**: `assign_spots_to_cells()` assigns spots to `cell_boundaries` (whole cell). Nuclear-vs-cytoplasmic spot classification (a second `within` predicate) is still a future refinement.
- **GeoJSON export** — **done**: `export_to_qupath()` exports `cell_boundaries` (cells) with cluster/community. Per-compartment measurements (nucleus/cytoplasm/cell) are still future.
- **Nuclei labels** — **done**: cluster/community are mapped back onto `stardist_boundaries` via the `nucleus_id` column recorded in `expand_nuclei` (survives `sopa.aggregate`'s re-indexing).

### Implementation plan

1. **Read physical pixel size** — **done**: `get_pixel_size(imagepath)` in `cellsurvey/segmentation.py` (returns µm/px or `None`).
2. **Expand nuclei (v3)** — **done**: `expand_nuclei(nuclei_gdf, expansion_um, pixel_size_um)` rasterizes nuclei and expands them via a marker-controlled watershed (`_watershed_expand`), producing non-overlapping cells. Records `nucleus_id` on each cell.
3. **`--cell-expansion` flag** — **done** (µm, default `5.0`).
4. **Build & store `cell_boundaries`** — **done**: inserted after Stardist, before aggregation; `stardist_boundaries` kept.
5. **Aggregate over cells** — **partially done**: `sopa.aggregate` now targets `cell_boundaries` → `table` (whole-cell). The separate `stardist_boundaries` → `nucleus` table is **pending**.
6. **Thread region through downstream** — **done**: `export.py`, `assign_spots_to_cells`, `network_analysis.py`, and `stability.py` now default to / are passed `cell_boundaries`; cluster/community are mapped onto `stardist_boundaries` via `nucleus_id`. `--cluster-regions` not added (single whole-cell table for now).

### Done: InstanSeg (whole-cell segmentation)
`--segmenter instanseg` (separate `instanseg` pixi env) segments nuclei **and** cells directly — InstanSeg run with `target='all_outputs'` through `sopa.segmentation.custom_staining_based` and our `_instanseg_paired_patch` (see the InstanSeg gotcha above for why not `sopa.segmentation.instanseg`) — no expansion approximation, and no whole-image watershed. This is the intended long-term replacement for Stardist + expansion; the Cellpose note below is superseded unless InstanSeg proves inadequate. A cytoplasm compartment (cell minus nucleus) is still not computed.

### Future: Cellpose (proper cytoplasm segmentation)
Integrate **Cellpose** for true whole-cell segmentation (cell + cytoplasm), replacing the nucleus-expansion approximation. Cellpose is the dominant pretrained whole-cell segmenter for fluorescent microscopy, and would also yield a genuine **cytoplasm** compartment (cell minus nucleus) — the piece this phase deliberately does not approximate. Note: Cellpose is PyTorch-based (see the TF/Torch co-existence note) — recommend a separate pixi environment or a deferred, opt-in `--segmenter cellpose` path. Recorded here for later; not in scope for the expansion work.

### Open questions
- **Primary `table` region**: whole cell (recommended) vs nucleus (current). Changing `table` to whole-cell alters existing results and makes prior `_seg.zarr` outputs non-comparable — acceptable if flagged as a breaking change.
- **Always-on vs opt-in**: default `--cell-expansion 5.0` (on) vs `0` (off, nuclear-only). A non-zero default is the point of the feature, but changes results vs today.
- **Fixed vs per-nucleus radius**: QuPath uses a fixed `cellExpansion`; a per-nucleus radius (scaled to nucleus area) is more accurate but adds a parameter.

### Risks
- ~~v1 overlaps double-count pixels and inflate aggregated intensity~~ — resolved: the v3 watershed produces non-overlapping cells (only the buffered fallback for nuclei with no watershed cell can overlap).
- **Watershed memory**: `_watershed_expand` rasterizes the whole image extent in one array (int32 labels + float64 EDT ≈ 12 bytes/px — ~30 GB for 50k×50k). Needs tiling for large images; InstanSeg (`--segmenter instanseg`) avoids it entirely.
- Two tables and a region-selection flag widen the downstream surface (clustering, network, export, sweeps) — easy to miss a hardcoded `'stardist_boundaries'`/`'table'` reference; grep thoroughly.
- ~~Missing OME physical pixel size breaks µm conversion~~ — resolved: the run stops at startup and asks for `--pixel-size`.

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

## Planned: Overview report

**Status**: planned (2026-11-21), not started.

**Goal**: every run writes a single, self-contained **overview report** (`<plot_dir>/report.html`) that a biologist can open in a browser, email, or print to PDF. It should say who is in the tissue (clusters), where they are (spatial maps over the tissue image), how the tissue is organised (communities), and which results to distrust (data-quality flags). Today the outputs are loose PNGs plus `summary.json`, with no explanation.

**Model**: the hand-written overview in [cellsurvey-analysis `explorations/cellsurvey_wholecell`](https://github.com/FrancisCrickInstitute/cellsurvey-analysis/tree/main/explorations/cellsurvey_wholecell):
- `overview_analysis.py` reads an exported `cells.csv`, plus the DAPI channel from a downsampled pyramid level (`s4`) of the `_seg.zarr`, and the sweep CSVs. It writes figures and a Markdown report.
- `make_report_html.py` converts that Markdown into one HTML file with the figures embedded as base64.
- `overview_analysis_filtered.py` repeats the analysis after removing the debris and collagen clusters. It re-clusters and rebuilds communities itself, and its docstring calls this an *approximation* of CellSurvey that should eventually move upstream.

### What to keep (generic) and what to leave out (dataset-specific)

From the example, kept and generalised:
- **Header** (new; not in the example): the source (file path, or OMERO host + Image ID + name), date, CellSurvey version/commit, segmenter and model, pixel size, cell count, channels, the key parameters (`--n-clusters`, `--cluster-on`, `--community-resolution`, `--max-edge-distance`, `--cell-expansion`), and counts from `summary.json` (cells dropped as unnucleated or dim, communities, edges).
- **"How to read this"**: the plain-language explanation of z-scores and of *clusters* (what cells express) versus *communities* (where cells sit). This text can be reused almost word for word.
- **Marker distributions**: one log-scale histogram per marker on a robust z scale (arcsinh, then median/IQR), with an illustrative "positive" line at z = +2, the percentage above it, and the warning that most markers are unimodal, so the cutoff is a choice.
- **Cluster profiles**: a heatmap of mean z per cluster (in addition to `cluster_intensity_heatmap.png`, which stays), plus a table of size, percentage, the top 5 markers and an automatic name built from the top 2 markers (e.g. `CD31/Fibronectin`).
- **Spatial maps over the tissue image** (DAPI in greyscale):
  - one panel per cluster;
  - one map of all clusters, with each hexbin coloured by its dominant cluster;
  - a map of "positive" cells for each marker. The example's hand-picked lineage list becomes all markers, capped at about 30 panels.
- **Communities**:
  - outlines over the tissue image, drawn as smoothed density contours per community;
  - a stacked bar of composition by cluster, grouped by dominant cluster;
  - a short auto-generated narrative: which clusters dominate the largest communities, how many communities are over 80% one cluster, and how many are dominated by flagged clusters;
  - a resolution → community-count table when `community_sweep_summary.csv` exists next to the Zarr.
- **Data-quality flags**: a generic version of the example's debris analysis.
  - The example flags debris with a fixed 4-marker signature (CD68 + FoxP3 + LamininA5 + H2AX). That doesn't transfer to other panels.
  - Use the example's other test instead: a cell is **non-specific** if it is bright (z > 2) in at least N markers, with N scaled to the panel (e.g. `max(3, 25%)` of the markers).
  - Flag clusters that are more than 50% non-specific. Say "non-specific / likely autofluorescence or debris" rather than asserting they are debris.
  - For each marker, report the share of "positive" cells that are non-specific: clean, suspect (> 25%) or dominated (> 50%).
  - Also list background or control channels found by name (`TRITC`, `Cy5`, `AF`, …) and unassigned cells or tiny communities.

New, not in the example:
- **Spots section** (only when the table has spot counts, i.e. after `--detect-blobs`):
  - a per-channel distribution of spot density (spots per 100 µm² of cell area, matching `build_clustering_features`);
  - the share of spots assigned to cells;
  - spot density per cluster (a heatmap like the marker one);
  - a map per spot channel over the tissue image (a hexbin of spot density).
  
  The section is shown whatever `--cluster-on` is set to; the header says which features the clusters were built from.

Left in cellsurvey-analysis (project-specific interpretation, not pipeline output): the thymus feasibility questions, the measurement-compartment table, CD68 and lineage attributions, and the CD31 community investigation. A later option could let users pass a panel file (marker → lineage, channels to treat as autofluorescent) to restore lineage names, but that's not part of v1.

### Design
- **New module `cellsurvey/report.py`**, which only reads the `_seg.zarr`, like the stability sweep. It's called at the end of `cli.main()` by default (`--no-report` to skip), and can also run standalone, `python -m cellsurvey.report --zarr <prefix>_seg.zarr -p <dir>`, so a report can be regenerated or restyled without re-running anything. It must not import TF, Torch or omero-py.
- **Inputs**, all taken from the `_seg.zarr`:
  - per-cell intensities via `utils.extract_feature_tables` (never `.X` directly);
  - `obs` `kmeans_cluster` / `community` / area;
  - centroids from the cell layer (`resolve_shape_keys`);
  - `table.uns['cellsurvey']` and `attrs['cellsurvey']` (pixel size, `cluster_on`, source);
  - optionally, the sweep CSVs next to the Zarr.
- **Tissue background**: read one channel from the image pyramid. Pick the level whose longest side is ≤ ~4000 px, rather than a hard-coded `s4`. Find DAPI by name (shared with the Normalisation plan's auto-detect), falling back to the first `--seg-channels` channel. If the image isn't readable, fall back to points on white, as the example does. To check: whether the `_seg.zarr` holds the image, or only `<prefix>.zarr` does (see the "Sopa writes into `<prefix>.zarr`" finding). Coordinates are pixels with y pointing down, so use `imshow(origin='upper')`; the example's `origin='lower'` reflects its CSV's frame.
- **Output**: build the HTML directly (Python string templates + `html.escape`; no Markdown step, no new dependency) with figures embedded as base64 PNG, light print-friendly CSS, and a contents list. Also keep the PNGs in `<plot_dir>/report/`. Record the report's settings (z cutoff, N for non-specific) in its header.
- **Scale** (617k cells on Image 602778): everything is vectorised (`histogram2d`, `hexbin`, `groupby`). Cap the panels: the top 20 communities by size, ≤ 30 markers. Rasterise scatter plots (subsample to ~200k points for the per-cluster panels) so figures stay small, and target an HTML file of a few MB. Time it on `omero_602778`; aim for under 2 min.
- **Robust z vs clustering features**: the report's arcsinh + median/IQR z is for display and differs from the plain z-score k-means uses (`build_clustering_features`). State this in the report. Revisit if the clustering moves to robust scores (see Known issues: k-means near-empty clusters, and the Normalisation plan).

### Dependencies and related work
- **Builds on Known issues to-do steps 1–2** (minimum community size, `-1` for unassigned, size-ordered labels, capped `cell_density.png`). Without them, the community sections are swamped by single-cell communities (1,723 on 602778). The report must treat `-1` as "unassigned", never as a community.
- **`--exclude-channels`** (Known issues, k-means): the report should leave the same channels out of the marker sections, as the example drops DAPI and the TRITC/Cy5 controls.
- **The filtered variant** (`overview_analysis_filtered.py`) is a request for CellSurvey itself: exclude flagged cells (and chosen clusters or channels) *before* clustering and communities, so the analysis repo needn't re-implement them. Later, the report's non-specific flag could feed an opt-in `--exclude-nonspecific` that labels those cells `-1` and re-runs clustering without them.
- **Export to OMERO** (Future section): `report.html` is a natural FileAnnotation to attach to the source image (option 1).

### Implementation plan
1. `report.py`: loading (table, centroids, metadata, background level) and the computations (robust z, cluster profiles, automatic names, non-specific flags, community composition), with no plotting. Check: the numbers match `overview_analysis.py` on the COMET whole-cell dataset it was written for (same cluster sizes, top markers and debris percentages; set N so its definition matches).
2. Figures, one function per figure, each a pure function of the computed tables.
3. HTML assembly, plus the header from `summary.json` / `uns`.
4. CLI: `--no-report` in `cli.py`, and the standalone `python -m cellsurvey.report` entry point.
5. Docs: this section, the flag table, `docs/outputs.md` (and a sample screenshot), `README.md`.

Check: `test_crop` (local TIFF, Stardist; with and without `--detect-blobs`), a resume of `omero_602778` (whole-slide scale and timing), a nuclei-only run (`--cell-expansion 0`), and the standalone command on an older `_seg.zarr`.

### Decisions (2026-11-21)
- **Alongside the existing plots**: `cluster_intensity_heatmap.png`, `cell_density.png` and the other PNGs stay as they are. The report is an addition.
- **Fixed thresholds for now**: positive at z > 2; non-specific if z > 2 in at least `max(3, ceil(25%))` of the markers; clusters flagged if more than 50% non-specific; markers suspect above 25% and dominated above 50%. These are module constants, not flags, and the report header states them.
- **A spots section**, as described above.
- **On by default**, with `--no-report` to skip.

### Open: the intensity transform (needs investigating before step 1 is final)
The example uses `arcsinh(x)` (cofactor 1), then median/IQR. Note that for x ≫ 1, `arcsinh(x) ≈ log(2x)`, and median/IQR removes any constant shift, so the result is then independent of per-channel gain (exposure). Inputs are always integer (Sopa casts float TIFFs to `uint8`, and the OMERO reader rejects floats). So the question is how the cofactor compares with each channel's background:
- for `uint16` with background in the hundreds, cofactor 1 is effectively `log`;
- for `uint8`, or dim channels near 0, the near-linear part of arcsinh (|x| ≲ cofactor) matters, and low-end noise may be stretched.

To investigate, on `test_crop`, the COMET whole-cell dataset and `omero_602778`: compare plain z (what k-means uses), `log1p`, arcsinh with cofactor 1 and 5 (the cytometry convention), and a per-channel cofactor (e.g. a low percentile of that channel, as an estimate of background). For each, look at:
- histogram shape (does a second peak become visible?);
- the per-marker share above z = 2;
- the non-specific cell count and the cluster flags;
- how stable these are when the image is converted to `uint8`.

The outcome should also inform the clustering features (Known issues: k-means near-empty clusters) and the Normalisation plan, so that the report and clustering eventually use the same transform. Until this is decided, use the example's transform, so step 1 can be checked against `overview_analysis.py`.

## Future: per-cell spatial analysis of protein and spot distribution

**Status**: idea only (2026-11-19); nothing planned in detail. The aim is to describe *where* signal sits within each cell, not just how much (today's per-cell mean intensity and spot count).

**Candidate per-cell features**:
- **Spots**:
  - nuclear vs cytoplasmic fraction (nucleus polygons already exist and are linked through `nucleus_id`);
  - radial position (distance to the nucleus boundary or centroid, normalised by cell size);
  - clustering within the cell (nearest-neighbour distances between spots);
  - polarity (spot centroid offset from the cell centroid).
- **Protein**:
  - compartment intensities and ratios: nucleus; cytoplasm (cell minus nucleus); membrane ring (a thin band inside the cell boundary);
  - intensity-weighted centroid offset (polarisation);
  - radial intensity profile;
  - robust positivity scores (see PIPEX: 90th percentile, fraction of pixels above Otsu).
- **Protein × spots**: marker intensity at spot locations vs the cell mean, i.e. spot–protein colocalisation per cell.

**Keeping it computationally feasible**:
- Whole slides have ~10⁵–10⁶ cells × 20+ channels. Avoid Python loops over cells and polygons, and pixel-level statistics that scale with cell count × channels × iterations, such as per-cell Ripley's K or Haralick texture on every channel.
- **Work per tile, vectorised**: reuse Sopa's image patches and read 1024²-chunked tiles from `<prefix>.zarr`. Rasterise cell and nucleus labels per tile (plus a membrane ring by erosion), then compute all cells' compartment statistics at once with `scipy.ndimage`/`skimage.measure.regionprops_table(intensity_image=…)`. Parallelise over tiles with dask. Cells crossing a tile edge belong to the tile containing their centroid, with tiles padded by the maximum cell radius.
- **Spots are cheap**: assignment is already a vectorised spatial join, and distances to nucleus/cell boundaries can use shapely 2's vectorised functions.
- **Limit the scope**: compute only for chosen channels (a flag) and store summary statistics only, as new `obs` columns or an `obsm` block. This could be an opt-in stage, or a separate read-only command like the stability sweep.
- **Fits the aim of owning aggregation** instead of relying on `sopa.aggregate`, which reduces Sopa coupling.

## Planned: OMERO input

**Status** (2026-11-19): **implemented on `feature/omero`** (`cellsurvey/omero_source.py`; plan steps 1–3 + reader). Tested locally against a fake server. The **first real-server test** on Image 602778 (below) got the metadata right but failed on pixel reads (one shared `RawPixelsStore`). The per-thread-store fix has so far only been checked against the fake server, so a **real-server re-test is pending** (also the first timed read of an uncached slide). User docs: `docs/omero.md`.

**Goal**: read the input image directly from an OMERO server instead of requiring a local OME-TIFF.

### Why it fits the current design
Stage 1 is the only place the source image is read. It loads the image and immediately writes `<prefix>.zarr`, and every later stage (including `--resume-from`) works from that local copy. So OMERO support only needs to replace image loading: produce a `(C, Y, X)` dask array, the channel names and the pixel size (µm/px), and build the SpatialData image from them. Nothing downstream needs to know the source.

### Proposed structure
- Add an **image-reader seam**: one function per source (local TIFF via `sopa.io.ome_tif`/`BioImage` as now, OMERO, and later OME-Zarr), each returning `(dask array, channel_names, pixel_size_um)`. This also covers the existing *broaden input support* aim (see the input-format gotcha).
- **Branch**: create a new feature branch from `feature/instanseg`, not `main`. The reader touches the same stage-1 code in `cli.py` that `feature/instanseg` changed (pixel-size handling). `feature/instanseg` will fast-forward into `main`, so the new branch will sit cleanly on `main` afterwards.

### Options
1. **Direct read from the server** via [omero-py](https://pypi.org/project/omero-py/) BlitzGateway or [ezomero](https://pypi.org/project/ezomero/), fetching tiles into dask on demand (e.g. `-i omero://host/Image:1234`).
   - Simplest for users, with no intermediate copy.
   - The main risk is the dependency: omero-py needs **ZeroC Ice 3.6**. There is no usable PyPI wheel (ezomero recommends Glencoe's wheels), and omero-py's Python 3.12 support was only "preliminary". conda-forge carries `omero-py`/`zeroc-ice`, which is the best route for pixi. It must be verified to install alongside TF (default env) and Torch (instanseg env). Put it in its own pixi feature.
   - Tile streaming a whole multiplexed slide may be slow.
   - Needs non-interactive auth on the HPC: reuse an `omero login` session or environment variables, **never a password on the command line**.
2. **Export, then read**: `omero zarr export Image:1234` ([omero-cli-zarr](https://pypi.org/project/omero-cli-zarr/)) → local `.ome.zarr` → CellSurvey reads OME-Zarr (SpatialData is Zarr-native; [bioio-ome-zarr](https://pypi.org/project/bioio-ome-zarr/) is an alternative).
   - The pipeline needs no OMERO dependencies (those stay in a small separate env), and OME-Zarr input is useful on its own.
   - Costs two steps and a temporary full-size copy on disk.
3. **Read server-published OME-Zarr over HTTP** (as IDR does), if the facility's OMERO exposes it. This would be the best of both, but it depends entirely on server setup. The OME-Zarr reader from option 2 would mostly cover it.

**Decision (2026-11)**: go with **option 1 — read pixels and metadata through the OMERO API**. OMERO has already decoded every format with Bio-Formats (VSI, CZI, ND2, multi-series files…) and holds the channel names and physical pixel sizes. Using it means CellSurvey gets every format for free instead of building its own readers. Option 2 (export → OME-Zarr) and reading in-place originals from Nemo are deprioritised; they remain fallbacks if API throughput proves too slow. Option 3 is worth adopting if the server ever publishes OME-Zarr.

Implementation notes for option 1 (not yet built):
- **Reader**: `omero.gateway.BlitzGateway` → `conn.getObject("Image", id)`. Pixels come from `image.getPrimaryPixels()` via **tiled reads** (`getTiles` / a `RawPixelsStore` per worker), wrapped as a lazy `dask.array` of `(C, Y, X)` at full resolution (Z=0, T=0). Never load a whole slide into memory (e.g. `ezomero.get_image` returns a full numpy array, which is unsuitable for whole slides).
- **Metadata**: channel names from `image.getChannelLabels()`; pixel size from `image.getPixelSizeX(units="MICROMETER")`/`getPixelSizeY` (feeds `pixel_size_um`; still stop at startup if absent and no `--pixel-size`); dtype from `getPixelsType()`.
- **CLI**: e.g. `-i omero://omero-prod.camp.thecrick.org/Image:1234` (or `--omero-image-id` + `--omero-host`), No group flag is needed: the reader looks the image up across all groups and then switches to its group (see Findings). Auth: reuse an existing `omero login` session or environment variables; never a password on the command line.
- **Environment**: put `omero-py` (conda-forge) in its own pixi feature and test co-installing it into both the `default` (TF) and `instanseg` (Torch) envs. Import it lazily, like TF/Torch, so `cellsurvey` imports without it.
- **Read once**: stage 1 writes `<prefix>.zarr` straight away, so OMERO is only contacted at the start. `--resume-from` never touches OMERO. Record the source (`host`, `Image ID`) in the Zarr/table `.uns` for provenance.

### Findings so far (Crick, 2026-11)
- **Server**: `omero-prod.camp.thecrick.org`, port **4064** (standard Ice port), reachable from the HPC. Sessions have a **100 min idle timeout**. A user's default group may differ from where images live (e.g. `STP-CALM`).
- **Client install works**: `omero-py` (+ `zeroc-ice`) from **conda-forge** installs and logs in under **Python 3.12** in a standalone pixi project (`pixi add "python=3.12" omero-py`). On first import it prints many harmless `SyntaxWarning: invalid escape sequence` lines from the generated Ice stubs' docstrings; silence with `PYTHONWARNINGS="ignore::SyntaxWarning"`. Still untested: co-installing with TF (default env) or Torch (instanseg env).
- **omero-py co-installs with both segmenter envs** (branch `feature/omero`, `[feature.omero]` in `pixi.toml`, added to `default`/`spatial`/`instanseg`): solves for linux-64 + win-64; installed and imported alongside TF (default) and Torch (instanseg) on **both win-64 and linux-64 (HPC)**, numpy 2.4.6. The HPC test ran on a **CPU node only**: GPU use with omero-py co-installed (TF `and-cuda`, Torch CUDA, `tf.config.list_physical_devices('GPU')` / `torch.cuda.is_available()`) is still to confirm on a GPU node. **Pitfall**: pixi locks conda packages *before* the PyPI solve, so omero-py's conda dependencies must stay compatible with the PyPI stack. Unpinned, conda chose `python-dateutil 4.0.0` (via omero-py → `matplotlib-base`), but every `botocore`/`aiobotocore` (s3fs ← fsspec[s3] ← ome-zarr ← spatialdata) needs `<3`. It also chose `numpy 2.5.x`, but numba/scipy need `<2.5`. The PyPI resolver then backtracked to `s3fs==0.3.5` and failed with the misleading `SafeConfigParser` build error. Fixed by pinning `python-dateutil >=2.8,<3` and `numpy >=2,<2.5` as conda deps in the omero feature, plus a PyPI floor `s3fs>=2024.1` so any future conflict is reported directly. omero-py's other conda pin, `urllib3<2`, is fine (nothing in the stack needs urllib3 2.x).
- **`omero-cli-zarr` from PyPI failed to resolve**: the resolver backtracked to `s3fs==0.3.5` (a 2019 sdist whose `versioneer` uses `configparser.SafeConfigParser`, removed in Python 3.12). Workaround (untested): install `s3fs` from conda-forge first (`pixi add s3fs`), or pin `--pypi "s3fs>=2023.1"`.
- **Crick uses in-place import from Nemo** ([Crick OMERO guide → Import](https://franciscrickinstitute.github.io/crick-omeroguide/import.html)): images are copied to `proj-omero/<id>` (`data2.thecrick.org`) and imported in place, so OMERO owns the files but they stay on Nemo, and users keep read access. **This may make an OMERO client unnecessary for Crick data**: CellSurvey could read the original files directly from `/nemo/...` on the HPC. OMERO would then only be needed to map an Image ID to its file path, or for images that weren't imported in place, or for external servers (e.g. IDR, which also offers bulk downloads). To check: `omero fs sets --with-transfer` (in-place transfer type and stored paths) and whether those paths are readable under `/nemo` from compute nodes.
- **Server holds many non-TIFF, multi-series files**: `omero hql "select i.id, i.name from Image i"` works over the saved session and shows e.g. Olympus `.vsi` files, where one file maps to several Images (overview + each scan, IDs 53–55). Consequences: (1) reading originals directly from Nemo only works for formats CellSurvey can read (OME-TIFF, like the COMET slides); `.vsi`/`.czi`/`.nd2` would need the separate *broaden input formats* work. (2) An Image ID can be **one series within a file**, so ID → file path is not enough (the series index is needed). (3) Reading pixels through the OMERO API avoids both, because the server has already decoded every format with Bio-Formats. This favours option 1 for anything beyond in-place OME-TIFFs. Image names are client-side import paths (often Windows `\` paths), not reliable file paths; use `omero fs sets --with-transfer` for stored paths.
- **Queries are group-scoped**: `omero hql` only searches the session's current group unless `--all` is passed (which sets `omero.group=-1`). An image in another group returns **no rows** rather than an error, while `omero obj get Image:<id>` still finds it. The reader must therefore look images up across all groups (`conn.SERVICE_OPTS.setOmeroGroup('-1')`, client-side only, for `getObject`), then read pixels with a call context set to the image's own group (`{'omero.group': str(image.getDetails().getGroup().getId())}` passed to `RawPixelsStore.setPixelsId` etc.). **Don't use `conn.setGroupForSession`**: with a reused `omero login` session it switches the group of the user's own CLI session too (the 2026-11 throughput test script did this; `omero sessions group <name>` switches it back). Users then only pass an Image ID; no `--omero-group` is needed. See the Implementation plan.
- **Example image metadata is usable** (Image 602778, an Akoya `.qptiff`, `[resolution #1]`, `series=0`): Pixels `44160 × 51840 × 8` channels, `uint16`, physical size **0.4987 µm/px with unit `MICROMETER`**. Pixel size is populated, so `--pixel-size` isn't needed for it. That's ≈ **36.6 GB** uncompressed, so tile-read throughput is the main risk: ~12 min at 50 MB/s, ~1 h at 10 MB/s. Parallel connections and tile size need testing. It is still a one-off cost per image because of the read-once design.
  - **Channel names are real marker names** (`select lc.name from Channel c join c.logicalChannel lc where c.pixels.image.id = …`): `DAPI, CNPASE, PSYN, GFAP, PTAU, IBA1, MAP2, AF`. They can be used directly as intensity column names; DAPI is channel 0 (works with `--seg-channels 0` and the planned DAPI normalisation's name lookup). `AF` (autofluorescence) is included by InstanSeg's all-channels default and is aggregated like a marker, so consider excluding it from segmentation/clustering.
  - **The pyramid is not split into separate Images.** Fileset 244467 holds 4 Images: series 0 `[resolution #1]` (the full-resolution scan, 602778; pyramid levels are internal to it), plus `[thumbnail]`, `[Overview]` (macro photo) and `[Label]` (slide-label photo). So one slide = one useful Image ID. The reader should check the fileset and **stop if the given Image is not the largest in its fileset**, suggesting the right ID ("did you mean Image:602778?"), so a user can't accidentally segment the Overview/Label photo.
  - The image lives in group **Lab-Gandhi** while the session's group was STP-CALM, which confirms the cross-group lookup is required.
- **Throughput is fine — option 1 is viable** (tile-read test on Image 602778 from an HPC compute node, 30 s per config, random channel/position tiles from the full-resolution level, one `RawPixelsStore` per worker on a shared session):

  | workers | tile | MB/s | median ms/tile | full 36.6 GB image |
  |---|---|---|---|---|
  | 1 | 512 | 11.5* | 36 | 53 min |
  | 4 | 512 | 62.0 | 23 | 9.8 min |
  | 8 | 512 | 196.1 | 14 | 3.1 min |
  | 1 | 1024 | 49.5 | 40 | 12.3 min |
  | 4 | 1024 | 180.9 | 42 | 3.4 min |
  | 8 | 1024 | 314.0 | 50 | 1.9 min |
  | 1 | 2048 | 58.2 | 146 | 10.5 min |
  | 4 | 2048 | 219.5 | 151 | 2.8 min |
  | 8 | 2048 | 319.9 | 201 | 1.9 min |

  \*First config run, likely connection warm-up. **Take-aways**: throughput scales with parallel readers, which is what matters most; ~300+ MB/s at 8 workers is likely near the network/server limit (16 workers untested). Use **1024 px tiles** (a multiple of the server's native 512 px tile; same throughput as 2048 at a quarter of the latency and memory per request) and **~8 parallel readers**. That's ~2 min for a 36.6 GB slide, negligible next to segmentation. Caveat: configs ran back-to-back on one file (~40 GB of reads), so later ones may have been helped by the server's file cache; a cold first read may be slower. Re-time a full read of an untouched slide once the reader exists.
- **Pyramid levels**: `RawPixelsStore` numbers resolution levels **smallest = 0 → largest = n−1** (602778: 6 levels, full resolution = level 5); `getResolutionDescriptions()` lists them largest first. The reader must `setResolutionLevel(n-1)` and should verify it by checking the level size equals `sizeX × sizeY`. Non-pyramidal images have a single level.
- **Design note**: because stage 1 copies the image to `<prefix>.zarr` immediately, OMERO is only contacted at the start of a run. The 100 min idle timeout only matters if the session is created long before the job starts (e.g. queue wait). Log in inside the job or refresh the session.

### Open questions (to test on the HPC / with the facility)
- ~~Does `omero-py` install on Python 3.12, alongside TF and Torch?~~ Yes, from conda-forge, in both envs on win-64 and linux-64 (CPU node), with the conda pins described in Findings. GPU confirmed on an A100 node (driver 580, 2026-11-18) with omero-py co-installed: Stardist runs on the GPU end to end (after the TF CUDA preload fix, see Key gotchas), and InstanSeg runs on the GPU end to end (CUDA 13.2 Torch wheels; on `test_crop`, 1663 of 1828 cells linked to a nucleus, 2 nuclei unlinked).
- **Auth from compute nodes**: login and session reuse work interactively on the HPC. Is the saved session in `~/.omero` visible from a compute node? What is the cleanest way to log in non-interactively inside a batch job (env vars vs. session key)?
- ~~**Throughput**~~: resolved — ~300 MB/s with 8 readers × 1024 px tiles (~2 min per 36.6 GB slide; see Findings). Remaining: confirm a cold full read, and whether 16 readers helps.
- ~~**Metadata**: are channel names and physical pixel sizes populated?~~ Yes on a PhenoImager `.qptiff` (see Findings). Still worth one check on a **COMET** OME-TIFF in OMERO before relying on it.
- (Lower priority now) Does the server publish OME-Zarr over HTTP (option 3)?

### Implementation plan (2026-11-18)

**Goal**: accept `-i omero://host/Image:<id>` wherever a TIFF path is accepted. The image is read once, in stage 1, into `<prefix>.zarr`. Everything after that, including `--resume-from`, works only from the local Zarr.

Already checked: `sopa.io.ome_tif` is a thin wrapper (dask array → `Image2DModel.parse` with `chunks=(1, 1024, 1024)`, `scale_factors=[2, 2, 2, 2]`, plus `attrs[SopaAttrs.CELL_SEGMENTATION]`, from `sopa.constants`). `SpatialData.attrs` are saved to the Zarr and kept through `read_zarr` and `write_element` (spatialdata 0.9.0). The installed omero-py API was inspected. `_PixelsWrapper.getTiles` opens a new `RawPixelsStore` on every call and decodes with `struct`, so it is too slow for whole slides.

1. **Reader layer, new `cellsurvey/image_source.py`** — **done (commit 1)**. Local TIFFs only, no behaviour change: on synthetic OME-TIFFs (named `uint16`, `float32` → Sopa's `uint8` cast, unnamed channels), the new path writes Zarrs identical to `sopa.io.ome_tif(...).write(...)` (5 scales, chunks, dtype, channel coords, data, transforms); the only addition is `attrs['cellsurvey']`, which survives `read_zarr`. The blob-detection branch still reads the TIFF itself, and resume still calls `get_pixel_size(imagepath)`; both change in commit 2. Details:
   - `ImageSource`: a lazy `(c, y, x)` dask array, channel names, `pixel_size_um`, provenance, and a reader-thread count.
   - `open_image_source(path)`: `omero://` paths go to the OMERO reader. Anything else goes through `sopa.io.ome_tif(..., as_image=True)` + `get_pixel_size`, so local TIFFs read exactly as before.
   - `to_spatialdata(pixel_size)`: copies `sopa.io.ome_tif` (same chunks, pyramid and Sopa attr) and adds `attrs['cellsurvey'] = {pixel_size_um, source}`.
   - `write()`: wraps `sdata.write` in `dask.config.set(scheduler='threads', num_workers=…)`.
2. **OMERO reader, new `cellsurvey/omero_source.py`** (omero-py imported only here, lazily) — **done (commit 3), real-server test pending**. As planned, plus: `ImageSource.closer`/`close()` (cli closes the source right after the stage-1 write, and before any `parser.error`); `open_image_source(path, readers=)`; `--omero-readers` (default 8); on resume with an old Zarr lacking a stored pixel size, the source is opened for metadata only and closed. Local tests (fake server mimicking omero-py objects):
   - URL parsing;
   - login order and the "run `omero login`" error, with a fake `BlitzGateway`;
   - cross-group lookup (`omero.group=-1`) and pixel-store context with the image's group;
   - full-resolution level `n-1`, with the size check;
   - big-endian decoding: pixel-identical Zarr after `to_spatialdata` + `write`;
   - a failed tile retried with a fresh store;
   - one store per thread, all closed, and a reused session closed with `hard=False`;
   - fileset stop (`[Label]`, 0.05% of the largest) and warn (a second scan);
   - `import cellsurvey` doesn't import omero.

   **First real-server test (2026-11-19, Image 602778, HPC)**:
   - **Metadata worked**: login via the saved `omero login` session; cross-group lookup (image in Lab-Gandhi, group 954); 44160 × 51840 px, 8 channels `DAPI, CNPASE, PSYN, GFAP, PTAU, IBA1, MAP2, AF`, uint16, 0.4987 µm/px.
   - **Pixel reads failed**, with `Hibernate session is not re-entrant. Either you have two threads operating on the same stateful object`. Cause: **`BlitzGateway.createRawPixelsStore()` does not create a store**; it returns the connection's single cached proxy (`self._proxies['rawpixels']`). So every dask thread shared one stateful service.
   - **Fix**: create one store per thread with `conn.c.sf.createRawPixelsStore()` (a new service each call), passing the context as a plain `dict`. The fake server now mimics both behaviours (cached `conn.createRawPixelsStore()`, and stores that raise when two threads use them at once), and the fixed reader passes: pixel-identical, one store per thread, all closed.

   **Pyramid risk resolved**: with a counting fake source, writing all 5 pyramid levels read each full-resolution tile exactly **once**, so no two-step write is needed. Details as planned:
   - **URL**: `omero://host[:port]/Image:123` or `omero://host/123`; port defaults to 4064.
   - **Login**, tried in order, never with a password on the command line:
     1. the `OMERO_SESSION_KEY` environment variable;
     2. the session saved by `omero login`, if it is for the same host (`omero.util.sessions.SessionsStore().get_current()`);
     3. the `OMERO_USER` and `OMERO_PASSWORD` environment variables.

     Otherwise stop and tell the user to run `omero login`. Keep-alive is switched on during the read.
   - **Cross-group lookup without changing the session's group**: set `conn.SERVICE_OPTS.setOmeroGroup('-1')` (client-side only) for `getObject`, then pass a context with the image's own group to `RawPixelsStore.setPixelsId`. This **replaces** the earlier `setGroupForSession` note: with a reused `omero login` session, that call would also switch the group for the user's own CLI session.
   - **Fileset check**: catches thumbnail, overview and label images (see decision c).
   - **Metadata**: `getChannelLabels()` (blank labels become their index); pixel size from `getPixelSizeX/Y(units="MICROMETER")` → `LengthI.getValue()`, with a warning for non-square pixels; pixel type from `getPrimaryPixels().getPixelsType().value` via `OMERO_NUMPY_TYPES`. If the image has more than one Z or T plane, warn and use Z=0, T=0.
   - **Pixels**:
     - one `RawPixelsStore` per dask thread (`threading.local`), set to full resolution (`setResolutionLevel(n-1)`, only if `getResolutionLevels() > 1`); that level is checked against `sizeX × sizeY` once at the start (`getResolutionDescriptions()` lists largest first);
     - 1024² tiles via `getTile(z, c, t, x, y, w, h)`, decoded with `np.frombuffer` (big-endian), wrapped with `dask.array.from_array` (fixed `name`, `meta`, `lock=False`);
     - up to 3 retries per tile with a fresh store;
     - 8 readers by default.
   - **Closing**: close all pixel stores, then `conn.close(hard=False)` if we reused a session (so `omero login` stays valid), or `hard=True` if we logged in ourselves.
3. **Stage 1 in `cli.py`** — **done (commit 2)**. Tested on a synthetic 4-channel OME-TIFF with Gaussian spots: the old flow (`imread` + `BioImage` names + `detect_blobs_tiled` on the TIFF) and the new flow (write the Zarr, `detect_spots` from it, `write_element`, read back) give identical spots (118; same x/y/sigma/gene). The stored pixel size survives `write_element`, and an out-of-range channel index gives a clear error. HPC-checked (2026-11-19, `test_crop`, A100): a new `--detect-blobs` run writes the Zarr, detects 12,518 spots from it, adds them with `write_element`, and segments on the GPU; spot assignment and the GeoJSON spot export work. Resuming (`_seg.zarr` removed) reads the stored pixel size and reuses the spots ("already has spots; skipping blob detection"). The resume-with-missing-spots path is untested. Details:
   - **New runs**: open the source once at startup. The pixel size comes from `--pixel-size`, or else the source. The startup check that stops the run without a pixel size stays.
   - **Resumed runs**: do not open the source. The pixel size comes from `--pixel-size`, or `attrs['cellsurvey']` in the resume Zarr, or, for older Zarrs only, the source's metadata. The existing `_seg.zarr` mismatch check stays.
   - **Write `<prefix>.zarr` first, then detect blobs** from the local copy (`sopa.utils.get_spatial_image`) and add spots with `write_element('spots')`, so the source is read once (decision d). Blob channel names come from `source.channel_names`. Channel indices are checked before any reading.
   - `-i` help mentions `omero://`; drop the `BioImage`/`imread` imports from `cli.py`.
4. **Docs**:
   - this section: mark the plan implemented, record test results;
   - `AGENTS.md` stage-1 description and flag table;
   - `CLAUDE.md` input bullet;
   - `docs/`: a "Reading from OMERO" section (URL form, login, which image ID, HPC batch-job session timeout), plus `getting-started`, `parameters`, `faq` and `installation`.
5. **Testing**:
   - **Locally, with a fake OMERO array that counts tile reads (no server)**: checks tile placement, byte order, and that each tile is read once. **Risk to check**: writing the 4-level pyramid may read the full-resolution tiles again for each lower level, which is cheap for a TIFF but 5× the reads from OMERO. If so, write full resolution first and build the pyramid from the local copy.
   - **On the HPC**:
     1. a local-TIFF run on `test_crop`, to confirm nothing changed;
     2. a small OMERO image (login, cross-group lookup, channels);
     3. Image 602778 as a full read, timing a first read the server hasn't cached;
     4. resuming from it, which should not contact OMERO;
     5. the fileset check, using the Overview or Label image ID.
6. **Commits**:
   1. the reader layer (local path, no behaviour change);
   2. the stage-1 reorder;
   3. the OMERO reader;
   4. docs.

**Decisions** (2026-11-18):
- a. **CLI**: the URL goes in `-i` (`omero://host[:port]/Image:<id>` or `omero://host/<id>`), with one new flag, `--omero-readers` (default 8). No separate host/ID flags.
- b. **Login order**: `OMERO_SESSION_KEY` → the saved `omero login` session for the same host → `OMERO_USER`/`OMERO_PASSWORD`. Otherwise stop with "run `omero login`". Never a password on the command line. Close with `hard=False` for reused sessions.
- c. **Fileset check**: stop if the image is under 1% of the largest image's area in its fileset (thumbnail, label or overview photo), suggesting the largest ID ("did you mean Image:602778?"). If it is smaller but not tiny (e.g. another scan in a multi-scan `.vsi`), warn, naming the largest, and continue.
- d. **Reorder blob detection for all inputs**: write `<prefix>.zarr` (image) first, then detect blobs from the local Zarr and add spots with `write_element('spots')`. The old order (detect, then a single write) dates from the March prototype (`78fbcc2`). It was built for Spotiflow, but Spotiflow (`Spotiflow.predict(img=...)`) takes an array, and that commit passed it a `dask_image.imread` dask array, so a Zarr-backed dask array works the same way. No other technical reason was found in the history, and `write_element` into an existing Zarr was verified on spatialdata 0.9.0. Structure detection as "takes a `(y, x)` array, returns spots", so the LoG detector and a future Spotiflow detector are interchangeable. Side effect: a crash during detection now leaves `<prefix>.zarr` without spots. So on `--resume-from` with `--detect-blobs`, if the Zarr has no `spots`, run detection from the Zarr. Needs an HPC regression check with `--detect-blobs` on a local TIFF.
- e. **Float pixel types**: reject with a clear error naming the type for now. Multiplexed slides are integer (the example qptiff is `uint16`). Revisit if a float image turns up.
- f. A small OMERO test image ID for the first test: *to do*.

### HPC test checklist (before merging `feature/omero` → `main`)

**Status**: **all tests 1–6 passed (2026-11-20)**. Test 6 ran interactively on a GPU node rather than as an `sbatch` job, so batch-job login is the only part not exercised. `feature/omero` is ready to fast-forward into `main`. History: on the first pass, test 1 failed only because `$OMERO` was empty in that shell, so the path became `/Image:602780` (now caught with a clear "Input image not found … must start with omero://" error). Test 3 exposed an old bug (zero blobs → float `gene` column → `sopa.aggregate` crash), now fixed (see stage 2). Second pass, after the fixes: all of 1–5 passed.

Minor, harmless log noise seen (possible follow-ups, not blockers):
- TF startup messages at the end of a resumed (no-segmentation) run come from `umap-learn` importing TensorFlow for its optional ParametricUMAP. Not our import; default env only.
- `python -m cellsurvey.stability` warns `'cellsurvey.stability' found in sys.modules…` because `cellsurvey/__init__.py` imports `stability` eagerly.
- On the first OMERO import in a fresh env, Ice's generated stubs print many `SyntaxWarning`s; they stop once Python has cached the compiled files. A `warnings` filter before the omero-py import in `omero_source` would hide them. Tick items off here as they pass and record anything surprising. `feature/omero` contains all of `feature/instanseg`, so passing tests 1–6 covers both branches; then `feature/omero` can fast-forward into `main` (no separate `feature/instanseg` merge needed). Paths in `<…>` are placeholders; use a small local TIFF (e.g. `test_crop`) for tests 2–5.

**Setup** (once):
```bash
cd <repo> && git pull                      # feature/omero
pixi install --locked && pixi install --locked -e instanseg
pixi run omero sessions list               # check the current group; the 2026-11 throughput script may have switched it:
pixi run omero sessions group STP-CALM     # or log in again: pixi run omero login -s omero-prod.camp.thecrick.org -u <user>
OMERO=omero://omero-prod.camp.thecrick.org   # set it in every new shell/job, or write the URL out in full
```

- [x] **1. Fileset check** (seconds; no pixels read) — **passed 2026-11-20** (interactive): logged in via the saved session, then stopped with "…`[Overview]` (2046 x 3950 px) is a small image… Did you mean omero://omero-prod.camp.thecrick.org/Image:602778?":
  `pixi run -e instanseg python run.py -i $OMERO/Image:602780 -o out/t_overview -p out/t_overview_plots --segmenter instanseg`
  ✅ Logs in via the saved session and finds the image in Lab-Gandhi, then stops with "…Did you mean omero://…/Image:602778?". No `out/t_overview.zarr` is created.
- [x] **2. Local TIFF, plain** (checks the default Stardist path for regressions; GPU node) — **passed 2026-11-20** on `test_data/test_crop.ome.tiff`: `Found 1 GPU(s)` (A100), 2109 cells, 6123 edges (≤ 3×), GeoJSON at `out/t_local.geojson`:
  `pixi run python run.py -i <test_crop.ome.tiff> -o out/t_local -p out/t_local_plots`
  ✅ Completes. The log shows `Found 1 GPU(s)` (the TF CUDA preload works). `out/t_local.zarr` has `attrs['cellsurvey']` (pixel size, source path). `summary.json` has `cluster_on: intensities`, and `n_edges` is ≤ ~3× the cell count.
- [x] **3. Local TIFF with `--detect-blobs`** (checks the reordered blob detection) — **passed 2026-11-20**. No spots + `both`: stopped straight after detection, before Stardist. No spots + `intensities`: warned and completed (2110 cells). Spots (`--channels 8,10 --thresholds 0.05,0.05 --cluster-on both`): correct log order; 5165 + 4888 spots, 98.4% assigned to cells; 21 intensity + 2 spot features = 23:
  `pixi run python run.py -i <test_crop.ome.tiff> -o out/t_blobs -p out/t_blobs_plots --detect-blobs --channels <i,j> --thresholds <a,b> --cluster-on both`
  ✅ Log order: "Saving Zarr", then "Finding blobs in channel …", then "Adding spots to out/t_blobs.zarr", then segmentation. The intensity and spot channels are listed separately, and the heatmap has `(spots)` columns. Any two channels will do; this only checks the plumbing.
  - **Also check the zero-spots path**: with thresholds that find nothing (e.g. `--thresholds 1,1` on `test_crop`), `--cluster-on both` must stop straight after detection with "No spots were detected in channels … Lower --thresholds, or use --cluster-on intensities". With `--cluster-on intensities` it should warn "no spots detected … continuing without spots" and complete. (The first pass found the bug this checks for.)
  - Then re-run with thresholds that find spots (try `0.05`–`0.1` on spot-like channels) for the main check above.
- [x] **4. Resume never opens the source** — **passed 2026-11-20**. Worked with `-i /nonexistent.tiff`, used the stored pixel size, skipped to clustering, and reproduced test 3's clusters, edges (6120) and communities (18) exactly:
  `pixi run python run.py -i /nonexistent.tiff -o out/t_blobs -p out/t_blobs_plots2 --resume-from out/t_blobs.zarr --detect-blobs --channels <i,j> --thresholds <a,b> --cluster-on both`
  ✅ Works despite the non-existent `-i`. Log shows `Physical pixel size (stored in out/t_blobs.zarr)`, "Segmented Zarr with table found, skipping to clustering", and no blob detection.
- [x] **5. Stability sweep uses the stored mode** — **passed 2026-11-20**: `Features: --cluster-on both (stored in the Zarr)`. Results on `test_crop` (2108 cells) are instructive. At resolution 0.1: 25 µm → 5963 edges, 51 communities; 50 µm → 6120 edges, 18; 100 µm → 6256 edges, 4. A 25 µm cutoff drops only ~160 edges against 50 µm but triples the community count (isolated cells become small communities); 50 µm looks like a reasonable default, but sweep per dataset:
  `pixi run python -m cellsurvey.stability --zarr out/t_blobs_seg.zarr --max-edge-distances 25,50,100`
  ✅ Prints `Features: --cluster-on both (stored in the Zarr)` and writes the CSVs.
- [x] **6. OMERO full run** — **passed 2026-11-19/20** (interactive GPU node, not `sbatch`); see Results below. (Originally specified as an `sbatch` job on a GPU node, covering the real-server re-test of the per-thread stores, the first read of an uncached slide, and batch-job login):
  `pixi run -e instanseg python run.py -i $OMERO/Image:602778 -o out/omero_602778 -p out/omero_602778_plots --segmenter instanseg`
  ✅ Log shows:
  - `Connected to OMERO … via saved `omero login` session` (batch-job login works);
  - the image in Lab-Gandhi, 8 named channels, 0.4987 µm/px;
  - `CUDA available`;
  - **no** "Hibernate session is not re-entrant" error.

  Record the time from "Saving Zarr to out/omero_602778.zarr" to "Done" (uncached read; ~2–5 min predicted; `date` before and after in the job script helps). InstanSeg afterwards takes ~1.5 h. Then spot-check visually: `out/omero_602778.zarr` in napari, or the GeoJSON in QuPath, against the image in OMERO. This catches tile placement or byte-order bugs the fake server could have missed.

**Results** (test 6, Image 602778, InstanSeg, A100; outputs in `out/omero_602778*`):
- **OMERO read**: logged in via the saved session; image found in Lab-Gandhi; 8 named channels; 0.4987 µm/px. **No "Hibernate session is not re-entrant" errors** (the per-thread `c.sf.createRawPixelsStore()` fix works on the real server).
- **Stage-1 copy of the uncached slide took ~5.6 min** for 36.6 GB (~110 MB/s), from Zarr file timestamps 21:23:59 → 21:29:38. That includes the 5-level pyramid, compression and writing to Nemo. The predicted range was 2–5 min.
- **Copy verified visually**:
  - a 1/16-scale render of DAPI/GFAP/MAP2 shows coherent tissue with no scrambled or missing tiles;
  - at a full-resolution junction of four 1024² tiles, nuclei run continuously across both edges (correct placement and byte order);
  - `attrs['cellsurvey']` records the pixel size and the OMERO source (host, image ID, name, group 954).

  MAP2 shows vertical banding, apparently in the source (acquisition or unmixing), not the copy.
- **Segmentation**: 621 patches.
  - Nuclei: 657,741 found → 622,853 after conflicts (~26 min).
  - Cells: 895,272 → 848,343 (~48 min including conflict resolution).
  - 621,844 cells linked to a nucleus (1,009 nuclei unlinked); **226,499 cells (27%) dropped with no nucleus**, against 9% on `test_crop`. That's plausible for brain (fragments of processes), but worth a visual check.
  - 4,317 cells were filtered as too dim, leaving 617,527 cells.
  - "No cell was returned" on ~14 patches is empty background, and harmless.
- **Clusters are biologically sensible**, mostly single-marker: MAP2 (neurons, the largest cells), DAPI-high small uniform cells, GFAP (astrocytes), CNPASE (oligodendrocytes), IBA1 (microglia), PSYN, PTAU (154 cells; rare cells or artefacts), and AF (autofluorescent; supports a way to exclude `AF` from clustering, see Known issues).
- **Timeline**:
  - stage-1 copy: 21:24 → 21:30;
  - nuclei: → 21:56;
  - cells: → 22:44;
  - GeoJSON: 22:52;
  - `_seg.zarr`: 22:57;
  - UMAP: 23:11;
  - `cell_density.png`: 23:16.

  About **1 h 53 min end to end**. Problems this whole-slide run exposed are in Known issues ("Whole-slide run findings").

## Future: export results back to OMERO

**Status**: idea only (2026-11-20). Goal: let people browse CellSurvey results (segmentation, k-means clusters, communities, per-cell measurements) in OMERO, next to the source image.

**Why it fits**: `<prefix>_seg.zarr` already records its OMERO source (`attrs['cellsurvey']['source']`: host, port, image ID, group). So a separate, read-only-on-the-Zarr command (e.g. `python -m cellsurvey.omero_export --zarr <output>_seg.zarr`) needs no other input, and can reuse `omero_source.connect` (same login order; and never `setGroupForSession`, see the OMERO notes).

**Options**, from cheapest to hardest:
1. **Attachments + key–value pairs** (easy):
   - `summary.json`, the plots and the run parameters attached to the source image as FileAnnotations;
   - a MapAnnotation with headline numbers: cells, clusters, communities, `--cluster-on`, segmenter, CellSurvey version or commit.

   Findable next to the image, but not spatially browsable.
2. **Label / cluster-map image** (medium; the best fit for browsing whole slides):
   - rasterise the results at the source resolution into a new multichannel image: cell labels, nucleus labels, k-means cluster ID, community ID;
   - write a pyramidal OME-TIFF and import it next to the source (e.g. a "CellSurvey results" dataset, with both images linked by a key–value or a comment);
   - OMERO.iviewer then shows it at full speed like any slide, with per-channel lookup tables, side by side with or overlaid on the original;
   - the file is large but very compressible (label rasters). Cluster and community maps fit uint8 or uint16; cell labels need uint32.
3. **ROIs** (hard at whole-slide scale): one polygon ROI per cell, with cluster and community, is the most native form. But iviewer and the API degrade badly with tens of thousands of ROIs per image, and whole slides have ~10⁵–10⁶ cells (617k for Image 602778). Workable variants:
   - ROIs for a subset only (a region or cluster of interest);
   - **one polygon per community** (the union or concave hull of its cells; ~50 for 602778), a compact, browsable map of tissue regions;
   - simplified polygons, uploaded in batches (`saveAndReturnArray`).
4. **Per-cell table** (medium; depends on the server): an OMERO.table (HDF5) on the image with cell ID, centroid, cluster, community and intensities. With ROIs, the OMERO.parade plugin can filter and colour cells by table values, if installed.
- **Related, no export needed**: QuPath's OMERO extension can open the image straight from OMERO, and users can load CellSurvey's GeoJSON onto it (the GeoJSON size issue, see Known issues, applies).

**To check first** (with the facility):
- **Write access**: can the user annotate or import into the image's group? That depends on the group's permission level (read-annotate vs read-only), and the image may belong to someone else.
- **Server features**: whether OMERO.tables and OMERO.parade are available, and the iviewer version.
- **ROI limit**: the recommended maximum ROI count per image.
- **Import route**: where derived images should go, and how imports are done from the HPC (`omero import` from the pixi env?).

**Suggested order**: 1 (attachments + key–values), then 2 (label/cluster-map image), then community-outline ROIs from 3. Per-cell ROIs or tables only if the server copes.

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

## Reference: MDV — Multi-Dimensional Viewer (Taylor-CCB-Group)

**Note (for future consideration):** [MDV](https://github.com/Taylor-CCB-Group/MDV) (Taylor-CCB-Group, University of Oxford; site [mdv.ndm.ox.ac.uk](https://mdv.ndm.ox.ac.uk/)) is an **open-source, web-based application for analyzing, annotating, and sharing multi-dimensional spatial/molecular datasets**. Not integrated into CellSurvey — retained as a potential **visualization/annotation alternative to QuPath** (CellSurvey's current GeoJSON export target).

### What it does
- Interactive, browser-based charts embeddable in a narrative "story": 2D/3D scatter, box/violin/heat/histogram/pie plots, annotation tools, and **interactive spatial-biology charts**.
- Handles large datasets — up to **~10 million data items** — via web workers, shared array buffers, and native arrays; lazy-loads raw bytes (≥4 GB RAM recommended).
- Multiple linked data tables (e.g. a gene-expression matrix referenced against cell/spatial coordinates); user-editable annotations; extensible data loaders.
- Runs in-browser with no installation for upload/viewing; Docker/devcontainer + Makefile for deployment.

### Relevance to CellSurvey
- **Overlaps the visualization stage, not the analysis stage**: CellSurvey ends at GeoJSON export → QuPath. MDV is a web-native alternative for exploring the *same* cell × marker / cell × community data CellSurvey already produces, with linked spatial + expression views rather than QuPath's image-centric paradigm.
- **Potential role**: a second export target (`export_to_qupath`-style writer emitting MDV's expected table/coordinate format) for interactive, shareable spatial results — useful for collaborators who want to explore communities/clusters without a pathology workstation.
- **Complementary to QuPath**: MDV is table/spatial-chart-centric and web-shareable; QuPath is whole-slide-image-centric with deep pathology tooling. They answer different "what does the tissue look like" vs. "what do the cell populations look like" questions.

### Caveats
- **No explicit segmentation or multiplexed-imaging ingest** documented in the README (no OME-TIFF/AnnData naming); data formats are not enumerated — the exact table/link schema for CellSurvey's output would need reverse-engineering or adapter work.
- GPL-3.0 (example data CC0); modest community traction (dozens of stars) but active (2000+ commits, CI, Playwright tests).
- Browser-based — not a headless/CLI pipeline stage; a separate deployment concern from CellSurvey's pixi/Docker compute env.

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

## Reference: PIPEX (CellProfiling)

**Note (for future consideration):** [PIPEX](https://github.com/CellProfiling/pipex) (CellProfiling; MIT) is a Python pipeline for **multiplexed imaging (built for CODEX)** that goes from images to per-cell tables, clusters, spatial analysis and viewer exports. It is a close sibling of CellSurvey (StarDist + expansion → per-cell measurement → k-means/Leiden → spatial analysis → QuPath/TissUUmaps export). It is not integrated and not worth adopting as a dependency. Several of its **per-cell measurement and annotation ideas** are worth borrowing.

### What it does
- **Input**: a flat folder of single-channel images named by marker (e.g. `DAPI.tif`), in any libvips-readable format, plus optional custom masks. It does *not* take a multichannel OME-TIFF. It is driven by a batch command file or a GUI.
- **Segmentation**: StarDist on DAPI (`nuclei_diameter`, `nuclei_definition`, `nuclei_closeness`, `nuclei_area_limit`) → non-overlapping expansion by `nuclei_expansion` **pixels** → optional **membrane-marker watershed** refinement (`membrane_diameter`, `membrane_compactness`, `membrane_keep`, `membrane_min_growth=0.05`). The refinement reshapes an expanded nucleus only if the watershed cell grows beyond the nucleus by more than that fraction. Watershed segments sharing a nucleus are merged.
- **Per-cell marker scores** (beyond the mean intensity in `cell_data.csv`): `_local_90` (90th-percentile intensity), `_ratio_pixels` (fraction of pixels above Otsu), `_otsu3`/`_triangle_score` (signed distance from a threshold), and `_gmm_prob` (2-component GMM posterior; NaN when `gmm_min_separation` < 0.5). `use_bin` clusters on these scores instead of raw intensities.
- **Image preprocessing**: threshold removal, top-hat background subtraction (`tophat_radius`), light-gradient reduction, tile balancing, stitch smoothing and exposure adjustment.
- **Normalisation / batch correction**: `log_norm`, `z_norm`, `minmax_norm`, `quantile_norm`, and ComBat (`batch_corr`). The README itself warns that ComBat may over-correct per-cell imaging data and produce negative values.
- **Clustering and annotation**: Leiden (`leiden_res=0.5`), k-means (`k_clusters=10`), and `k_estimation` (elbow/inertia/silhouette for k=1–19). Rule-based **cell-type annotation** comes from `cell_types.csv`: high/medium/low marker rules are matched against each cluster's ranked markers (`rank_filter`, `min_confidence`).
- **Spatial analysis**: **k-NN neighbourhood composition** (fraction of each cell type among the k nearest neighbours, `neigh_k_values=1,5,10`) and **DBSCAN per cell type** (scattered vs clustered, `neigh_density_threshold=0.05`).
- **Outputs**: label mask TIFF, `cell_data*.csv`, `anndata.h5ad` (+ a TissUUmaps variant), GeoJSON plus compressed **`.pbf`** for QuPath/TissUUmaps, a TissUUmaps web export, and Leica **laser-microdissection** XML.
- **Scale**: `max_res` downscaling and temporary swap space; quoted at ~16 GB RAM for 12k-pixel images. Large whole-slide images are not its target.

### Relevance to CellSurvey
- **Highest value — per-cell marker scores**: CellSurvey clusters only on `sopa.aggregate` mean intensity. With InstanSeg whole cells, a cell-wide mean dilutes membrane/punctate markers (CD3, CD8, PD-1), so `_local_90`, `_ratio_pixels` and `_gmm_prob` are more robust positivity measures and better clustering inputs. Computing them means owning aggregation instead of calling `sopa.aggregate`, which fits the **reduce-Sopa-coupling** aim. `_gmm_prob` is the same per-marker GMM idea as CELESTA's first step.
- **Normalisation menu** feeds the Planned Normalisation (DAPI) section and any future multi-sample mode (with the ComBat over-correction caveat).
- **Choosing k**: the cluster sweep (`run_cluster_sweep`) already records inertia. A **silhouette score on a random subsample** (O(n²) on all ~270k cells is infeasible) would cheaply mirror `k_estimation`.
- **Rule-based cluster naming** via a marker-rule CSV is a lightweight step towards cell-type labels (much simpler than CELESTA; no spatial prior).
- **k-NN neighbourhood composition** is the standard CODEX "cellular neighbourhoods" method and a complementary alternative to our Delaunay + Louvain communities. It is cheap with `scipy.spatial.cKDTree`.
- **Compressed `.pbf` export** addresses GeoJSON size at whole-slide scale (a 270k-cell `qupath_export.geojson`).
- **Low value**: the membrane-guided watershed refines StarDist + expansion, which InstanSeg already supersedes (and Stardist is being retired). The image-preprocessing steps (tile balancing, stitch smoothing) might matter only if stitched COMET images show tile artefacts. LMD export is niche.

### Caveats
- Small project (≈19 stars, no releases) with no paper cited in the README.
- Single-channel-folder input and pixel-unit parameters: not interoperable with CellSurvey's OME-TIFF/SpatialData/µm conventions. Borrow algorithms, not code paths.
- TensorFlow/StarDist stack. Its memory guidance targets ≤12k-pixel images, not whole slides.

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
│   ├── blob_detection.py             # detect_spots ((c,y,x) DataArray → spots PointsModel) → detect_blobs_tiled → detect_blobs_in_tile
│   ├── network_analysis.py           # run_network_analysis (µm → px) + compute_louvain_communities + delaunay_edges (scipy Delaunay + networkx Louvain, expr-similarity weights, summary.json)
│   ├── stability.py                  # run_stability_sweep (community) + run_cluster_sweep (k-means), read-only → CSV
│   ├── image_source.py               # ImageSource + open_image_source (stage-1 input: local TIFF or omero://) → to_spatialdata / write / close; stored_pixel_size
│   ├── omero_source.py               # OMERO reader: parse_omero_url, connect (login order), fileset check, per-thread RawPixelsStore tiles → lazy dask array
│   ├── segmentation.py               # run_stardist / run_instanseg, link_nuclei_to_cells, get_pixel_size, expand_nuclei (watershed)
│   ├── export.py                     # export_to_qupath
│   └── utils.py                      # shape-layer keys + resolve_shape_keys, remove_channel_suffix, cluster_data, assign_spots_to_cells, get_colors_for_communities
├── docs/                             # MkDocs Material site (ReadTheDocs-hosted)
│   ├── index.md                      # Landing / overview + pipeline diagram
│   ├── installation.md               # pixi + Docker setup
│   ├── getting-started.md            # First-run walkthrough
│   ├── pipeline.md                   # 10-stage walkthrough
│   ├── parameters.md                 # Full CLI reference + tuning guidance
│   ├── outputs.md                    # Every output file explained
│   ├── visualization.md              # Odon / QuPath / TissUUmaps / napari
│   └── faq.md                        # Question-phrased how-to + troubleshooting
├── Dockerfile                        # Ubuntu 24.04 + pixi + GPU-ready container (build arg PIXI_ENV)
├── pixi.toml                         # Pixi config: linux-64 + win-64; envs default (Stardist/TF) and instanseg (Torch)
├── mkdocs.yml                        # MkDocs + Material config
├── .readthedocs.yaml                 # ReadTheDocs build config
├── requirements.txt                  # Minimal pip requirements (sopa)
├── requirements-docs.txt             # Docs build deps (mkdocs-material)
├── README.md                         # User-facing installation and usage docs
├── LICENSE                           # License
├── .gitignore / .gitattributes       # Git config
└── .pixi/                            # Pixi environment directory (gitignored)
```
