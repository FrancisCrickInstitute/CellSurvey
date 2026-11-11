# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**`AGENTS.md` is the canonical, detailed project reference.** It covers the full CLI flag table, every pipeline stage, gotchas, planned work (stability sweep phases 2–3, normalisation, Cellpose) and notes on related tools. Read the relevant section before making non-trivial changes. When a change alters behaviour, flags, layer names or plans, update `AGENTS.md` too (and `README.md` / `docs/` if users will see the change).

## What this is

CellSurvey is a spatial-biology pipeline built on Sopa / SpatialData. It takes a multichannel OME-TIFF, optionally detects RNA spots, segments it (Stardist or InstanSeg), aggregates per-cell channel intensities, then runs k-means clustering, a Delaunay network with Louvain communities, and Scanpy neighbourhood/UMAP/Leiden plots. It exports a QuPath GeoJSON file. `run.py` is a thin shim; all logic lives in `cellsurvey/`, and `cli.main()` orchestrates the stages in order.

## Commands

Environments are managed by pixi (`linux-64` with GPU, `win-64` CPU-only; macOS is not supported). Each segmenter gets its own environment, so TensorFlow and PyTorch are never installed together:

```bash
pixi install                  # default env: Stardist + TensorFlow
pixi install -e instanseg     # instanseg env: PyTorch, no TF

pixi run python run.py -i <input.tiff> -o <output_prefix> -p <plot_dir> [--detect-blobs]
pixi run -e instanseg python run.py --segmenter instanseg -i <input.tiff> -o <output_prefix> -p <plot_dir>

# Resume from an existing Zarr (skips loading/spot detection; reuses <prefix>_seg.zarr if valid)
pixi run python run.py -i <input.tiff> -o <output_prefix> -p <plot_dir> --resume-from <output_prefix>.zarr

# Read-only parameter sweeps over an existing segmented Zarr (writes CSVs next to it)
pixi run python -m cellsurvey.stability --zarr <output>_seg.zarr --resolutions 0.1,0.05 --max-edge-distances 25,50,100   # µm
pixi run python -m cellsurvey.stability --zarr <output>_seg.zarr --mode clusters

# Docker (build arg PIXI_ENV selects the env)
docker build -t cellsurvey .
docker build --build-arg PIXI_ENV=instanseg -t cellsurvey-instanseg .

# Docs (MkDocs Material, hosted on ReadTheDocs)
pip install -r requirements-docs.txt && mkdocs serve
```

There are no tests, no linter config and no CI. To check a change, run the pipeline end to end on a small TIFF (`test-input/` is a local, untracked folder of sample data). Use `--resume-from` to iterate on later stages without re-segmenting.

## Architecture essentials

- **Outputs**: `<prefix>.zarr` is the raw image plus optional `spots`. `<prefix>_seg.zarr` adds shapes, the AnnData `table` and labels; it is written to a temp directory and then atomically renamed. Plots and `summary.json` go to `-p`, and the GeoJSON goes to `--geojson-path`.
- **Shape layer keys vary by segmenter, so never hard-code one.** Nuclei are `stardist_boundaries` (Stardist) or `nucleus_boundaries` (InstanSeg). Cells are `cell_boundaries`, which comes from the Stardist watershed expansion (`--cell-expansion`, in µm) or from InstanSeg's cell pass. Use `utils.resolve_shape_keys(sdata)` to get `(nuclei_key, cells_key)`; `cells_key` falls back to the nuclei layer for nuclei-only runs. Aggregation, clustering, the network, spot assignment, export and the sweeps all target the cell layer. Labels are mapped back onto nuclei through `nucleus_id`, and unlinked nuclei get `-1`.
- **Lazy framework imports**: TensorFlow is imported only inside `segmentation.run_stardist`, and torch/InstanSeg only inside `run_instanseg`. Importing `cellsurvey` must not pull in either one, or the other env breaks. The long-term plan is to drop Stardist/TF, so keep Stardist-specific code isolated there.
- **Sopa coupling should shrink over time**: Sopa has a single maintainer and is a risk. Prefer owning the logic over adding new Sopa dependencies. Do not use `sopa.shapes.expand.remove_overlap` or `sopa.aggregate(no_overlap=True)`: they empty cells on dense data (see the untracked `SOPA_REMOVE_OVERLAP_FIX.md`). The project uses its own watershed expansion instead.
- **Distances are in µm, coordinates in pixels.** `--max-edge-distance`, `--radius-min/max` and `--cell-expansion` are µm and are converted at the point of use (`value_um / pixel_size_um`). The pixel size comes from OME metadata (bioio) or `--pixel-size`; if neither exists the run stops at startup. It is stored in `table.uns['cellsurvey']['pixel_size_um']` (`utils.store_pixel_size`), and `utils.resolve_pixel_size` reads it back, so the stability sweep works without `--pixel-size` (except on older Zarrs). InstanSeg also uses it to rescale to its training resolution, so a wrong value can make it find few cells.
- **Input** must be a 3D/4D multichannel TIFF (read via `sopa.io.ome_tif` / dask_image). ND2, CZI and similar formats are not supported yet.

## Gotchas worth remembering

- `run.py` sets `TF_USE_LEGACY_KERAS=1` before any imports (Stardist needs `tf_keras` on TF ≥2.16). It also preloads the env's `libstdc++.so.6`. That preload is duplicated at the top of `cli.py`, so change both together.
- `matplotlib.use('Agg')` must come before `pyplot` imports, and Scanpy plots need `show=False`.
- AnnData `.X` may be sparse or dense; check with `hasattr(X, 'toarray')`.
- Pandas ArrowStringArray breaks Zarr writes. Keep the `pd.option_context('future.infer_string', False)` around `sopa.aggregate` and `obs.index.astype(str)`.
- Seeds are fixed at 42 (k-means, Louvain). Logging uses `print`.
