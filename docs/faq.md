# FAQ

Question-phrased how-to and troubleshooting, in one place. If you hit a problem
that isn't covered here, check the error message and the sections below.

## What are the inputs and outputs of CellSurvey?

- **Input**: one multichannel TIFF image (OME-TIFF recommended).
- **Output**: a segmented SpatialData Zarr store (`*_seg.zarr`), a QuPath
  GeoJSON file, `summary.json`, and a directory of plots.

See [Getting started](getting-started.md) and [Outputs](outputs.md) for the full
picture.

## What's the difference between clusters and communities?

- **Clusters** (`kmeans_cluster`) group cells by **expression** — cells with a
  similar marker profile, wherever they are in the tissue.
- **Communities** (`community`) group cells by **space** — cells that sit next to
  each other, regardless of marker similarity.

A single cell type therefore usually spans several communities (e.g. CD31⁺ blood
vessels are separate communities), while a cluster labels that cell type as one
group across the whole tissue. See [Pipeline](pipeline.md) for details.

## Should I use Stardist or InstanSeg?

- **Stardist** (`--segmenter stardist`, the default) segments nuclei on one
  channel, then approximates whole cells by expanding each nucleus by
  `--cell-expansion` µm.
- **InstanSeg** (`--segmenter instanseg`) segments nuclei *and* whole cells
  directly, using all channels by default. It runs in its own pixi environment:

  ```bash
  pixi install -e instanseg
  pixi run -e instanseg python run.py -i <image> -o <output> -p <plots> --segmenter instanseg
  ```

The project plans to migrate fully to InstanSeg and eventually drop Stardist.
See [Pipeline](pipeline.md) for how the two differ.

## How do I get whole cells instead of nuclei?

Both segmenters produce a `cell_boundaries` layer by default, and the
measurements, clusters, and communities are computed over it. With Stardist,
cell size is set by `--cell-expansion` (µm, default `5.0`); use
`--cell-expansion 0` to analyse nuclei only. Expansion is in µm, so it uses the
image's pixel size — if the OME metadata lacks it, pass `--pixel-size`.

## Do I need a GPU?

Strictly, no — but practically, yes. Segmentation is the slowest stage, and
running it on CPU is dramatically slower. A GPU (an A100 or similar) is
strongly recommended. See [Installation](installation.md).

## GPU isn't detected — what now?

1. Confirm the framework sees the GPU. For Stardist (TensorFlow):

   ```bash
   pixi run python -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"
   ```

   For InstanSeg (PyTorch):

   ```bash
   pixi run -e instanseg python -c "import torch; print(torch.cuda.is_available())"
   ```

2. If a GPU is present but not listed, try requesting it with `--use-gpu`
   (Stardist). InstanSeg can only use the GPU if PyTorch reports CUDA as
   available.

3. If there is genuinely no GPU, the pipeline still runs — segmentation just
   takes much longer on CPU. (On Windows the environments are CPU-only.)

## What does "Autotuner could not find any supported configs" mean?

This is a known TensorFlow / cuDNN failure when TensorFlow ≥2.16 defaults to
Keras 3. CellSurvey avoids it by:

- setting `TF_USE_LEGACY_KERAS=1` before imports (in `run.py`), and
- requiring the `tf_keras` pip package.

If you see this error, you are likely running without that flag set (for
example, invoking a module directly in an environment missing `tf_keras`). Make
sure you use `run.py` (or `python -m cellsurvey.cli`) from the pixi environment.

## Why do I get `ModuleNotFoundError: No module named 'tf_keras'`?

`tf_keras` is a required dependency for the legacy-Keras path. Re-run
`pixi install` to ensure the environment is complete.

## How do I resume a crashed run?

Use `--resume-from` with the path to an existing Zarr:

```bash
pixi run python run.py -i <image> -o <output> -p <plots> --resume-from <output>.zarr
```

This skips image loading and spot detection. If a segmented Zarr (`*_seg.zarr`)
already exists, the pipeline also skips the stages it can reuse. Pass the same
`--segmenter` (and use the same pixi environment) as the original run. See
[Pipeline](pipeline.md) for the exact reuse rules.

## Which parameters should I change from their defaults?

For a first run, none — the defaults are sensible. When you want to tune:

- `--n-clusters` — the k-means resolution (there is no automatic selection).
- `--community-resolution` — higher → more, smaller communities.
- `--max-edge-distance` — how dense the spatial network is (µm, default `50`).
- `--thresholds` (with `--detect-blobs`) — spot-detection sensitivity.

See [Parameters](parameters.md) for full guidance.

## Can I process formats other than OME-TIFF?

Not currently. Only multi-channel TIFF is supported (`.ome.tiff`, `.tiff`,
`.tif`). Formats such as ND2, CZI, LIF, and DV are not accepted yet.

## What does "No cell was returned by the segmentation" mean?

The segmenter found no cells, typically because the patches were empty or only
contained image edges with no nuclei. Check that:

- the input image actually contains nuclei on the channel being segmented
  (`--seg-channels`),
- the image isn't blank or mislabelled, and
- for InstanSeg, the pixel size is correct — from the OME metadata or
  `--pixel-size`. InstanSeg rescales the image to its training resolution with
  it, so a wrong value can make it find few or no cells.

## Why do I get an error saying InstanSeg needs the `instanseg` environment?

InstanSeg (PyTorch) is only installed in the `instanseg` pixi environment, and
Stardist (TensorFlow) only in the `default` one. Run `--segmenter instanseg`
with `pixi run -e instanseg ...`, and Stardist with plain `pixi run ...`.

## My GeoJSON ended up in the wrong place — how do I control that?

The QuPath GeoJSON defaults to `./qupath_export.geojson`. Change it with
`--geojson-path`:

```bash
pixi run python run.py -i <image> -o <output> -p <plots> --geojson-path path/to/export.geojson
```

---

Back to [Home](index.md) or [Pipeline](pipeline.md) for the full stage
walkthrough.
