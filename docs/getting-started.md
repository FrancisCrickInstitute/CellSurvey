# Getting started

This guide takes you through your first full run of CellSurvey, from a raw
image to labelled cells and communities.

!!! tip "Before you start"
    Make sure you have [installed](installation.md) the environment and confirmed
    your GPU is detected.

## The basic command

CellSurvey runs through `run.py` with three required arguments:

```bash
pixi run python run.py \
  -i path/to/image.ome.tiff \
  -o path/to/output \
  -p path/to/plots
```

| Argument | Meaning |
|---|---|
| `-i`, `--input_file` | Path to the input multichannel image (OME-TIFF) |
| `-o`, `--output_file` | Prefix for the output Zarr file (`.zarr` is appended automatically) |
| `-p`, `--plot_dir` | Directory where summary plots are written |

!!! note "Input formats"
    The input must be a **multi-channel TIFF** — OME-TIFF (`.ome.tiff`) is
    recommended. Its embedded channel names become the marker column headers in
    your results. Plain `.tiff`/`.tif` also works, but channels then fall back to
    numeric names. Non-TIFF formats such as ND2, CZI, LIF, or DV are **not**
    currently supported.

## A realistic example

```bash
pixi run python run.py \
  -i ~/data/sample.ome.tiff \
  -o ~/results/output \
  -p ~/results/plots
```

This reads `sample.ome.tiff`, writes two Zarr stores in `~/results/`, and puts
all plots in `~/results/plots/`.

## What happens

When you run the command, CellSurvey works through its stages in order:

1. **Load** the image and detect channels.
2. **Segment** cell nuclei with Stardist (GPU-accelerated).
3. **Measure** the intensity of each marker channel inside every cell.
4. **Cluster** cells by their expression profile (k-means).
5. **Build a spatial network** and find communities (Delaunay + Louvain).
6. **Export** results and **plot** summaries.

Progress is printed to the terminal as each stage completes. You'll see lines
like:

```
Found 1 GPU(s), using GPU backend for Stardist
Run stardist...
Aggregating...
Performing k-means clustering...
Found 10 communities
```

!!! note "Segmentation can take a while"
    Stardist segmentation is the slowest stage. For large images it can take
    minutes to tens of minutes. This is expected.

## What you get

When the run finishes, you'll have:

| Output | Location | What it is |
|---|---|---|
| Initial Zarr | `~/results/output.zarr` | The image converted to SpatialData format |
| Segmented Zarr | `~/results/output_seg.zarr` | Cells, boundaries, and per-cell measurements |
| QuPath GeoJSON | `./qupath_export.geojson` | Cell boundaries + clusters + communities |
| Summary | `~/results/plots/summary.json` | Counts of cells, clusters, communities |
| Plots | `~/results/plots/*.png` | Density maps, UMAP, heatmaps |

The **segmented Zarr** (`*_seg.zarr`) is the main result — it contains
everything downstream analysis and visualisation need.

!!! tip "Interpreting the outputs"
    See [Outputs](outputs.md) for a full description of every file, and
    [Visualising results](visualization.md) to open them in a viewer.

## Enabling RNA spot detection

By default, spot detection is **off**. If your image contains RNA spots you want
to assign to cells, add `--detect-blobs`:

```bash
pixi run python run.py \
  -i ~/data/sample.ome.tiff \
  -o ~/results/output \
  -p ~/results/plots \
  --detect-blobs
```

Spot detection adds a `spots` layer to the output, which is exported as
detection objects in the QuPath GeoJSON. See [Parameters](parameters.md) for the
channel/threshold options.

---

Next: [Pipeline](pipeline.md) to understand each stage in detail.
