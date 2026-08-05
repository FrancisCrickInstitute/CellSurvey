[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/) [![Built with Pixi](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json)](https://pixi.sh) ![Commit activity](https://img.shields.io/github/commit-activity/y/FrancisCrickInstitute/Spatial-Biology-Pipeline?style=plastic) ![GitHub](https://img.shields.io/github/license/FrancisCrickInstitute/Spatial-Biology-Pipeline?color=green&style=plastic)

# Overview

CellSurvey is a Python pipeline for image-based spatial biology/omics analysis built on [Sopa](https://gustaveroussy.github.io/sopa/). It processes multichannel microscopy images (OME-TIFF) into spatial data objects, segments nuclei with Stardist, detects RNA spots via blob detection, clusters cells with k-means, builds Delaunay networks via scipy, detects Louvain communities via networkx, and exports GeoJSON for QuPath visualization.

* Blampey, Q., Mulder, K., Gardet, M. et al. Sopa: a technology-invariant pipeline for analyses of image-based spatial omics. _Nat Commun_ 15, 4981 (2024).

<img width="350" height="350" alt="cell_type_to_cell_type" src="https://github.com/user-attachments/assets/9807689f-f471-49a5-b1ef-d701cb2db1c8" />
<img width="466" height="350" alt="umap_leiden" src="https://github.com/user-attachments/assets/ae31a25c-889b-41e6-ac7d-4af4df775a6b" />

# Installation

CellSurvey uses **pixi** for environment management targeting Linux (64-bit) with GPU support.

> [!NOTE]
> CellSurvey depends on TensorFlow and while TensorFlow will run on all operating systems, support for GPU processing is generally only supported on Linux ,  see [here](https://www.tensorflow.org/install) for more information. Windows and macOS are not supported via pixi.

## Pixi (recommended)

[Pixi](https://pixi.sh/latest/) manages the full environment including Python, CUDA, cuDNN, and all Python packages in a single command.

First, [install pixi](https://pixi.sh/latest/#installation) if you haven't already:

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

From the repository root:

```bash
pixi install
```

That's it ,  pixi reads `pixi.toml` and sets up everything. To run:

```bash
pixi run python run.py -i <input_tiff> -o <output_zarr> -p <plot_dir>
```

## Docker

```bash
docker build -t cellsurvey .
docker run --gpus all -v /path/to/data:/data cellsurvey -i /data/input.tiff -o /data/output -p /data/plots
```

# Usage

## Basic Usage

```bash
pixi run python run.py -i <path_to_input_file> -o <path_to_output_zarr> -p <path_to_output_plots_directory>
```

## Parameters

### Required
* `-i`, `--input_file`: Path to the input image. While only TIFF files have been tested, most common microscopy formats should work.
* `-o`, `--output_file`: Path for the output Zarr file (`.zarr` suffix appended automatically if missing). The input image is converted to a [SpatialData](https://www.nature.com/articles/s41592-024-02212-x) Zarr object.
* `-p`, `--plot_dir`: Directory where all output plots will be saved.

### Blob detection
* `--detect-blobs`: Enable RNA spot blob detection on the specified channels. Without this flag, blob detection is skipped entirely (default: off).
* `--channels`: Comma-separated channel indices to process (default: `9,10,11,12`)
* `--thresholds`: Comma-separated detection thresholds per channel (default: `0.01,0.1,0.1,0.1`)
* `--tile-size`: Tile size in pixels for tiled processing (default: `2048`)
* `--overlap`: Overlap between tiles in pixels (default: `50`)
* `--workers`: Number of worker threads for parallel detection (default: `14`)
* `--min-sigma`: Minimum blob radius for spot detection (default: `2`)
* `--max-sigma`: Maximum blob radius for spot detection (default: `5`)
* `--num-sigma`: Number of sigma steps for blob detection (default: `5`)

### Clustering and network analysis
* `--n-clusters`: Number of k-means clusters (default: `10`)
* `--community-resolution`: Louvain community detection resolution (default: `0.1`)
* `--max-edge-distance`: Maximum edge distance for Delaunay network (default: `1000`)

### Spatial analysis
* `--radius-min`: Minimum radius for spatial neighbors graph (default: `0`)
* `--radius-max`: Maximum radius for spatial neighbors graph (default: `1000`)

### GPU and Stardist segmentation
* `--use-gpu`: Force GPU usage for Stardist segmentation. Without this flag, GPU is auto-detected and used if available. Useful when auto-detection fails (default: off).

### Output and visualization
* `--geojson-path`: Output path for QuPath GeoJSON (default: `./qupath_export.geojson`)
* `--fig-size`: Figure size for plots (default: `20`)
* `--font-size`: Font size for plots (default: `20`)
* `--axes-linewidth`: Axes line width for plots (default: `3`)

### Crash recovery
* `--resume-from`: Path to an existing Zarr file to resume from. Skips image loading and spot detection, resuming directly at Stardist segmentation.

## Full Example

```bash
pixi run python run.py -i ~/data/sample.tiff -o ~/results/output -p ~/results/plots/
```

# Visualising Results

## 1. Odon (recommended)

[Odon](https://github.com/alexcoulton/odon) is a lightweight, GPU-accelerated viewer for SpatialData Zarr files. It supports multiscale image viewing with shape overlays colored by metadata columns.

To visualise CellSurvey output in Odon:

1. [Install Odon](https://github.com/alexcoulton/odon#installation)
2. Open the segmented Zarr file (`*_seg.zarr`)
3. The `stardist_boundaries` shape layer can be colored by cluster or community using the **"Color by"** dropdown

**Note:** Only categorical columns with ≤24 distinct values appear in the dropdown. `kmeans_cluster` (10 values) and `community` (up to 24) will be visible. Channel intensity columns are available as continuous properties but are not shown as color-by options.

## 2. QuPath

For traditional pathology workflows, CellSurvey exports a GeoJSON file compatible with [QuPath](https://qupath.github.io/). Open `qupath_export.geojson` via File → Import → GeoJSON. Cell boundaries appear as annotations colored by community, and RNA spot detections (when `--detect-blobs` is enabled) appear as detection objects.

## 3. TissUUmaps

[TissUUmaps](https://github.com/TissUUmaps/TissUUmaps4) is a GPU-accelerated browser-based viewer for spatial biology data. It runs entirely in-browser ,  no installation required ,  and supports SpatialData Zarr files natively via the OME-Zarr + SpatialData plugin.

1. Open [TissUUmaps live](https://tissuumaps.github.io/TissUUmaps4/live/) or download the [latest release](https://github.com/TissUUmaps/TissUUmaps4/releases)
2. Load the segmented Zarr file (`*_seg.zarr`)
3. Cell boundaries appear as a shapes layer; link to metadata columns for color-by-cluster or color-by-community

**Note:** TissUUmaps 4 is under active development. The stable release (v3) may not include SpatialData Zarr support ,  use the v4 development builds. Shapes require matching coordinate systems and an ID column for metadata linkage. A modern browser with WebGL 2 and File System API is required.

## 4. napari + napari-spatialdata

[napari](https://napari.org/) is a multi-dimensional image viewer for Python, and [napari-spatialdata](https://github.com/scverse/napari-spatialdata) adds native SpatialData Zarr support.

```bash
pip install "napari-spatialdata[all]"
```

Open the segmented Zarr in napari:

```python
from napari_spatialdata import Interactive
from spatialdata import SpatialData

sdata = SpatialData.read("path/to/output_seg.zarr")
Interactive(sdata).run()
```

Select a coordinate system, click `stardist_boundaries` to load cell shapes, then use the **View** widget (Plugins → napari-spatialdata → View) to color cells by cluster, community, or any channel intensity column. Double-click any `obs` column to apply it as the face color.

**Performance notes for large datasets:**
- Shape loading is slow above ~50K polygons due to triangulation. Use the `bermuda` backend for faster loading: `pip install "napari-spatialdata[all,bermuda]"`
- Polygons are simplified when the shape count exceeds 100 (configurable via `napari_spatialdata.constants.config.POLYGON_THRESHOLD`)
- If boundaries appear too simplified at high zoom, increase the threshold:<br/>
  `from napari_spatialdata.constants import config; config.POLYGON_THRESHOLD = 50000`
