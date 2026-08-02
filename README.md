[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/) [![Built with Pixi](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json)](https://pixi.sh) ![Commit activity](https://img.shields.io/github/commit-activity/y/FrancisCrickInstitute/Spatial-Biology-Pipeline?style=plastic) ![GitHub](https://img.shields.io/github/license/FrancisCrickInstitute/Spatial-Biology-Pipeline?color=green&style=plastic)

# Overview

CellSurvey is a Python pipeline for image-based spatial biology/omics analysis built on [Sopa](https://gustaveroussy.github.io/sopa/). It processes multichannel microscopy images (OME-TIFF) into spatial data objects, segments nuclei with Stardist, detects RNA spots via blob detection, clusters cells with k-means, builds Delaunay networks via scipy, detects Louvain communities via networkx, and exports GeoJSON for QuPath visualization.

* Blampey, Q., Mulder, K., Gardet, M. et al. Sopa: a technology-invariant pipeline for analyses of image-based spatial omics. _Nat Commun_ 15, 4981 (2024).

<img width="350" height="350" alt="cell_type_to_cell_type" src="https://github.com/user-attachments/assets/9807689f-f471-49a5-b1ef-d701cb2db1c8" />
<img width="466" height="350" alt="umap_leiden" src="https://github.com/user-attachments/assets/ae31a25c-889b-41e6-ac7d-4af4df775a6b" />

# Installation

CellSurvey uses **pixi** for environment management targeting Linux (64-bit) with GPU support.

> [!NOTE]
> CellSurvey depends on TensorFlow and while TensorFlow will run on all operating systems, support for GPU processing is generally only supported on Linux — see [here](https://www.tensorflow.org/install) for more information. Windows and macOS are not supported via pixi.

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

That's it — pixi reads `pixi.toml` and sets up everything. To run:

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

### Clustering and network analysis
* `--n-clusters`: Number of k-means clusters (default: `10`)
* `--community-resolution`: Louvain community detection resolution (default: `0.1`)
* `--max-edge-distance`: Maximum edge distance for Delaunay network (default: `1000`)

### Spatial analysis
* `--radius-min`: Minimum radius for spatial neighbors graph (default: `0`)
* `--radius-max`: Maximum radius for spatial neighbors graph (default: `1000`)

### GPU and Stardist segmentation
* `--use-gpu`: Force GPU usage for Stardist segmentation. Without this flag, GPU is auto-detected and used if available. Useful when auto-detection fails (default: off).

### Crash recovery
* `--resume-from`: Path to an existing Zarr file to resume from. Skips image loading and spot detection, resuming directly at Stardist segmentation.

## Full Example

```bash
pixi run python run.py -i ~/data/sample.tiff -o ~/results/output -p ~/results/plots/
```
