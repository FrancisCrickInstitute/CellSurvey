[![Documentation Status](https://readthedocs.org/projects/cell-survey/badge/?version=latest)](https://cell-survey.readthedocs.io/en/latest/) [![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/) [![Built with Pixi](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json)](https://pixi.sh) ![Commit activity](https://img.shields.io/github/commit-activity/y/FrancisCrickInstitute/CellSurvey?style=plastic) ![GitHub](https://img.shields.io/github/license/FrancisCrickInstitute/CellSurvey?color=green&style=plastic)

# Overview

CellSurvey is a Python pipeline for image-based spatial biology/omics analysis built on [Sopa](https://gustaveroussy.github.io/sopa/). It processes multichannel microscopy images (OME-TIFF) into spatial data objects, segments nuclei with Stardist, detects RNA spots via blob detection, clusters cells with k-means, builds Delaunay networks via scipy, detects Louvain communities via networkx, and exports GeoJSON for QuPath visualization.

* Blampey, Q., Mulder, K., Gardet, M. et al. Sopa: a technology-invariant pipeline for analyses of image-based spatial omics. _Nat Commun_ 15, 4981 (2024).

<img width="350" height="350" alt="cell_type_to_cell_type" src="https://github.com/user-attachments/assets/9807689f-f471-49a5-b1ef-d701cb2db1c8" />
<img width="466" height="350" alt="umap_leiden" src="https://github.com/user-attachments/assets/ae31a25c-889b-41e6-ac7d-4af4df775a6b" />

# Documentation

Full documentation is available at **[cell-survey.readthedocs.io](https://cell-survey.readthedocs.io/en/latest/)** — installation, a pipeline walkthrough, every parameter, the outputs reference, the visualisation guide, and an FAQ.

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

A full reference of every flag, with tuning guidance, is in the
[Parameters](https://cell-survey.readthedocs.io/en/latest/parameters/) page. The
three required arguments are:

* `-i`, `--input_file`: Path to the input multichannel TIFF image (OME-TIFF recommended).
* `-o`, `--output_file`: Path for the output Zarr (`.zarr` appended automatically).
* `-p`, `--plot_dir`: Directory where output plots are saved.

## Full Example

```bash
pixi run python run.py -i ~/data/sample.tiff -o ~/results/output -p ~/results/plots/
```

## Stability sweep

To explore how clustering responds to its parameters without re-running the
expensive segmentation, use the read-only sweep over an existing segmented Zarr.

Sweep the Louvain resolution (and Delaunay edge distance) for communities:

```bash
pixi run python -m cellsurvey.stability \
  --zarr <output>_seg.zarr \
  --resolutions 0.1,0.05,0.02,0.01 \
  --max-edge-distances 1000
```

Or sweep the k-means `k`:

```bash
pixi run python -m cellsurvey.stability \
  --zarr <output>_seg.zarr \
  --n-clusters 5,8,10,12,15,20
```

Each writes per-cell labels (plus a summary CSV with the community count, or the
k-means inertia) and leaves the Zarr untouched.

# Visualising Results

For the full visualisation guide — Odon (recommended), QuPath, TissUUmaps, and
napari + napari-spatialdata, with performance notes — see the
[Visualising results](https://cell-survey.readthedocs.io/en/latest/visualization/)
page.
