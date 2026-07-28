[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](https://www.python.org/downloads/release/python-3100/) [![Built with Pixi](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json)](https://pixi.sh) ![Commit activity](https://img.shields.io/github/commit-activity/y/FrancisCrickInstitute/Spatial-Biology-Pipeline?style=plastic) ![GitHub](https://img.shields.io/github/license/FrancisCrickInstitute/Spatial-Biology-Pipeline?color=green&style=plastic)

# Overview

SopaSpan is a Python library for the analysis of spatial biology/omics data. It is heavily based on Sopa and MuSpan, combining elements of both into a single generic workflow:
* Blampey, Q., Mulder, K., Gardet, M. et al. Sopa: a technology-invariant pipeline for analyses of image-based spatial omics. _Nat Commun_ 15, 4981 (2024).
* Bull, J. A., Moore, J. W., Corry S. M., el al. MuSpAn: A Toolbox for Multiscale Spatial Analysis. _bioRxiv_ 2024.12.06.627195

<img width="350" height="350" alt="cell_type_to_cell_type" src="https://github.com/user-attachments/assets/9807689f-f471-49a5-b1ef-d701cb2db1c8" />
<img width="466" height="350" alt="umap_leiden" src="https://github.com/user-attachments/assets/ae31a25c-889b-41e6-ac7d-4af4df775a6b" />

# Installation

SopaSpan can be installed using either pixi (recommended) or conda.

> [!NOTE]
> SopaSpan depends on Tensorflow and while Tensorflow will run on all operating systems, support for GPU processing is generally only supported on Linux - see [here](https://www.tensorflow.org/install) for more information.

## Option A: Pixi (recommended)

[Pixi](https://pixi.sh/latest/) manages the full environment including Python, CUDA, cuDNN, and all Python packages in a single command.

First, [install pixi](https://pixi.sh/latest/#installation) if you haven't already:

```bash
# Linux / macOS
curl -fsSL https://pixi.sh/install.sh | bash

# Windows (PowerShell)
iwr -useb https://pixi.sh/install.ps1 | iex
```

Then, from the repository root:

```bash
pixi install
```

That's it — pixi reads `pixi.toml` and sets up everything. To run:

```bash
pixi run python sopaspan.py -i <input_tiff> -o <output_zarr> -p <plot_dir>
```

## Option B: Conda + pip

### Step 1: Install a Python Distribution

We recommend using conda as it's relatively straightforward and makes the management of different Python environments simple. You can install conda from [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation) (miniconda will suffice).

## Step 2: Create Environment and Install Pip Dependencies

### 2.1: Create an environment

Once conda is installed, open a terminal (Mac) or Anaconda Prompt (Windows) and run the following series of commands:

```bash
conda create --name spatial-bio python=3.12
conda activate spatial-bio
```

You will be presented with a list of packages to be downloaded and installed. The following prompt will appear:
```bash
Proceed ([y]/n)?
```
Hit Enter and all necessary packages will be downloaded and installed - this may take some time.

### 2.2: Install Tensorflow

SopaSpan depends on [Stardist](https://github.com/stardist/stardist) to segment cell nuclei, which in turn depends on Tensorflow.

>[!NOTE]
>Tensorflow can be run on CPUs, but this can be quite slow. To speed things up, a GPU-compatible installation is recommended. In order to enable this, you need to have the [the necessary CUDA drivers](https://developer.nvidia.com/cuda/toolkit). On linux, you can typically load the drivers ([CUDA](https://developer.nvidia.com/cuda) and [cuDNN](https://developer.nvidia.com/cudnn)) with commands such as:
>```shell
>ml CUDA/12.5.1 
>ml cuDNN/9.3.0.75-CUDA-12.5.1
>```

Then, install tensorflow as follows:

```bash
python -m pip install tensorflow[and-cuda]
```

On any other operating system, or for a CPU-only installation, use the following:

```bash
python -m pip install tensorflow
```

### 2.3: Install Sopa and bioio

Install [Sopa](https://gustaveroussy.github.io/sopa/) with support for stardist and wsi (whole slide imaging), plus [bioio](https://github.com/bioio-devs/bioio) for image metadata:

```bash
python -m pip install 'sopa[stardist,wsi]' bioio bioio-ome-tiff
```

### 2.4: Install MuSpan

Unfortunately, at this time, MuSpan requires a username and password to install. You can obtain these by completing the form [here](https://www.muspan.co.uk/get-the-code). Once you receive a response by email, MuSpan can be installed as follows:

```bash
python -m pip install https://docs.muspan.co.uk/code/latest.zip
```

You will then be prompted to enter the login credentials you received by email and the installation should proceed.

## Step 3: Get the code for this repository

To get the necessary python code to run SopaSpan, the recommended approach is to [clone this repository using Git](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository). Alternatively, you can download a Zip file of the repo by clicking on the green code button above and then clicking "Download Zip":

<img width="472" height="369" alt="image" src="https://github.com/user-attachments/assets/ee52fdf5-1574-4342-aa85-77f623d60709" />

Unzip the contents of the zip file once downloaded - the contents should contain a file called `sopaspan.py` and a `sopaspan/` directory of supporting modules.

## Installation Complete

That's it - your set up is complete. You can deactivate the environment you have created with the following command.

```bash
conda deactivate
```

# Usage

## Basic Usage

To run SopaSpan, use the following:

```bash
conda activate spatial-bio
python <path_to_sopaspan.py> -i <path_to_input_file> -o <path_to_output_zarr> -p <path_to_output_plots_directory>
```

## Parameters

### Required
* `-i`, `--input_file`: Path to the input image. While only TIFF files have been tested, most common microscopy formats should work.
* `-o`, `--output_file`: Path for the output Zarr file (`.zarr` suffix appended automatically if missing). The input image is converted to a [SpatialData](https://www.nature.com/articles/s41592-024-02212-x) Zarr object.
* `-p`, `--plot_dir`: Directory where all output plots will be saved.

### Blob detection
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

### Crash recovery
* `--resume-from`: Path to an existing Zarr file to resume from. Skips image loading and spot detection, resuming directly at Stardist segmentation.

## Full Example

```bash
python ~/Downloads/SopaSpan/sopaspan.py -i ~/data/sample.tiff -o ~/results/output.zarr -p ~/results/plots/
```
