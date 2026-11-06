# Installation

CellSurvey runs on **Linux (64-bit)** with an NVIDIA GPU. It is installed and
run through [pixi](https://pixi.sh/), which manages Python, TensorFlow or
PyTorch, CUDA, and every other dependency in self-contained environments.

!!! note "Supported platforms"
    **Linux (64-bit)** is the supported platform, with GPU acceleration.
    **Windows (64-bit)** environments also resolve, but are CPU-only and meant
    for local testing. macOS is **not** supported. See the
    [TensorFlow install guide](https://www.tensorflow.org/install) for details
    of GPU support.

!!! info "Two environments, one per segmenter"
    pixi defines a separate environment for each segmenter, so TensorFlow and
    PyTorch are never installed together:

    | Environment | Segmenter | Framework |
    |---|---|---|
    | `default` (alias `spatial`) | Stardist (`--segmenter stardist`, the default) | TensorFlow |
    | `instanseg` | InstanSeg (`--segmenter instanseg`) | PyTorch |

    You only need the environment for the segmenter you plan to use. The
    project plans to migrate fully to InstanSeg and eventually drop Stardist.

## Prerequisites

- **Linux 64-bit** with an NVIDIA GPU (an A100 or similar is ideal).
- **Git**, to clone the repository.
- A **pixi** installation (see below).

## 1. Install pixi

pixi is a single binary. Install it from the official script:

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

Then reload your shell so `pixi` is on your `PATH`:

```bash
source ~/.bashrc
```

Check it works:

```bash
pixi --version
```

## 2. Get the code

```bash
git clone https://github.com/FrancisCrickInstitute/CellSurvey.git
cd CellSurvey
```

## 3. Install the environment

From the repository root, run:

```bash
pixi install
```

pixi reads `pixi.toml`, resolves the pinned dependency set, and creates the
environment in `.pixi/`. This step downloads TensorFlow (with bundled CUDA 12 /
cuDNN 9), StarDist, Sopa, and the analysis libraries — it can take a few minutes
the first time.

This installs the `default` (Stardist) environment. To use InstanSeg instead,
install the `instanseg` environment, which has PyTorch and InstanSeg in place of
TensorFlow and StarDist:

```bash
pixi install -e instanseg
```

!!! tip "No separate CUDA install"
    You do **not** need to install CUDA or cuDNN yourself, and you should **not**
    load `CUDA`/`cuDNN` modules on a cluster. TensorFlow 2.18 bundles its own
    CUDA 12 and cuDNN 9 libraries inside the pixi environment, and the Linux
    PyTorch wheels in the `instanseg` environment bundle CUDA likewise.

## 4. Check the environment

Confirm Python and GPU detection work. For the `default` (Stardist)
environment:

```bash
pixi run python -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"
```

If all is well, you should see a `PhysicalDevice` entry for your GPU.

For the `instanseg` environment:

```bash
pixi run -e instanseg python -c "import torch; print(torch.cuda.is_available())"
```

This should print `True`.

!!! warning "If no GPU is listed"
    See [FAQ — "GPU isn't detected"](faq.md) for what to try. CellSurvey can
    still run on CPU, but segmentation will be very slow.

## Docker (alternative)

A Docker image is also available for containerised runs:

```bash
docker build -t cellsurvey .
docker run --gpus all -v /path/to/data:/data cellsurvey \
  -i /data/input.ome.tiff -o /data/output -p /data/plots
```

The image installs pixi, resolves the environment, copies in the pipeline code,
and sets the same `TF_USE_LEGACY_KERAS=1` flag used by the local entry point.

The default build uses the Stardist environment. To build an InstanSeg image,
pass the `PIXI_ENV` build argument, then add `--segmenter instanseg` when running:

```bash
docker build --build-arg PIXI_ENV=instanseg -t cellsurvey-instanseg .
docker run --gpus all -v /path/to/data:/data cellsurvey-instanseg \
  -i /data/input.ome.tiff -o /data/output -p /data/plots --segmenter instanseg
```

---

Next: [Getting started](getting-started.md) to run your first analysis.
