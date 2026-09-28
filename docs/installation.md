# Installation

CellSurvey runs on **Linux (64-bit)** with an NVIDIA GPU. It is installed and
run through [pixi](https://pixi.sh/), which manages Python, TensorFlow, CUDA,
and every other dependency in a single self-contained environment.

!!! note "Linux only"
    Windows and macOS are **not** supported. GPU acceleration for TensorFlow is
    only available on Linux — see the
    [TensorFlow install guide](https://www.tensorflow.org/install) for details.

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

!!! tip "No separate CUDA install"
    You do **not** need to install CUDA or cuDNN yourself, and you should **not**
    load `CUDA`/`cuDNN` modules on a cluster. TensorFlow 2.18 bundles its own
    CUDA 12 and cuDNN 9 libraries inside the pixi environment.

## 4. Check the environment

Confirm Python and GPU detection work:

```bash
pixi run python -c "import tensorflow as tf; print(tf.config.list_physical_devices('GPU'))"
```

If all is well, you should see a `PhysicalDevice` entry for your GPU.

!!! warning "If no GPU is listed"
    See [FAQ — "GPU isn't detected"](faq.md) for what to try. CellSurvey can
    still run on CPU, but Stardist segmentation will be very slow.

## Docker (alternative)

A Docker image is also available for containerised runs:

```bash
docker build -t cellsurvey .
docker run --gpus all -v /path/to/data:/data cellsurvey \
  -i /data/input.ome.tiff -o /data/output -p /data/plots
```

The image installs pixi, resolves the environment, and sets the same
`TF_USE_LEGACY_KERAS=1` flag used by the local entry point.

---

Next: [Getting started](getting-started.md) to run your first analysis.
