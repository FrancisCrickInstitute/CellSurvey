# Parameters

CellSurvey is controlled entirely through command-line flags. For most runs you
only need the three required arguments — the rest have sensible defaults.

```bash
pixi run python run.py -i <image> -o <output> -p <plots> [options]
```

## Required arguments

| Flag | Meaning |
|---|---|
| `-i`, `--input_file` | Path to the input multichannel TIFF (OME-TIFF recommended). |
| `-o`, `--output_file` | Prefix for the output Zarr (`.zarr` is appended automatically). |
| `-p`, `--plot_dir` | Directory for the summary plots (default `.`). |

## Spot / blob detection

Spot detection is **off by default**. It runs only when `--detect-blobs` is
passed, and detects RNA spots on the channels you specify.

| Flag | Default | What it controls |
|---|---|---|
| `--detect-blobs` | off | Enable RNA spot detection on the specified channels. |
| `--channels` | `9,10,11,12` | Comma-separated channel indices to process. |
| `--thresholds` | `0.01,0.1,0.1,0.1` | Detection threshold per channel (one per `--channels`). |
| `--tile-size` | `2048` | Tile size in pixels for tiled processing. |
| `--overlap` | `50` | Overlap between tiles in pixels. |
| `--workers` | `14` | Worker threads for parallel detection. |
| `--min-sigma` | `2` | Smallest blob radius to look for. |
| `--max-sigma` | `5` | Largest blob radius to look for. |
| `--num-sigma` | `5` | Number of scale steps between `min` and `max`. |

!!! warning "Channels and thresholds must match"
    The number of entries in `--channels` and `--thresholds` must be equal, or
    the run stops with an error.

!!! tip "Tuning spot detection"
    - `--thresholds` is the main **sensitivity** dial — lower values detect more
      (and noisier) spots, higher values are stricter.
    - `--min-sigma` / `--max-sigma` / `--num-sigma` set the **range of blob sizes**
      searched (the Laplacian-of-Gaussian scale space).
    - `--tile-size` / `--overlap` / `--workers` are **performance** knobs. The
      overlap prevents double-counting at tile seams.

## Clustering and network analysis

| Flag | Default | What it controls |
|---|---|---|
| `--n-clusters` | `10` | Number of k-means clusters. |
| `--community-resolution` | `0.1` | Louvain community resolution. |
| `--max-edge-distance` | `1000` | Max edge length for the Delaunay network. |

!!! tip "Tuning clustering and communities"
    - `--n-clusters` is the k-means **k** — there is no automatic selection, so
      this is a value you must choose.
    - `--community-resolution` is the Louvain **resolution**: higher values give
      more, smaller communities; lower values give fewer, larger ones.
    - `--max-edge-distance` controls **network density**: larger values connect
      more distant cells (a denser graph), smaller values keep only close
      neighbours.

## Spatial analysis

| Flag | Default | What it controls |
|---|---|---|
| `--radius-min` | `0` | Inner radius of the spatial-neighbours band. |
| `--radius-max` | `1000` | Outer radius of the spatial-neighbours band. |

These set the radius band used for the mean-hop-distance heatmap between
clusters (`cell_type_to_cell_type.png`).

## Segmentation

| Flag | Default | What it controls |
|---|---|---|
| `--segmenter` | `stardist` | Segmentation method: `stardist` or `instanseg`. |
| `--seg-channels` | `0` (Stardist) / all channels (InstanSeg) | Comma-separated channel indices to segment on. Stardist takes exactly one. |
| `--cell-expansion` | `5.0` | **Stardist only.** Radius in µm to grow each nucleus into a whole cell; `0` = nuclei only. Ignored for InstanSeg. |
| `--instanseg-model` | `fluorescence_nuclei_and_cells` | **InstanSeg only.** Pretrained model name, or a path to a model. |
| `--pixel-size` | *(from image metadata)* | Pixel size in µm/px, overriding the OME metadata. |
| `--use-gpu` | off | Request GPU segmentation (both segmenters). The GPU is auto-detected and used when available. |

!!! info "Stardist vs InstanSeg"
    - **Stardist** (default) segments **nuclei** only, on a single channel. Whole
      cells are then approximated by expanding each nucleus outward by
      `--cell-expansion` µm with a watershed, so neighbouring cells stop where
      they meet (as in QuPath's cell detection).
    - **InstanSeg** segments **nuclei and whole cells** directly. It is
      channel-invariant, so by default it uses every channel. It needs the
      separate `instanseg` pixi environment — see
      [Installation](installation.md):

    ```bash
    pixi run -e instanseg python run.py -i <image> -o <output> -p <plots> \
      --segmenter instanseg
    ```

    The project plans to migrate fully to InstanSeg and eventually drop Stardist.

!!! warning "Pixel size matters"
    The pixel size is read from the OME-TIFF metadata. If the image has none,
    pass `--pixel-size`. Without it, Stardist cannot expand nuclei into cells
    (the run continues with nuclei only), and InstanSeg cannot rescale the image
    to its training resolution, so it may find few or no cells.

By default the GPU is auto-detected. With Stardist, `--use-gpu` tells the
pipeline to attempt the GPU backend when detection fails but a GPU is present.
InstanSeg always uses the GPU when PyTorch can see CUDA. If it can't, `--use-gpu`
only adds a warning. Without a GPU, segmentation runs on CPU and is very slow.

## Crash recovery

| Flag | Default | What it controls |
|---|---|---|
| `--resume-from` | *(none)* | Path to an existing Zarr to resume from. |

`--resume-from` skips image loading and spot detection, resuming directly at
segmentation. Combined with the segmented-Zarr reuse (see
[Pipeline](pipeline.md)), it provides two levels of checkpoint restart.

## Output and visualisation

| Flag | Default | What it controls |
|---|---|---|
| `--geojson-path` | `./qupath_export.geojson` | Output path for the QuPath GeoJSON. |
| `--font-size` | `20` | Font size for plots. |
| `--axes-linewidth` | `3` | Axes line width for plots. |

!!! note "Where things go"
    Plots and `summary.json` are written to `--plot_dir`; the GeoJSON goes to
    `--geojson-path` (independent of `-o` and `-p`).

## Parameter sweeps

You can explore how clustering responds to its parameters without re-running the
expensive segmentation, using the read-only sweep over an existing segmented
Zarr. Results are written to CSVs in the same directory as the Zarr; the Zarr
itself is never modified.

### Community sweep

Sweep the Louvain resolution (and Delaunay edge distance):

```bash
pixi run python -m cellsurvey.stability \
  --zarr <output>_seg.zarr \
  --resolutions 0.1,0.05,0.02,0.01 \
  --max-edge-distances 1000
```

This writes `community_sweep.csv` (one column per combination,
`community_r<res>_d<dist>`, plus `cell_id`/`x`/`y`) and
`community_sweep_summary.csv` (community count per combination).

### Cluster sweep

Sweep the k-means `k` (defaults to `5,8,10,12,15,20`):

```bash
pixi run python -m cellsurvey.stability \
  --zarr <output>_seg.zarr \
  --mode clusters
```

This writes `cluster_sweep.csv` (a `cluster_k<k>` column per value) and
`cluster_sweep_summary.csv` (the k-means inertia per `k`, for an elbow plot).

### Which shapes are swept

Both modes use the `cell_boundaries` layer when the Zarr has one, otherwise the
nuclei layer (`stardist_boundaries` or `nucleus_boundaries`), so nuclei-only
Zarrs work too. Override it with `--cell-boundaries <layer>`.

!!! note "Clusters are independent of communities"
    `--n-clusters` does **not** affect the communities: Louvain weights its graph
    from the raw intensity matrix, not the k-means labels. Sweep them separately
    — `--n-clusters` for the expression clusters, `--community-resolution` /
    `--max-edge-distance` for the spatial communities.

---

Next: [Outputs](outputs.md) for a description of every file the pipeline
produces.
