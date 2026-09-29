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

## GPU and Stardist segmentation

| Flag | Default | What it controls |
|---|---|---|
| `--use-gpu` | off | Force GPU for Stardist even if auto-detection fails. |

By default the GPU is auto-detected. Use `--use-gpu` when detection fails but a
GPU is present. Without a GPU, Stardist runs on CPU and is very slow.

## Crash recovery

| Flag | Default | What it controls |
|---|---|---|
| `--resume-from` | *(none)* | Path to an existing Zarr to resume from. |

`--resume-from` skips image loading and spot detection, resuming directly at
Stardist segmentation. Combined with the segmented-Zarr reuse (see
[Pipeline](pipeline.md)), it provides two levels of checkpoint restart.

## Output and visualisation

| Flag | Default | What it controls |
|---|---|---|
| `--geojson-path` | `./qupath_export.geojson` | Output path for the QuPath GeoJSON. |
| `--fig-size` | `20` | Figure size for plots. |
| `--font-size` | `20` | Font size for plots. |
| `--axes-linewidth` | `3` | Axes line width for plots. |

!!! note "Where things go"
    Plots and `summary.json` are written to `--plot_dir`; the GeoJSON goes to
    `--geojson-path` (independent of `-o` and `-p`).

## Parameter sweeps

To see how many communities result at different `--community-resolution` (and
`--max-edge-distance`) values — without re-running segmentation — use the
read-only sweep over an existing segmented Zarr:

```bash
pixi run python -m cellsurvey.stability \
  --zarr <output>_seg.zarr \
  --resolutions 0.1,0.05,0.02,0.01 \
  --max-edge-distances 1000 \
  --output community_sweep.csv \
  --summary community_sweep_summary.csv
```

It writes `community_sweep.csv` (one column per parameter combination,
`community_r<res>_d<dist>`, plus `cell_id`/`x`/`y`) and
`community_sweep_summary.csv` (community count per combination). The Zarr is
read but never modified.

!!! note "Clusters are independent of communities"
    `--n-clusters` does **not** affect the communities: Louvain weights its graph
    from the raw intensity matrix, not the k-means labels. The sweep therefore
    varies only `--community-resolution` and `--max-edge-distance`.

---

Next: [Outputs](outputs.md) for a description of every file the pipeline
produces.
