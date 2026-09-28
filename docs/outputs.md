# Outputs

A CellSurvey run produces a set of files organised around two Zarr stores, a
GeoJSON file, and a directory of plots. This page lists every artifact, what it
contains, and how to open it.

!!! tip "Where to look"
    - Zarr stores go next to `--output_file` (the segmented store replaces
      `.zarr` with `_seg.zarr`).
    - Plots and `summary.json` go in `--plot_dir`.
    - The GeoJSON goes to `--geojson-path` (default `./qupath_export.geojson`).

## Zarr stores

| File | What it is |
|---|---|
| `<output>.zarr` | The image converted to SpatialData format (the initial checkpoint, written before segmentation). |
| `<output>_seg.zarr` | The **main result** — cells, boundaries, and per-cell measurements. |

Both are [SpatialData](https://www.nature.com/articles/s41592-024-02212-x) Zarr
stores. The segmented store contains:

- `stardist_boundaries` — cell polygon shapes, with `kmeans_cluster` and
  `community` labels attached.
- `table` — an AnnData table of cells × markers (per-cell mean intensities).
- `spots` — detected RNA spots (only if `--detect-blobs` was used).

Open them with any SpatialData-aware viewer (see
[Visualising results](visualization.md)).

## `summary.json`

Written to `--plot_dir`. A small JSON summary of the network analysis:

```json
{
  "n_cells": 12345,
  "n_clusters": 10,
  "n_communities": 7,
  "n_edges": 48210,
  "max_edge_distance": 1000,
  "community_resolution": 0.1,
  "cluster_sizes": { "0": 1204, "1": 980, "...": "..." },
  "community_sizes": { "0": 2100, "1": 1550, "...": "..." }
}
```

Useful for a quick sanity check of how many cells, clusters, and communities
were found.

## Plots

All plots are PNG files written to `--plot_dir`.

| File | What it shows |
|---|---|
| `cell_type_to_cell_type.png` | Mean hop distance between k-means clusters in the spatial network. |
| `umap_kmeans_cluster.png` | UMAP embedding of cells coloured by k-means cluster. |
| `umap_leiden.png` | UMAP embedding coloured by Leiden cluster. |
| `cell_density.png` | Per-cluster and per-community spatial density maps. |
| `cluster_intensity_heatmap.png` | Mean (z-scored) channel intensity per k-means cluster. |
| `morphology_by_cluster.png` | Cell area by cluster (only present if an `area` column exists). |

## QuPath GeoJSON

| File | What it is |
|---|---|
| `qupath_export.geojson` (or `--geojson-path`) | Cell boundaries + spots for QuPath. |

The GeoJSON contains:

- **Cell boundaries** as annotations, coloured by community, with cluster and
  per-marker mean-intensity measurements attached.
- **RNA spots** as detection objects (only when `--detect-blobs` was used),
  linked to their parent cell where assignment succeeded.

Open it in QuPath via **File → Import → GeoJSON**.

---

Next: [Visualising results](visualization.md) to open these in a viewer.
