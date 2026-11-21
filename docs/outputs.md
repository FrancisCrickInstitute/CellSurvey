# Outputs

A CellSurvey run produces a set of files organised around two Zarr stores, a
GeoJSON file, and a directory of plots. This page lists every artifact, what it
contains, and how to open it.

!!! tip "Where to look"
    - Zarr stores go next to `--output_file` (the segmented store replaces
      `.zarr` with `_seg.zarr`).
    - Plots and `summary.json` go in `--plot_dir`.
    - The GeoJSON goes next to the Zarr stores as `<output>.geojson` (or to
      `--geojson-path`).

## Zarr stores

| File | What it is |
|---|---|
| `<output>.zarr` | The image converted to SpatialData format (the initial checkpoint, written before segmentation). |
| `<output>_seg.zarr` | The **main result** — cells, boundaries, and per-cell measurements. |

Both are [SpatialData](https://www.nature.com/articles/s41592-024-02212-x) Zarr
stores. The segmented store contains:

- `cell_boundaries` — whole-cell polygons, with `kmeans_cluster` and
  `community` labels and a `nucleus_id` column linking each cell to its nucleus.
  Cells in communities smaller than `--min-community-size` have community `-1`
  (unassigned).
  With Stardist these are nuclei expanded by `--cell-expansion`; with InstanSeg
  they are segmented directly.
- A nuclei layer — `stardist_boundaries` (Stardist) or `nucleus_boundaries`
  (InstanSeg) — with the same `kmeans_cluster` and `community` labels, copied
  from each nucleus's cell. Nuclei not linked to a cell get `-1`.
- `table` — an AnnData table of cells, measured over `cell_boundaries`.
  Without spots, `.X` holds the per-cell mean intensities (cells × markers).
  With `--detect-blobs`, `.X` holds the **spot counts** per spot channel and the
  intensities are in `obsm['intensities']`. That is how Sopa stores them.
- `spots` — detected RNA spots (only if `--detect-blobs` was used).

!!! note "Nuclei-only runs"
    When there is no `cell_boundaries` layer (Stardist with
    `--cell-expansion 0`, or a Zarr from an older nuclei-only run), the nuclei layer is the
    cell layer: the table, labels, network, and GeoJSON all use it instead.

Open them with any SpatialData-aware viewer (see
[Visualising results](visualization.md)).

## `summary.json`

Written to `--plot_dir`. A small JSON summary of the network analysis:

```json
{
  "n_cells": 12345,
  "n_clusters": 10,
  "n_communities": 7,
  "n_communities_raw": 31,
  "min_community_size": 10,
  "n_unassigned_cells": 58,
  "n_edges": 36950,
  "max_edge_distance_um": 50,
  "max_edge_distance_px": 178.57,
  "pixel_size_um": 0.28,
  "community_resolution": 0.1,
  "cluster_on": "intensities",
  "n_features": 21,
  "cluster_sizes": { "0": 1204, "1": 980, "...": "..." },
  "community_sizes": { "0": 2100, "1": 1550, "...": "..." }
}
```

Useful for a quick sanity check of how many cells, clusters, and communities
were found. `n_communities` counts communities of at least
`--min-community-size` cells; `n_communities_raw` is Louvain's count before
that minimum, and `n_unassigned_cells` the cells in the smaller communities,
which are labelled `-1`. Communities are numbered by size (0 = largest), so
`community_sizes` lists them largest first. `n_edges` counts unique Delaunay edges after the distance filter (at
most about 3× `n_cells`). The pixel size is also stored in the segmented Zarr
(`table.uns['cellsurvey']['pixel_size_um']`) so the stability sweep can convert
µm distances.

## Plots

All plots are PNG files written to `--plot_dir`.

| File | What it shows |
|---|---|
| `cell_type_to_cell_type.png` | Mean hop distance between k-means clusters in the spatial network. |
| `umap_kmeans_cluster.png` | UMAP embedding of cells coloured by k-means cluster. |
| `umap_leiden.png` | UMAP embedding coloured by Leiden cluster. |
| `cell_density.png` | Spatial density maps: one panel per k-means cluster, then one for each of the 20 largest communities. The title says how many smaller communities and unassigned cells are not shown. Every panel spans the whole tissue and is oriented like the image (y increasing downwards, as in QuPath). |
| `cluster_intensity_heatmap.png` | The clustering features (`--cluster-on`) per k-means cluster: channel intensities and/or log spot densities (labelled `(spots)`), each z-scored across all cells (0 = average cell; red = higher, blue = lower). The colour scale is capped at ±3, so more extreme values show as the darkest colour. |
| `morphology_by_cluster.png` | Cell area (µm²) by cluster (only present if an `area` column exists). |

## Overview report

`report.html` in `--plot_dir` is a single, self-contained web page summarising
the run for a biologist. It has all figures embedded, so it can be emailed,
opened in any browser, or printed to PDF. It contains:

- **Header**: the source image (file path, or OMERO server and Image ID),
  segmentation method, pixel size, number of cells, which channels count as
  markers (nuclear and control channels such as `DAPI`, `TRITC`, `Cy5` or `AF`
  are recognised by name), the clustering and community parameters, and the
  CellSurvey code version.
- **How to read this**: z-scores, and the difference between *clusters* (what
  cells express) and *communities* (where cells sit).
- **Cell populations**: a heatmap of each k-means cluster's average marker
  levels, and a table of cluster sizes and distinguishing markers. Clusters get
  automatic names from the markers they have raised (e.g. `GFAP`,
  `PSYN/PTAU`), `low X`, or `no distinct marker`. These describe marker
  profiles, not confirmed cell types.
- **Where the clusters are**: each cluster drawn over the tissue image (the
  nuclear channel), and one map coloured by the most common cluster in each
  area.
- **Communities**: outlines of the largest communities over the tissue,
  their composition by cluster, a short generated summary, and the resolution
  sweep if one was run for this result (see [Parameters](parameters.md)).
- **Data-quality flags**: cells that are bright in many unrelated markers at
  once ("non-specific": likely autofluorescence or debris), clusters made
  mostly of them, markers whose "positive" cells are largely non-specific, and
  control channels.

Marker levels in the report use a robust z-score (arcsinh, then median and
interquartile range), which resists outliers. It is for display and differs
slightly from the z-score used for clustering. PNG copies of the figures are
saved in `--plot_dir/report/`.

The report is written at the end of every run (`--no-report` skips it). To
re-make it from an existing result without re-running anything:

```bash
pixi run python -m cellsurvey.report --zarr <output>_seg.zarr -p <plot_dir>
```

## QuPath GeoJSON

| File | What it is |
|---|---|
| `<output>.geojson` (or `--geojson-path`) | Cell boundaries + spots for QuPath. |

The GeoJSON contains:

- **Cell boundaries** (from `cell_boundaries`, or the nuclei layer in a
  nuclei-only run) as annotations, coloured by community (unassigned cells are
  grey, classified `Unassigned`), with cluster and
  per-marker mean-intensity measurements (`Cell: <channel> mean`) attached, and,
  with `--detect-blobs`, the raw spot count per spot channel
  (`Cell: <channel> spot count`).
- **RNA spots** as detection objects (only when `--detect-blobs` was used),
  linked to their parent cell where assignment succeeded.

Open it in QuPath via **File → Import → GeoJSON**.

---

Next: [Visualising results](visualization.md) to open these in a viewer.
