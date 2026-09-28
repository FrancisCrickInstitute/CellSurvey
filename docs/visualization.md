# Visualising results

The segmented Zarr store (`*_seg.zarr`) is the main result, and the QuPath
GeoJSON is the companion export. This page covers the viewers you can open them
with, in order of recommendation.

## 1. Odon (recommended)

[Odon](https://github.com/alexcoulton/odon) is a lightweight, GPU-accelerated
viewer for SpatialData Zarr files. It supports multiscale image viewing with
shape overlays coloured by metadata columns.

1. [Install Odon](https://github.com/alexcoulton/odon#installation).
2. Open the segmented Zarr file (`*_seg.zarr`).
3. The `stardist_boundaries` shape layer can be coloured by cluster or community
   using the **"Color by"** dropdown.

!!! note "Which columns appear"
    Only categorical columns with ≤24 distinct values appear in the dropdown.
    `kmeans_cluster` (10 values) and `community` (up to 24) will be visible.
    Channel intensity columns are available as continuous properties but are not
    shown as colour-by options.

## 2. QuPath

For traditional pathology workflows, open the exported GeoJSON in
[QuPath](https://qupath.github.io/):

1. **File → Import → GeoJSON**, then select `qupath_export.geojson`.
2. Cell boundaries appear as annotations coloured by community.
3. RNA spots (when `--detect-blobs` was used) appear as detection objects.

Each cell annotation carries its cluster, community, and per-marker mean
intensity as measurements.

## 3. TissUUmaps

[TissUUmaps](https://github.com/TissUUmaps/TissUUmaps4) is a GPU-accelerated,
browser-based viewer. It runs entirely in-browser — no installation required —
and supports SpatialData Zarr files natively via the OME-Zarr + SpatialData
plugin.

1. Open [TissUUmaps live](https://tissuumaps.github.io/TissUUmaps4/live/) or
   download the [latest release](https://github.com/TissUUmaps/TissUUmaps4/releases).
2. Load the segmented Zarr file (`*_seg.zarr`).
3. Cell boundaries appear as a shapes layer; link to metadata columns for
   colour-by-cluster or colour-by-community.

!!! note "TissUUmaps version"
    TissUUmaps 4 is under active development. The stable release (v3) may not
    include SpatialData Zarr support — use the v4 development builds. Shapes
    require matching coordinate systems and an ID column for metadata linkage. A
    modern browser with WebGL 2 and the File System API is required.

## 4. napari + napari-spatialdata

[napari](https://napari.org/) is a multi-dimensional image viewer for Python;
[napari-spatialdata](https://github.com/scverse/napari-spatialdata) adds native
SpatialData Zarr support.

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

1. Select a coordinate system.
2. Click `stardist_boundaries` to load the cell shapes.
3. Use the **View** widget (Plugins → napari-spatialdata → View) to colour cells
   by cluster, community, or any channel intensity column. Double-click any
   `obs` column to apply it as the face colour.

### Performance notes for large datasets

- Shape loading is slow above ~50K polygons due to triangulation. Use the
  `bermuda` backend for faster loading:
  `pip install "napari-spatialdata[all,bermuda]"`.
- Polygons are simplified when the shape count exceeds 100 (configurable via
  `napari_spatialdata.constants.config.POLYGON_THRESHOLD`).
- If boundaries appear too simplified at high zoom, increase the threshold:

  ```python
  from napari_spatialdata.constants import config
  config.POLYGON_THRESHOLD = 50000
  ```

---

Next: [FAQ](faq.md) for common questions and troubleshooting.
