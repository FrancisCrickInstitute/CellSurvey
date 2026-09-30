"""Cell segmentation: nucleus expansion into whole-cell boundaries.

Reads the physical pixel size (µm/px) from OME metadata and expands nuclei into
approximate whole-cell boundaries (QuPath-style, non-overlapping).
"""

import numpy as np
import geopandas as gpd
import shapely
from shapely.geometry import MultiPoint

from bioio import BioImage


def get_pixel_size(imagepath):
    """Return the physical pixel size (µm/px) of an image.

    Reads the in-plane (Y, X) physical pixel size from the OME-TIFF metadata via
    bioio. Assumes isotropic pixels; if Y and X differ materially, a warning is
    printed and the X value is used.

    Parameters
    ----------
    imagepath : str or Path
        Path to the input image.

    Returns
    -------
    float or None
        The pixel size in µm/px, or ``None`` if the image carries no physical
        size metadata (bioio returns ``PhysicalPixelSizes(None, None, None)``).
    """
    pps = BioImage(imagepath).physical_pixel_sizes
    y, x = pps.Y, pps.X

    # Missing metadata → None for both in-plane axes.
    if x is None and y is None:
        return None

    # Take whichever in-plane axis is present.
    if x is None:
        x = y
    if y is None:
        y = x

    # Treat zero/non-positive as missing (some readers may report 0 instead of None).
    if x <= 0 or y <= 0:
        return None

    # Assume isotropic pixels; warn if the axes differ by more than 1%.
    if abs(y - x) > 0.01 * max(y, x):
        print(f"WARNING: anisotropic pixels (Y={y:.4f}, X={x:.4f} µm/px); using X={x:.4f}")

    return float(x)


def _resolve_overlap(cell_gdf):
    """Clip overlapping polygons to their Voronoi cells (non-overlapping).

    A robust replacement for ``sopa.shapes.expand.remove_overlap``, which empties
    most cells on dense data. Each polygon is intersected with the Voronoi cell
    of its centroid, so adjacent cells stop where they meet (QuPath-style).
    """
    centroids = cell_gdf.geometry.centroid
    coords = np.column_stack([centroids.x.values, centroids.y.values])

    bbox = shapely.box(*cell_gdf.total_bounds)
    voronoi = shapely.voronoi_polygons(MultiPoint(coords), extend_to=bbox)
    tree = shapely.STRtree(list(voronoi.geoms))

    new_geometries = []
    for i in range(len(cell_gdf)):
        centroid = centroids.iloc[i]
        matches = tree.query(centroid, predicate="contains")
        if len(matches) == 0:
            matches = tree.query(centroid, predicate="intersects")
        region = voronoi.geoms[int(matches[0])]
        new_geometries.append(cell_gdf.geometry.iloc[i].intersection(region))

    result = cell_gdf.copy()
    result.geometry = new_geometries
    return result


def expand_nuclei(nuclei_gdf, expansion_um, pixel_size_um):
    """Expand nucleus polygons into non-overlapping whole-cell boundaries (v2).

    Buffers each nucleus geometry outward by ``expansion_um`` (converted to
    pixels via ``pixel_size_um``), then clips each to its Voronoi cell so
    adjacent cells stop where they meet (QuPath-style).

    Parameters
    ----------
    nuclei_gdf : geopandas.GeoDataFrame
        Nucleus polygons (e.g. ``sdata.shapes['stardist_boundaries']``). The
        index is preserved and treated as the cell id.
    expansion_um : float
        Expansion radius in microns.
    pixel_size_um : float
        Physical pixel size in µm/px (from :func:`get_pixel_size`).

    Returns
    -------
    geopandas.GeoDataFrame
        A copy of ``nuclei_gdf`` with each geometry buffered outward and overlaps
        removed. Index and columns are preserved.
    """
    expansion_px = expansion_um / pixel_size_um
    cell_gdf = nuclei_gdf.copy()
    cell_gdf["geometry"] = nuclei_gdf.geometry.buffer(expansion_px)
    return _resolve_overlap(cell_gdf)
