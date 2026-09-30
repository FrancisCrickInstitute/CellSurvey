"""Cell segmentation: nucleus expansion into whole-cell boundaries.

Step 1 of the cell-segmentation plan: read the physical pixel size (µm/px) from
the input image's OME metadata so that cell expansion can be specified in microns
rather than pixels.
"""

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


def expand_nuclei(nuclei_gdf, expansion_um, pixel_size_um):
    """Expand nucleus polygons into approximate whole-cell boundaries (v1).

    Buffers each nucleus geometry outward by ``expansion_um`` (converted to
    pixels via ``pixel_size_um``). This is the simple, overlap-permitting v1 of
    the cell-expansion plan; upgrade to non-overlapping expansion later.

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
        A copy of ``nuclei_gdf`` with each geometry buffered outward by
        ``expansion_um / pixel_size_um`` pixels. Index and columns are preserved.
    """
    expansion_px = expansion_um / pixel_size_um
    cell_gdf = nuclei_gdf.copy()
    cell_gdf["geometry"] = nuclei_gdf.geometry.buffer(expansion_px)
    return cell_gdf
