"""Cell segmentation: segmenter runners and nucleus/cell boundary handling.

Two segmenters are supported, each in its own pixi environment:

- **Stardist** (``default`` env, TensorFlow): nuclei only, written to
  ``stardist_boundaries``. Whole cells are approximated by expanding nuclei via
  a marker-controlled watershed (QuPath-style, non-overlapping).
- **InstanSeg** (``instanseg`` env, PyTorch): nuclei and whole cells, written to
  ``nucleus_boundaries`` and ``cell_boundaries``.

TensorFlow and InstanSeg are imported lazily inside their runners so that each
environment only needs its own framework.
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import shapely
import sopa
from scipy.ndimage import distance_transform_edt
from skimage.segmentation import watershed
from skimage.draw import polygon as _draw_polygon
from skimage.measure import find_contours, regionprops

from bioio import BioImage

from cellsurvey.utils import STARDIST_NUCLEI_KEY, NUCLEI_KEY, CELLS_KEY


def run_stardist(sdata, channel, use_gpu=False):
    """Segment nuclei with Stardist (``2D_versatile_fluo``) into ``stardist_boundaries``."""
    try:
        import tensorflow as tf
    except ImportError as e:
        raise ImportError("Stardist needs TensorFlow, which is only in the default pixi "
                          "environment: `pixi run python run.py ...`") from e

    gpus = tf.config.list_physical_devices('GPU')
    if gpus:
        print(f"Found {len(gpus)} GPU(s), using GPU backend for Stardist")
    elif use_gpu:
        print("No GPU detected but --use-gpu specified, attempting GPU backend")
    else:
        print("WARNING: No GPU detected. Stardist segmentation will run on CPU and may be very slow.")
    sopa.settings.parallelization_backend = None

    print(f"Run stardist on channel {channel}...")
    sopa.segmentation.stardist(sdata, model_type='2D_versatile_fluo', channels=channel)


def run_instanseg(sdata, channels, pixel_size_um, model_type='fluorescence_nuclei_and_cells',
                  use_gpu=False):
    """Segment nuclei and whole cells with InstanSeg.

    Sopa's wrapper returns one mask per call, so the model is run twice
    (``target='nuclei'`` then ``target='cells'``), writing ``nucleus_boundaries``
    and ``cell_boundaries``. Each cell is then linked to the nucleus it contains
    via a ``nucleus_id`` column (see :func:`link_nuclei_to_cells`).
    """
    if not hasattr(sopa.segmentation, 'instanseg'):
        raise ImportError("InstanSeg needs sopa>=2.3.0 and instanseg-torch, which are only in the "
                          "instanseg pixi environment: `pixi run -e instanseg python run.py ...`")
    try:
        import torch
        has_cuda = torch.cuda.is_available()
    except ImportError as e:
        raise ImportError("InstanSeg needs PyTorch: `pixi run -e instanseg python run.py ...`") from e

    if has_cuda:
        print("CUDA available, running InstanSeg on GPU")
    elif use_gpu:
        print("WARNING: --use-gpu specified but CUDA is not available to PyTorch; running on CPU")
    else:
        print("WARNING: No GPU detected. InstanSeg will run on CPU and may be slow.")
    sopa.settings.parallelization_backend = None

    if pixel_size_um is None:
        print("WARNING: no pixel size for InstanSeg; the image is not rescaled to the model's "
              "training resolution and it may find few or no cells. Pass --pixel-size.")

    common = dict(model_type=model_type, channels=channels, pixel_size=pixel_size_um,
                  device='cuda' if has_cuda else 'cpu')
    print("Run InstanSeg (nuclei)...")
    sopa.segmentation.instanseg(sdata, target='nuclei', key_added=NUCLEI_KEY, **common)
    print("Run InstanSeg (cells)...")
    sopa.segmentation.instanseg(sdata, target='cells', key_added=CELLS_KEY, **common)

    linked = link_nuclei_to_cells(sdata.shapes[NUCLEI_KEY], sdata.shapes[CELLS_KEY])
    sdata.shapes[CELLS_KEY]['nucleus_id'] = linked['nucleus_id'].values


def link_nuclei_to_cells(nuclei_gdf, cells_gdf):
    """Return a copy of ``cells_gdf`` with a ``nucleus_id`` column.

    Each nucleus is assigned to the cell containing its representative point.
    A cell containing several nuclei keeps the largest; a cell with none gets
    NaN. Nuclei and cells come from separate InstanSeg passes, so their ids do
    not correspond and must be matched spatially.
    """
    points = gpd.GeoDataFrame(
        {'_nucleus_id': nuclei_gdf.index, '_area': nuclei_gdf.geometry.area.values},
        geometry=nuclei_gdf.geometry.representative_point().values, crs=nuclei_gdf.crs)
    cells = gpd.GeoDataFrame({'_cell_id': cells_gdf.index},
                             geometry=cells_gdf.geometry.values, crs=cells_gdf.crs)
    joined = gpd.sjoin(points, cells, how='inner', predicate='within')
    primary = (joined.sort_values('_area', ascending=False)
               .drop_duplicates(subset='_cell_id')
               .set_index('_cell_id')['_nucleus_id'])

    result = cells_gdf.copy()
    result['nucleus_id'] = primary.reindex(result.index).values
    n_linked = int(result['nucleus_id'].notna().sum())
    print(f"Linked {n_linked} of {len(result)} cells to a nucleus "
          f"({len(nuclei_gdf) - n_linked} nuclei unlinked)")
    return result


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


def _watershed_expand(nuclei_gdf, expansion_px):
    """Rasterize nuclei and expand them via a marker-controlled watershed.

    Mirrors QuPath's ``WatershedCellDetection``: each nucleus is a watershed
    marker, the elevation is the distance to the nearest nucleus, and the flood is
    capped at ``expansion_px`` so cells grow outward but stop where they meet.
    """
    minx, miny, maxx, maxy = nuclei_gdf.total_bounds
    pad = expansion_px + 1.0
    ox, oy = minx - pad, miny - pad
    height = int(np.ceil(maxy - miny + 2 * pad)) + 1
    width = int(np.ceil(maxx - minx + 2 * pad)) + 1

    # Rasterize nuclei into a labeled array (1..n; 0 = background).
    labeled = np.zeros((height, width), dtype=np.int32)
    for i, geom in enumerate(nuclei_gdf.geometry, start=1):
        if geom is None or geom.is_empty:
            continue
        parts = geom.geoms if geom.geom_type == "MultiPolygon" else [geom]
        for part in parts:
            xs, ys = np.asarray(part.exterior.coords.xy)
            rr, cc = _draw_polygon(ys - oy, xs - ox, shape=(height, width))
            labeled[rr, cc] = i

    # Distance to the nearest nucleus (0 inside a nucleus, growing outward).
    dist = distance_transform_edt(labeled == 0)

    # Marker-controlled watershed, capped at expansion_px from any nucleus.
    cells = watershed(dist, markers=labeled, mask=(dist <= expansion_px), connectivity=1)

    # Vectorize back to polygons, preserving the cell id (label -> index).
    label_to_cell = {i: cell_id for i, cell_id in enumerate(nuclei_gdf.index, start=1)}

    geometries = []
    cell_ids = []
    for region in regionprops(cells):
        label_id = region.label
        minr, minc = region.bbox[0], region.bbox[1]
        # Pad with a background border so find_contours also works when the
        # region exactly fills its bounding box (no interior background).
        padded = np.pad(region.image, 1, constant_values=False)
        contours = find_contours(padded, level=0.5)
        if not contours:
            continue
        contour = max(contours, key=len)
        coords = [(ox + minc + c - 1, oy + minr + r - 1) for r, c in contour]
        if len(coords) < 4:
            continue
        poly = shapely.Polygon(coords)
        if not poly.is_valid:
            poly = poly.buffer(0)
        if poly.is_empty:
            continue
        geometries.append(poly)
        cell_ids.append(label_to_cell[label_id])

    result = nuclei_gdf.loc[cell_ids].copy()
    result.geometry = gpd.GeoSeries(geometries, index=cell_ids, crs=nuclei_gdf.crs)

    # Fallback: any nucleus that produced no watershed cell gets a buffered copy,
    # so the cell count always matches the input nuclei.
    missing = nuclei_gdf.index.difference(result.index)
    if len(missing) > 0:
        print(f"WARNING: {len(missing)} nuclei produced no watershed cell; using buffered fallback")
        fallback = nuclei_gdf.loc[missing].copy()
        fallback.geometry = nuclei_gdf.loc[missing].geometry.buffer(expansion_px)
        result = pd.concat([result, fallback])

    # Record the original nucleus id so it survives sopa.aggregate's re-indexing
    # (which assigns new cell ids to the shapes).
    result["nucleus_id"] = result.index

    return result


def expand_nuclei(nuclei_gdf, expansion_um, pixel_size_um):
    """Expand nucleus polygons into non-overlapping whole-cell boundaries.

    Expands each nucleus outward by ``expansion_um`` (converted to pixels via
    ``pixel_size_um``) using a marker-controlled watershed, so adjacent cells
    stop where they meet (QuPath-style).

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
        A copy of ``nuclei_gdf`` with each geometry expanded into a
        non-overlapping cell boundary. Index and columns are preserved.
    """
    expansion_px = expansion_um / pixel_size_um
    return _watershed_expand(nuclei_gdf, expansion_px)
