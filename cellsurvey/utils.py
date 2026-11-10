import re
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans


# Shape layer names. Stardist writes nuclei to ``stardist_boundaries``;
# InstanSeg writes nuclei to ``nucleus_boundaries``. Whole cells always live in
# ``cell_boundaries`` (InstanSeg cells, or Stardist nuclei expanded by watershed).
STARDIST_NUCLEI_KEY = 'stardist_boundaries'
NUCLEI_KEY = 'nucleus_boundaries'
CELLS_KEY = 'cell_boundaries'
_NUCLEI_KEYS = (NUCLEI_KEY, STARDIST_NUCLEI_KEY)


def resolve_shape_keys(sdata):
    """Return ``(nuclei_key, cells_key)`` for the segmentation layers in ``sdata``.

    ``cells_key`` is ``cell_boundaries`` when present, otherwise it falls back to
    the nuclei layer (nuclei-only runs: ``--cell-expansion 0``, missing pixel
    size, or a ``_seg.zarr`` written before whole-cell support). ``nuclei_key``
    is ``None`` if there is no nuclei layer.
    """
    nuclei_key = next((k for k in _NUCLEI_KEYS if k in sdata.shapes), None)
    cells_key = CELLS_KEY if CELLS_KEY in sdata.shapes else nuclei_key
    if cells_key is None:
        raise KeyError(f"No segmentation boundaries found; expected one of "
                       f"{[CELLS_KEY, *_NUCLEI_KEYS]}, got {list(sdata.shapes.keys())}")
    return nuclei_key, cells_key


# Pixel size is stored in the table's .uns so read-only consumers of _seg.zarr
# (e.g. the stability sweep) can convert µm parameters to pixel coordinates.
_UNS_KEY = 'cellsurvey'


def store_pixel_size(sdata, pixel_size_um, table_key='table'):
    """Record ``pixel_size_um`` (µm/px) in ``sdata.tables[table_key].uns``."""
    uns = sdata.tables[table_key].uns
    meta = dict(uns.get(_UNS_KEY, {}))
    meta['pixel_size_um'] = float(pixel_size_um)
    uns[_UNS_KEY] = meta


def resolve_pixel_size(sdata, pixel_size_um=None, table_key='table'):
    """Return ``pixel_size_um`` if given, else the value stored by :func:`store_pixel_size`.

    Raises ``ValueError`` if neither is available. Cell coordinates are in
    pixels, so µm distance parameters cannot be applied without it.
    """
    if pixel_size_um is not None:
        return float(pixel_size_um)
    stored = None
    if table_key in sdata.tables:
        stored = sdata.tables[table_key].uns.get(_UNS_KEY, {}).get('pixel_size_um')
    if stored is None:
        raise ValueError("No pixel size available: this _seg.zarr predates stored pixel sizes. "
                         "Pass the image's pixel size (µm/px), e.g. --pixel-size 0.28.")
    return float(stored)


def remove_channel_suffix(name):
    return re.sub(r'_ch_\d+$', '', name)


def cluster_data(data, n_clusters=10, random_seed=42):
    print("\nPerforming k-means clustering...")

    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(data)

    kmeans = KMeans(n_clusters=n_clusters, random_state=random_seed)
    cluster_labels = kmeans.fit_predict(scaled_data)

    distances = kmeans.transform(scaled_data)

    print(f"Clustering complete. Found {n_clusters} clusters.")
    print("Cluster distribution:")
    for i in range(n_clusters):
        count = np.sum(cluster_labels == i)
        print(f"  Cluster {i}: {count} cells ({100 * count / len(cluster_labels):.1f}%)")

    return cluster_labels, distances


def assign_spots_to_cells(spatial_data, spots_key='spots', cell_boundaries=None):
    if spots_key not in spatial_data.points:
        print("No spots found, skipping spot-to-cell assignment")
        return None

    if cell_boundaries is None:
        cell_boundaries = resolve_shape_keys(spatial_data)[1]

    print("Computing spot-to-cell assignments...")

    spots_df = spatial_data.points[spots_key].compute()
    spots_gdf = gpd.GeoDataFrame(
        spots_df,
        geometry=gpd.points_from_xy(spots_df['x'], spots_df['y']),
        crs=spatial_data.shapes[cell_boundaries].crs
    )

    spots_with_cells = gpd.sjoin(
        spots_gdf,
        spatial_data.shapes[cell_boundaries][['geometry']],
        how='left',
        predicate='within'
    )

    spots_with_cells['cell_id'] = spots_with_cells['cell_id'].where(
        spots_with_cells['cell_id'].notna(), other=None
    )

    assigned = spots_with_cells['cell_id'].notna().sum()
    total = len(spots_with_cells)
    print(f"Assigned {assigned} of {total} spots to cells ({100 * assigned / total:.1f}%)")
    print(f"Unassigned spots: {total - assigned} ({100 * (total - assigned) / total:.1f}%)")

    return spots_with_cells


def get_colors_for_communities(n_communities):
    if n_communities <= 10:
        cmap = plt.get_cmap('tab10')
    elif n_communities <= 20:
        cmap = plt.get_cmap('tab20')
    else:
        cmap = plt.get_cmap('hsv')

    colors = []
    for i in range(n_communities):
        rgba = cmap(i / max(n_communities - 1, 1))
        rgb = [int(rgba[0] * 255), int(rgba[1] * 255), int(rgba[2] * 255)]
        colors.append(rgb)

    return colors
