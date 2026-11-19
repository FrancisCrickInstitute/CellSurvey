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
    the nuclei layer (nuclei-only runs: ``--cell-expansion 0``, or a ``_seg.zarr`` written before whole-cell support). ``nuclei_key``
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


def store_table_meta(sdata, table_key='table', **values):
    """Record ``values`` in ``sdata.tables[table_key].uns['cellsurvey']``."""
    uns = sdata.tables[table_key].uns
    meta = dict(uns.get(_UNS_KEY, {}))
    meta.update(values)
    uns[_UNS_KEY] = meta


def store_pixel_size(sdata, pixel_size_um, table_key='table'):
    """Record ``pixel_size_um`` (µm/px) in ``sdata.tables[table_key].uns``."""
    store_table_meta(sdata, table_key, pixel_size_um=float(pixel_size_um))


def resolve_pixel_size(sdata, pixel_size_um=None, table_key='table'):
    """Return ``pixel_size_um`` if given, else the value stored by :func:`store_pixel_size`.

    Raises ``ValueError`` if neither is available, or if the value is not a
    positive number. Cell coordinates are in pixels, so µm distance parameters
    cannot be applied without it.
    """
    if pixel_size_um is not None:
        value, source = pixel_size_um, "--pixel-size"
    else:
        stored = None
        if table_key in sdata.tables:
            stored = sdata.tables[table_key].uns.get(_UNS_KEY, {}).get('pixel_size_um')
        if stored is None:
            raise ValueError("No pixel size available: this _seg.zarr predates stored pixel sizes. "
                             "Pass the image's pixel size (µm/px), e.g. --pixel-size 0.28.")
        value, source = stored, "stored pixel size"
    value = float(value)
    # `not value > 0` also rejects NaN.
    if not value > 0:
        raise ValueError(f"{source} must be a positive number of µm/px, got {value}.")
    return value


def remove_channel_suffix(name):
    return re.sub(r'_ch_\d+$', '', name)


# What clustering, the network edge weights and UMAP/Leiden are based on.
CLUSTER_ON_CHOICES = ('intensities', 'spots', 'both')
# Suffix for spot-count features, whose names would otherwise match the
# intensity channel they were detected in.
SPOTS_SUFFIX = ' (spots)'
# sopa.aggregate puts channel intensities here when it also counts spots.
_INTENSITIES_OBSM = 'intensities'


def _clean_channel_columns(df):
    df = df.copy()
    df.columns = [remove_channel_suffix(str(c)) for c in df.columns]
    return df.loc[:, ~df.columns.duplicated(keep='first')]


def extract_feature_tables(table):
    """Return ``(intensity_df, counts_df)`` from an aggregated AnnData table.

    Without spots, ``sopa.aggregate`` puts the per-cell mean channel intensities
    in ``.X``. With spots (``aggregate_genes=True``), ``.X`` holds the spot
    counts per cell and the intensities move to ``obsm['intensities']``.
    ``counts_df`` is ``None`` when there are no spot counts. Both are indexed by
    cell (``table.obs.index``), with ``_ch_<i>`` suffixes removed and duplicate
    channel names dropped.
    """
    X = table.X.toarray() if hasattr(table.X, 'toarray') else np.asarray(table.X)
    x_df = pd.DataFrame(X, index=table.obs.index, columns=table.var.index)
    if _INTENSITIES_OBSM not in table.obsm:
        return _clean_channel_columns(x_df), None
    intensities = table.obsm[_INTENSITIES_OBSM]
    if not isinstance(intensities, pd.DataFrame):
        intensities = pd.DataFrame(np.asarray(intensities))
    intensities = pd.DataFrame(intensities.to_numpy(), index=table.obs.index, columns=intensities.columns)
    return _clean_channel_columns(intensities), _clean_channel_columns(x_df)


def _zscore(df):
    """Z-score each column across cells; constant columns become 0 (as StandardScaler)."""
    std = df.std(ddof=0).replace(0, np.nan)
    return ((df - df.mean()) / std).fillna(0)


def build_clustering_features(intensity_df, counts_df, mode, cell_area_um2=None, weight_groups=True):
    """Return the scaled per-cell feature matrix used for clustering, the network and UMAP.

    - ``intensities``: each channel's mean intensity, z-scored across cells.
    - ``spots``: spot density per channel (spots per 100 µm² of cell area, so
      clusters don't just track cell size), ``log1p``, then z-scored. Columns
      get the ``' (spots)'`` suffix.
    - ``both``: both groups; with ``weight_groups``, each group is divided by
      √(its number of columns) so the intensity channels and the spot channels
      contribute equally, however many each has.

    ``cell_area_um2`` (a Series aligned to ``counts_df``) is needed for the spot
    features. Raises ``ValueError`` if a spot mode is requested without counts.
    """
    if mode not in CLUSTER_ON_CHOICES:
        raise ValueError(f"Unknown --cluster-on mode {mode!r}; expected one of {CLUSTER_ON_CHOICES}")
    groups = []
    if mode in ('intensities', 'both'):
        groups.append(_zscore(intensity_df))
    if mode in ('spots', 'both'):
        if counts_df is None:
            raise ValueError(f"--cluster-on {mode} needs spot counts, but this table has none "
                             f"(run with --detect-blobs)")
        if cell_area_um2 is None:
            raise ValueError("Spot features need the cell area (µm²)")
        area = pd.Series(cell_area_um2, index=counts_df.index).where(lambda a: a > 0)
        density = counts_df.div(area, axis=0).mul(100).fillna(0)   # spots per 100 µm²
        spots = _zscore(np.log1p(density))
        spots.columns = [f"{c}{SPOTS_SUFFIX}" for c in spots.columns]
        groups.append(spots)
    if len(groups) == 1:
        return groups[0]
    if weight_groups:
        groups = [g / np.sqrt(g.shape[1]) for g in groups]
    return pd.concat(groups, axis=1)


def cluster_data(data, n_clusters=10, random_seed=42, scale=True):
    """k-means on ``data``. With ``scale=False``, ``data`` is used as given
    (e.g. already scaled by :func:`build_clustering_features`)."""
    print("\nPerforming k-means clustering...")

    scaled_data = StandardScaler().fit_transform(data) if scale else np.asarray(data)

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
