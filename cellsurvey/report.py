"""Overview report: one self-contained HTML page summarising a CellSurvey run.

Reads only the segmented Zarr (``<prefix>_seg.zarr``), like the stability
sweep, and writes ``<plot_dir>/report.html`` (figures embedded as base64 PNG,
so the file can be emailed or printed) plus the PNGs in ``<plot_dir>/report/``.

Sections (v1, "core"): header (source, parameters, counts); how to read the
report; cluster profiles with automatic names; spatial maps over the tissue
image; communities (outlines, composition, narrative, resolution sweep); and
data-quality flags (non-specific cells, flagged clusters, suspect markers,
control channels). Marker histograms, per-marker "positive" maps and a spots
section are planned (see AGENTS.md, "Planned: Overview report").

Display scale: each marker is shown as a robust z-score of arcsinh intensity
(median/IQR), following the cellsurvey-analysis overview. This is for display
only and differs from the plain z-score k-means uses (see
:func:`cellsurvey.utils.build_clustering_features`).

Run standalone (re-generate without re-running the pipeline):
    python -m cellsurvey.report --zarr <prefix>_seg.zarr -p <plot_dir>

Must not import TensorFlow, PyTorch or omero-py.
"""

import argparse
import base64
import datetime
import html
import io
import json
import math
import os
import re
import subprocess

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import patheffects
from scipy.ndimage import gaussian_filter

from cellsurvey.utils import extract_feature_tables, remove_channel_suffix

# --- Fixed report settings (stated in the report header) ---------------------
POSITIVE_Z = 2.0                  # illustrative "positive" cutoff on the robust z scale
NONSPECIFIC_MIN_MARKERS = 3       # non-specific: z > POSITIVE_Z in at least max(this,
NONSPECIFIC_MIN_FRACTION = 0.25   #   ceil(this × number of markers)) markers
FLAG_CLUSTER_FRACTION = 0.5       # flag clusters that are more than this non-specific
MARKER_SUSPECT = 0.25             # share of a marker's positives that are non-specific
MARKER_DOMINATED = 0.5
MIN_POSITIVES = 100               # fewer positive cells than this: too few to judge a marker
PURE_COMMUNITY = 0.8              # "mostly one cluster" in the community narrative
MAX_COMMUNITY_PANELS = 20         # communities drawn in the outline / composition figures
MAX_BACKGROUND_SIDE = 4000        # pick the largest pyramid level no bigger than this
MAX_POINTS_PER_PANEL = 60_000     # subsample scatter panels (rasterised)
MAX_MARKER_PANELS = 30            # per-marker "positive" maps

# Channel-name tokens (first word of the name, upper-case).
NUCLEAR_TOKENS = {'DAPI', 'HOECHST', 'DRAQ5', 'SYTOX', 'NUCLEI', 'NUCLEAR'}
CONTROL_TOKENS = {'TRITC', 'CY3', 'CY5', 'CY7', 'FITC', 'AF', 'AUTOFLUORESCENCE', 'AUTOFLUOR',
                  'BLANK', 'EMPTY', 'BACKGROUND', 'BG'}

# No greys in the cluster palette: grey is reserved for flagged clusters, and the
# tissue background is greyscale.
_TAB10 = plt.get_cmap('tab10').colors
CLUSTER_COLOURS = tuple(c for i, c in enumerate(_TAB10) if i != 7) + tuple(
    c for i, c in enumerate(plt.get_cmap('tab20b').colors) if i % 4 in (0, 2))
FLAGGED_COLOUR = (0.62, 0.62, 0.62)
NAME_MIN_Z = 0.5                  # a marker must average at least this (robust z) to name a cluster
NAME_LOW_Z = -0.75                # ... or at most this to name it "low <marker>"
MIN_CELLS_PER_TILE = 3            # overview map: emptier tiles (isolated debris) are left blank


# =============================================================================
# Loading
# =============================================================================

def _channel_token(name):
    parts = re.split(r'[\s_()\-]+', str(name).strip())
    return parts[0].upper() if parts and parts[0] else ''


def classify_channel(name):
    """Return ``'nuclear'``, ``'control'`` or ``'marker'`` for a channel name."""
    token = _channel_token(name)
    if token in NUCLEAR_TOKENS:
        return 'nuclear'
    if token in CONTROL_TOKENS:
        return 'control'
    return 'marker'


def code_version():
    """Short git commit of the CellSurvey code (with ``-dirty`` if modified), or ``None``."""
    repo = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    try:
        sha = subprocess.run(['git', '-C', repo, 'rev-parse', '--short', 'HEAD'], capture_output=True,
                             text=True, timeout=10, check=True).stdout.strip()
        dirty = subprocess.run(['git', '-C', repo, 'status', '--porcelain', '--untracked-files=no'],
                               capture_output=True, text=True, timeout=10).stdout.strip()
        return sha + ('-dirty' if dirty else '') if sha else None
    except Exception:
        return None


def _read_sdata(zarr_path):
    import spatialdata
    try:
        return spatialdata.read_zarr(zarr_path, selection=('images', 'tables'))
    except TypeError:   # older spatialdata without `selection`
        return spatialdata.read_zarr(zarr_path)


def _shape_layers(zarr_path):
    try:
        return sorted(n for n in os.listdir(os.path.join(zarr_path, 'shapes')) if not n.startswith(('.', 'zarr')))
    except OSError:
        return []


def _load_background(sdata):
    """Nuclear channel from the largest pyramid level whose longest side <= MAX_BACKGROUND_SIDE.

    Returns ``(image, (sx, sy), channel_name)`` with ``sx, sy`` = full-resolution
    pixels per background pixel, or ``(None, None, None)``.
    """
    if not len(sdata.images):
        print("Report: the Zarr has no image; drawing points on white")
        return None, None, None
    try:
        name = sdata.attrs.get('cell_segmentation_image') or next(iter(sdata.images))
        img = sdata.images[name]
        if hasattr(img, 'children') and len(img.children):
            levels = [img[s].ds['image'] for s in img.children]
        else:
            levels = [img]
        full = levels[0]
        level = levels[-1]
        for lv in levels:   # largest level that fits
            if max(lv.shape[-2:]) <= MAX_BACKGROUND_SIDE:
                level = lv
                break
        channels = [remove_channel_suffix(str(c)) for c in level.coords['c'].values]
        idx = next((i for i, c in enumerate(channels) if classify_channel(c) == 'nuclear'), 0)
        arr = np.asarray(level.isel(c=idx).data, dtype='float32')
        sy = full.shape[-2] / arr.shape[0]
        sx = full.shape[-1] / arr.shape[1]
        return arr, (sx, sy), channels[idx]
    except Exception as e:
        print(f"Report: no tissue background ({type(e).__name__}: {e}); drawing points on white")
        return None, None, None


def _load_matching_sweep(zarr_path, n_cells, mode, folder=None):
    """This run's sweep summary (``mode``: ``communities`` or ``clusters``), or ``None``.

    Looks in ``folder`` (default: the Zarr's directory) for the run-named files
    ``python -m cellsurvey.stability`` writes (``<run>_<community|cluster>_sweep*.csv``),
    then for the old fixed names (``community_sweep*.csv``), which runs sharing a
    folder overwrite. Either is accepted only if the per-cell CSV has one row
    per cell of this run.
    """
    from cellsurvey.stability import sweep_paths
    folder = folder or os.path.dirname(os.path.abspath(zarr_path))
    stem = {'communities': 'community', 'clusters': 'cluster', 'niches': 'niche'}[mode]
    candidates = [sweep_paths(zarr_path, mode, folder),
                  (os.path.join(folder, f'{stem}_sweep.csv'), os.path.join(folder, f'{stem}_sweep_summary.csv'))]
    for cells_path, summary_path in candidates:
        if not (os.path.exists(cells_path) and os.path.exists(summary_path)):
            continue
        try:
            with open(cells_path, 'rb') as f:
                n_rows = sum(1 for _ in f) - 1
            if n_rows != n_cells:
                print(f"Report: ignoring {summary_path}: it's from another run ({n_rows:,} cells, not {n_cells:,})")
                continue
            return pd.read_csv(summary_path)
        except Exception:
            continue
    return None


def load_report_data(zarr_path, summary_path=None, sweep_dir=None):
    """Read everything the report needs from ``<prefix>_seg.zarr`` (and summary/sweep files)."""
    sdata = _read_sdata(zarr_path)
    table = sdata.tables['table']
    intensity_df, counts_df = extract_feature_tables(table)
    obs = table.obs
    if 'kmeans_cluster' not in obs:
        raise ValueError(f"{zarr_path} has no k-means clusters; was it written by a full CellSurvey run?")
    if 'spatial' not in table.obsm:
        raise ValueError(f"{zarr_path} has no cell centroids (table.obsm['spatial'])")

    meta = dict(table.uns.get('cellsurvey', {}) or {})
    attrs = dict(sdata.attrs.get('cellsurvey', {}) or {})
    summary = {}
    if summary_path and os.path.exists(summary_path):
        with open(summary_path) as f:
            summary = json.load(f)
    sweep = _load_matching_sweep(zarr_path, table.n_obs, 'communities', sweep_dir)
    cluster_sweep = _load_matching_sweep(zarr_path, table.n_obs, 'clusters', sweep_dir)
    niche_sweep = _load_matching_sweep(zarr_path, table.n_obs, 'niches', sweep_dir)

    background, scale, background_channel = _load_background(sdata)
    xy = np.asarray(table.obsm['spatial'], dtype=float)[:, :2]
    community = (obs['community'].astype(int).to_numpy() if 'community' in obs else None)
    area_px = obs['area'].astype(float).to_numpy() if 'area' in obs else None
    return {
        'zarr_path': os.path.abspath(zarr_path),
        'intensity_df': intensity_df,
        'counts_df': counts_df,
        'cluster': obs['kmeans_cluster'].astype(int).to_numpy(),
        'community': community,
        'area_px': area_px,
        'x': xy[:, 0], 'y': xy[:, 1],
        'meta': meta, 'attrs': attrs, 'summary': summary, 'sweep': sweep, 'cluster_sweep': cluster_sweep,
        'niche_sweep': niche_sweep,
        'niche': obs['niche'].astype(int).to_numpy() if 'niche' in obs else None,
        'niche_region': obs['niche_region'].astype(int).to_numpy() if 'niche_region' in obs else None,
        'shape_layers': _shape_layers(zarr_path),
        'background': background, 'scale': scale, 'background_channel': background_channel,
    }


# =============================================================================
# Computations (no plotting)
# =============================================================================

def robust_z_params(values):
    """``(median, spread)`` of arcsinh intensity, as used by :func:`robust_z` (IQR; 0 → std, then 1)."""
    a = np.arcsinh(np.asarray(values, dtype='float64'))
    med = np.median(a)
    q75, q25 = np.percentile(a, [75, 25])
    spread = q75 - q25
    if not spread > 0:
        spread = a.std() if a.std() > 0 else 1.0
    return float(med), float(spread)


def robust_z(values):
    """Robust z-score of arcsinh intensity: (asinh(x) - median) / IQR (IQR 0 → std, then 1)."""
    med, spread = robust_z_params(values)
    return (np.arcsinh(np.asarray(values, dtype='float64')) - med) / spread


def raw_at_z(values, z):
    """The raw intensity at robust z = ``z`` for this channel (``z = 0``: the median cell)."""
    med, spread = robust_z_params(values)
    return float(np.sinh(med + z * spread))


def nonspecific_threshold(n_markers):
    """Number of markers a cell must be bright in to count as non-specific (None if too few markers)."""
    if n_markers < 2:
        return None
    return min(n_markers, max(NONSPECIFIC_MIN_MARKERS, math.ceil(NONSPECIFIC_MIN_FRACTION * n_markers)))


def cluster_name(row):
    """Automatic name from a cluster's mean robust-z profile (a Series indexed by marker).

    The (up to) two markers raised by at least ``NAME_MIN_Z``, e.g. ``GFAP`` or
    ``PSYN/PTAU``. If none is raised: ``low X/Y`` for markers at or below
    ``NAME_LOW_Z``, else ``no distinct marker``. So weakly-raised markers never
    end up in a name.
    """
    high = row[row >= NAME_MIN_Z].sort_values(ascending=False)
    if len(high):
        return '/'.join(high.index[:2])
    low = row[row <= NAME_LOW_Z].sort_values()
    if len(low):
        return 'low ' + '/'.join(low.index[:2])
    return 'no distinct marker'


def compute(data):
    """All tables the figures and text use."""
    df = data['intensity_df']
    classes = {c: classify_channel(c) for c in df.columns}
    markers = [c for c in df.columns if classes[c] == 'marker']
    if not markers:   # e.g. only DAPI + controls: fall back to every channel
        markers = list(df.columns)
    Z = pd.DataFrame({m: robust_z(df[m].to_numpy()) for m in markers}, index=df.index)
    cluster = data['cluster']
    clusters = sorted(np.unique(cluster).tolist())

    n_needed = nonspecific_threshold(len(markers))
    if n_needed is not None:
        nonspecific = ((Z.to_numpy() > POSITIVE_Z).sum(axis=1) >= n_needed)
    else:
        nonspecific = np.zeros(len(Z), dtype=bool)

    profile = Z.groupby(cluster).mean()
    sizes = pd.Series(cluster).value_counts().reindex(clusters).fillna(0).astype(int)
    ns_frac = pd.Series(nonspecific).groupby(cluster).mean().reindex(clusters).fillna(0)
    flagged = [c for c in clusters if ns_frac[c] > FLAG_CLUSTER_FRACTION]
    # Most distinctive markers: largest |mean robust z|, so clearly *low* markers count too.
    top = {c: list(profile.loc[c].reindex(profile.loc[c].abs().sort_values(ascending=False).index)
                   .head(5).items()) for c in clusters}
    names = {c: ('non-specific' if c in flagged else cluster_name(profile.loc[c])) for c in clusters}

    marker_rows = []
    for m in markers:
        pos = Z[m].to_numpy() > POSITIVE_Z
        n_pos = int(pos.sum())
        share = float((pos & nonspecific).sum() / n_pos) if n_pos else 0.0
        if n_pos < MIN_POSITIVES:
            status = 'too few'
        else:
            status = ('dominated' if share > MARKER_DOMINATED else 'suspect' if share > MARKER_SUSPECT
                      else 'clean')
        marker_rows.append({'marker': m, 'n_positive': n_pos, 'pct_positive': 100 * pos.mean(),
                            'pct_nonspecific': 100 * share, 'status': status})

    # Raw (un-normalised) levels: the "positive" line and the typical cell back on the intensity scale,
    # for the markers and, as a background reference, the control channels.
    raw_rows = []
    for ch in markers + [c for c in df.columns if classes[c] == 'control' and c not in markers]:
        v = df[ch].to_numpy(dtype=float)
        v = v[np.isfinite(v)]
        median, line = raw_at_z(v, 0), raw_at_z(v, POSITIVE_Z)
        raw_rows.append({'channel': ch, 'kind': classes[ch] if ch not in markers else 'marker',
                         'median': median, 'positive_line': line,
                         'ratio': line / median if median > 0 else np.nan,
                         'p99': float(np.percentile(v, 99)) if v.size else np.nan})

    comm = None
    if data['community'] is not None:
        community = data['community']
        valid = community >= 0
        if valid.any():
            ct = pd.crosstab(community[valid], cluster[valid])
            frac = ct.div(ct.sum(axis=1), axis=0)
            comm_sizes = ct.sum(axis=1).sort_values(ascending=False)
            dom = frac.idxmax(axis=1)
            comm = {
                'sizes': comm_sizes,
                'frac': frac,
                'dominant': dom,
                'dominant_share': frac.max(axis=1),
                'n_unassigned': int((~valid).sum()),
            }
    niche = None
    if data.get('niche') is not None:
        nl = np.asarray(data['niche'])
        valid = nl >= 0
        if valid.any():
            ct = pd.crosstab(nl[valid], cluster[valid])
            frac = ct.div(ct.sum(axis=1), axis=0)
            ids = sorted(ct.index.tolist())
            niche = {'ids': ids, 'frac': frac, 'sizes': ct.sum(axis=1),
                     'names': {n: niche_name(frac.loc[n], names) for n in ids},
                     'n_unassigned': int((~valid).sum())}
    region = None
    if niche is not None and data.get('niche_region') is not None:
        rl = np.asarray(data['niche_region'])
        nl = np.asarray(data['niche'])
        valid = (rl >= 0) & (nl >= 0)
        if valid.any():
            ct = pd.crosstab(rl[valid], nl[valid])
            sizes_r = pd.Series(np.bincount(rl[rl >= 0]), name='cells')
            sizes_r = sizes_r[sizes_r > 0].sort_values(ascending=False, kind='stable')
            region = {
                'sizes': sizes_r,
                'niche': ct.idxmax(axis=1),                       # majority niche of each region
                'purity': ct.max(axis=1) / ct.sum(axis=1),        # share of its cells in that niche
                'n_unassigned': int((rl < 0).sum()),
                'per_niche': ct.idxmax(axis=1).value_counts().reindex(niche['ids']).fillna(0).astype(int),
            }
    return {
        'markers': markers, 'classes': classes, 'Z': Z, 'clusters': clusters, 'sizes': sizes,
        'profile': profile, 'top': top, 'names': names, 'flagged': flagged, 'ns_frac': ns_frac,
        'nonspecific': nonspecific, 'n_needed': n_needed, 'marker_table': pd.DataFrame(marker_rows),
        'raw_table': pd.DataFrame(raw_rows), 'comm': comm, 'niche': niche, 'region': region,
    }


def niche_name(frac_row, cluster_names):
    """Automatic niche name from its cluster mix, e.g. ``mostly 3 (PSYN/CNPASE)`` or ``1 (…) + 4 (…)``."""
    top = frac_row.sort_values(ascending=False)
    c1, s1 = top.index[0], top.iloc[0]
    if s1 >= 0.5 or len(top) == 1:
        return f'mostly {c1} ({cluster_names.get(c1, c1)}), {100 * s1:.0f}%'
    c2, s2 = top.index[1], top.iloc[1]
    return (f'{c1} ({cluster_names.get(c1, c1)}) {100 * s1:.0f}% + '
            f'{c2} ({cluster_names.get(c2, c2)}) {100 * s2:.0f}%')


# =============================================================================
# Figures (each a function of the loaded data + computed tables)
# =============================================================================

def _cluster_colour(c, res):
    if c in res['flagged']:
        return FLAGGED_COLOUR
    return CLUSTER_COLOURS[res['clusters'].index(c) % len(CLUSTER_COLOURS)]


def _display_xy(data):
    """Cell coordinates in the background image's pixel frame (or full-res if no background)."""
    if data['background'] is None:
        return data['x'], data['y']
    sx, sy = data['scale']
    return data['x'] / sx, data['y'] / sy


def _tissue_axes(ax, data, xs=None, ys=None, dim=0.45):
    """Draw the tissue background (image coordinates: y down) and set equal aspect.

    The background is dimmed (shown at ``dim`` of full brightness) so coloured
    points and outlines on top stay visible.
    """
    bg = data['background']
    if bg is not None:
        finite = bg[np.isfinite(bg)]
        lo, hi = (np.percentile(finite, [1, 99.8]) if finite.size else (0, 1))
        shown = np.clip((bg - lo) / max(hi - lo, 1e-9), 0, 1) * dim
        ax.imshow(shown, cmap='gray', origin='upper', vmin=0, vmax=1, interpolation='nearest')
        ax.set_xlim(0, bg.shape[1])
        ax.set_ylim(bg.shape[0], 0)
    else:
        ax.set_facecolor('white')
        if xs is not None and len(xs):
            ax.set_xlim(np.nanmin(xs), np.nanmax(xs))
            ax.set_ylim(np.nanmax(ys), np.nanmin(ys))
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])


def fig_cluster_profile(res):
    prof = res['profile'].loc[res['clusters'], res['markers']]
    n_c, n_m = prof.shape
    fig, ax = plt.subplots(figsize=(max(6, 0.55 * n_m + 4), max(3.5, 0.45 * n_c + 1.5)), layout='constrained')
    im = ax.imshow(prof.to_numpy(), aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
    ax.set_xticks(range(n_m))
    ax.set_xticklabels(res['markers'], rotation=60, ha='right', fontsize=8)
    ax.set_yticks(range(n_c))
    ax.set_yticklabels([f"{c}: {res['names'][c]}" + (' ⚠' if c in res['flagged'] else '')
                        for c in res['clusters']], fontsize=8)
    ax.set_title('Mean marker level per k-means cluster (robust z)', fontsize=11)
    fig.colorbar(im, ax=ax, shrink=0.8, label='mean robust z (capped at ±2)')
    return fig


def fig_cluster_maps(data, res, rng):
    xs, ys = _display_xy(data)
    clusters = res['clusters']
    ncols = 4 if len(clusters) <= 8 else 5
    nrows = math.ceil(len(clusters) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.6 * ncols, 3.8 * nrows), layout='constrained', squeeze=False)
    for ax in axes.ravel():
        ax.axis('off')
    cluster = data['cluster']
    for ax, c in zip(axes.ravel(), clusters):
        ax.axis('on')
        _tissue_axes(ax, data, xs, ys)
        idx = np.flatnonzero(cluster == c)
        n = len(idx)
        if n > MAX_POINTS_PER_PANEL:
            idx = rng.choice(idx, MAX_POINTS_PER_PANEL, replace=False)
        # Small clusters get bigger, more opaque points so a few hundred cells still show.
        size, alpha = (0.3, 0.45) if n > 20_000 else (1.5, 0.7) if n > 2_000 else (6.0, 0.9)
        ax.scatter(xs[idx], ys[idx], s=size, alpha=alpha, color=_cluster_colour(c, res), linewidths=0,
                   rasterized=True)
        ax.set_title(f"{c}: {res['names'][c]} ({int(res['sizes'][c]):,})", fontsize=8)
    fig.suptitle('Where each k-means cluster sits (over the tissue image)', fontsize=12)
    return fig


def _dominant_tile_map(data, labels, ids, colours, legend_labels, title, bins=220):
    """Each square tile coloured by its most common label (``ids``), over the tissue image.

    Tiles with fewer than ``MIN_CELLS_PER_TILE`` labelled cells stay blank;
    negative labels are ignored.
    """
    xs, ys = _display_xy(data)
    keep = np.isfinite(xs) & np.isfinite(ys) & (np.asarray(labels) >= 0)
    xs_f, ys_f, lab = xs[keep], ys[keep], np.asarray(labels)[keep]
    if data['background'] is not None:
        h, w = data['background'].shape
        rng_x, rng_y = (0, w), (0, h)
    else:
        f = np.isfinite(xs) & np.isfinite(ys)
        rng_x, rng_y = (xs[f].min(), xs[f].max()), (ys[f].min(), ys[f].max())
    long_side = max(rng_x[1] - rng_x[0], rng_y[1] - rng_y[0])
    nx = max(1, int(bins * (rng_x[1] - rng_x[0]) / long_side))
    ny = max(1, int(bins * (rng_y[1] - rng_y[0]) / long_side))
    counts = np.stack([np.histogram2d(ys_f[lab == c], xs_f[lab == c], bins=(ny, nx), range=(rng_y, rng_x))[0]
                       for c in ids])
    dominant = counts.argmax(axis=0)
    occupied = counts.sum(axis=0) >= MIN_CELLS_PER_TILE
    rgba = np.zeros((ny, nx, 4))
    rgba[..., :3] = np.array(colours)[dominant]
    rgba[..., 3] = np.where(occupied, 0.85, 0.0)

    fig, ax = plt.subplots(figsize=(9, 9 * (rng_y[1] - rng_y[0]) / max(rng_x[1] - rng_x[0], 1)))
    _tissue_axes(ax, data, xs, ys)
    ax.imshow(rgba, origin='upper', extent=(rng_x[0], rng_x[1], rng_y[1], rng_y[0]), interpolation='nearest')
    handles = [mpatches.Patch(color=colours[i], label=legend_labels[i]) for i in range(len(ids))]
    ax.legend(handles=handles, fontsize=7, loc='upper left', bbox_to_anchor=(1.01, 1), frameon=False)
    ax.set_title(title, fontsize=11)
    return fig


def fig_cluster_overview_map(data, res, bins=220):
    """Each square bin coloured by its most common cluster, over the tissue image."""
    return _dominant_tile_map(data, data['cluster'], res['clusters'],
                              [_cluster_colour(c, res) for c in res['clusters']],
                              [f"{c}: {res['names'][c]}" for c in res['clusters']],
                              'All clusters at once: each tile shows its most common cluster', bins)


NICHE_COLOURS = plt.get_cmap('Set1').colors + plt.get_cmap('Set3').colors


def fig_niche_map(data, res):
    """Each tile coloured by its most common niche, over the tissue image."""
    niche = res['niche']
    ids = niche['ids']
    return _dominant_tile_map(data, data['niche'], ids,
                              [NICHE_COLOURS[i % len(NICHE_COLOURS)] for i in range(len(ids))],
                              [f"{n}: {niche['names'][n]}" for n in ids],
                              'Niches: each tile shows its most common niche')


def fig_niche_composition(res):
    """Heatmap: the share of each k-means cluster among the cells of each niche."""
    niche = res['niche']
    frac = niche['frac'].reindex(index=niche['ids'], columns=res['clusters']).fillna(0)
    fig, ax = plt.subplots(figsize=(max(6, 0.6 * len(res['clusters']) + 3.5), max(3, 0.45 * len(frac) + 1.5)),
                           layout='constrained')
    im = ax.imshow(frac.to_numpy(), aspect='auto', cmap='viridis', vmin=0, vmax=1)
    for i in range(frac.shape[0]):
        for j in range(frac.shape[1]):
            v = frac.iat[i, j]
            if v >= 0.05:
                ax.text(j, i, f'{100 * v:.0f}', ha='center', va='center', fontsize=7,
                        color='white' if v < 0.6 else 'black')
    ax.set_xticks(range(len(res['clusters'])))
    ax.set_xticklabels([f"{c}: {res['names'][c]}" for c in res['clusters']], rotation=45, ha='right', fontsize=8)
    ax.set_yticks(range(len(frac)))
    ax.set_yticklabels([f"niche {n} ({int(niche['sizes'][n]):,})" for n in frac.index], fontsize=8)
    ax.set_title('What each niche is made of (% of its cells in each k-means cluster)', fontsize=11)
    fig.colorbar(im, ax=ax, shrink=0.8, label='fraction of the niche\'s cells')
    return fig


MAX_REGION_LABELS = 30


def fig_region_map(data, res, bins=260):
    """Niche regions: tiles coloured by their region's niche, black lines between regions, the largest numbered."""
    niche, region = res['niche'], res['region']
    xs, ys = _display_xy(data)
    rl = np.asarray(data['niche_region'])
    keep = np.isfinite(xs) & np.isfinite(ys) & (rl >= 0)
    if data['background'] is not None:
        h, w = data['background'].shape
        rng_x, rng_y = (0, w), (0, h)
    else:
        f = np.isfinite(xs) & np.isfinite(ys)
        rng_x, rng_y = (xs[f].min(), xs[f].max()), (ys[f].min(), ys[f].max())
    long_side = max(rng_x[1] - rng_x[0], rng_y[1] - rng_y[0])
    nx_ = max(1, int(bins * (rng_x[1] - rng_x[0]) / long_side))
    ny_ = max(1, int(bins * (rng_y[1] - rng_y[0]) / long_side))
    tx = np.clip(((xs[keep] - rng_x[0]) / max(rng_x[1] - rng_x[0], 1e-9) * nx_).astype(int), 0, nx_ - 1)
    ty = np.clip(((ys[keep] - rng_y[0]) / max(rng_y[1] - rng_y[0], 1e-9) * ny_).astype(int), 0, ny_ - 1)
    # Most common region per tile (sparse tile × region counts).
    from scipy.sparse import coo_matrix
    counts = coo_matrix((np.ones(len(tx)), (ty * nx_ + tx, rl[keep])),
                        shape=(ny_ * nx_, int(rl.max()) + 1)).tocsr()
    occupied = (np.asarray(counts.sum(axis=1)).ravel() >= MIN_CELLS_PER_TILE).reshape(ny_, nx_)
    dominant = np.asarray(counts.argmax(axis=1)).ravel().reshape(ny_, nx_)
    niche_of = region['niche']
    colours = [NICHE_COLOURS[i % len(NICHE_COLOURS)] for i in range(len(niche['ids']))]
    lut = np.zeros((int(rl.max()) + 1, 3))
    for r, n in niche_of.items():
        lut[r] = colours[niche['ids'].index(n)]
    rgba = np.zeros((ny_, nx_, 4))
    rgba[..., :3] = lut[dominant]
    rgba[..., 3] = np.where(occupied, 0.8, 0.0)
    # Region borders: occupied tiles whose right or lower neighbour belongs to another region.
    border = np.zeros((ny_, nx_), dtype=bool)
    diff_x = (dominant[:, 1:] != dominant[:, :-1]) & occupied[:, 1:] & occupied[:, :-1]
    diff_y = (dominant[1:, :] != dominant[:-1, :]) & occupied[1:, :] & occupied[:-1, :]
    border[:, 1:] |= diff_x
    border[1:, :] |= diff_y
    rgba[border] = (0, 0, 0, 0.9)

    fig, ax = plt.subplots(figsize=(9, 9 * (rng_y[1] - rng_y[0]) / max(rng_x[1] - rng_x[0], 1)))
    _tissue_axes(ax, data, xs, ys)
    ax.imshow(rgba, origin='upper', extent=(rng_x[0], rng_x[1], rng_y[1], rng_y[0]), interpolation='nearest')
    halo = [patheffects.withStroke(linewidth=2.5, foreground='white')]
    for r in region['sizes'].index[:MAX_REGION_LABELS]:
        m = keep & (rl == r)
        if not m.any():
            continue
        # Label on the region's cell nearest its median point (a curved band's median can lie outside it).
        mx, my = np.median(xs[m]), np.median(ys[m])
        i = np.argmin((xs[m] - mx) ** 2 + (ys[m] - my) ** 2)
        ax.text(xs[m][i], ys[m][i], str(r), fontsize=8, ha='center', va='center',
                weight='bold', color='black', path_effects=halo)
    handles = [mpatches.Patch(color=colours[i], label=f"niche {n}: {niche['names'][n]}")
               for i, n in enumerate(niche['ids'])]
    ax.legend(handles=handles, fontsize=7, loc='upper left', bbox_to_anchor=(1.01, 1), frameon=False)
    n_reg = len(region['sizes'])
    ax.set_title(f"{n_reg} niche regions, coloured by niche"
                 + (f" (the {MAX_REGION_LABELS} largest numbered)" if n_reg > MAX_REGION_LABELS else ''),
                 fontsize=11)
    return fig


def fig_niche_sweep(sweep, radius_used, k_used):
    """Silhouette and seed stability against niche count, one line per radius."""
    sw = sweep.sort_values(['niche_radius_um', 'n_niches'])
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.6), layout='constrained')
    for radius, g in sw.groupby('niche_radius_um'):
        style = dict(lw=2.4) if radius_used is not None and np.isclose(radius, radius_used) else dict(lw=1.2)
        axes[0].plot(g['n_niches'], g['silhouette'], 'o-', label=f'{radius:g} µm', **style)
        axes[1].plot(g['n_niches'], g['seed_ari'], 'o-', label=f'{radius:g} µm', **style)
    axes[0].set_ylabel('silhouette (higher = better separated)')
    axes[0].set_title('Separation (rises with radius: can\'t pick the radius)', fontsize=10)
    axes[1].set_ylabel('seed stability (ARI; 1 = identical)')
    axes[1].set_ylim(0, 1.05)
    axes[1].set_title('Stability across random seeds', fontsize=10)
    for a in axes:
        a.set_xlabel('number of niches')
        a.set_xticks(sorted(sw['n_niches'].unique()))
        if k_used:
            a.axvline(k_used, color='#C1622D', ls='--', lw=1)
        a.legend(fontsize=7, title='radius', title_fontsize=7)
    return fig


def fig_community_outlines(data, res, grid=320):
    comm = res['comm']
    xs, ys = _display_xy(data)
    community = data['community']
    top = comm['sizes'].index[:MAX_COMMUNITY_PANELS].tolist()
    if data['background'] is not None:
        h, w = data['background'].shape
        rng_x, rng_y = (0, w), (0, h)
    else:
        f = np.isfinite(xs) & np.isfinite(ys)
        rng_x, rng_y = (xs[f].min(), xs[f].max()), (ys[f].min(), ys[f].max())
    fig, ax = plt.subplots(figsize=(9, 9 * (rng_y[1] - rng_y[0]) / max(rng_x[1] - rng_x[0], 1)))
    _tissue_axes(ax, data, xs, ys)
    cmap = plt.get_cmap('tab20')
    halo = [patheffects.withStroke(linewidth=2.5, foreground='black')]
    for i, c in enumerate(top):
        m = community == c
        H, ye, xe = np.histogram2d(ys[m], xs[m], bins=grid, range=(rng_y, rng_x))
        H = gaussian_filter(H, sigma=1.5)
        if H.max() <= 0:
            continue
        xc, yc = (xe[:-1] + xe[1:]) / 2, (ye[:-1] + ye[1:]) / 2
        ax.contour(xc, yc, H, levels=[max(H.max() * 0.05, 1e-9)], colors=[cmap(i % 20)], linewidths=1.3)
        ax.text(np.nanmedian(xs[m]), np.nanmedian(ys[m]), str(c), color=cmap(i % 20), fontsize=9,
                ha='center', va='center', weight='bold', path_effects=halo)
    shown = len(top)
    total = len(comm['sizes'])
    ax.set_title(f"Outlines of the {shown} largest communities"
                 + (f" (of {total})" if total > shown else ""), fontsize=11)
    return fig


def fig_community_composition(res):
    comm = res['comm']
    top = comm['sizes'].index[:max(MAX_COMMUNITY_PANELS, 30)].tolist()
    dom = comm['dominant'].loc[top]
    order = sorted(top, key=lambda c: (res['clusters'].index(dom[c]), -comm['sizes'][c]))
    frac = comm['frac'].reindex(index=order, columns=res['clusters']).fillna(0)
    fig, ax = plt.subplots(figsize=(max(8, 0.32 * len(order) + 4), 4.8), layout='constrained')
    bottom = np.zeros(len(order))
    for c in res['clusters']:
        vals = frac[c].to_numpy()
        ax.bar(range(len(order)), vals, bottom=bottom, width=0.9, color=_cluster_colour(c, res),
               label=f"{c}: {res['names'][c]}")
        bottom += vals
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([str(c) for c in order], rotation=90, fontsize=7)
    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.set_ylim(0, 1)
    ax.set_xlabel('community (grouped by its most common cluster, then by size)')
    ax.set_ylabel('fraction of cells')
    ax.set_title(f'What each of the {len(order)} largest communities is made of', fontsize=11)
    ax.legend(fontsize=7, loc='upper left', bbox_to_anchor=(1.01, 1), frameon=False)
    return fig


def fig_marker_distributions(res):
    """One log-scale histogram per marker on the robust z scale, with the "positive" line."""
    markers = res['markers']
    ncols = min(4, len(markers))
    nrows = math.ceil(len(markers) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.4 * ncols, 2.5 * nrows), squeeze=False,
                             layout='constrained')
    for ax in axes.ravel():
        ax.axis('off')
    for ax, m in zip(axes.ravel(), markers):
        ax.axis('on')
        z = res['Z'][m].to_numpy()
        lo, hi = np.percentile(z, [0.01, 99.99])
        lo, hi = min(lo, -3), max(hi, POSITIVE_Z + 2)
        ax.hist(np.clip(z, lo, hi), bins=120, range=(lo, hi), color='#9A9FA8', log=True)
        ax.axvline(POSITIVE_Z, color='#C1622D', lw=1.5)
        ax.axvline(0, color='#6C5B9E', lw=0.8, ls=':')
        pct = 100 * (z > POSITIVE_Z).mean()
        ax.set_title(f'{m}  ({pct:.1f}% "positive")', fontsize=9)
        ax.tick_params(labelsize=7)
    fig.suptitle(f'Marker distributions (robust z; log count). Orange = "positive" line (z = {POSITIVE_Z:g}); '
                 'dotted = typical cell (z = 0)', fontsize=11)
    return fig


def fig_raw_distributions(data, res):
    """Raw (un-normalised) per-cell intensities of every marker, plus control channels for reference.

    All panels share one log intensity axis, so absolute levels can be compared:
    a marker whose cells all sit near the background (and near the control
    channels) is probably absent, and its "positive" cells are then just the
    brightest noise. The "positive" line and the typical cell are the same
    cutoffs as on the robust z scale, converted back to intensity.
    """
    df = data['intensity_df']
    raw = res['raw_table']
    channels = raw['channel'].tolist()
    values = {ch: np.clip(df[ch].to_numpy(dtype=float), 0, None) + 1 for ch in channels}
    values = {ch: v[np.isfinite(v)] for ch, v in values.items()}
    lo = max(1.0, min(np.percentile(v, 0.01) for v in values.values() if v.size))
    hi = max(np.percentile(v, 99.99) for v in values.values() if v.size)
    hi = max(hi, raw['positive_line'].max() + 1, lo * 10)
    bins = np.logspace(np.log10(lo), np.log10(hi), 120)
    ncols = min(4, len(channels))
    nrows = math.ceil(len(channels) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.4 * ncols, 2.5 * nrows), squeeze=False,
                             layout='constrained', sharex=True)
    for ax in axes.ravel():
        ax.axis('off')
    for ax, r in zip(axes.ravel(), raw.itertuples()):
        ax.axis('on')
        control = r.kind != 'marker'
        ax.hist(np.clip(values[r.channel], lo, hi), bins=bins, log=True,
                color='#C9CCD1' if control else '#9A9FA8')
        ax.axvline(r.positive_line + 1, color='#C1622D', lw=1.5)
        ax.axvline(r.median + 1, color='#6C5B9E', lw=0.8, ls=':')
        ax.set_xscale('log')
        ax.set_title(f'{r.channel}' + (' (control)' if control else '')
                     + f'\nline at {r.positive_line:,.0f} = {r.ratio:.1f}× typical', fontsize=8.5)
        ax.tick_params(labelsize=7)
    for ax in axes[-1]:
        if ax.axison:
            ax.set_xlabel('mean intensity per cell + 1 (log)', fontsize=7)
    fig.suptitle('Raw marker intensities (same axis for every panel; log count). Orange = the same "positive" line; '
                 'dotted = typical cell', fontsize=11)
    return fig


def fig_positive_maps(data, res, rng):
    """Cells above the "positive" line for each marker, over the tissue; non-specific cells in grey."""
    xs, ys = _display_xy(data)
    markers = res['markers'][:MAX_MARKER_PANELS]
    ncols = 4 if len(markers) <= 8 else 5
    nrows = math.ceil(len(markers) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.6 * ncols, 3.8 * nrows), layout='constrained', squeeze=False)
    for ax in axes.ravel():
        ax.axis('off')
    ns = res['nonspecific']
    for ax, m in zip(axes.ravel(), markers):
        ax.axis('on')
        _tissue_axes(ax, data, xs, ys)
        pos = res['Z'][m].to_numpy() > POSITIVE_Z
        n = int(pos.sum())
        size, alpha = (0.4, 0.5) if n > 20_000 else (1.5, 0.75) if n > 2_000 else (5.0, 0.9)
        for sel, colour in ((pos & ns, FLAGGED_COLOUR), (pos & ~ns, '#F28E2B')):
            idx = np.flatnonzero(sel)
            if len(idx) > MAX_POINTS_PER_PANEL:
                idx = rng.choice(idx, MAX_POINTS_PER_PANEL, replace=False)
            ax.scatter(xs[idx], ys[idx], s=size, alpha=alpha, color=colour, linewidths=0, rasterized=True)
        ax.set_title(f'{m} ({n:,} positive)', fontsize=8)
    extra = (f' (first {MAX_MARKER_PANELS} of {len(res["markers"])} markers)'
             if len(res['markers']) > MAX_MARKER_PANELS else '')
    fig.suptitle(f'"Positive" cells (robust z > {POSITIVE_Z:g}) for each marker{extra}: orange; '
                 'grey = non-specific cells', fontsize=12)
    return fig


def fig_cluster_sweep(sweep, k_used):
    """Inertia and silhouette against k, with the run's k marked."""
    sw = sweep.sort_values('n_clusters')
    has_sil = 'silhouette' in sw.columns and sw['silhouette'].notna().any()
    fig, axes = plt.subplots(1, 2 if has_sil else 1, figsize=(10 if has_sil else 5.5, 3.6), squeeze=False,
                             layout='constrained')
    ax = axes[0, 0]
    ax.plot(sw['n_clusters'], sw['inertia'], 'o-', color='#4C72B0')
    ax.set_xlabel('number of clusters (k)')
    ax.set_ylabel('inertia (within-cluster spread)')
    ax.set_title('Elbow: where adding clusters stops helping', fontsize=10)
    if has_sil:
        ax2 = axes[0, 1]
        ax2.plot(sw['n_clusters'], sw['silhouette'], 'o-', color='#55A868')
        ax2.set_xlabel('number of clusters (k)')
        ax2.set_ylabel('mean silhouette (higher = better separated)')
        ax2.set_title('Silhouette: how distinct the clusters are', fontsize=10)
    for a in axes.ravel():
        a.axvline(k_used, color='#C1622D', ls='--', lw=1.2, label=f'this run (k = {k_used})')
        a.set_xticks(sw['n_clusters'])
        a.legend(fontsize=8)
    return fig


def fig_quality(res):
    clusters = res['clusters']
    vals = [100 * res['ns_frac'][c] for c in clusters]
    fig, ax = plt.subplots(figsize=(max(6, 0.6 * len(clusters) + 2), 3.6), layout='constrained')
    ax.bar(range(len(clusters)), vals, color=[_cluster_colour(c, res) for c in clusters])
    ax.axhline(100 * FLAG_CLUSTER_FRACTION, color='#C1622D', lw=1.2, ls='--')
    ax.set_xticks(range(len(clusters)))
    ax.set_xticklabels([f"{c}: {res['names'][c]}" for c in clusters], rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('% non-specific cells')
    ax.set_ylim(0, 100)
    ax.set_title('Share of each cluster that is non-specific (dashed line = flag threshold)', fontsize=11)
    return fig


# =============================================================================
# HTML
# =============================================================================

CSS = """
body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;max-width:1100px;margin:2em auto;padding:0 1.2em;
color:#222;line-height:1.5}
h1{font-size:1.7em;margin-bottom:.2em}h2{border-bottom:1px solid #ddd;padding-bottom:.2em;margin-top:2em}
h3{margin-top:1.4em}table{border-collapse:collapse;margin:.8em 0;font-size:.9em}
th,td{border:1px solid #ddd;padding:.3em .6em;text-align:left;vertical-align:top}th{background:#f5f5f5}
td.num{text-align:right;font-variant-numeric:tabular-nums}img{max-width:100%;height:auto;margin:.6em 0}
.note{background:#f7f7f2;border-left:4px solid #c9b458;padding:.6em 1em;margin:1em 0}
.warn{background:#fbf1ee;border-left:4px solid #C1622D;padding:.6em 1em;margin:1em 0}
.muted{color:#666;font-size:.9em}nav ol{columns:2}code{background:#f3f3f3;padding:0 .25em}
@media print{body{max-width:none;margin:0}h2{page-break-before:always}h2:first-of-type{page-break-before:auto}
img{page-break-inside:avoid}}
"""


def _esc(x):
    return html.escape(str(x))


def _table(headers, rows, num_cols=()):
    out = ['<table><tr>' + ''.join(f'<th>{_esc(h)}</th>' for h in headers) + '</tr>']
    for r in rows:
        cells = []
        for i, v in enumerate(r):
            cls = ' class="num"' if i in num_cols else ''
            cells.append(f'<td{cls}>{v}</td>')
        out.append('<tr>' + ''.join(cells) + '</tr>')
    out.append('</table>')
    return '\n'.join(out)


def _fmt(v, nd=None):
    if v is None or (isinstance(v, (float, np.floating)) and not np.isfinite(v)):
        return '—'
    if nd is not None:
        return f'{v:.{nd}f}'
    if isinstance(v, (float, np.floating)) and float(v).is_integer():
        v = int(v)   # whole numbers read back from CSV as floats (51.0 → 51)
    if isinstance(v, (int, np.integer)):
        return f'{v:,}'
    if isinstance(v, (float, np.floating)):
        return f'{v:g}'
    return str(v)


class _Figures:
    """Saves figures as PNG files and returns base64 <img> tags."""

    def __init__(self, out_dir):
        self.out_dir = out_dir
        os.makedirs(out_dir, exist_ok=True)

    def tag(self, fig, name, alt, dpi=130):
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight', pad_inches=0.1)
        plt.close(fig)
        data = buf.getvalue()
        with open(os.path.join(self.out_dir, f'{name}.png'), 'wb') as f:
            f.write(data)
        return (f'<img alt="{_esc(alt)}" src="data:image/png;base64,'
                f'{base64.b64encode(data).decode("ascii")}">')


def _source_text(attrs, zarr_path):
    src = attrs.get('source') or {}
    if src.get('type') == 'omero':
        return (f"OMERO {src.get('host')} — Image:{src.get('image_id')} "
                f"({src.get('image_name', '')}), group {src.get('group_id', '?')}")
    if src.get('type') == 'file':
        return src.get('path', '')
    return f"(not recorded) — results in {zarr_path}"


def _segmenter_text(meta, shape_layers):
    seg = meta.get('segmenter')
    if not seg:
        if 'nucleus_boundaries' in shape_layers:
            seg = 'instanseg'
        elif 'stardist_boundaries' in shape_layers:
            seg = 'stardist'
    if seg == 'instanseg':
        model = meta.get('instanseg_model')
        return 'InstanSeg' + (f' ({model})' if model else '') + ' — nuclei + whole cells'
    if seg == 'stardist':
        exp = meta.get('cell_expansion_um')
        cells = ('nuclei expanded by ' + f'{exp:g} µm' if exp else
                 'nuclei, expanded to cells' if 'cell_boundaries' in shape_layers else 'nuclei only')
        return f'Stardist (2D_versatile_fluo) — {cells}'
    return '(not recorded)'


def build_html(data, res, figs):
    meta, attrs, summary = data['meta'], data['attrs'], data['summary']
    get = lambda k: meta.get(k, summary.get(k))   # noqa: E731 — stored in the Zarr, else summary.json
    n_cells = len(data['cluster'])
    markers = res['markers']
    nuclear = [c for c, k in res['classes'].items() if k == 'nuclear']
    controls = [c for c, k in res['classes'].items() if k == 'control']
    comm = res['comm']
    px = meta.get('pixel_size_um') or attrs.get('pixel_size_um')

    parts = [f'<!doctype html><html lang="en"><head><meta charset="utf-8">'
             f'<meta name="viewport" content="width=device-width,initial-scale=1">'
             f'<title>CellSurvey overview</title><style>{CSS}</style></head><body>']
    parts.append('<h1>CellSurvey overview report</h1>')
    parts.append(f'<p class="muted">Generated {datetime.datetime.now():%Y-%m-%d %H:%M} from '
                 f'<code>{_esc(data["zarr_path"])}</code></p>')

    # --- Header -------------------------------------------------------------
    area_note = ''
    if data['area_px'] is not None and px:
        med_area = float(np.nanmedian(data['area_px'])) * px ** 2
        area_note = f' (median cell area {med_area:,.0f} µm²)'
    rows = [
        ('Source image', _esc(_source_text(attrs, data['zarr_path']))),
        ('Segmentation', _esc(_segmenter_text(meta, data['shape_layers']))),
        ('Pixel size', f'{px:.4f} µm/px' if px else '—'),
        ('Cells analysed', f'{n_cells:,}{area_note}'),
        ('Marker channels', _esc(', '.join(markers)) + f' ({len(markers)})'),
        ('Nuclear / control channels', _esc(', '.join(nuclear + controls) or 'none')
         + ' — shown as background or excluded from the marker sections'),
        ('Clustering', f"k-means, {len(res['clusters'])} clusters, on "
                       f"<code>{_esc(get('cluster_on') or 'intensities')}</code> features"),
        ('Communities', (f"Louvain at resolution {_fmt(get('community_resolution'))}, Delaunay edges up to "
                         f"{_fmt(get('max_edge_distance_um'))} µm, minimum {_fmt(get('min_community_size'))} "
                         f"cells: <b>{len(comm['sizes'])} communities</b>, "
                         f"{comm['n_unassigned']:,} cells unassigned" if comm else 'not available')),
        ('Niches', (f"<b>{len(res['niche']['ids'])} niches</b> from each cell's mix of the "
                    f"{_fmt(get('niche_n_clusters') or len(res['clusters']))} k-means clusters within "
                    f"{_fmt(get('niche_radius_um'))} µm"
                    + (f"; <b>{len(res['region']['sizes'])} niche regions</b> of at least "
                       f"{_fmt(get('min_region_size'))} cells" if res['region'] else '')
                    if res['niche'] else 'not computed')),
        ('CellSurvey code', _esc(f"run: {meta.get('code_version') or 'not recorded'}; "
                                 f"report: {code_version() or 'unknown'}")),
        ('Report settings', f'"positive" = robust z &gt; {POSITIVE_Z:g}; non-specific = positive in ≥ '
                            f'{res["n_needed"] if res["n_needed"] else "—"} of {len(markers)} markers; '
                            f'clusters flagged if &gt; {100 * FLAG_CLUSTER_FRACTION:.0f}% non-specific'),
    ]
    parts.append(_table(['', ''], [(f'<b>{a}</b>', b) for a, b in rows]))

    parts.append('<nav><b>Contents</b><ol><li><a href="#read">How to read this</a></li>'
                 '<li><a href="#markers">Markers</a></li>'
                 '<li><a href="#clusters">Cell populations (clusters)</a></li>'
                 '<li><a href="#maps">Where things are</a></li>'
                 '<li><a href="#niches">Niches and niche regions</a></li>'
                 '<li><a href="#communities">Communities</a></li>'
                 '<li><a href="#quality">Data-quality flags</a></li></ol></nav>')

    # --- How to read ----------------------------------------------------------
    parts.append('<h2 id="read">How to read this</h2>')
    parts.append(
        '<p><b>"z-score"</b> puts every marker on the same ruler. Raw brightness can\'t be compared between '
        'markers (each is imaged at its own settings, and the slide is lit unevenly), so for each marker a '
        'typical cell is set to 0, and 1 is one "spread" of the data. Here the spread is the interquartile '
        'range of the arcsinh-transformed intensity, which resists outliers. A z-score of +2 means '
        '"noticeably brighter than the typical cell for that marker"; it does <i>not</i> mean the cell is '
        'definitively positive. That line is a separate judgement call.</p>')
    parts.append(
        '<p><b>Different kinds of grouping</b>, easy to confuse:</p><ul>'
        '<li><b>Clusters</b> (k-means) group cells by <i>what they express</i>: cells with similar marker '
        'profiles go together, wherever they sit in the tissue.</li>'
        '<li><b>Niches</b> group cells by <i>what surrounds them</i>: cells whose neighbourhoods have a similar mix '
        'of clusters go together, so a niche is a kind of place (e.g. "mostly neurons", or a tissue layer), which can '
        'occur in several parts of the tissue. <b>Niche regions</b> split each niche into its separate, connected '
        'pieces, so each region is one place in the tissue.</li>'
        '<li><b>Communities</b> group cells by <i>where they are</i>: neighbouring cells go together, '
        'whatever they express. A cell type that is physically scattered shows up in many communities. '
        'Communities and niche regions both divide the tissue into places, by different methods; both are shown '
        'while we learn which works better.</li>'
        '</ul><p>All of them come from the same data; they answer different questions.</p>')
    parts.append(f'<p class="muted">Cluster names are generated automatically from the (up to two) markers a '
                 f'cluster has raised (mean robust z ≥ {NAME_MIN_Z:g}); "low X" means no marker is raised but X '
                 f'is clearly below typical, and "no distinct marker" means neither. Names describe marker '
                 'profiles, not confirmed cell types, and can repeat when clusters differ mainly in overall '
                 'brightness; the cluster number always identifies them. The display scale here (robust z) '
                 'differs slightly from the plain z-score used for clustering.</p>')

    # --- Markers ----------------------------------------------------------------
    parts.append('<h2 id="markers">1. Markers, and what "positive" means</h2>')
    mt = res['marker_table']
    unimodal_note = ('Most markers show one broad pile with a tail rather than two clear peaks ("off" and '
                     '"on"), so ')
    raw = res['raw_table']
    controls_note = (' Control channels (pale grey) show what background and autofluorescence look like on '
                     'this slide.' if (raw['kind'] != 'marker').any() else '')
    parts.append('<h3>Raw intensities</h3>')
    parts.append(
        '<p>First, each marker\'s <b>raw</b> intensity (the mean over each cell), with every panel on the '
        'same log axis. This shows what the normalised view below hides: how bright each marker actually is. '
        f'The orange "positive" line and the dotted typical cell are the same cutoffs used below, converted back '
        f'to intensity.{controls_note}</p>')
    parts.append(figs['raw_distributions'])
    parts.append(
        '<div class="note">Raw levels aren\'t strictly comparable between markers (each is imaged at its own '
        'exposure and gain), but they show whether a marker has real signal. <b>If a marker\'s cells all sit at '
        'a low level, close to the control channels, with no bright tail, the marker is probably absent or '
        'failed</b>: its "positive" cells are just the brightest background noise and artefacts, and the robust '
        'z scale below stretches that noise to look like any other marker. The ratio of the "positive" line to '
        'the typical cell is another hint: a line only slightly above typical means "positive" cells are barely '
        'brighter than the rest.</div>')
    parts.append(_table(
        ['Channel', 'Typical cell (median)', '"Positive" line', 'Line ÷ typical', '99th percentile'],
        [(_esc(r.channel) + (' <span class="muted">(control)</span>' if r.kind != 'marker' else ''),
          f'{r.median:,.1f}', f'{r.positive_line:,.1f}', f'{r.ratio:.1f}×', f'{r.p99:,.1f}')
         for r in raw.itertuples()], num_cols=(1, 2, 3, 4)))
    parts.append('<h3>Normalised (robust z)</h3>')
    parts.append(f'<p>The same distributions on the robust z scale, which puts every marker on the same ruler. '
                 f'The <b>orange line</b> is the illustrative "positive" cutoff used throughout this report '
                 f'(z = {POSITIVE_Z:g}); the dotted line is the typical cell (z = 0). {unimodal_note}the positive '
                 'line is a <i>choice</i>, not a property of the data: move it and different cells flip sides. '
                 'Where the biology depends on it, choose the line deliberately for each marker, or treat the '
                 'marker as a continuous quantity. Because every marker is rescaled to the same spread, this view '
                 'makes absent markers look as "real" as bright ones; check the raw view above.</p>')
    parts.append(figs['marker_distributions'])
    parts.append(_table(['Marker', '"Positive" cells', '%'],
                        [(_esc(r.marker), f'{r.n_positive:,}', f'{r.pct_positive:.1f}%') for r in mt.itertuples()],
                        num_cols=(1, 2)))

    # --- Clusters ---------------------------------------------------------------
    parts.append('<h2 id="clusters">2. Cell populations (k-means clusters)</h2>')
    parts.append(f'<p>The cells fall into <b>{len(res["clusters"])} clusters</b>. The heatmap shows each '
                 'cluster\'s average level of every marker; the table lists its size and the five markers '
                 'that distinguish it most, i.e. furthest from typical in either direction (+ = higher, '
                 '− = lower than the typical cell).</p>')
    parts.append(figs['profile'])
    rows = []
    for c in res['clusters']:
        tops = ', '.join(f'{_esc(m)} ({z:+.2f})' for m, z in res['top'][c])
        flag = ' ⚠ non-specific' if c in res['flagged'] else ''
        rows.append((str(c), _esc(res['names'][c]) + flag, f"{int(res['sizes'][c]):,}",
                     f"{100 * res['sizes'][c] / n_cells:.1f}%", f"{100 * res['ns_frac'][c]:.0f}%", tops))
    parts.append(_table(['Cluster', 'Name', 'Cells', '%', 'Non-specific',
                         'Most distinctive markers (mean robust z, largest |z| first)'],
                        rows, num_cols=(2, 3, 4)))
    if res['flagged']:
        parts.append(f'<div class="warn"><b>Cluster(s) {", ".join(map(str, res["flagged"]))}</b> are mostly '
                     'non-specific cells (bright in many markers at once; likely autofluorescence, debris or '
                     'edge artefacts). Don\'t read them as real marker-positive populations.</div>')

    k_used = len(res['clusters'])
    csw = data['cluster_sweep']
    parts.append('<h3>Choosing the number of clusters</h3>')
    if csw is not None and 'cluster_sweep' in figs:
        parts.append('<p>k-means needs the number of clusters (k) chosen in advance. These are the results of '
                     're-running it with other values of k on the same features. <b>Inertia</b> always falls as '
                     'k grows; look for the "elbow" where it stops falling quickly. <b>Silhouette</b> (computed '
                     'on a random sample of cells) measures how well separated the clusters are: higher is '
                     'better, and values below about 0.25 mean the clusters overlap heavily (common for cell '
                     'data, where many populations grade into each other).</p>')
        parts.append(figs['cluster_sweep'])
        sw = csw.sort_values('n_clusters')
        cols = [c for c in ('n_clusters', 'silhouette', 'min_cluster_size', 'min_cluster_pct', 'inertia')
                if c in sw.columns]
        labels = {'n_clusters': 'k', 'silhouette': 'silhouette', 'min_cluster_size': 'smallest cluster (cells)',
                  'min_cluster_pct': 'smallest cluster (%)', 'inertia': 'inertia'}
        fmt = {'silhouette': lambda v: _fmt(v, 3), 'min_cluster_pct': lambda v: _fmt(v, 2),
               'inertia': lambda v: f'{v:,.0f}'}
        rows = []
        for _, r in sw.iterrows():
            row = [fmt.get(c, _fmt)(r[c]) for c in cols]
            if int(r['n_clusters']) == k_used:
                row[0] = f'<b>{row[0]} (this run)</b>'
            rows.append(row)
        parts.append(_table([labels[c] for c in cols], rows, num_cols=tuple(range(len(cols)))))
        if 'silhouette' in sw.columns and sw['silhouette'].notna().any():
            best = sw.loc[sw['silhouette'].idxmax()]
            here = sw[sw['n_clusters'] == k_used]
            text = f'The best-separated choice in this sweep is k = {int(best["n_clusters"])} ' \
                   f'(silhouette {best["silhouette"]:.3f})'
            if int(best['n_clusters']) == k_used:
                text += ', which is what this run used'
            elif len(here):
                text += f'; this run used k = {k_used} (silhouette {float(here["silhouette"].iloc[0]):.3f})'
            text += ('. Silhouette is only one guide: a higher k can still be useful if it separates a rare, '
                     'biologically meaningful population; check the cluster profiles above.')
            parts.append(f'<p>{text}</p>')
    else:
        parts.append('<p class="muted">No cluster sweep was found for this result. To see how the clusters change '
                     'with k, run <code>python -m cellsurvey.stability --zarr &lt;output&gt;_seg.zarr --mode '
                     'clusters</code>, then re-make this report.</p>')

    # --- Maps -------------------------------------------------------------------
    parts.append('<h2 id="maps">3. Where things are</h2>')
    bg = data['background_channel']
    parts.append('<p>Each cluster\'s cells drawn over the tissue image'
                 + (f' ({_esc(bg)}, greyscale)' if bg else ' (no image available: points only)')
                 + '. Large clusters are subsampled for drawing.</p>')
    parts.append(figs['cluster_maps'])
    parts.append('<p>The same, all at once: each tile is coloured by its most common cluster '
                 '(grey = flagged non-specific clusters).</p>')
    parts.append(figs['overview_map'])
    parts.append('<p>Cells above the "positive" line for each marker, independent of the clustering. '
                 'Positive cells that are also non-specific (bright in many markers) are drawn in grey, so '
                 'artefacts are easy to tell apart from real marker patterns.</p>')
    parts.append(figs['positive_maps'])

    # --- Communities (built here, then moved after the niches below) -------------
    comm_start = len(parts)
    parts.append('<h2 id="communities">5. Communities</h2>')
    if comm is None:
        parts.append('<p>No community labels in this result.</p>')
    else:
        parts.append(
            '<p>Communities are built in two steps: a network joins neighbouring cells (Delaunay '
            'triangulation, edges up to the maximum distance, each weighted by how similar the two cells\' '
            'features are), then <b>Louvain</b> community detection splits it at a <b>resolution</b> setting. '
            'Higher resolution gives more, smaller communities. Communities smaller than the minimum size are '
            'labelled "unassigned" (mostly isolated cells at tissue edges or gaps) and left out below. '
            'Communities are another way of dividing the tissue into places, alongside the niche regions in '
            'section 4; compare the two maps.</p>')
        sizes = comm['sizes']
        top10 = sizes.index[:10]
        dom_counts = comm['dominant'].loc[top10].value_counts()
        dom_text = ', '.join(f'{_esc(res["names"][c])} (cluster {c}: {n})' for c, n in dom_counts.items())
        n_pure = int((comm['dominant_share'] > PURE_COMMUNITY).sum())
        n_flag_dom = int(comm['dominant'].isin(res['flagged']).sum())
        cover = sizes.iloc[:MAX_COMMUNITY_PANELS].sum() / max(sizes.sum(), 1)
        bullets = [
            f'<b>{len(sizes)} communities</b> of at least {_fmt(get("min_community_size"))} cells; the '
            f'{min(MAX_COMMUNITY_PANELS, len(sizes))} largest hold {100 * cover:.0f}% of assigned cells.',
            f'The 10 largest are most often dominated by: {dom_text}.',
            f'{n_pure} of {len(sizes)} communities are more than {100 * PURE_COMMUNITY:.0f}% one cluster'
            + (' (the tissue is organised into fairly uniform blocks).' if n_pure > len(sizes) / 2
               else ' (most communities are mixtures of cell types).'),
        ]
        if res['flagged']:
            bullets.append(f'{n_flag_dom} communit{"y is" if n_flag_dom == 1 else "ies are"} dominated by '
                           'flagged non-specific clusters.')
        parts.append('<ul>' + ''.join(f'<li>{b}</li>' for b in bullets) + '</ul>')
        if len(markers) < 12:
            parts.append(f'<div class="note">With only {len(markers)} marker channels, the similarity weights on '
                         'the network carry little information, so communities mostly reflect distance (compact '
                         'patches of tissue) rather than tissue domains. Interpret them with care.</div>')
        parts.append(figs['outlines'])
        parts.append(figs['composition'])
        sweep = data['sweep']
        if sweep is not None and {'resolution', 'n_communities'} <= set(sweep.columns):
            cols = [c for c in ('resolution', 'max_edge_distance_um', 'n_communities', 'n_communities_raw',
                                'n_unassigned_cells') if c in sweep.columns]
            parts.append('<h3>Resolution sweep</h3><p>From the community sweep run on this result: how the number '
                         'of communities changes with the parameters.</p>')
            parts.append(_table(cols, [[_fmt(r[c]) for c in cols] for _, r in sweep.iterrows()],
                                num_cols=tuple(range(len(cols)))))

    # --- Niches -------------------------------------------------------------------
    niche_start = len(parts)
    parts.append('<h2 id="niches">4. Niches and niche regions</h2>')
    niche = res['niche']
    if niche is None:
        parts.append('<p>No niches in this result (run with <code>--n-niches</code> greater than 0).</p>')
    else:
        radius = get('niche_radius_um')
        parts.append(
            f'<p><b>Niches</b> group cells by the <i>mix of cell types around them</i>: for each cell, the share of '
            f'each k-means cluster among the cells within {_fmt(radius)} µm, then cells with similar mixes are '
            f'grouped. So a niche is a kind of <i>place</i> (e.g. "mostly cluster 3" or "half 1, half 4"), and '
            'unlike communities the same niche can occur in several parts of the tissue. They often follow '
            'anatomical regions and layers better than communities, which tend to be compact patches. Niches are '
            'numbered by size (0 = largest) and named after their most common clusters.</p>')
        parts.append(f'<div class="note">This result uses <b>{len(niche["ids"])} niches at {_fmt(radius)} µm</b>, '
                     f'built on the {_fmt(get("niche_n_clusters") or len(res["clusters"]))} k-means clusters above '
                     '(changing the clusters changes the niches). Suitable values depend heavily on the tissue, '
                     'panel, cell density and the structures of interest: treat them as a starting point and check '
                     'them with a niche sweep (<code>python -m cellsurvey.stability --zarr &lt;output&gt;_seg.zarr '
                     '--mode niches</code>).</div>')
        parts.append(figs['niche_map'])
        parts.append(figs['niche_composition'])
        rows = [(str(n), _esc(niche['names'][n]), f"{int(niche['sizes'][n]):,}",
                 f"{100 * niche['sizes'][n] / n_cells:.1f}%") for n in niche['ids']]
        parts.append(_table(['Niche', 'Made of (most common clusters)', 'Cells', '%'], rows, num_cols=(2, 3)))

        region = res['region']
        parts.append('<h3>Niche regions</h3>')
        if region is None:
            parts.append('<p class="muted">No niche regions in this result (they were added in a later CellSurvey '
                         'version; re-run with <code>--resume-from</code> to compute them).</p>')
        else:
            min_size = get('min_region_size')
            n_reg = len(region['sizes'])
            parts.append(
                '<p>A niche can occur in several places; <b>niche regions</b> split each niche into its separate, '
                'connected pieces, so each region is one place in the tissue (e.g. "this stretch of white matter"). '
                f'Pieces smaller than {_fmt(min_size)} cells are merged into the neighbouring region they border '
                'most, so a region can contain a few cells of other niches. Unlike communities, region borders '
                'follow the niche borders. Regions are numbered by size (0 = largest) and coloured by their niche; '
                'black lines mark the borders between regions.</p>')
            parts.append(f'<div class="note">This result has <b>{n_reg} regions</b> of at least {_fmt(min_size)} '
                         f'cells ({region["n_unassigned"]:,} cells unassigned). The minimum size '
                         '(<code>--min-region-size</code>) is a starting point: lower values give more, smaller '
                         'regions, higher values merge more. Regions inherit every niche setting above.</div>')
            parts.append(figs['region_map'])
            top_r = region['sizes'].index[:MAX_REGION_LABELS]
            rows = []
            for r in top_r:
                n = region['niche'][r]
                rows.append((str(r), f"{n}: {_esc(niche['names'][n])}", f"{int(region['sizes'][r]):,}",
                             f"{100 * region['sizes'][r] / n_cells:.1f}%", f"{100 * region['purity'][r]:.0f}%"))
            parts.append(_table(['Region', 'Niche', 'Cells', '% of cells', 'Cells in that niche'], rows,
                                num_cols=(2, 3, 4)))
            if n_reg > len(top_r):
                parts.append(f'<p class="muted">Showing the {len(top_r)} largest of {n_reg} regions.</p>')
            per = region['per_niche']
            parts.append('<p class="muted">Regions per niche: ' + ', '.join(
                f'niche {n}: {int(per[n])}' for n in per.index) + '.</p>')

        nsw = data['niche_sweep']
        parts.append('<h3>Choosing the radius and number of niches</h3>')
        if nsw is not None and 'niche_sweep' in figs:
            parts.append('<p>From the niche sweep run on this result. <b>Silhouette</b> (separation) rises with the '
                         'radius at every niche count, because bigger neighbourhoods average more cells and give '
                         'smoother compositions, so it can compare niche counts at one radius but <i>can\'t</i> '
                         'choose the radius. That is a choice of scale: larger radii give smoother, blockier regions; '
                         'smaller ones keep thin structures but get noisy. <b>Seed stability</b> shows whether the '
                         'same niches come back from different random starts (1 = identical); prefer settings near '
                         '1.</p>')
            parts.append(figs['niche_sweep'])
            sw = nsw.sort_values(['niche_radius_um', 'n_niches'])
            rows = []
            for _, r in sw.iterrows():
                row = [_fmt(r['niche_radius_um']), _fmt(r['n_niches']), _fmt(r.get('median_neighbours')),
                       _fmt(r['silhouette'], 3), _fmt(r['seed_ari'], 2), _fmt(r['min_niche_size'])]
                if radius is not None and np.isclose(r['niche_radius_um'], radius) \
                        and int(r['n_niches']) == len(niche['ids']):
                    row = [f'<b>{v}</b>' for v in row]
                    row[1] += ' <b>(this run)</b>'
                rows.append(row)
            parts.append(_table(['radius (µm)', 'niches', 'median neighbours', 'silhouette', 'seed stability',
                                 'smallest niche (cells)'], rows, num_cols=tuple(range(6))))
        else:
            parts.append('<p class="muted">No niche sweep was found for this result. Run <code>python -m '
                         'cellsurvey.stability --zarr &lt;output&gt;_seg.zarr --mode niches</code>, then re-make this '
                         'report.</p>')

    # Niches and regions (section 4) come before communities (section 5).
    parts[comm_start:] = parts[niche_start:] + parts[comm_start:niche_start]

    # --- Quality ------------------------------------------------------------------
    parts.append('<h2 id="quality">6. Data-quality flags</h2>')
    n_ns = int(res['nonspecific'].sum())
    mt = res['marker_table']
    if res['n_needed']:
        parts.append(f'<p>A cell is counted as <b>non-specific</b> if it is "positive" (robust z &gt; '
                     f'{POSITIVE_Z:g}) in at least {res["n_needed"]} of the {len(markers)} markers at once. Real '
                     'cells rarely express that many unrelated markers, so these are likely autofluorescent '
                     f'cells, debris or edge artefacts. <b>{n_ns:,} cells ({100 * n_ns / n_cells:.1f}%)</b> are '
                     'non-specific.</p>')
    else:
        parts.append('<p>Too few marker channels to assess non-specific cells.</p>')
    parts.append(figs['quality'])
    mt = res['marker_table']
    status_text = {'clean': 'clean', 'suspect': '⚠ suspect', 'dominated': '⚠⚠ dominated',
                   'too few': f'too few positive cells to judge (&lt; {MIN_POSITIVES})'}
    rows = [(_esc(r.marker), f'{r.n_positive:,}', f'{r.pct_positive:.1f}%',
             f'{r.pct_nonspecific:.0f}%' if r.n_positive else '—', status_text[r.status])
            for r in mt.itertuples()]
    parts.append('<p>For each marker: how many cells are above the "positive" line, and how many of those '
                 f'positives are non-specific cells (suspect above {100 * MARKER_SUSPECT:.0f}%, dominated above '
                 f'{100 * MARKER_DOMINATED:.0f}%). A suspect marker\'s "positive" cells shouldn\'t be counted '
                 'naively.</p>')
    parts.append(_table(['Marker', 'Positive cells', '%', 'of which non-specific', 'Status'], rows,
                        num_cols=(1, 2, 3)))
    notes = []
    if controls:
        notes.append(f'Control / autofluorescence channels found by name: {_esc(", ".join(controls))}. They are '
                     'left out of the marker sections above, but the clustering currently uses every channel, '
                     'so they can still shape the clusters (a channel-exclusion option is planned).')
    if comm is not None and comm['n_unassigned']:
        notes.append(f'{comm["n_unassigned"]:,} cells ({100 * comm["n_unassigned"] / n_cells:.1f}%) are in '
                     'communities below the minimum size ("unassigned").')
    if notes:
        parts.append('<ul>' + ''.join(f'<li>{n}</li>' for n in notes) + '</ul>')

    parts.append('<p class="muted">Produced by CellSurvey (<code>python -m cellsurvey.report</code>). '
                 'PNG copies of every figure are saved next to this file in <code>report/</code>.</p>')
    parts.append('</body></html>')
    return '\n'.join(parts)


# =============================================================================
# Entry points
# =============================================================================

def write_report(zarr_path, plot_dir, summary_path=None, filename='report.html', seed=42, sweep_dir=None):
    """Build ``<plot_dir>/report.html`` (and ``<plot_dir>/report/*.png``) from ``zarr_path``. Returns the path.

    ``sweep_dir``: where to look for this run's sweep CSVs (default: the Zarr's directory).
    """
    t0 = datetime.datetime.now()
    os.makedirs(plot_dir, exist_ok=True)
    if summary_path is None:
        summary_path = os.path.join(plot_dir, 'summary.json')
    print(f"Report: reading {zarr_path} ...")
    data = load_report_data(zarr_path, summary_path, sweep_dir)
    res = compute(data)
    rng = np.random.default_rng(seed)
    out = _Figures(os.path.join(plot_dir, 'report'))
    figs = {
        'raw_distributions': out.tag(fig_raw_distributions(data, res), 'raw_marker_distributions',
                                     'Raw marker intensities'),
        'marker_distributions': out.tag(fig_marker_distributions(res), 'marker_distributions',
                                        'Marker distributions'),
        'profile': out.tag(fig_cluster_profile(res), 'cluster_profile', 'Cluster marker profiles'),
        'cluster_maps': out.tag(fig_cluster_maps(data, res, rng), 'cluster_maps', 'Cluster maps', dpi=110),
        'overview_map': out.tag(fig_cluster_overview_map(data, res), 'cluster_overview_map', 'All clusters'),
        'positive_maps': out.tag(fig_positive_maps(data, res, rng), 'positive_maps', 'Positive cells per marker',
                                 dpi=110),
        'quality': out.tag(fig_quality(res), 'quality_clusters', 'Non-specific share per cluster'),
    }
    if data['cluster_sweep'] is not None:
        figs['cluster_sweep'] = out.tag(fig_cluster_sweep(data['cluster_sweep'], len(res['clusters'])),
                                        'cluster_sweep', 'Cluster sweep')
    if res['niche'] is not None:
        figs['niche_map'] = out.tag(fig_niche_map(data, res), 'niche_map', 'Niche map')
        figs['niche_composition'] = out.tag(fig_niche_composition(res), 'niche_composition', 'Niche composition')
        if res['region'] is not None:
            figs['region_map'] = out.tag(fig_region_map(data, res), 'niche_region_map', 'Niche regions')
        if data['niche_sweep'] is not None:
            radius = data['meta'].get('niche_radius_um', data['summary'].get('niche_radius_um'))
            figs['niche_sweep'] = out.tag(fig_niche_sweep(data['niche_sweep'], radius, len(res['niche']['ids'])),
                                          'niche_sweep', 'Niche sweep')
    if res['comm'] is not None:
        figs['outlines'] = out.tag(fig_community_outlines(data, res), 'community_outlines', 'Community outlines')
        figs['composition'] = out.tag(fig_community_composition(res), 'community_composition',
                                      'Community composition')
    path = os.path.join(plot_dir, filename)
    with open(path, 'w', encoding='utf-8') as f:
        f.write(build_html(data, res, figs))
    secs = (datetime.datetime.now() - t0).total_seconds()
    print(f"Report written to {path} ({os.path.getsize(path) / 1e6:.1f} MB, {secs:.0f} s)")
    return path


def main():
    parser = argparse.ArgumentParser(description="Write a self-contained HTML overview report from a "
                                                 "segmented Zarr (read-only).")
    parser.add_argument('--zarr', required=True, help='Path to the segmented Zarr (<prefix>_seg.zarr).')
    parser.add_argument('-p', '--plot-dir', default=None,
                        help='Output directory (default: the Zarr\'s directory). Its summary.json, if present, '
                             'supplies run parameters that older Zarrs don\'t store.')
    parser.add_argument('--summary', default=None,
                        help='Path to the run\'s summary.json (default: <plot-dir>/summary.json).')
    parser.add_argument('--output', default='report.html', help='Report file name (default: report.html).')
    parser.add_argument('--sweep-dir', default=None,
                        help='Directory holding this run\'s sweep CSVs from python -m cellsurvey.stability '
                             '(default: the Zarr\'s directory).')
    args = parser.parse_args()
    zarr_path = args.zarr.rstrip('/\\')
    if not os.path.exists(zarr_path):
        parser.error(f"Zarr not found: {zarr_path}")
    plot_dir = args.plot_dir or (os.path.dirname(os.path.abspath(zarr_path)) or '.')
    try:
        write_report(zarr_path, plot_dir, summary_path=args.summary, filename=args.output,
                     sweep_dir=args.sweep_dir)
    except ValueError as e:
        parser.error(str(e))


if __name__ == '__main__':
    main()
