import json
import numpy as np
import pandas as pd
from shapely.geometry import mapping, Point
from cellsurvey.utils import get_colors_for_communities, resolve_shape_keys


def export_to_qupath(cell_ids, community_labels, cluster_labels, output_path, sdata, intensity_df,
                     spots_with_cells=None, cell_boundaries=None, counts_df=None, niche_labels=None,
                     region_labels=None):
    """Write cells (and spots) as QuPath GeoJSON.

    Each cell gets ``Cell: <channel> mean`` measurements from ``intensity_df``
    and, if ``counts_df`` is given (runs with spots), ``Cell: <channel> spot
    count`` measurements (raw counts). ``niche_labels`` and ``region_labels``
    (dicts: cell id as a string → label) add ``Niche ID`` and ``Region ID``
    measurements (-1 if the cell has none).

    Classification: by niche region (``Region_<id>``, coloured by the region's
    niche, as in the report's region map) when ``region_labels`` has any
    region; otherwise by community (``Community_<id>``). Cells without one are
    ``Unassigned`` (grey). ``Community ID`` is always a measurement.
    """
    if cell_boundaries is None:
        cell_boundaries = resolve_shape_keys(sdata)[1]
    print("Exporting to QuPath GeoJSON format...")

    cell_to_community = dict(zip(cell_ids, community_labels))
    cell_to_cluster = dict(zip(cell_ids, cluster_labels))

    # Communities are numbered 0..n-1; -1 (unassigned) is grey below.
    community_labels = np.asarray(community_labels)
    n_communities = int(community_labels.max()) + 1 if (community_labels >= 0).any() else 0
    community_colors = get_colors_for_communities(n_communities)

    by_region = region_labels is not None and any(v >= 0 for v in region_labels.values())
    region_colours = _region_colours(region_labels, niche_labels) if by_region else {}

    def classify(cell_id):
        """(classification name, RGB) for a cell id (as stored in the shapes / spots)."""
        if by_region:
            r = int(region_labels.get(str(cell_id), -1))
            if r >= 0:
                return f"Region_{r}", region_colours[r]
        else:
            community = cell_to_community.get(cell_id, -1)
            if 0 <= community < len(community_colors):
                return f"Community_{int(community)}", community_colors[int(community)]
        return "Unassigned", [128, 128, 128]

    features = []

    print(f"Exporting cell boundaries, classified by {'niche region' if by_region else 'community'}...")
    boundaries = sdata.shapes[cell_boundaries]
    for idx, row in boundaries.iterrows():
        community = cell_to_community.get(idx, -1)
        cluster_id = cell_to_cluster.get(idx, -1)
        classification_name, color = classify(idx)

        measurements_dict = {
            "Community ID": float(community),
            "Cluster ID": float(cluster_id)
        }
        if niche_labels is not None:
            measurements_dict["Niche ID"] = float(niche_labels.get(str(idx), -1))
        if region_labels is not None:
            measurements_dict["Region ID"] = float(region_labels.get(str(idx), -1))

        if idx in intensity_df.index:
            for col in intensity_df.columns:
                if col != 'Cluster':
                    try:
                        val = intensity_df.loc[idx, col]
                        if pd.notna(val):
                            measurements_dict[f"Cell: {col} mean"] = float(val)
                    except (ValueError, TypeError, KeyError):
                        pass

        if counts_df is not None and idx in counts_df.index:
            for col, val in counts_df.loc[idx].items():
                if pd.notna(val):
                    measurements_dict[f"Cell: {col} spot count"] = float(val)

        feature = {
            "type": "Feature",
            "id": str(idx),
            "geometry": mapping(row.geometry),
            "properties": {
                "objectType": "cell",
                "classification": {
                    "name": classification_name,
                    "color": color
                },
                "measurements": measurements_dict
            }
        }
        features.append(feature)

    print(f"Exported {len(boundaries)} cell boundaries")

    if spots_with_cells is not None:
        print("Exporting spot detections with cell assignments...")
        for idx, row in spots_with_cells.iterrows():
            cell_id_val = row['cell_id']
            sigma = float(row['sigma'])

            spot_classification, spot_color = classify(cell_id_val)

            circle = Point(float(row['x']), float(row['y'])).buffer(sigma)

            feature = {
                "type": "Feature",
                "geometry": mapping(circle),
                "properties": {
                    "objectType": "detection",
                    "parent_id": str(cell_id_val) if cell_id_val is not None else None,
                    "channel_name": str(row['gene']),
                    "classification": {
                        "name": spot_classification,
                        "color": spot_color
                    },
                    "measurements": {
                        "sigma": sigma,
                    }
                }
            }
            features.append(feature)

        assigned = spots_with_cells['cell_id'].notna().sum()
        print(f"Exported {len(spots_with_cells)} spots ({assigned} assigned to cells, "
              f"{len(spots_with_cells) - assigned} unassigned)")

    geojson = {
        "type": "FeatureCollection",
        "features": features
    }

    with open(output_path, 'w') as f:
        json.dump(geojson, f, indent=2)

    print(f"\n\u2713 Exported {len(features)} total features")
    print(f"Output file: {output_path}")


def _region_colours(region_labels, niche_labels):
    """RGB (0\u2013255) per region: the colour of its majority niche, as in the report's region map."""
    from cellsurvey.report import NICHE_COLOURS
    counts = {}
    for cid, r in region_labels.items():
        n = int(niche_labels.get(cid, -1)) if niche_labels is not None else -1
        if r >= 0 and n >= 0:
            counts.setdefault(int(r), {}).setdefault(n, 0)
            counts[int(r)][n] += 1
    colours = {}
    for r in {int(v) for v in region_labels.values() if v >= 0}:
        niche = max(counts[r], key=counts[r].get) if r in counts else 0
        rgb = NICHE_COLOURS[niche % len(NICHE_COLOURS)]
        colours[r] = [int(round(255 * c)) for c in rgb[:3]]
    return colours
