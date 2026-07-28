import json
import numpy as np
import pandas as pd
from shapely.geometry import mapping, Point
from sopaspan.muspan_workflow import get_colors_for_communities


def export_to_qupath(domain, communities, clusters, output_path, sdata, intensity_df,
                     cell_id='cell_id', spots_with_cells=None):
    print("Exporting to QuPath GeoJSON format...")

    cell_ids = domain.labels[cell_id]['labels']
    communities_labels = domain.labels[communities]['labels']
    cluster_labels_from_domain = domain.labels[clusters]['labels']

    cell_to_community = dict(zip(cell_ids, communities_labels))
    cell_to_cluster = dict(zip(cell_ids, cluster_labels_from_domain))

    unique_communities = np.unique(communities_labels)
    community_colors = get_colors_for_communities(len(unique_communities))

    features = []

    print("Exporting cell boundaries with community labels...")
    boundaries = sdata.shapes['stardist_boundaries']
    for idx, row in boundaries.iterrows():
        community = cell_to_community.get(idx, -1)
        cluster_id = cell_to_cluster.get(idx, -1)

        if 0 <= community < len(community_colors):
            color = community_colors[int(community)]
            classification_name = f"Community_{int(community)}"
        else:
            color = [128, 128, 128]
            classification_name = "Unassigned"

        measurements_dict = {
            "Community ID": float(community),
            "Cluster ID": float(cluster_id)
        }

        if idx in intensity_df.index:
            for col in intensity_df.columns:
                if col != 'Cluster':
                    try:
                        val = intensity_df.loc[idx, col]
                        if pd.notna(val):
                            measurements_dict[f"Cell: {col} mean"] = float(val)
                    except (ValueError, TypeError, KeyError):
                        pass

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

            community = cell_to_community.get(cell_id_val, -1)
            if 0 <= community < len(community_colors):
                spot_color = community_colors[int(community)]
                spot_classification = f"Community_{int(community)}"
            else:
                spot_color = [128, 128, 128]
                spot_classification = "Unassigned"

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
