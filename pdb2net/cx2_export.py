"""Headless Cytoscape CX2 export."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

from .visual_style import VISUAL_TUNING, get_network_visual_profile


def _cx_attr_value(value: Any) -> str | int | float | bool:
    """Convert dataframe values to CX2-friendly scalar attributes."""
    if value is None:
        return ""
    try:
        if value != value:
            return ""
    except Exception:
        pass
    if isinstance(value, bool):
        return value
    if isinstance(value, int) and not isinstance(value, bool):
        return value
    if isinstance(value, float):
        return value
    return str(value)


def _cx_attr_type(value: Any) -> str:
    if isinstance(value, bool):
        return "boolean"
    if isinstance(value, int) and not isinstance(value, bool):
        return "integer"
    if isinstance(value, float):
        return "double"
    return "string"


def export_cx2_headless(
    network_title,
    run_output_path,
    nodes_df,
    edges_df_for_export,
    color_map,
    positions,
) -> None:
    """Write a portable CX2 file without requiring a running Cytoscape session."""
    os.makedirs(run_output_path, exist_ok=True)
    out_path = Path(run_output_path) / f"{network_title}.cx2"
    profile = get_network_visual_profile(network_title)

    node_ids = list(nodes_df["id"].astype(str))
    nid_map = {nid: index + 1 for index, nid in enumerate(node_ids)}

    nodes_aspect = []
    node_attributes_aspect = []
    cartesian_layout_aspect = []
    for _, row in nodes_df.iterrows():
        nid = str(row["id"])
        raw_node_id = nid_map[nid]
        pos = positions.get(nid, {"x": 0.0, "y": 0.0})
        border_color = (
            str(row.get("linked_identity_border_color", row.get("uniprot_border_color", profile["node_border_color"])))
            if profile["is_linked_identity_network"]
            else profile["node_border_color"]
        )
        border_transparency = (
            int(row.get("linked_identity_border_transparency", profile["node_border_transparency_headless"]))
            if profile["is_linked_identity_network"]
            else profile["node_border_transparency_headless"]
        )

        node_values = {
            "name": str(row.get("name", nid)),
            "tooltip": str(row.get("tooltip", "")),
            "color_group": str(row.get("color_group", "Unknown")),
            "uniprot_border_color": str(row.get("uniprot_border_color", border_color)),
            "linked_identity_border_color": border_color,
            "linked_identity_border_transparency": border_transparency,
        }
        for column in nodes_df.columns:
            if column == "id" or column in node_values:
                continue
            node_values[column] = _cx_attr_value(row.get(column))

        nodes_aspect.append(
            {
                "id": raw_node_id,
                "x": float(pos["x"]),
                "y": float(pos["y"]),
                "v": node_values,
            }
        )
        cartesian_layout_aspect.append(
            {"node": raw_node_id, "x": float(pos["x"]), "y": float(pos["y"])}
        )
        for attr_name, attr_value in nodes_aspect[-1]["v"].items():
            node_attributes_aspect.append({"po": raw_node_id, "n": attr_name, "v": attr_value})

    edges_aspect = []
    edge_attributes_aspect = []
    for index, row in edges_df_for_export.iterrows():
        source = nid_map.get(str(row["source"]))
        target = nid_map.get(str(row["target"]))
        if source is None or target is None or source == target:
            continue
        edge_id = index + 1
        edge_values = {
            "interaction": str(row.get("interaction", "interacts_with")),
            "all_atoms_count": int(row.get("all_atoms_count", 1)),
        }
        for column in edges_df_for_export.columns:
            if column in {"source", "target"} or column in edge_values:
                continue
            edge_values[column] = _cx_attr_value(row.get(column))

        edges_aspect.append(
            {
                "id": edge_id,
                "s": source,
                "t": target,
                "v": edge_values,
            }
        )
        for attr_name, attr_value in edges_aspect[-1]["v"].items():
            edge_attributes_aspect.append({"po": edge_id, "n": attr_name, "v": attr_value})

    discrete_map = [{"v": key, "vp": value} for key, value in color_map.items()]

    visual_props = {
        "visualProperties": [
            {
                "default": {
                    "network": {
                        "NETWORK_BACKGROUND_COLOR": VISUAL_TUNING["network_background_color"]
                    },
                    "node": {
                        "NODE_SHAPE": VISUAL_TUNING["node_shape"],
                        "NODE_WIDTH": VISUAL_TUNING["node_width"],
                        "NODE_HEIGHT": VISUAL_TUNING["node_height"],
                        "NODE_BORDER_WIDTH": profile["node_border_width"],
                        "NODE_BORDER_TRANSPARENCY": profile["node_border_transparency_headless"],
                        "NODE_BORDER_COLOR": profile["node_border_color"],
                    },
                    "edge": {
                        "EDGE_LINE_COLOR": profile["default_edge_color"],
                        "EDGE_OPACITY": profile["default_edge_opacity_headless"],
                        "EDGE_CURVED": True,
                        "EDGE_LINE_STYLE": profile["default_edge_style_headless"],
                        "EDGE_WIDTH": profile["default_edge_width"],
                    },
                },
                "nodeMapping": {
                    "NODE_LABEL": {
                        "type": "PASSTHROUGH",
                        "definition": {"attribute": "name", "type": "string"},
                    },
                    "NODE_TOOLTIP": {
                        "type": "PASSTHROUGH",
                        "definition": {"attribute": "tooltip", "type": "string"},
                    },
                    "NODE_BACKGROUND_COLOR": {
                        "type": "DISCRETE",
                        "definition": {
                            "attribute": "color_group",
                            "type": "string",
                            "map": discrete_map,
                        },
                    },
                },
                "edgeMapping": {},
            }
        ]
    }

    if profile["is_linked_identity_network"]:
        visual_props["visualProperties"][0]["nodeMapping"]["NODE_BORDER_COLOR"] = {
            "type": "PASSTHROUGH",
            "definition": {"attribute": "linked_identity_border_color", "type": "string"},
        }
        visual_props["visualProperties"][0]["nodeMapping"]["NODE_BORDER_TRANSPARENCY"] = {
            "type": "PASSTHROUGH",
            "definition": {"attribute": "linked_identity_border_transparency", "type": "integer"},
        }
        visual_props["visualProperties"][0]["edgeMapping"]["EDGE_LINE_STYLE"] = {
            "type": "DISCRETE",
            "definition": {
                "attribute": "interaction",
                "type": "string",
                "map": [
                    {"v": "identity", "vp": profile["identity_edge_style_headless"]},
                    {"v": "interacts_with", "vp": profile["default_edge_style_headless"]},
                ],
            },
        }
        visual_props["visualProperties"][0]["edgeMapping"]["EDGE_LINE_COLOR"] = {
            "type": "DISCRETE",
            "definition": {
                "attribute": "interaction",
                "type": "string",
                "map": [
                    {"v": "identity", "vp": profile["identity_edge_color"]},
                    {"v": "interacts_with", "vp": profile["default_edge_color"]},
                ],
            },
        }
        visual_props["visualProperties"][0]["edgeMapping"]["EDGE_WIDTH"] = {
            "type": "DISCRETE",
            "definition": {
                "attribute": "interaction",
                "type": "string",
                "map": [
                    {"v": "identity", "vp": profile["identity_edge_width"]},
                    {"v": "interacts_with", "vp": profile["default_edge_width"]},
                ],
            },
        }
        visual_props["visualProperties"][0]["edgeMapping"]["EDGE_OPACITY"] = {
            "type": "DISCRETE",
            "definition": {
                "attribute": "interaction",
                "type": "string",
                "map": [
                    {"v": "identity", "vp": profile["identity_edge_opacity_headless"]},
                    {"v": "interacts_with", "vp": profile["default_edge_opacity_headless"]},
                ],
            },
        }

    node_attr_types = {
        "name": "string",
        "tooltip": "string",
        "color_group": "string",
        "uniprot_border_color": "string",
        "linked_identity_border_color": "string",
        "linked_identity_border_transparency": "integer",
    }
    for column in nodes_df.columns:
        if column == "id" or column in node_attr_types:
            continue
        sample = ""
        for value in nodes_df[column].tolist():
            converted = _cx_attr_value(value)
            if converted != "":
                sample = converted
                break
        node_attr_types[column] = _cx_attr_type(sample)

    edge_attr_types = {
        "interaction": "string",
        "all_atoms_count": "integer",
    }
    for column in edges_df_for_export.columns:
        if column in {"source", "target"} or column in edge_attr_types:
            continue
        sample = ""
        for value in edges_df_for_export[column].tolist():
            converted = _cx_attr_value(value)
            if converted != "":
                sample = converted
                break
        edge_attr_types[column] = _cx_attr_type(sample)

    attr_decls = {
        "attributeDeclarations": [
            {
                "networkAttributes": {"name": {"d": "string"}},
                "nodes": {name: {"d": attr_type} for name, attr_type in node_attr_types.items()},
                "edges": {name: {"d": attr_type} for name, attr_type in edge_attr_types.items()},
            }
        ]
    }

    meta = {
        "metaData": [
            {"name": "attributeDeclarations", "elementCount": 1},
            {"name": "networkAttributes", "elementCount": 1},
            {"name": "nodes", "elementCount": len(nodes_aspect)},
            {"name": "edges", "elementCount": len(edges_aspect)},
            {"name": "nodeAttributes", "elementCount": len(node_attributes_aspect)},
            {"name": "edgeAttributes", "elementCount": len(edge_attributes_aspect)},
            {"name": "cartesianLayout", "elementCount": len(cartesian_layout_aspect)},
            {
                "name": "visualProperties",
                "elementCount": len(visual_props.get("visualProperties", [])),
            },
        ]
    }

    cx = [
        {"CXVersion": "2.0", "hasFragments": False},
        meta,
        attr_decls,
        {"networkAttributes": [{"name": network_title}]},
        {"nodes": nodes_aspect},
        {"edges": edges_aspect},
        {"nodeAttributes": node_attributes_aspect},
        {"edgeAttributes": edge_attributes_aspect},
        {"cartesianLayout": cartesian_layout_aspect},
        visual_props,
        {"status": [{"error": "", "success": True}]},
    ]

    with out_path.open("w", encoding="utf-8") as handle:
        json.dump(cx, handle, ensure_ascii=False, indent=2)
