"""Deterministic, specification-native Cytoscape CX2 export."""

from __future__ import annotations

import json
import math
import numbers
import os
from pathlib import Path
from typing import Any, Mapping

from .artifact_names import portable_artifact_stem
from .visual_style import VISUAL_TUNING, get_network_visual_profile


def _missing_scalar(value: Any) -> bool:
    return value is None or value.__class__.__name__ == "NAType"


def _cx_attr_value(value: Any) -> str | int | float | bool:
    """Convert a dataframe value to one finite CX2 scalar."""
    if _missing_scalar(value):
        return ""
    if isinstance(value, bool):
        return value
    if isinstance(value, numbers.Integral):
        return int(value)
    if isinstance(value, numbers.Real):
        converted = float(value)
        if not math.isfinite(converted):
            raise ValueError("CX2 numeric attributes must be finite.")
        return converted
    return str(value)


def _cx_attr_type(value: Any) -> str:
    if isinstance(value, bool):
        return "boolean"
    if isinstance(value, numbers.Integral):
        return "integer"
    if isinstance(value, numbers.Real):
        return "double"
    return "string"


def _finite_coordinate(value: Any, *, node_id: str, axis: str) -> float:
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise ValueError(f"CX2 coordinate {axis} for node {node_id!r} must be numeric.")
    converted = float(value)
    if not math.isfinite(converted):
        raise ValueError(f"CX2 coordinate {axis} for node {node_id!r} must be finite.")
    return converted


def _records(frame: Any) -> list[dict[str, Any]]:
    if len(set(str(column) for column in frame.columns)) != len(frame.columns):
        raise ValueError("CX2 dataframes must not contain duplicate column names.")
    return [dict(record) for record in frame.to_dict(orient="records")]


def _attribute_types(records: list[Mapping[str, Any]], fields: list[str]) -> dict[str, str]:
    declarations: dict[str, str] = {}
    for field in fields:
        observed: set[str] = set()
        for record in records:
            converted = _cx_attr_value(record.get(field))
            if converted != "":
                observed.add(_cx_attr_type(converted))
        if len(observed) > 1:
            raise ValueError(
                f"CX2 attribute {field!r} contains incompatible scalar types: "
                f"{', '.join(sorted(observed))}."
            )
        declarations[field] = next(iter(observed), "string")
    return declarations


def _visual_properties(
    *,
    color_map: Mapping[str, str],
    profile: Mapping[str, Any],
) -> dict[str, Any]:
    discrete_map = [
        {"v": str(key), "vp": str(color_map[key])}
        for key in sorted(color_map, key=str)
    ]
    visual_props: dict[str, Any] = {
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
        item = visual_props["visualProperties"][0]
        item["nodeMapping"]["NODE_BORDER_COLOR"] = {
            "type": "PASSTHROUGH",
            "definition": {"attribute": "linked_identity_border_color", "type": "string"},
        }
        item["nodeMapping"]["NODE_BORDER_TRANSPARENCY"] = {
            "type": "PASSTHROUGH",
            "definition": {
                "attribute": "linked_identity_border_transparency",
                "type": "integer",
            },
        }
        mappings = {
            "EDGE_LINE_STYLE": (
                profile["identity_edge_style_headless"],
                profile["default_edge_style_headless"],
            ),
            "EDGE_LINE_COLOR": (
                profile["identity_edge_color"],
                profile["default_edge_color"],
            ),
            "EDGE_WIDTH": (
                profile["identity_edge_width"],
                profile["default_edge_width"],
            ),
            "EDGE_OPACITY": (
                profile["identity_edge_opacity_headless"],
                profile["default_edge_opacity_headless"],
            ),
        }
        for visual_name, (identity_value, interaction_value) in mappings.items():
            item["edgeMapping"][visual_name] = {
                "type": "DISCRETE",
                "definition": {
                    "attribute": "interaction",
                    "type": "string",
                    "map": [
                        {"v": "identity", "vp": identity_value},
                        {"v": "interacts_with", "vp": interaction_value},
                    ],
                },
            }
    return visual_props


def export_cx2_headless(
    network_title: str,
    run_output_path: str,
    nodes_df: Any,
    edges_df_for_export: Any,
    color_map: Mapping[str, str],
    positions: Mapping[str, Mapping[str, Any]],
) -> None:
    """Write one deterministic CX2 document with inline attributes/layout."""
    os.makedirs(run_output_path, exist_ok=True)
    out_path = Path(run_output_path) / f"{portable_artifact_stem(network_title)}.cx2"
    profile = get_network_visual_profile(network_title)

    raw_nodes = _records(nodes_df)
    if any("id" not in record for record in raw_nodes):
        raise ValueError("Every CX2 node must have an id.")
    raw_nodes.sort(key=lambda record: str(record["id"]))
    node_labels = [str(record["id"]) for record in raw_nodes]
    if len(set(node_labels)) != len(node_labels):
        raise ValueError("CX2 node IDs must be unique.")
    node_id_map = {label: index for index, label in enumerate(node_labels)}

    node_records_for_types: list[dict[str, Any]] = []
    nodes_aspect: list[dict[str, Any]] = []
    base_node_fields = [
        "name",
        "tooltip",
        "color_group",
        "uniprot_border_color",
        "linked_identity_border_color",
        "linked_identity_border_transparency",
    ]
    node_extra_fields = sorted(
        {
            str(column)
            for column in nodes_df.columns
            if str(column) != "id" and str(column) not in base_node_fields
        }
    )
    node_fields = base_node_fields + node_extra_fields
    for raw in raw_nodes:
        label = str(raw["id"])
        position = positions.get(label, {"x": 0.0, "y": 0.0})
        if not isinstance(position, Mapping):
            raise ValueError(f"CX2 layout position for node {label!r} must be an object.")
        border_color = (
            str(
                raw.get(
                    "linked_identity_border_color",
                    raw.get("uniprot_border_color", profile["node_border_color"]),
                )
            )
            if profile["is_linked_identity_network"]
            else str(profile["node_border_color"])
        )
        border_transparency = (
            int(
                raw.get(
                    "linked_identity_border_transparency",
                    profile["node_border_transparency_headless"],
                )
            )
            if profile["is_linked_identity_network"]
            else int(profile["node_border_transparency_headless"])
        )
        values: dict[str, Any] = {
            "name": str(raw.get("name", label)),
            "tooltip": str(raw.get("tooltip", "")),
            "color_group": str(raw.get("color_group", "Unknown")),
            "uniprot_border_color": str(raw.get("uniprot_border_color", border_color)),
            "linked_identity_border_color": border_color,
            "linked_identity_border_transparency": border_transparency,
        }
        for field in node_extra_fields:
            converted = _cx_attr_value(raw.get(field))
            if converted != "":
                values[field] = converted
        node_records_for_types.append(values)
        nodes_aspect.append(
            {
                "id": node_id_map[label],
                "x": _finite_coordinate(position.get("x", 0.0), node_id=label, axis="x"),
                "y": _finite_coordinate(position.get("y", 0.0), node_id=label, axis="y"),
                "v": values,
            }
        )

    raw_edges = _records(edges_df_for_export)
    edge_extra_fields = sorted(
        {
            str(column)
            for column in edges_df_for_export.columns
            if str(column) not in {"source", "target", "interaction", "all_atoms_count"}
        }
    )
    if "id" in edge_extra_fields:
        raise ValueError("CX2 edge attribute objects must not contain the reserved id field.")
    normalized_edges: list[tuple[tuple[str, str, str, str], dict[str, Any]]] = []
    seen_edges: set[tuple[str, str, str]] = set()
    for raw in raw_edges:
        source_label = str(raw.get("source", ""))
        target_label = str(raw.get("target", ""))
        if source_label not in node_id_map or target_label not in node_id_map:
            raise ValueError("CX2 edge endpoints must reference exported nodes.")
        if source_label == target_label:
            raise ValueError("CX2 self-loop edges are not supported by this exporter.")
        interaction = str(raw.get("interaction", "interacts_with"))
        duplicate_key = (*sorted((source_label, target_label)), interaction)
        if duplicate_key in seen_edges:
            raise ValueError("CX2 edges must not contain duplicate endpoint/interaction pairs.")
        seen_edges.add(duplicate_key)
        values: dict[str, Any] = {
            "interaction": interaction,
            "all_atoms_count": _cx_attr_value(raw.get("all_atoms_count", 1)),
        }
        if not isinstance(values["all_atoms_count"], int) or isinstance(
            values["all_atoms_count"], bool
        ):
            raise ValueError("CX2 all_atoms_count must be an integer.")
        for field in edge_extra_fields:
            converted = _cx_attr_value(raw.get(field))
            if converted != "":
                values[field] = converted
        sort_payload = json.dumps(values, ensure_ascii=False, sort_keys=True, allow_nan=False)
        normalized_edges.append(
            ((duplicate_key[0], duplicate_key[1], interaction, sort_payload), {
                "s": node_id_map[source_label],
                "t": node_id_map[target_label],
                "v": values,
            })
        )
    normalized_edges.sort(key=lambda item: item[0])
    edges_aspect = [
        {"id": index, **record}
        for index, (_key, record) in enumerate(normalized_edges)
    ]

    node_attr_types = _attribute_types(node_records_for_types, node_fields)
    edge_fields = ["interaction", "all_atoms_count"] + edge_extra_fields
    edge_attr_types = _attribute_types(
        [record["v"] for _key, record in normalized_edges], edge_fields
    )
    declarations = {
        "attributeDeclarations": [
            {
                "networkAttributes": {"name": {"d": "string"}},
                "nodes": {name: {"d": node_attr_types[name]} for name in node_fields},
                "edges": {name: {"d": edge_attr_types[name]} for name in edge_fields},
            }
        ]
    }
    visual_props = _visual_properties(color_map=color_map, profile=profile)
    metadata = {
        "metaData": [
            {"name": "attributeDeclarations", "elementCount": 1},
            {"name": "networkAttributes", "elementCount": 1},
            {"name": "nodes", "elementCount": len(nodes_aspect)},
            {"name": "edges", "elementCount": len(edges_aspect)},
            {
                "name": "visualProperties",
                "elementCount": len(visual_props["visualProperties"]),
            },
        ]
    }
    cx = [
        {"CXVersion": "2.0", "hasFragments": False},
        metadata,
        declarations,
        {"networkAttributes": [{"name": str(network_title)}]},
        {"nodes": nodes_aspect},
        {"edges": edges_aspect},
        visual_props,
        {"status": [{"error": "", "success": True}]},
    ]
    with out_path.open("x", encoding="utf-8") as handle:
        json.dump(cx, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write("\n")


__all__ = ["export_cx2_headless"]
