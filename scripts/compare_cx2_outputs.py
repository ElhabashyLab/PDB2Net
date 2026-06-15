"""Inspect and compare Cytoscape CX2 output files semantically.

The goal is to compare exported network meaning rather than JSON bytes. CX2
files can differ in aspect order, object order, metadata, timestamps, or export
IDs while still representing the same PDB2Net network.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import defaultdict, deque
from pathlib import Path
from typing import Any


ANNOTATION_FIELDS = [
    "name",
    "label",
    "molecule_type",
    "node_type",
    "uniprot_id",
    "color_group",
    "node_color",
]

EDGE_VALUE_FIELDS = [
    "interaction",
    "interaction_type",
    "weight",
    "all_atoms_count",
    "ca_count",
    "nn_count",
    "count",
]

STYLE_DEFAULT_PROPERTIES = {
    "network": [
        "NETWORK_BACKGROUND_COLOR",
    ],
    "node": [
        "NODE_BACKGROUND_COLOR",
        "NODE_FILL_COLOR",
        "NODE_BORDER_COLOR",
        "NODE_BORDER_PAINT",
        "NODE_BORDER_WIDTH",
        "NODE_BORDER_OPACITY",
        "NODE_BORDER_TRANSPARENCY",
        "NODE_SHAPE",
        "NODE_WIDTH",
        "NODE_HEIGHT",
    ],
    "edge": [
        "EDGE_LINE_COLOR",
        "EDGE_STROKE_UNSELECTED_PAINT",
        "EDGE_WIDTH",
        "EDGE_LINE_WIDTH",
        "EDGE_OPACITY",
        "EDGE_TRANSPARENCY",
        "EDGE_LINE_STYLE",
        "EDGE_CURVED",
    ],
}

STYLE_MAPPING_PROPERTIES = {
    "node": [
        "NODE_BACKGROUND_COLOR",
        "NODE_FILL_COLOR",
        "NODE_BORDER_COLOR",
        "NODE_BORDER_PAINT",
        "NODE_BORDER_TRANSPARENCY",
        "NODE_LABEL",
        "NODE_TOOLTIP",
    ],
    "edge": [
        "EDGE_LINE_COLOR",
        "EDGE_STROKE_UNSELECTED_PAINT",
        "EDGE_WIDTH",
        "EDGE_LINE_WIDTH",
        "EDGE_OPACITY",
        "EDGE_TRANSPARENCY",
        "EDGE_LINE_STYLE",
    ],
}

NODE_COLOR_FIELDS = [
    "node_color",
    "color",
    "fill_color",
    "NODE_BACKGROUND_COLOR",
    "NODE_FILL_COLOR",
]

IGNORED_ASPECTS = {
    "CXVersion",
    "metaData",
    "status",
}


def load_cx2(path: str | Path) -> Any:
    """Load a CX2 JSON file."""
    with Path(path).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _as_aspect_map(cx2_data: Any) -> dict[str, list[Any]]:
    """Return a map from CX2 aspect name to aspect payload entries."""
    aspects: dict[str, list[Any]] = defaultdict(list)

    if isinstance(cx2_data, dict):
        iterable = [{key: value} for key, value in cx2_data.items()]
    elif isinstance(cx2_data, list):
        iterable = cx2_data
    else:
        return {}

    for item in iterable:
        if not isinstance(item, dict):
            continue
        for key, value in item.items():
            if isinstance(value, list):
                aspects[key].extend(value)
            else:
                aspects[key].append(value)

    return dict(aspects)


def _stringify(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text if text else None


def _clean_value(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, float):
        return round(value, 6)
    if isinstance(value, (str, int, bool)):
        return value
    return json.dumps(value, sort_keys=True, ensure_ascii=False)


def _attrs_from_record(record: dict[str, Any]) -> dict[str, Any]:
    attrs: dict[str, Any] = {}
    embedded = record.get("v")
    if isinstance(embedded, dict):
        attrs.update(embedded)

    for key, value in record.items():
        if key not in {"id", "@id", "s", "t", "x", "y", "z", "v"}:
            attrs.setdefault(key, value)
    return attrs


def _apply_attribute_aspect(
    target: dict[str, dict[str, Any]],
    records: list[Any],
) -> None:
    """Apply separate CX/CX2 attribute aspect rows to node or edge attrs."""
    for record in records:
        if not isinstance(record, dict):
            continue
        owner = (
            record.get("po")
            or record.get("propertyOf")
            or record.get("id")
            or record.get("@id")
        )
        name = record.get("n") or record.get("name")
        if owner is None or name is None:
            continue
        value = record.get("v", record.get("value"))
        target.setdefault(str(owner), {})[str(name)] = value


def _extract_type_from_tooltip(tooltip: str | None) -> str | None:
    if not tooltip:
        return None
    for line in tooltip.splitlines():
        if line.strip().lower().startswith("type:"):
            return line.split(":", 1)[1].strip() or None
    return None


def _extract_uniprot_from_tooltip(tooltip: str | None) -> str | None:
    if not tooltip:
        return None
    for line in tooltip.splitlines():
        if line.strip().lower().startswith("uniprot:"):
            return line.split(":", 1)[1].strip() or None
    return None


def _node_color(attrs: dict[str, Any]) -> Any:
    for key in NODE_COLOR_FIELDS:
        if key in attrs:
            return attrs[key]
    return None


def _node_stable_id(raw_id: str, attrs: dict[str, Any]) -> str:
    """Pick a semantic node identifier, preferring visible stable labels."""
    for key in ("name", "label", "shared_name"):
        value = _stringify(attrs.get(key))
        if value:
            return value
    return raw_id


def _edge_interaction(attrs: dict[str, Any]) -> str:
    return str(
        attrs.get("interaction")
        or attrs.get("interaction_type")
        or attrs.get("type")
        or "interacts_with"
    )


def _edge_values(attrs: dict[str, Any]) -> dict[str, Any]:
    values = {}
    for key in EDGE_VALUE_FIELDS:
        if key in attrs:
            values[key] = _clean_value(attrs[key])
    return values


def _normalise_style_value(value: Any) -> Any:
    """Canonicalize visual style values while preserving useful semantics."""
    if isinstance(value, dict):
        return {
            key: _normalise_style_value(val)
            for key, val in sorted(value.items())
            if key not in {"timestamp", "created", "modified"}
        }
    if isinstance(value, list):
        return [_normalise_style_value(item) for item in value]
    return _clean_value(value)


def _extract_visual_properties(aspects: dict[str, list[Any]]) -> list[Any]:
    styles = []
    for value in aspects.get("visualProperties", []):
        styles.append(_normalise_style_value(value))
    return sorted(styles, key=lambda item: json.dumps(item, sort_keys=True))


def _mapping_attribute(mapping: dict[str, Any]) -> str | None:
    definition = mapping.get("definition")
    if isinstance(definition, dict):
        return _stringify(definition.get("attribute") or definition.get("mappingColumn"))
    return _stringify(mapping.get("mappingColumn") or mapping.get("attribute"))


def _mapping_rows(mapping: dict[str, Any]) -> list[Any]:
    definition = mapping.get("definition")
    if isinstance(definition, dict) and isinstance(definition.get("map"), list):
        return definition["map"]
    if isinstance(mapping.get("map"), list):
        return mapping["map"]
    return []


def _extract_discrete_mapping(mapping: Any) -> dict[str, Any] | None:
    """Return a canonical discrete mapping, or None for passthrough/unsupported mappings."""
    if not isinstance(mapping, dict):
        return None

    mapping_kind = _stringify(mapping.get("type") or mapping.get("mappingType"))
    if mapping_kind and mapping_kind.lower() not in {"discrete"}:
        return None

    rows = _mapping_rows(mapping)
    if not rows:
        return None

    values: dict[str, Any] = {}
    for row in rows:
        if not isinstance(row, dict):
            continue
        key = row.get("v", row.get("key"))
        value = row.get("vp", row.get("value"))
        if key is None:
            continue
        values[str(key)] = _clean_value(value)

    return {
        "attribute": _mapping_attribute(mapping),
        "values": dict(sorted(values.items())),
    }


def _extract_style_semantics(aspects: dict[str, list[Any]]) -> dict[str, Any]:
    """Extract style defaults and mappings that carry visual semantics."""
    semantics: dict[str, Any] = {
        "defaults": {"network": {}, "node": {}, "edge": {}},
        "mappings": {"node": {}, "edge": {}},
    }

    for style in aspects.get("visualProperties", []):
        if not isinstance(style, dict):
            continue

        defaults = style.get("default")
        if isinstance(defaults, dict):
            for scope, properties in STYLE_DEFAULT_PROPERTIES.items():
                scoped_defaults = defaults.get(scope)
                if not isinstance(scoped_defaults, dict):
                    continue
                for prop in properties:
                    if prop in scoped_defaults:
                        semantics["defaults"][scope][prop] = _clean_value(scoped_defaults[prop])

        for scope in ("node", "edge"):
            mapping_block = style.get(f"{scope}Mapping")
            if not isinstance(mapping_block, dict):
                continue
            for prop in STYLE_MAPPING_PROPERTIES[scope]:
                mapping = _extract_discrete_mapping(mapping_block.get(prop))
                if mapping is not None:
                    semantics["mappings"][scope][prop] = mapping

    return semantics


def _apply_layout_aspects(
    positions: dict[str, dict[str, float]],
    aspects: dict[str, list[Any]],
) -> None:
    """Apply common CX/CX2 cartesian layout aspect variants if present."""
    for aspect_name in ("cartesianLayout", "cyVisualProperties"):
        for record in aspects.get(aspect_name, []):
            if not isinstance(record, dict):
                continue
            raw_id = (
                record.get("node")
                or record.get("nodeId")
                or record.get("id")
                or record.get("@id")
            )
            if raw_id is None:
                continue
            if "x" in record and "y" in record:
                try:
                    positions[str(raw_id)] = {
                        "x": float(record["x"]),
                        "y": float(record["y"]),
                    }
                except (TypeError, ValueError):
                    pass


def extract_canonical_network(cx2_data: Any) -> dict[str, Any]:
    """Extract a canonical semantic representation from CX2 JSON data."""
    aspects = _as_aspect_map(cx2_data)

    node_attrs_by_raw_id: dict[str, dict[str, Any]] = {}
    edge_attrs_by_raw_id: dict[str, dict[str, Any]] = {}
    _apply_attribute_aspect(node_attrs_by_raw_id, aspects.get("nodeAttributes", []))
    _apply_attribute_aspect(edge_attrs_by_raw_id, aspects.get("edgeAttributes", []))

    raw_to_stable: dict[str, str] = {}
    positions_by_raw_id: dict[str, dict[str, float]] = {}
    nodes: dict[str, dict[str, Any]] = {}

    for record in aspects.get("nodes", []):
        if not isinstance(record, dict):
            continue
        raw_id = str(record.get("id", record.get("@id", "")))
        if not raw_id:
            continue
        attrs = _attrs_from_record(record)
        attrs.update(node_attrs_by_raw_id.get(raw_id, {}))

        if "x" in record and "y" in record:
            try:
                positions_by_raw_id[raw_id] = {
                    "x": float(record["x"]),
                    "y": float(record["y"]),
                }
            except (TypeError, ValueError):
                pass

        tooltip = _stringify(attrs.get("tooltip"))
        stable_id = _node_stable_id(raw_id, attrs)
        raw_to_stable[raw_id] = stable_id
        nodes[stable_id] = {
            "id": stable_id,
            "raw_id": raw_id,
            "name": _stringify(attrs.get("name")),
            "label": _stringify(attrs.get("label") or attrs.get("shared_name")),
            "molecule_type": _stringify(attrs.get("molecule_type"))
            or _extract_type_from_tooltip(tooltip),
            "node_type": _stringify(attrs.get("node_type") or attrs.get("type")),
            "uniprot_id": _stringify(attrs.get("uniprot_id"))
            or _extract_uniprot_from_tooltip(tooltip),
            "color_group": _stringify(attrs.get("color_group")),
            "node_color": _clean_value(_node_color(attrs)),
            "attributes": {
                str(key): _clean_value(value)
                for key, value in sorted(attrs.items())
                if key not in {"tooltip"}
            },
        }

    _apply_layout_aspects(positions_by_raw_id, aspects)

    positions = {
        raw_to_stable[raw_id]: pos
        for raw_id, pos in positions_by_raw_id.items()
        if raw_id in raw_to_stable
    }

    edges: dict[str, dict[str, Any]] = {}
    for record in aspects.get("edges", []):
        if not isinstance(record, dict):
            continue
        raw_id = str(record.get("id", record.get("@id", "")))
        attrs = _attrs_from_record(record)
        if raw_id:
            attrs.update(edge_attrs_by_raw_id.get(raw_id, {}))

        source_raw = record.get("s", record.get("source"))
        target_raw = record.get("t", record.get("target"))
        if source_raw is None or target_raw is None:
            continue
        source = raw_to_stable.get(str(source_raw), str(source_raw))
        target = raw_to_stable.get(str(target_raw), str(target_raw))
        interaction = _edge_interaction(attrs)
        a, b = sorted([source, target])
        key = f"{a} -- {b} [{interaction}]"
        edges[key] = {
            "key": key,
            "source": a,
            "target": b,
            "interaction": interaction,
            "values": _edge_values(attrs),
            "attributes": {
                str(k): _clean_value(v)
                for k, v in sorted(attrs.items())
            },
        }

    network_name = None
    for record in aspects.get("networkAttributes", []):
        if isinstance(record, dict) and record.get("name"):
            network_name = str(record["name"])
            break

    return {
        "network_name": network_name,
        "nodes": dict(sorted(nodes.items())),
        "edges": dict(sorted(edges.items())),
        "positions": dict(sorted(positions.items())),
        "styles": _extract_visual_properties(aspects),
        "style_semantics": _extract_style_semantics(aspects),
    }


def compute_layout_metrics(canonical_network: dict[str, Any], overlap_threshold: float = 20.0) -> dict[str, Any]:
    """Compute basic graph and layout diagnostics for a canonical network."""
    nodes = canonical_network["nodes"]
    edges = canonical_network["edges"]
    positions = canonical_network["positions"]

    adjacency: dict[str, set[str]] = {node_id: set() for node_id in nodes}
    for edge in edges.values():
        source = edge["source"]
        target = edge["target"]
        adjacency.setdefault(source, set()).add(target)
        adjacency.setdefault(target, set()).add(source)

    seen = set()
    components = 0
    for node_id in adjacency:
        if node_id in seen:
            continue
        components += 1
        queue = deque([node_id])
        seen.add(node_id)
        while queue:
            current = queue.popleft()
            for neighbor in adjacency.get(current, set()):
                if neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)

    xs = [pos["x"] for pos in positions.values()]
    ys = [pos["y"] for pos in positions.values()]
    bbox = None
    if xs and ys:
        bbox = {
            "min_x": min(xs),
            "max_x": max(xs),
            "min_y": min(ys),
            "max_y": max(ys),
            "width": max(xs) - min(xs),
            "height": max(ys) - min(ys),
        }

    min_distance = None
    overlap_pairs: list[dict[str, Any]] = []
    positioned = sorted(positions.items())
    for i, (left_id, left) in enumerate(positioned):
        for right_id, right in positioned[i + 1 :]:
            distance = math.hypot(left["x"] - right["x"], left["y"] - right["y"])
            if min_distance is None or distance < min_distance:
                min_distance = distance
            if distance < overlap_threshold:
                overlap_pairs.append({
                    "source": left_id,
                    "target": right_id,
                    "distance": round(distance, 3),
                })

    return {
        "node_count": len(nodes),
        "edge_count": len(edges),
        "positioned_node_count": len(positions),
        "connected_components": components,
        "bounding_box": bbox,
        "minimum_node_distance": round(min_distance, 3) if min_distance is not None else None,
        "overlap_threshold": overlap_threshold,
        "possible_node_overlap_count": len(overlap_pairs),
        "possible_node_overlaps": overlap_pairs[:25],
    }


def _field_status(actual_value: Any, expected_value: Any) -> str:
    if actual_value == expected_value:
        return "PASS"
    if actual_value is None or expected_value is None:
        return "WARN"
    return "FAIL"


def _format_style_path(path: tuple[str, ...]) -> str:
    return ".".join(path)


def _compare_style_values(
    actual_value: Any,
    expected_value: Any,
    path: tuple[str, ...],
    differences: list[dict[str, Any]],
) -> None:
    """Compare nested style semantics and record every meaningful difference as FAIL."""
    if isinstance(actual_value, dict) and isinstance(expected_value, dict):
        actual_keys = set(actual_value)
        expected_keys = set(expected_value)
        for key in sorted(expected_keys - actual_keys):
            style_path = path + (str(key),)
            differences.append({
                "status": "FAIL",
                "path": _format_style_path(style_path),
                "actual": None,
                "expected": expected_value[key],
                "message": f"Missing expected style semantic {_format_style_path(style_path)!r}",
            })
        for key in sorted(actual_keys - expected_keys):
            style_path = path + (str(key),)
            differences.append({
                "status": "FAIL",
                "path": _format_style_path(style_path),
                "actual": actual_value[key],
                "expected": None,
                "message": f"Extra actual style semantic {_format_style_path(style_path)!r}",
            })
        for key in sorted(actual_keys & expected_keys):
            _compare_style_values(
                actual_value[key],
                expected_value[key],
                path + (str(key),),
                differences,
            )
        return

    if actual_value != expected_value:
        style_path = _format_style_path(path)
        differences.append({
            "status": "FAIL",
            "path": style_path,
            "actual": actual_value,
            "expected": expected_value,
            "message": f"Style semantic {style_path!r} changed",
        })


def _compare_style_semantics(
    actual_styles: dict[str, Any],
    expected_styles: dict[str, Any],
) -> list[dict[str, Any]]:
    """Compare visual semantics such as color mappings, widths, and line styles."""
    differences: list[dict[str, Any]] = []
    _compare_style_values(actual_styles, expected_styles, ("style_semantics",), differences)
    return differences


def compare_networks(actual: dict[str, Any], expected: dict[str, Any]) -> dict[str, Any]:
    """Compare two canonical networks and classify semantic differences."""
    result: dict[str, Any] = {
        "status": "PASS",
        "failures": [],
        "warnings": [],
        "node_differences": [],
        "edge_differences": [],
        "style_differences": [],
        "layout_warnings": [],
    }

    actual_nodes = set(actual["nodes"])
    expected_nodes = set(expected["nodes"])
    missing_nodes = sorted(expected_nodes - actual_nodes)
    extra_nodes = sorted(actual_nodes - expected_nodes)
    if missing_nodes:
        result["failures"].append(f"Missing expected nodes: {', '.join(missing_nodes)}")
    if extra_nodes:
        result["failures"].append(f"Extra actual nodes: {', '.join(extra_nodes)}")

    for node_id in sorted(actual_nodes & expected_nodes):
        actual_node = actual["nodes"][node_id]
        expected_node = expected["nodes"][node_id]
        for field in ANNOTATION_FIELDS:
            status = _field_status(actual_node.get(field), expected_node.get(field))
            if status == "PASS":
                continue
            diff = {
                "node": node_id,
                "field": field,
                "actual": actual_node.get(field),
                "expected": expected_node.get(field),
                "status": status,
            }
            result["node_differences"].append(diff)
            if status == "FAIL":
                result["failures"].append(f"Node {node_id!r} changed {field!r}")
            else:
                result["warnings"].append(f"Node {node_id!r} has one-sided {field!r}")

    actual_edges = set(actual["edges"])
    expected_edges = set(expected["edges"])
    missing_edges = sorted(expected_edges - actual_edges)
    extra_edges = sorted(actual_edges - expected_edges)
    if missing_edges:
        result["failures"].append(f"Missing expected edges: {', '.join(missing_edges)}")
    if extra_edges:
        result["failures"].append(f"Extra actual edges: {', '.join(extra_edges)}")

    for edge_key in sorted(actual_edges & expected_edges):
        actual_values = actual["edges"][edge_key]["values"]
        expected_values = expected["edges"][edge_key]["values"]
        for field in sorted(set(actual_values) | set(expected_values)):
            status = _field_status(actual_values.get(field), expected_values.get(field))
            if status == "PASS":
                continue
            diff = {
                "edge": edge_key,
                "field": field,
                "actual": actual_values.get(field),
                "expected": expected_values.get(field),
                "status": status,
            }
            result["edge_differences"].append(diff)
            if status == "FAIL":
                result["failures"].append(f"Edge {edge_key!r} changed {field!r}")
            else:
                result["warnings"].append(f"Edge {edge_key!r} has one-sided {field!r}")

    style_differences = _compare_style_semantics(
        actual.get("style_semantics", {}),
        expected.get("style_semantics", {}),
    )
    for diff in style_differences:
        result["style_differences"].append(diff)
        result["failures"].append(diff["message"])

    if actual.get("styles") != expected.get("styles") and not style_differences:
        result["style_differences"].append({
            "status": "WARN",
            "message": "Visual properties differ",
        })
        result["warnings"].append("Visual properties differ")

    actual_positions = set(actual["positions"])
    expected_positions = set(expected["positions"])
    if actual_positions != expected_positions:
        missing_positions = sorted(expected_positions - actual_positions)
        extra_positions = sorted(actual_positions - expected_positions)
        result["layout_warnings"].append({
            "status": "WARN",
            "message": "Positioned node set differs",
            "missing_positions": missing_positions,
            "extra_positions": extra_positions,
        })
        result["warnings"].append("Positioned node set differs")
    elif actual_positions:
        changed = []
        for node_id in sorted(actual_positions):
            if actual["positions"][node_id] != expected["positions"][node_id]:
                changed.append(node_id)
        if changed:
            result["layout_warnings"].append({
                "status": "WARN",
                "message": "Node positions differ",
                "changed_nodes": changed[:50],
                "changed_node_count": len(changed),
            })
            result["warnings"].append(f"Node positions differ for {len(changed)} nodes")

    if result["failures"]:
        result["status"] = "FAIL"
    elif result["warnings"]:
        result["status"] = "WARN"
    return result


def _find_cx2_files(root: Path) -> dict[Path, Path]:
    return {path.relative_to(root): path for path in sorted(root.rglob("*.cx2"))}


def _match_directory_files(actual_root: Path, expected_root: Path) -> dict[str, Any]:
    actual_by_rel = _find_cx2_files(actual_root)
    expected_by_rel = _find_cx2_files(expected_root)

    matched: list[tuple[str, Path, Path]] = []
    unmatched_actual = dict(actual_by_rel)
    unmatched_expected = dict(expected_by_rel)

    for rel_path in sorted(set(actual_by_rel) & set(expected_by_rel)):
        matched.append((str(rel_path), actual_by_rel[rel_path], expected_by_rel[rel_path]))
        unmatched_actual.pop(rel_path, None)
        unmatched_expected.pop(rel_path, None)

    actual_by_name: dict[str, list[tuple[Path, Path]]] = defaultdict(list)
    for rel_path, path in unmatched_actual.items():
        actual_by_name[path.name].append((rel_path, path))
    expected_by_name: dict[str, list[tuple[Path, Path]]] = defaultdict(list)
    for rel_path, path in unmatched_expected.items():
        expected_by_name[path.name].append((rel_path, path))

    for filename in sorted(set(actual_by_name) & set(expected_by_name)):
        actual_candidates = actual_by_name[filename]
        expected_candidates = expected_by_name[filename]
        while actual_candidates and expected_candidates:
            actual_rel, actual_path = actual_candidates.pop(0)
            expected_rel, expected_path = expected_candidates.pop(0)
            label = f"{actual_rel} ↔ {expected_rel}"
            matched.append((label, actual_path, expected_path))
            unmatched_actual.pop(actual_rel, None)
            unmatched_expected.pop(expected_rel, None)

    return {
        "matched": matched,
        "missing_expected_files": [str(path) for path in sorted(unmatched_expected)],
        "extra_actual_files": [str(path) for path in sorted(unmatched_actual)],
    }


def inspect_file(path: Path, overlap_threshold: float) -> dict[str, Any]:
    cx2_data = load_cx2(path)
    canonical = extract_canonical_network(cx2_data)
    metrics = compute_layout_metrics(canonical, overlap_threshold)
    return {
        "file": str(path),
        "network_name": canonical.get("network_name"),
        "metrics": metrics,
        "canonical": canonical,
    }


def compare_directories(
    actual_root: Path,
    expected_root: Path,
    overlap_threshold: float,
) -> dict[str, Any]:
    matches = _match_directory_files(actual_root, expected_root)
    report: dict[str, Any] = {
        "actual_root": str(actual_root),
        "expected_root": str(expected_root),
        "status": "PASS",
        "missing_expected_files": matches["missing_expected_files"],
        "extra_actual_files": matches["extra_actual_files"],
        "files": [],
    }

    if report["missing_expected_files"]:
        report["status"] = "FAIL"
    elif report["extra_actual_files"]:
        report["status"] = "WARN"

    for label, actual_path, expected_path in matches["matched"]:
        actual = extract_canonical_network(load_cx2(actual_path))
        expected = extract_canonical_network(load_cx2(expected_path))
        comparison = compare_networks(actual, expected)
        actual_metrics = compute_layout_metrics(actual, overlap_threshold)
        expected_metrics = compute_layout_metrics(expected, overlap_threshold)
        file_result = {
            "label": label,
            "actual_file": str(actual_path),
            "expected_file": str(expected_path),
            "status": comparison["status"],
            "comparison": comparison,
            "actual_metrics": actual_metrics,
            "expected_metrics": expected_metrics,
        }
        report["files"].append(file_result)

        if comparison["status"] == "FAIL":
            report["status"] = "FAIL"
        elif comparison["status"] == "WARN" and report["status"] == "PASS":
            report["status"] = "WARN"

    return report


def _markdown_metrics(metrics: dict[str, Any]) -> list[str]:
    bbox = metrics.get("bounding_box")
    lines = [
        f"- Nodes: {metrics['node_count']}",
        f"- Edges: {metrics['edge_count']}",
        f"- Positioned nodes: {metrics['positioned_node_count']}",
        f"- Connected components: {metrics['connected_components']}",
    ]
    if bbox:
        lines.append(f"- Bounding box: {bbox['width']:.3f} x {bbox['height']:.3f}")
    else:
        lines.append("- Bounding box: n/a")
    lines.extend([
        f"- Minimum node distance: {metrics['minimum_node_distance']}",
        (
            "- Possible overlaps "
            f"(< {metrics['overlap_threshold']} px): {metrics['possible_node_overlap_count']}"
        ),
    ])
    return lines


def write_markdown_report(report: dict[str, Any], path: str | Path) -> None:
    """Write a human-readable Markdown report for inspect or compare output."""
    lines: list[str] = []
    mode = report.get("mode")
    if mode == "inspect":
        inspected = report["inspection"]
        lines.append("# CX2 Inspect Report")
        lines.append("")
        lines.append(f"- File: `{inspected['file']}`")
        lines.append(f"- Network: `{inspected.get('network_name') or 'unknown'}`")
        lines.append("")
        lines.append("## Layout Metrics")
        lines.extend(_markdown_metrics(inspected["metrics"]))
        lines.append("")
        lines.append("## Sample Nodes")
        nodes = inspected["canonical"]["nodes"]
        for node_id, node in list(nodes.items())[:20]:
            lines.append(
                f"- `{node_id}` type={node.get('molecule_type')} "
                f"uniprot={node.get('uniprot_id')} color_group={node.get('color_group')}"
            )
        lines.append("")
        lines.append("## Sample Edges")
        edges = inspected["canonical"]["edges"]
        for edge_key, edge in list(edges.items())[:20]:
            lines.append(f"- `{edge_key}` values={edge.get('values')}")
        lines.append("")
        lines.append("## Style Semantics")
        style_semantics = inspected["canonical"].get("style_semantics", {})
        node_color_mapping = (
            style_semantics
            .get("mappings", {})
            .get("node", {})
            .get("NODE_BACKGROUND_COLOR")
            or style_semantics
            .get("mappings", {})
            .get("node", {})
            .get("NODE_FILL_COLOR")
        )
        edge_style_mappings = style_semantics.get("mappings", {}).get("edge", {})
        node_defaults = style_semantics.get("defaults", {}).get("node", {})
        edge_defaults = style_semantics.get("defaults", {}).get("edge", {})
        if node_color_mapping:
            lines.append("- Node color mapping:")
            for key, value in node_color_mapping.get("values", {}).items():
                lines.append(f"  - `{key}` -> `{value}`")
        elif node_defaults:
            lines.append(f"- Node defaults: `{node_defaults}`")
        else:
            lines.append("- Node color mapping: n/a")
        if edge_style_mappings:
            lines.append("- Edge style mappings:")
            for prop, mapping in sorted(edge_style_mappings.items()):
                lines.append(f"  - `{prop}` via `{mapping.get('attribute')}`: `{mapping.get('values')}`")
        elif edge_defaults:
            lines.append(f"- Edge defaults: `{edge_defaults}`")
        else:
            lines.append("- Edge style mappings: n/a")
    else:
        lines.append("# CX2 Compare Report")
        lines.append("")
        lines.append(f"- Status: **{report['status']}**")
        lines.append(f"- Actual root: `{report['actual_root']}`")
        lines.append(f"- Expected root: `{report['expected_root']}`")
        lines.append(f"- Compared files: {len(report['files'])}")
        lines.append("")

        if report["missing_expected_files"]:
            lines.append("## Missing Expected Files")
            for rel_path in report["missing_expected_files"]:
                lines.append(f"- `{rel_path}`")
            lines.append("")
        if report["extra_actual_files"]:
            lines.append("## Extra Actual Files")
            for rel_path in report["extra_actual_files"]:
                lines.append(f"- `{rel_path}`")
            lines.append("")

        lines.append("## File Results")
        for file_result in report["files"]:
            lines.append(f"### {file_result['label']}")
            lines.append(f"- Status: **{file_result['status']}**")
            lines.append("- Actual metrics:")
            lines.extend(f"  {line}" for line in _markdown_metrics(file_result["actual_metrics"]))
            lines.append("- Expected metrics:")
            lines.extend(f"  {line}" for line in _markdown_metrics(file_result["expected_metrics"]))

            comparison = file_result["comparison"]
            for title, key in [
                ("Failures", "failures"),
                ("Warnings", "warnings"),
            ]:
                if comparison[key]:
                    lines.append(f"- {title}:")
                    for item in comparison[key][:25]:
                        lines.append(f"  - {item}")
            lines.append("")

    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_json_report(report: dict[str, Any], path: str | Path) -> None:
    """Write a JSON report."""
    with Path(path).open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True, ensure_ascii=False)
        handle.write("\n")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--overlap-threshold",
        type=float,
        default=20.0,
        help="Distance in pixels below which positioned nodes are reported as possible overlaps.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    inspect_parser = subparsers.add_parser("inspect", help="Inspect one CX2 file.")
    inspect_parser.add_argument("--file", required=True, type=Path)
    inspect_parser.add_argument("--out", required=True, type=Path)
    inspect_parser.add_argument("--json-out", type=Path)

    compare_parser = subparsers.add_parser("compare", help="Compare actual and expected CX2 directories.")
    compare_parser.add_argument("--actual", required=True, type=Path)
    compare_parser.add_argument("--expected", required=True, type=Path)
    compare_parser.add_argument("--out", required=True, type=Path)
    compare_parser.add_argument("--json-out", type=Path)

    return parser


def main() -> int:
    parser = _build_parser()
    args = parser.parse_args()

    if args.command == "inspect":
        report = {
            "mode": "inspect",
            "inspection": inspect_file(args.file, args.overlap_threshold),
        }
        write_markdown_report(report, args.out)
        if args.json_out:
            write_json_report(report, args.json_out)
        return 0

    report = compare_directories(args.actual, args.expected, args.overlap_threshold)
    report["mode"] = "compare"
    write_markdown_report(report, args.out)
    if args.json_out:
        write_json_report(report, args.json_out)
    return 1 if report["status"] == "FAIL" else 0


if __name__ == "__main__":
    raise SystemExit(main())
