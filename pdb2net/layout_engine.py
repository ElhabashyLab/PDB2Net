"""Layout backends for PDB2Net network exports.

This module is intentionally independent from Cytoscape export. It receives the
prepared node/edge data frames and returns final node coordinates that can be
written into CX2 by the export layer.
"""

from __future__ import annotations

import math
import zlib
from typing import Dict


PositionMap = Dict[str, Dict[str, float]]


def _node_ids(nodes_df) -> list[str]:
    if nodes_df is None or len(nodes_df) == 0 or "id" not in nodes_df.columns:
        return []
    return [str(node_id) for node_id in nodes_df["id"]]


def _validate_positions(node_ids: list[str], positions: PositionMap | None) -> PositionMap:
    """Ensure every node has a numeric x/y position."""
    final: PositionMap = {}
    positions = positions or {}
    for node_id in node_ids:
        raw = positions.get(node_id, {})
        try:
            x = float(raw.get("x", 0.0))
            y = float(raw.get("y", 0.0))
        except (AttributeError, TypeError, ValueError):
            x, y = 0.0, 0.0
        final[node_id] = {"x": x, "y": y}
    return final


def calculate_python_fast_layout(
    nodes_df,
    edges_df,
    network_title: str,
    scale: float = 1000.0,
) -> PositionMap:
    """Calculate deterministic NetworkX-based coordinates.

    This is the previous headless layout logic, moved out of `cytoscape.py`.
    Combined networks are laid out per connected component and packed into a
    simple grid; standard networks retain the original spring-layout scaling.
    """
    import networkx as nx

    node_ids = _node_ids(nodes_df)
    node_count = len(node_ids)

    if node_count == 0:
        return {}
    if node_count == 1:
        return {node_ids[0]: {"x": 0.0, "y": 0.0}}

    graph = nx.Graph()
    graph.add_nodes_from(node_ids)
    if edges_df is not None and len(edges_df) > 0:
        for _, edge in edges_df.iterrows():
            source = str(edge["source"])
            target = str(edge["target"])
            if source == target:
                continue
            weight = float(edge.get("all_atoms_count", 1.0))
            weight = 1.0 + math.log10(max(1.0, weight))
            graph.add_edge(source, target, weight=weight)

    seed = zlib.adler32(str(network_title).encode("utf-8")) & 0xFFFFFFFF
    is_combined = "combined" in str(network_title).lower()

    if is_combined:
        components = list(nx.connected_components(graph))
        components.sort(key=len, reverse=True)
        grid_cols = math.ceil(math.sqrt(len(components))) if components else 1
        final_positions: PositionMap = {}

        for index, component in enumerate(components):
            sub_graph = graph.subgraph(component)
            sub_node_count = len(component)
            sub_scale = 100.0 + (math.sqrt(sub_node_count) * 60.0)
            sub_seed = (seed + index) & 0xFFFFFFFF
            sub_positions = nx.spring_layout(
                sub_graph,
                seed=sub_seed,
                weight="weight",
                scale=sub_scale,
                center=(0, 0),
                k=None,
            )

            row = index // grid_cols
            col = index % grid_cols
            stride = 600.0 + (math.sqrt(node_count / grid_cols) * 50.0)
            offset_x = col * stride
            offset_y = row * stride

            for node_id, coords in sub_positions.items():
                final_positions[str(node_id)] = {
                    "x": float(coords[0] + offset_x),
                    "y": float(coords[1] + offset_y),
                }

        return _validate_positions(node_ids, final_positions)

    if node_count == 2:
        distance = scale * 0.08
        return {
            node_ids[0]: {"x": -distance, "y": 0.0},
            node_ids[1]: {"x": distance, "y": 0.0},
        }

    iterations = max(100, min(250, 8 * node_count))
    scale_used = scale * 0.8
    if node_count < 40:
        scale_used *= max(0.12, (node_count / 40.0))

    k = 1.0 / math.sqrt(max(node_count, 100))
    if node_count < 40:
        k *= 0.75

    positions = nx.spring_layout(
        graph,
        seed=seed,
        weight="weight",
        iterations=iterations,
        dim=2,
        center=(0.0, 0.0),
        scale=scale_used,
        k=k,
    )
    return _validate_positions(
        node_ids,
        {str(node_id): {"x": float(x), "y": float(y)} for node_id, (x, y) in positions.items()},
    )


def calculate_positions(
    nodes_df,
    edges_df,
    network_title: str,
    layout_mode: str,
) -> PositionMap:
    """Calculate final positions for all nodes using the selected layout backend."""
    if not isinstance(layout_mode, str) or layout_mode.strip().lower() != "python_fast":
        raise ValueError(f"Unsupported layout_mode: {layout_mode!r}")
    return calculate_python_fast_layout(nodes_df, edges_df, network_title)
