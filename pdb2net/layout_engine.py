"""Layout backends for PDB2Net network exports.

This module is intentionally independent from Cytoscape export. It receives the
prepared node/edge data frames and returns final node coordinates that can be
written into CX2 by the export layer.
"""

from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
import tempfile
import zlib
from pathlib import Path
from typing import Any, Dict


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


def _prefuse_layout_job(nodes_df, edges_df, network_title: str) -> dict[str, Any]:
    edges = []
    if edges_df is not None and len(edges_df) > 0:
        for _, edge in edges_df.iterrows():
            edges.append({
                "source": str(edge["source"]),
                "target": str(edge["target"]),
            })

    return {
        "network_title": str(network_title),
        "nodes": [{"id": node_id} for node_id in _node_ids(nodes_df)],
        "edges": edges,
        "layout": {
            "algorithm": "prefuse_force_directed",
            "iterations": 100,
            "spring_coefficient": 0.0001,
            "spring_length": 50,
            "node_mass": 3,
            "use_edge_weights": False,
            "deterministic": True,
        },
    }


def _engine_command(engine: Path) -> list[str]:
    if engine.suffix.lower() == ".jar":
        java = shutil.which("java")
        if not java:
            raise RuntimeError("java was not found in PATH")
        return [java, "-jar", str(engine)]
    return [str(engine)]


def calculate_prefuse_headless_layout(
    nodes_df,
    edges_df,
    network_title: str,
    config: dict[str, Any] | None = None,
) -> PositionMap:
    """Placeholder hook for a future external Java/Prefuse layout engine.

    If `layout_engine_path` points to an executable, the layout job is sent as
    JSON input/output files are exchanged through a temporary directory. If the
    engine is not configured or fails, this falls back to the Python layout.
    """
    config = config or {}
    node_ids = _node_ids(nodes_df)
    engine_path = str(config.get("layout_engine_path") or "").strip()
    if not engine_path:
        print("[layout] prefuse_headless requested but layout_engine_path is not set; using python_fast.")
        return calculate_python_fast_layout(nodes_df, edges_df, network_title)

    engine = Path(os.path.expanduser(os.path.expandvars(engine_path)))
    if not engine.exists():
        print(f"[layout] prefuse_headless engine not found at {engine}; using python_fast.")
        return calculate_python_fast_layout(nodes_df, edges_df, network_title)

    job = _prefuse_layout_job(nodes_df, edges_df, network_title)
    temp_dir: tempfile.TemporaryDirectory[str] | None = None
    try:
        command = _engine_command(engine)
        if config.get("layout_keep_temp_files", False):
            work_dir = Path(tempfile.mkdtemp(prefix="pdb2net-layout-"))
        else:
            temp_dir = tempfile.TemporaryDirectory(prefix="pdb2net-layout-")
            work_dir = Path(temp_dir.name)

        input_path = work_dir / "layout_job.json"
        output_path = work_dir / "positions.json"
        with input_path.open("w", encoding="utf-8") as handle:
            json.dump(job, handle, ensure_ascii=False, indent=2)

        subprocess.run(
            command + ["--input", str(input_path), "--output", str(output_path)],
            capture_output=True,
            text=True,
            check=True,
        )

        with output_path.open("r", encoding="utf-8") as handle:
            response = json.load(handle)

        if config.get("layout_keep_temp_files", False):
            print(f"[layout] Kept prefuse_headless temp files in {work_dir}")

        positions = response.get("positions")
        if not isinstance(positions, dict):
            raise ValueError("layout engine response did not contain a positions object")
        return _validate_positions(node_ids, positions)
    except Exception as exc:
        print(f"[layout] prefuse_headless failed ({exc}); using python_fast.")
        return calculate_python_fast_layout(nodes_df, edges_df, network_title)
    finally:
        if temp_dir is not None:
            temp_dir.cleanup()


def calculate_positions(
    nodes_df,
    edges_df,
    network_title: str,
    layout_mode: str | None,
    config: dict[str, Any] | None = None,
) -> PositionMap:
    """Calculate final positions for all nodes using the selected layout backend."""
    node_ids = _node_ids(nodes_df)
    mode = (layout_mode or "python_fast").strip().lower()
    if mode == "python_fast":
        positions = calculate_python_fast_layout(nodes_df, edges_df, network_title)
    elif mode == "prefuse_headless":
        positions = calculate_prefuse_headless_layout(nodes_df, edges_df, network_title, config=config)
    elif mode == "cytoscape_live":
        print("[layout] cytoscape_live is handled by the live Cytoscape path; using python_fast for CX2 preset positions.")
        positions = calculate_python_fast_layout(nodes_df, edges_df, network_title)
    else:
        print(f"[layout] Unknown layout_mode {layout_mode!r}; using python_fast.")
        positions = calculate_python_fast_layout(nodes_df, edges_df, network_title)

    return _validate_positions(node_ids, positions)
