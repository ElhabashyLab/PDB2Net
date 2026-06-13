"""Smoke-test the prefuse_headless layout backend.

This is intentionally a tiny diagnostic script, not a pytest suite. It verifies
that the Python layout interface returns complete numeric positions, that
repeated calls are deterministic, and that missing engines fall back safely.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any

import pandas as pd

from pdb2net.layout_engine import calculate_positions


def _small_network() -> tuple[pd.DataFrame, pd.DataFrame]:
    nodes_df = pd.DataFrame({"id": ["A", "B", "C", "D", "E"]})
    edges_df = pd.DataFrame([
        {"source": "A", "target": "B", "all_atoms_count": 10},
        {"source": "B", "target": "C", "all_atoms_count": 15},
        {"source": "D", "target": "E", "all_atoms_count": 5},
    ])
    return nodes_df, edges_df


def _assert_positions(label: str, positions: dict[str, dict[str, Any]], expected_nodes: set[str]) -> None:
    if set(positions) != expected_nodes:
        raise AssertionError(f"{label}: expected nodes {sorted(expected_nodes)}, got {sorted(positions)}")
    for node_id, pos in positions.items():
        if not isinstance(pos, dict):
            raise AssertionError(f"{label}: position for {node_id} is not an object")
        for axis in ("x", "y"):
            value = pos.get(axis)
            if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
                raise AssertionError(f"{label}: invalid {axis} for {node_id}: {value!r}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine-path", type=Path, help="Optional layout engine jar or executable.")
    args = parser.parse_args()

    nodes_df, edges_df = _small_network()
    expected_nodes = set(nodes_df["id"].astype(str))

    fallback_one = calculate_positions(
        nodes_df,
        edges_df,
        "Smoke_Network",
        "prefuse_headless",
        config={"layout_engine_path": ""},
    )
    fallback_two = calculate_positions(
        nodes_df,
        edges_df,
        "Smoke_Network",
        "prefuse_headless",
        config={"layout_engine_path": "/missing/layout-engine.jar"},
    )
    _assert_positions("fallback without engine", fallback_one, expected_nodes)
    _assert_positions("fallback missing engine", fallback_two, expected_nodes)

    print("fallback_without_engine: ok")
    print("fallback_missing_engine: ok")

    if args.engine_path:
        engine_config = {
            "layout_engine_path": str(args.engine_path),
            "layout_keep_temp_files": False,
        }
        engine_one = calculate_positions(nodes_df, edges_df, "Smoke_Network", "prefuse_headless", config=engine_config)
        engine_two = calculate_positions(nodes_df, edges_df, "Smoke_Network", "prefuse_headless", config=engine_config)
        _assert_positions("engine run", engine_one, expected_nodes)
        _assert_positions("engine repeat", engine_two, expected_nodes)
        if engine_one != engine_two:
            raise AssertionError("engine output is not deterministic for repeated calls")
        print("engine_positions: ok")
        print("engine_deterministic: ok")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
