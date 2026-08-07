"""Helpers for optional size limits on combined network exports."""

from __future__ import annotations

from typing import Any, Mapping


LIMIT_KEYS = ("max_nodes", "max_edges")


def normalize_combined_graph_limits(raw: Mapping[str, Any] | None) -> dict[str, int | None]:
    """Return validated per-network combined graph limits.

    ``None`` means unlimited. Positive integer values are hard caps. The core
    defaults remain unlimited so existing local runs keep their prior behavior.
    """
    values = raw if isinstance(raw, Mapping) else {}
    normalized: dict[str, int | None] = {}
    for key in LIMIT_KEYS:
        value = values.get(key)
        if value in (None, ""):
            normalized[key] = None
            continue
        if isinstance(value, bool):
            raise ValueError(f"{key} must be a positive integer or null")
        try:
            parsed = int(value)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"{key} must be a positive integer or null") from exc
        if parsed <= 0:
            raise ValueError(f"{key} must be a positive integer or null")
        normalized[key] = parsed
    return normalized


def combined_graph_skip(
    *,
    network_kind: str,
    name: str,
    node_count: int,
    edge_count: int,
    limits: Mapping[str, int | None],
) -> dict[str, Any] | None:
    """Describe a skipped combined graph, or return ``None`` when within caps."""
    max_nodes = limits.get("max_nodes")
    max_edges = limits.get("max_edges")
    exceeded = []
    if max_nodes is not None and node_count > max_nodes:
        exceeded.append("max_nodes")
    if max_edges is not None and edge_count > max_edges:
        exceeded.append("max_edges")
    if not exceeded:
        return None

    return {
        "output_type": "network",
        "network_kind": network_kind,
        "name": name,
        "reason": "combined_graph_limit_exceeded",
        "exceeded_limits": exceeded,
        "counts": {"nodes": int(node_count), "edges": int(edge_count)},
        "limits": {"max_nodes": max_nodes, "max_edges": max_edges},
    }
