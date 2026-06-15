"""Shared helpers for linked chain components across PDB files."""

from __future__ import annotations

import hashlib
from typing import Any, Dict, Iterable, List, Mapping, Set, Tuple


def sanitize_filename_part(text: object) -> str:
    """Return a filesystem-safe identifier fragment."""
    allowed = []
    for char in str(text):
        if char.isalnum() or char in {"-", "_"}:
            allowed.append(char)
        else:
            allowed.append("_")
    sanitized = "".join(allowed).strip("_")
    return sanitized or "Unknown"


def make_component_title(prefix: str, identifiers: Iterable[object], *, max_preview: int = 5) -> str:
    """Build a deterministic component title from stable identifiers."""
    sanitized_ids = [sanitize_filename_part(value) for value in sorted(str(value) for value in identifiers if value)]
    if not sanitized_ids:
        return f"{prefix}_Unknown"

    preview = "_".join(sanitized_ids[:max_preview])
    if len(sanitized_ids) > max_preview:
        digest_source = "|".join(sanitized_ids).encode("utf-8")
        digest = hashlib.md5(digest_source).hexdigest()[:8]
        return f"{prefix}_{preview}__{digest}"
    return f"{prefix}_{preview}"


def build_identity_edges(
    uniprot_to_chains: Mapping[str, Iterable[str]],
    chain_to_pdb: Mapping[str, str],
) -> List[Dict[str, Any]]:
    """Create cross-PDB identity bridge edges between chains with the same UniProt ID."""
    identity_edges: List[Dict[str, Any]] = []
    for chains in uniprot_to_chains.values():
        sorted_chains = sorted(chain for chain in chains if chain in chain_to_pdb)
        if len(sorted_chains) < 2:
            continue

        for index in range(len(sorted_chains) - 1):
            chain_a = sorted_chains[index]
            chain_b = sorted_chains[index + 1]
            if chain_to_pdb.get(chain_a) == chain_to_pdb.get(chain_b):
                continue

            identity_edges.append({
                "chain_a": chain_a,
                "chain_b": chain_b,
                "all_atoms_count": 1000,
                "interaction_type": "identity",
            })
    return identity_edges


def find_linked_components(
    interaction_edges: Iterable[Mapping[str, Any]],
    identity_edges: Iterable[Mapping[str, Any]],
    valid_nodes: Iterable[str] | None = None,
) -> List[Set[str]]:
    """Return components connected by interactions and cross-PDB identity edges."""
    valid = set(valid_nodes) if valid_nodes is not None else None
    adjacency: Dict[str, Set[str]] = {}
    identity_nodes: Set[str] = set()

    def add_edge(left: object, right: object, *, identity: bool = False) -> None:
        left_id = str(left)
        right_id = str(right)
        if valid is not None and (left_id not in valid or right_id not in valid):
            return
        adjacency.setdefault(left_id, set()).add(right_id)
        adjacency.setdefault(right_id, set()).add(left_id)
        if identity:
            identity_nodes.add(left_id)
            identity_nodes.add(right_id)

    for edge in interaction_edges:
        add_edge(edge.get("chain_a", ""), edge.get("chain_b", ""))
    for edge in identity_edges:
        add_edge(edge.get("chain_a", ""), edge.get("chain_b", ""), identity=True)

    visited: Set[str] = set()
    components: List[Set[str]] = []
    for start_node in sorted(identity_nodes):
        if start_node in visited:
            continue

        stack = [start_node]
        component: Set[str] = set()
        while stack:
            node = stack.pop()
            if node in visited:
                continue

            visited.add(node)
            component.add(node)
            for neighbor in adjacency.get(node, set()):
                if neighbor not in visited:
                    stack.append(neighbor)

        if component:
            components.append(component)

    return components


def edge_tuple(edge: Mapping[str, Any]) -> Tuple[str, str]:
    """Return the chain-pair tuple for an edge-like mapping."""
    return str(edge.get("chain_a", "")), str(edge.get("chain_b", ""))
