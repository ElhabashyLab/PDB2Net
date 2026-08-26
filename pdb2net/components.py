"""Shared helpers for linked chain components across PDB files."""

from __future__ import annotations

import hashlib
from typing import Any, Dict, Iterable, List, Mapping, Set, Tuple

COMBINED_COMPONENT_SEMANTICS = "pdb2net-cross-pdb-uniprot-linked-components-v1"


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
    """Create a deterministic cross-PDB spanning tree for each UniProt group.

    Every chain participates when the accession occurs in at least two PDBs,
    while same-PDB identity edges remain excluded.  The spanning tree avoids a
    quadratic complete multipartite graph without fragmenting components.
    """
    identity_edges: List[Dict[str, Any]] = []
    for chains in uniprot_to_chains.values():
        by_pdb: Dict[str, List[str]] = {}
        for chain in sorted(set(chains)):
            pdb_id = chain_to_pdb.get(chain)
            if pdb_id:
                by_pdb.setdefault(pdb_id, []).append(chain)
        ordered_groups = [by_pdb[pdb_id] for pdb_id in sorted(by_pdb)]
        if len(ordered_groups) < 2:
            continue

        first_anchor = ordered_groups[0][0]
        second_anchor = ordered_groups[1][0]
        pairs = [
            (first_anchor, chain)
            for group in ordered_groups[1:]
            for chain in group
        ]
        pairs.extend(
            (second_anchor, chain)
            for chain in ordered_groups[0][1:]
        )
        for chain_a, chain_b in pairs:
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
