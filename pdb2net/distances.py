"""Distance-based interaction detection using cKDTree.

This module computes chain–chain interactions by counting nearby atom pairs
within configurable distance thresholds. It supports:
- Protein–Protein interactions: first filter by Cα proximity, then confirm with
  all-atom proximity.
- Nucleic-acid interactions (DNA/RNA/DNA-RNA/Nucleic Acid): use all-atom
  proximity directly.
- Fine-grained interaction labels (e.g., Protein-RNA, DNA-RNA, etc.).

Caches
------
- `coords_cache`: memoizes extracted coordinate arrays per chain object/extraction
- `tree_cache`: memoizes cKDTree instances per chain object/extraction

Configuration
-------------
- Distance thresholds are read from `config["distance_thresholds"]`:
  * "ca_radius"         → Cα cutoff for protein–protein prefilter
  * "all_atoms_radius"  → all-atom cutoff for confirmation / NA interactions
- Minimum contact filters are read from `config["interaction_filters"]`.
"""

from __future__ import annotations

import itertools
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.spatial import cKDTree

from .config_loader import config
from .contact_atoms import ContactView, contact_positions, count_contacts
from .input_contract import InputValidationError

INTERACTION_PIPELINE_VERSION = "pdb2net-distance-classification-v3"
ALLOWED_INTERACTION_TYPES = frozenset(
    {
        "Protein-Protein",
        "Protein-DNA",
        "Protein-RNA",
        "Protein-DNA/RNA",
        "Protein-Nucleic Acid",
        "DNA-DNA",
        "DNA-RNA",
        "DNA-DNA/RNA",
        "RNA-RNA",
        "RNA-DNA/RNA",
        "DNA/RNA-DNA/RNA",
        "Nucleic Acid-Nucleic Acid",
    }
)
DISTANCE_THRESHOLD_RULES = {
    "ca_radius": {"type": "number", "default": 12.0, "minimum": 2.0, "maximum": 30.0},
    "all_atoms_radius": {"type": "number", "default": 5.0, "minimum": 1.0, "maximum": 15.0},
}
INTERACTION_FILTER_RULES = {
    "protein_protein_min_ca_neighbors": {"type": "integer", "default": 10, "minimum": 1, "maximum": 100_000},
    "protein_protein_min_all_atom_contacts": {"type": "integer", "default": 1, "minimum": 1, "maximum": 100_000},
    "protein_nucleic_acid_min_all_atom_contacts": {"type": "integer", "default": 1, "minimum": 1, "maximum": 100_000},
    "nucleic_acid_min_all_atom_contacts": {"type": "integer", "default": 1, "minimum": 1, "maximum": 100_000},
}

# Caches for KD-Trees and extracted coordinates. Include object identity so
# repeated chain IDs from different structures in one process never share trees.
CacheKey = Tuple[int, str, str]
tree_cache: Dict[CacheKey, cKDTree] = {}
coords_cache: Dict[CacheKey, np.ndarray] = {}


def _cache_key(chain: Dict[str, Any], extraction_type: str) -> CacheKey:
    return (id(chain), chain["unique_chain_id"], extraction_type)


def extract_coordinates(chain: Dict[str, Any], extraction_type: str) -> np.ndarray:
    """Extract and cache atom coordinates for a chain.

    Parameters
    ----------
    chain : dict
        Chain dictionary with a 'residues' list where each residue contains an
        'atoms' list of dicts: {"atom_name": str, "coordinates": [x, y, z]}.
        Must contain a unique identifier under 'unique_chain_id'.
    extraction_type : {"ca", "all_atoms"}
        - "ca":    only Cα atoms (atom_name == "CA")
        - "all_atoms": all heavy atoms (coordinates as provided upstream)

    Returns
    -------
    numpy.ndarray, shape (N, 3)
        Extracted coordinates (may be empty with shape (0, 3)).
    """
    key = _cache_key(chain, extraction_type)
    if key in coords_cache:
        return coords_cache[key]

    if extraction_type in {"ca", "all_atoms"}:
        coords = np.array([position.coordinates for position in contact_positions(chain, extraction_type)])
    else:
        coords = np.empty((0, 3))

    # Normalize empty arrays to a consistent shape
    coords_cache[key] = coords if coords.size else np.empty((0, 3))
    return coords_cache[key]


def get_or_create_tree(chain: Dict[str, Any], extraction_type: str) -> Optional[cKDTree]:
    """Create or retrieve a cached cKDTree for the given chain/extraction.

    Parameters
    ----------
    chain : dict
        Chain dictionary; see `extract_coordinates` for expected structure.
    extraction_type : {"ca", "all_atoms"}
        Coordinate set to build a KD-tree on.

    Returns
    -------
    scipy.spatial.cKDTree | None
        KD-tree instance if there are any coordinates, otherwise None.
    """
    key = _cache_key(chain, extraction_type)
    if key not in tree_cache:
        points = extract_coordinates(chain, extraction_type)
        if points.size:
            tree_cache[key] = cKDTree(points)
    return tree_cache.get(key)


def count_nearby_atoms(tree_a: Optional[cKDTree], tree_b: Optional[cKDTree], radius: float) -> int:
    """Count nearby atom pairs within a distance threshold.

    Parameters
    ----------
    tree_a : cKDTree | None
        KD-tree for chain A (query target).
    tree_b : cKDTree | None
        KD-tree for chain B (query points).
    radius : float
        Distance cutoff for neighborhood queries.

    Returns
    -------
    int
        Total number of atom–atom pairs within `radius`.
    """
    if not tree_a or not tree_b:
        return 0
    # count_neighbors returns the same pair count without materializing every
    # matching atom index, which can otherwise require quadratic Python memory.
    return int(tree_a.count_neighbors(tree_b, radius))


def determine_interaction_type(mol_type_a: Optional[str], mol_type_b: Optional[str]) -> Optional[str]:
    """Return a fine-grained interaction label for a pair of molecule types.

    Rules
    -----
    - If either type is 'Unknown' → return None (skip).
    - Protein–Protein → "Protein-Protein"
    - Protein–(DNA|RNA|DNA/RNA|Nucleic Acid) → "Protein-<type>" or "Protein-Nucleic Acid"
    - Nucleic-acid pairs:
        DNA-DNA, RNA-RNA, DNA-RNA, DNA-DNA/RNA, RNA-DNA/RNA, DNA/RNA-DNA/RNA
      If either side is the generic 'Nucleic Acid', fall back to:
        "Nucleic Acid-Nucleic Acid"

    Parameters
    ----------
    mol_type_a, mol_type_b : str | None
        Molecule types as assigned upstream.

    Returns
    -------
    str | None
        Canonical label or None if no meaningful interaction type applies.
    """
    a = (mol_type_a or "Unknown").strip()
    b = (mol_type_b or "Unknown").strip()

    if "Unknown" in (a, b):
        return None

    specific_na = {"DNA", "RNA", "DNA/RNA"}
    generic_na = "Nucleic Acid"

    # Protein—Protein
    if a == "Protein" and b == "Protein":
        return "Protein-Protein"

    # Protein—Nucleic Acid variants
    if "Protein" in (a, b):
        other = b if a == "Protein" else a
        if other in specific_na:
            return f"Protein-{other}"
        if other == generic_na:
            return "Protein-Nucleic Acid"
        return None  # unsupported combination

    # Nucleic Acid—Nucleic Acid variants
    if a in specific_na and b in specific_na:
        # Sort for canonical labels independent of order
        order = {"DNA": 0, "RNA": 1, "DNA/RNA": 2}
        pair = "-".join(sorted([a, b], key=lambda x: order.get(x, 3)))
        return pair

    # If either side is only the generic NA, keep it generic
    if generic_na in (a, b) and (a in specific_na | {generic_na}) and (b in specific_na | {generic_na}):
        return "Nucleic Acid-Nucleic Acid"

    return None


def _interaction_filters() -> Dict[str, int]:
    """Return configured minimum contact filters with backward-compatible defaults."""
    filters = config.get("interaction_filters", {})
    if not isinstance(filters, dict):
        filters = {}
    return {
        "protein_protein_min_ca_neighbors": int(filters.get("protein_protein_min_ca_neighbors", 10)),
        "protein_protein_min_all_atom_contacts": int(filters.get("protein_protein_min_all_atom_contacts", 1)),
        "protein_nucleic_acid_min_all_atom_contacts": int(
            filters.get("protein_nucleic_acid_min_all_atom_contacts", 1)
        ),
        "nucleic_acid_min_all_atom_contacts": int(filters.get("nucleic_acid_min_all_atom_contacts", 1)),
    }


def _distance_thresholds() -> Tuple[float, float]:
    thresholds = config["distance_thresholds"]
    return float(thresholds["ca_radius"]), float(thresholds["all_atoms_radius"])


def _required_all_atom_contacts(interaction_type: str, filters: Dict[str, int]) -> int:
    if interaction_type == "Protein-Protein":
        return filters["protein_protein_min_all_atom_contacts"]
    if interaction_type.startswith("Protein-"):
        return filters["protein_nucleic_acid_min_all_atom_contacts"]
    return filters["nucleic_acid_min_all_atom_contacts"]


# Backward-compatible module constants for callers that imported them directly.
RADIUS_CA: float = float(config["distance_thresholds"]["ca_radius"])
RADIUS_ALL_ATOMS: float = float(config["distance_thresholds"]["all_atoms_radius"])


def chain_pairs(chains):
    """Yield each pair of distinct author chains within one model once."""
    unique = {}
    for chain in chains:
        previous = unique.setdefault(chain["unique_chain_id"], chain)
        if (previous.get("model_index", 1), previous["chain_id"]) != (chain.get("model_index", 1), chain["chain_id"]):
            raise InputValidationError("DUPLICATE_CHAIN_IDENTITY", "One chain identity describes different chains or models.")
    for left, right in itertools.combinations(unique.values(), 2):
        if left.get("model_index", 1) != right.get("model_index", 1) or left["chain_id"] == right["chain_id"]:
            continue
        yield left, right


def select_interaction(structure, left, right, all_atoms_count, ca_neighbors, *, filters=None):
    """Apply the same annotation-dependent rules to measured or prepared counts."""
    filters = _interaction_filters() if filters is None else filters
    kind = determine_interaction_type(left.get("molecule_type"), right.get("molecule_type"))
    if kind is None or all_atoms_count < _required_all_atom_contacts(kind, filters):
        return None
    if kind == "Protein-Protein" and ca_neighbors < filters["protein_protein_min_ca_neighbors"]:
        return None
    result = {
        "file_path": structure["file_path"],
        "structure_key": str(structure.get("structure_identity", {}).get("key") or left.get("structure_key") or ""),
        "model_index": int(left.get("model_index", 1)),
        "chain_a": left["unique_chain_id"], "chain_b": right["unique_chain_id"],
        "all_atoms_count": int(all_atoms_count), "interaction_type": kind,
    }
    if kind == "Protein-Protein":
        result["ca_neighbors"] = int(ca_neighbors)
    return result


def calculate_distances_with_ckdtree(combined_data: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Measure contacts, retaining the raw path's type and C-alpha prefilters."""
    results = []
    radius_ca, radius_all = _distance_thresholds()
    filters = _interaction_filters()
    for structure in combined_data:
        chains = structure["atom_data"]
        views = {(c["unique_chain_id"], t): ContactView(c, t) for c in chains for t in ("ca", "all_atoms")}
        for left, right in chain_pairs(chains):
            kind = determine_interaction_type(left.get("molecule_type"), right.get("molecule_type"))
            if kind is None:
                continue
            a, b = left["unique_chain_id"], right["unique_chain_id"]
            ca_count = 0
            if kind == "Protein-Protein":
                ca_count = count_contacts(views[a, "ca"], views[b, "ca"], radius_ca)
                if ca_count < filters["protein_protein_min_ca_neighbors"]:
                    continue
            all_count = count_contacts(views[a, "all_atoms"], views[b, "all_atoms"], radius_all)
            edge = select_interaction(structure, left, right, all_count, ca_count, filters=filters)
            if edge is not None:
                results.append(edge)
    return results
