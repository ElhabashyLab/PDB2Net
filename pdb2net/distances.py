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
- `coords_cache`: memoizes extracted coordinate arrays per (chain, extraction_type)
- `tree_cache`: memoizes cKDTree instances per (chain, extraction_type)

Configuration
-------------
- Distance thresholds are read from `config["distance_thresholds"]`:
  * "ca_radius"         → Cα cutoff for protein–protein prefilter
  * "all_atoms_radius"  → all-atom cutoff for confirmation / NA interactions
- Minimum contact filters are read from `config["interaction_filters"]`.
"""

from __future__ import annotations

import itertools
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
from scipy.spatial import cKDTree

from .config_loader import config

# Caches for KD-Trees and extracted coordinates
tree_cache: Dict[Tuple[str, str], cKDTree] = {}
coords_cache: Dict[Tuple[str, str], np.ndarray] = {}


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
    key = (chain["unique_chain_id"], extraction_type)
    if key in coords_cache:
        return coords_cache[key]

    atoms_iter = (atom for residue in chain["residues"] for atom in residue["atoms"])
    if extraction_type == "ca":
        coords = np.array([atom["coordinates"] for atom in atoms_iter if atom["atom_name"] == "CA"])
    elif extraction_type == "all_atoms":
        coords = np.array([atom["coordinates"] for atom in atoms_iter])
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
    key = (chain["unique_chain_id"], extraction_type)
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
    # query_ball_point returns a list of neighbor indices per query point
    return sum(len(neigh) for neigh in tree_a.query_ball_point(tree_b.data, r=radius))


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


def calculate_distances_with_ckdtree(combined_data: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Compute inter-chain interactions using cKDTree.

    For each structure:
      1) Build or reuse KD-trees for "ca" and "all_atoms" per chain.
      2) Iterate over unique chain pairs (unordered) and infer interaction type.
      3) For nucleic acids, count all-atom neighbors at the configured radius.
      4) For proteins, require the configured Cα-neighbor minimum first,
         then confirm with the configured all-atom contact minimum.
      5) Emit a result for pairs that meet the criteria.

    Parameters
    ----------
    combined_data : list[dict]
        Items produced by earlier stages, each with:
          - "file_path" (str)
          - "atom_data" (list of chain dicts, each with "unique_chain_id",
            "chain_id", "molecule_type", and "residues" including atoms)

    Returns
    -------
    list[dict]
        Each result dict contains:
          - "file_path": str
          - "chain_a": str  (unique_chain_id)
          - "chain_b": str  (unique_chain_id)
          - "all_atoms_count": int
          - "interaction_type": str
    """
    results: List[Dict[str, Any]] = []

    nucleic_acid_types = {"Nucleic Acid", "DNA", "RNA", "DNA/RNA"}
    radius_ca, radius_all_atoms = _distance_thresholds()
    filters = _interaction_filters()

    for file_data in combined_data:
        file_path = file_data["file_path"]
        atom_data = file_data["atom_data"]

        # Prepare KD-trees for both extraction types per chain
        local_trees: Dict[Tuple[str, str], Optional[cKDTree]] = {
            (c["unique_chain_id"], t): get_or_create_tree(c, t)
            for c in atom_data for t in ("ca", "all_atoms")
        }

        seen_pairs: set[Tuple[str, str]] = set()  # deduplicate unordered chain pairs

        for chain_a, chain_b in itertools.combinations(atom_data, 2):
            # Skip explicit same-chain interactions (should not occur with combinations)
            if chain_a["chain_id"] == chain_b["chain_id"]:
                continue

            # Unordered pair key based on unique chain ids
            key = tuple(sorted([chain_a["unique_chain_id"], chain_b["unique_chain_id"]]))
            if key in seen_pairs:
                continue
            seen_pairs.add(key)

            interaction_type = determine_interaction_type(
                chain_a.get("molecule_type", "Unknown"),
                chain_b.get("molecule_type", "Unknown"),
            )
            if not interaction_type:
                continue

            # Any Nucleic Acid involvement → use all-atom proximity directly
            if (chain_a.get("molecule_type") in nucleic_acid_types or
                chain_b.get("molecule_type") in nucleic_acid_types):

                all_atoms_count = count_nearby_atoms(
                    local_trees.get((chain_a["unique_chain_id"], "all_atoms")),
                    local_trees.get((chain_b["unique_chain_id"], "all_atoms")),
                    radius=radius_all_atoms,
                )
                if all_atoms_count >= _required_all_atom_contacts(interaction_type, filters):
                    results.append({
                        "file_path": file_path,
                        "chain_a": chain_a["unique_chain_id"],
                        "chain_b": chain_b["unique_chain_id"],
                        "all_atoms_count": int(all_atoms_count),
                        "interaction_type": interaction_type,
                    })
                continue

            # Protein–Protein: require a minimum number of nearby Cα atoms first
            ca_neighbors = count_nearby_atoms(
                local_trees.get((chain_a["unique_chain_id"], "ca")),
                local_trees.get((chain_b["unique_chain_id"], "ca")),
                radius=radius_ca,
            )
            if ca_neighbors >= filters["protein_protein_min_ca_neighbors"]:
                all_atoms_count = count_nearby_atoms(
                    local_trees.get((chain_a["unique_chain_id"], "all_atoms")),
                    local_trees.get((chain_b["unique_chain_id"], "all_atoms")),
                    radius=radius_all_atoms,
                )
                if all_atoms_count >= filters["protein_protein_min_all_atom_contacts"]:
                    results.append({
                        "file_path": file_path,
                        "chain_a": chain_a["unique_chain_id"],
                        "chain_b": chain_b["unique_chain_id"],
                        "all_atoms_count": int(all_atoms_count),
                        "interaction_type": interaction_type,
                    })

    return results
