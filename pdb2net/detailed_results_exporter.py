"""Detailed interaction export utilities.

This module writes a CSV with atom-level interaction details between chains:
- per-atom pairs within a fixed cutoff (currently 5.0 Å),
- residue identifiers and atom names on both sides,
- Euclidean distance of the interacting atoms,
- optional UniProt accessions carried along from the chain dictionaries.

Notes
-----
- Input `structure_data` is expected to contain a parsed `atom_data` list as produced
  by earlier stages (see data_processor.process_structure).
- Atom coordinates are prepacked into NumPy arrays and queried with a KD-tree for speed.
- No filtering beyond the provided interactions list is performed here.
"""

from __future__ import annotations

import os
from typing import Any, Dict, List

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


def export_detailed_interactions(
    structure_data: Dict[str, Any],
    interactions: List[Dict[str, Any]],
    run_output_path: str,
) -> None:
    pdb_id = structure_data["pdb_id"]
    atom_data = structure_data["atom_data"]

    # NEW: use the same cutoff as the interaction detection
    from config_loader import config
    radius = float(config["distance_thresholds"]["all_atoms_radius"])

    residues_atoms_lookup: Dict[str, Dict[str, Any]] = {}
    for chain in atom_data:
        chain_id = chain["chain_id"]
        atoms: List[Dict[str, str]] = []
        atom_coords: List[List[float]] = []

        for residue in chain.get("residues", []):
            res_id = residue.get("residue_name", "?")
            res_number = residue.get("residue_number", "?")
            for atom in residue.get("atoms", []):
                atoms.append(
                    {"residue": f"{res_number}:{res_id}", "atom_name": atom.get("atom_name", "?")}
                )
                atom_coords.append(atom.get("coordinates"))

        coords_arr = np.array(atom_coords, dtype=float) if atom_coords else np.empty((0, 3), dtype=float)
        residues_atoms_lookup[chain_id] = {"atoms": atoms, "coords": coords_arr}

    detailed_interactions: List[Dict[str, Any]] = []

    for interaction in interactions:
        chain_a_raw = interaction.get("chain_a", "")
        chain_b_raw = interaction.get("chain_b", "")
        chain_a_id = chain_a_raw.split(":")[1] if ":" in chain_a_raw else chain_a_raw
        chain_b_id = chain_b_raw.split(":")[1] if ":" in chain_b_raw else chain_b_raw

        data_a = residues_atoms_lookup.get(chain_a_id, {"atoms": [], "coords": np.empty((0, 3))})
        data_b = residues_atoms_lookup.get(chain_b_id, {"atoms": [], "coords": np.empty((0, 3))})

        if data_a["coords"].size == 0 or data_b["coords"].size == 0:
            continue

        chain_a_data = next((ch for ch in atom_data if ch.get("chain_id") == chain_a_id), None)
        chain_b_data = next((ch for ch in atom_data if ch.get("chain_id") == chain_b_id), None)
        if not chain_a_data or not chain_b_data:
            continue

        uniprot_a = chain_a_data.get("uniprot_id", "UNKNOWN")
        uniprot_b = chain_b_data.get("uniprot_id", "UNKNOWN")
        interaction_type = interaction.get("interaction_type", "AA")

        tree_a = cKDTree(data_a["coords"])
        # CHANGED: r=radius (from config) instead of r=5.0
        pairs = tree_a.query_ball_point(data_b["coords"], r=radius)

        for idx_b, idx_list_a in enumerate(pairs):
            for idx_a in idx_list_a:
                atom_a = data_a["atoms"][idx_a]
                atom_b = data_b["atoms"][idx_b]
                distance = float(np.linalg.norm(data_a["coords"][idx_a] - data_b["coords"][idx_b]))

                detailed_interactions.append(
                    {
                        "PDB_ID": pdb_id,
                        "Chain_A": chain_a_id,
                        "Residue_A": atom_a["residue"],
                        "Atom_A": atom_a["atom_name"],
                        "Chain_B": chain_b_id,
                        "Residue_B": atom_b["residue"],
                        "Atom_B": atom_b["atom_name"],
                        "Distance": round(distance, 2),
                        "UniProt_A": uniprot_a,
                        "UniProt_B": uniprot_b,
                        "Interaction_Type": interaction_type,
                    }
                )

    df = pd.DataFrame(detailed_interactions)
    os.makedirs(run_output_path, exist_ok=True)
    output_file = os.path.join(run_output_path, f"{pdb_id}_detailed_interactions.csv")
    df.to_csv(output_file, index=False)
