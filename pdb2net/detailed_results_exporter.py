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
    """Export atom-level interaction details to CSV.

    Parameters
    ----------
    structure_data : dict
        Dictionary containing structure context:
          - "pdb_id" (str): PDB identifier.
          - "atom_data" (list[dict]): Chain dictionaries with residue/atom info.
            Each chain dict should have:
              - "chain_id" (str)
              - "residues" (list[dict] with:
                    "residue_name", "residue_number", "atoms" (list with
                    {"atom_name": str, "coordinates": [x, y, z]})
                )
              - optional "uniprot_id" (str)
    interactions : list[dict]
        Detected chain-level interactions. Each entry should provide:
          - "chain_a" (str, typically "PDBID:Chain")
          - "chain_b" (str, typically "PDBID:Chain")
          - "all_atoms_count" (int)  # may be present from previous stage
          - "interaction_type" (str)  # carried through into the CSV
    run_output_path : str
        Base directory where the CSV file is written. Output will be placed under:
        {run_output_path}/{pdb_id}/{pdb_id}_detailed_interactions.csv

    Returns
    -------
    None
        Writes a CSV file to disk; does not return a value.

    Notes
    -----
    - Atom-level pairs are gathered using a 5.0 Å radius search via cKDTree.
    - Chains without coordinates are skipped silently.
    """
    pdb_id = structure_data["pdb_id"]
    atom_data = structure_data["atom_data"]

    # Precompute residue/atom lookup per chain for efficient spatial queries.
    residues_atoms_lookup: Dict[str, Dict[str, Any]] = {}
    for chain in atom_data:
        chain_id = chain["chain_id"]
        atoms: List[Dict[str, str]] = []
        atom_coords: List[List[float]] = []

        for residue in chain["residues"]:
            res_id = residue["residue_name"]
            res_number = residue.get("residue_number", "?")  # Fallback if the number is missing
            for atom in residue["atoms"]:
                atoms.append(
                    {
                        "residue": f"{res_number}:{res_id}",
                        "atom_name": atom["atom_name"],
                    }
                )
                atom_coords.append(atom["coordinates"])

        residues_atoms_lookup[chain_id] = {
            "atoms": atoms,
            "coords": np.array(atom_coords) if atom_coords else np.empty((0, 3)),
        }

    detailed_interactions: List[Dict[str, Any]] = []

    # For each chain-level interaction, enumerate atom-atom pairs within 5.0 Å.
    for interaction in interactions:
        chain_a_id = interaction["chain_a"].split(":")[1]
        chain_b_id = interaction["chain_b"].split(":")[1]

        data_a = residues_atoms_lookup.get(chain_a_id, {"atoms": [], "coords": np.empty((0, 3))})
        data_b = residues_atoms_lookup.get(chain_b_id, {"atoms": [], "coords": np.empty((0, 3))})

        # Skip if coordinates are missing in either chain.
        if data_a["coords"].size == 0 or data_b["coords"].size == 0:
            continue

        chain_a_data = next((ch for ch in atom_data if ch["chain_id"] == chain_a_id), None)
        chain_b_data = next((ch for ch in atom_data if ch["chain_id"] == chain_b_id), None)
        if not chain_a_data or not chain_b_data:
            continue

        uniprot_a = chain_a_data.get("uniprot_id", "UNKNOWN")
        uniprot_b = chain_b_data.get("uniprot_id", "UNKNOWN")
        interaction_type = interaction["interaction_type"]

        # Build KD-tree for chain A; query all atoms of chain B within 5.0 Å.
        tree_a = cKDTree(data_a["coords"])
        pairs = tree_a.query_ball_point(data_b["coords"], r=5.0)

        for idx_b, idx_list_a in enumerate(pairs):
            for idx_a in idx_list_a:
                atom_a = data_a["atoms"][idx_a]
                atom_b = data_b["atoms"][idx_b]

                distance = np.linalg.norm(data_a["coords"][idx_a] - data_b["coords"][idx_b])

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

    # Save interaction details as a CSV file (always creates the directory).
    df = pd.DataFrame(detailed_interactions)

    output_dir = os.path.join(run_output_path, pdb_id)
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{pdb_id}_detailed_interactions.csv")
    df.to_csv(output_file, index=False)
