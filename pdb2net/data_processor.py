"""Lightweight structure processing using Gemmi.

This module extracts per-chain residue and atom information from a Gemmi
Structure object produced earlier in the pipeline. It keeps only non-hydrogen
atoms (excludes H and D) and flags chains as valid if they contain at least one
residue from a curated allow-list covering standard protein amino acids and
canonical nucleic-acid residues.

Notes
-----
- No molecule typing or naming is inferred here; fields like `molecule_name`,
  `molecule_type`, and `sequence` are intentionally left as placeholders for
  later stages (e.g., UniProt matching).
- The function preserves the original data flow and returns a dict with
  `file_path`, `pdb_id`, and `atom_data` (list of chain dicts).
"""

from __future__ import annotations

import gemmi
from typing import Any, Dict, List

from .config_loader import config
from .residue_types import AMINO_ACIDS, DNA_RESIDUES, RNA_RESIDUES, normalize_residue_name


# Set of allowed residue names for proteins and nucleic acids.
# This is used to decide whether a chain is considered "valid".
ALLOWED_RESIDUES = AMINO_ACIDS | RNA_RESIDUES | DNA_RESIDUES


def process_structure(structure_data: Dict[str, Any]) -> Dict[str, Any]:
    """Extract atom and residue information from a parsed structure.

    Iterates over all models and chains in the Gemmi structure and collects:
    - residues with their residue name, number, and non-hydrogen atom coordinates
    - valid chains (contain at least one residue from `ALLOWED_RESIDUES`)

    Parameters
    ----------
    structure_data : dict
        Dictionary with keys:
        - 'file_path' (str): original input file path
        - 'pdb_id' (str): PDB identifier (lower/upper case as provided upstream)
        - 'structure' (gemmi.Structure): parsed structure object

    Returns
    -------
    dict
        {
          "file_path": str,
          "pdb_id": str,
          "atom_data": [
            {
              "chain_id": str,
              "unique_chain_id": "PDBID:Chain",
              "molecule_name": "UNKNOWN",
              "molecule_type": "UNKNOWN",
              "sequence": "",
              "is_hetatm": False,
              "residues": [
                {
                  "residue_name": str,
                  "residue_number": int,
                  "atoms": [
                    {"atom_name": str, "coordinates": [x, y, z]}, ...
                  ]
                }, ...
              ]
            }, ...
          ]
        }
    """
    file_path = structure_data["file_path"]
    pdb_id = structure_data["pdb_id"]
    structure: gemmi.Structure = structure_data["structure"]

    atom_data: List[Dict[str, Any]] = []

    model_policy = str(config.get("structure_model_policy", "first")).strip().lower()
    models = list(structure)
    if model_policy != "all":
        models = models[:1]

    # Traverse the selected model(s) and chains in the structure.
    for model_index, model in enumerate(models, start=1):
        for chain in model:
            chain_id = chain.name.strip()
            unique_chain_id = f"{pdb_id}:{chain_id}"
            if model_policy == "all":
                unique_chain_id = f"{pdb_id}:model{model_index}:{chain_id}"
            residues: List[Dict[str, Any]] = []
            is_valid_chain = False

            # Collect residues and their heavy-atom coordinates.
            for res in chain:
                original_res_name = res.name.upper()
                res_name = normalize_residue_name(original_res_name)

                atoms = []
                for atom in res:
                    # Exclude hydrogens (H) and deuterium (D) to keep heavy atoms only.
                    if atom.element.name not in ["H", "D"]:
                        atoms.append(
                            {
                                "atom_name": atom.name,
                                "coordinates": list(atom.pos),
                            }
                        )

                if atoms:
                    residue_record = {
                        "residue_name": res_name,
                        "residue_number": res.seqid.num,
                        "atoms": atoms,
                    }
                    if original_res_name != res_name:
                        residue_record["original_residue_name"] = original_res_name
                    residues.append(residue_record)
                    # Mark chain as valid if it contains at least one allowed residue.
                    if res_name in ALLOWED_RESIDUES:
                        is_valid_chain = True

            # Emit chain only if it passed the allow-list criterion.
            if is_valid_chain:
                atom_data.append(
                    {
                        "chain_id": chain_id,
                        "unique_chain_id": unique_chain_id,
                        "model_index": model_index,
                        "molecule_name": "UNKNOWN",    # filled downstream
                        "molecule_type": "UNKNOWN",    # filled downstream
                        "sequence": "",                # filled downstream (e.g., from PDB FASTA)
                        "is_hetatm": False,            # reserved flag; not inferred here
                        "residues": residues,
                    }
                )

    return {
        "file_path": file_path,
        "pdb_id": pdb_id,
        "atom_data": atom_data,
    }
