"""Assign molecule names/types and UniProt IDs using SIFTS + FASTA sources.

This module:
- Loads a SIFTS TSV mapping (PDB-chain → UniProt).
- Lazily loads UniProt FASTA (UniProt ID → protein name).
- Parses `pdb_seqres.txt` FASTA to obtain per-chain names/types and sequences.
- Derives molecule info for each chain:
    * If SIFTS provides a UniProt ID, use it; prefer UniProt name over PDB FASTA name.
    * Otherwise, infer molecule type from the PDB FASTA header (Protein / DNA / RNA / DNA/RNA / Nucleic Acid / Unknown).

Notes
-----
- All mapping files are loaded once (lazy-loading).
- The functions mutate the provided `combined_data` structure in place.
"""

from __future__ import annotations

import csv
import os
import re
from typing import Any, Dict, List, Optional, Tuple

from Bio import SeqIO

from config_loader import config

# --- Paths from configuration ---
PDB_FASTA_PATH: str = config["pdb_fasta_path"]
UNIPROT_FASTA_PATH: str = config["uniprot_fasta_path"]
SIFTS_TSV_PATH: str = config["sifts_tsv_path"]

# --- In-memory lookup tables ---
pdb_to_uniprot: Dict[str, str] = {}   # "pdbid_CHAIN" → UniProt ID
uniprot_dict: Dict[str, str] = {}     # UniProt ID → Protein display name

# --- Lazy-loading guards ---
_sifts_loaded: bool = False
_uniprot_loaded: bool = False


def load_sifts_mapping(tsv_path: str) -> None:
    """Load SIFTS PDB→UniProt chain mapping from a TSV file (once)."""
    global pdb_to_uniprot, _sifts_loaded
    if _sifts_loaded:
        return

    with open(tsv_path, "r", encoding="utf-8", errors="ignore") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or len(row) < 3:
                continue
            pdb_id = row[0].strip().lower()
            chain = row[1].strip().upper()
            uniprot_id = row[2].strip()
            key = f"{pdb_id}_{chain}"
            pdb_to_uniprot[key] = uniprot_id

    _sifts_loaded = True


def load_uniprot_fasta(fasta_path: str) -> None:
    """Load UniProt FASTA to map UniProt ID → protein name (once)."""
    global uniprot_dict, _uniprot_loaded
    if _uniprot_loaded:
        return

    for record in SeqIO.parse(fasta_path, "fasta"):
        # FASTA record.id often like: sp|P12345|...
        parts = record.id.split("|")
        if len(parts) >= 2:
            uniprot_id = parts[1]
        else:
            # Fallback: if no pipe fields, use the whole id
            uniprot_id = record.id

        # Description: take everything after the first space as "name"
        protein_name = record.description.split(" ", 1)[1] if " " in record.description else record.description
        uniprot_dict[uniprot_id] = protein_name

    _uniprot_loaded = True


def load_pdb_fasta(pdb_fasta_path: str) -> Dict[str, Dict[str, str]]:
    """Parse pdb_seqres.txt-like FASTA and return per-chain info/sequences.

    Returns
    -------
    dict
        {
          "pdbid_CHAIN": {"info": <header tail>, "sequence": <seq>},
          ...
        }
    """
    pdb_sequences: Dict[str, Dict[str, str]] = {}
    with open(pdb_fasta_path, "r", encoding="utf-8", errors="ignore") as f:
        current_key: Optional[str] = None
        current_seq: List[str] = []
        for line in f:
            if line.startswith(">"):
                # Flush previous record
                if current_key and current_seq:
                    pdb_sequences[current_key]["sequence"] = "".join(current_seq)

                parts = line.split()
                fasta_header = parts[0][1:]  # strip '>'
                if "_" in fasta_header:
                    pdb_id, chain_id = fasta_header.split("_", 1)
                    formatted_key = f"{pdb_id.lower()}_{chain_id.upper()}"
                    pdb_sequences[formatted_key] = {"info": " ".join(parts[1:]), "sequence": ""}
                    current_key = formatted_key
                    current_seq = []
                else:
                    current_key = None
                    current_seq = []
            else:
                if current_key is not None:
                    current_seq.append(line.strip())

        # Flush last record
        if current_key and current_seq:
            pdb_sequences[current_key]["sequence"] = "".join(current_seq)

    return pdb_sequences


def determine_from_fasta(search_key: str, pdb_fasta: Dict[str, Dict[str, str]]) -> Tuple[str, str, Optional[str]]:
    """Determine display name and molecule type from PDB FASTA header.

    Parameters
    ----------
    search_key : str
        Key formatted as "pdbid_CHAIN".
    pdb_fasta : dict
        Output of `load_pdb_fasta`.

    Returns
    -------
    tuple[str, str, str | None]
        (cleaned_name, molecule_type, None)
    """
    if search_key in pdb_fasta:
        fasta_info = pdb_fasta[search_key]["info"]
        # sequence = pdb_fasta[search_key]["sequence"]  # not used here

        # Remove mol:... and length:... tokens from the info field to keep it readable
        cleaned_info = re.sub(r"mol:\w+\s*", "", fasta_info)
        cleaned_info = re.sub(r"length:\d+\s*", "", cleaned_info).strip()

        lower = fasta_info.lower()
        if "mol:protein" in lower:
            molecule_type = "Protein"
        elif "mol:na" in lower:
            # Try to extract subcategory when available
            if "dna/rna" in lower:
                molecule_type = "DNA/RNA"
            elif "dna" in lower:
                molecule_type = "DNA"
            elif "rna" in lower:
                molecule_type = "RNA"
            elif "polyribonucleotide" in lower:
                molecule_type = "RNA"
            elif "polydeoxyribonucleotide" in lower:
                molecule_type = "DNA"
            else:
                molecule_type = "Nucleic Acid"
        else:
            molecule_type = "Unknown"

        return cleaned_info, molecule_type, None
    return "Unknown", "Unknown", None


def determine_molecule_info(pdb_id: str, chain_id: str, pdb_fasta: Dict[str, Dict[str, str]]) -> Tuple[str, str, Optional[str]]:
    """Combine SIFTS and PDB FASTA to determine name/type/UniProt for a chain."""
    search_key = f"{pdb_id.lower()}_{chain_id.upper()}"
    fasta_name, mol_type, _ = determine_from_fasta(search_key, pdb_fasta)

    name = fasta_name
    uniprot_id: Optional[str] = None

    # Only trust the 4-char PDB convention for SIFTS lookup
    if len(pdb_id) == 4:
        uniprot_id = pdb_to_uniprot.get(search_key)
        if uniprot_id:
            better_name = uniprot_dict.get(uniprot_id)
            if better_name and better_name != "Unknown Protein":
                name = better_name
            mol_type = "Protein"

    return name, mol_type, uniprot_id


def process_molecule_info(combined_data: List[Dict[str, Any]]) -> None:
    """Assign molecule names and types to all chains in the dataset (in place).

    This function lazy-loads:
      - SIFTS mapping (TSV)
      - UniProt FASTA (for display names)
      - PDB FASTA (per-run)

    Parameters
    ----------
    combined_data : list[dict]
        Per-structure dictionaries with keys:
          - "pdb_id": str
          - "atom_data": list[dict] of chain dicts with at least "chain_id".
    """
    load_sifts_mapping(SIFTS_TSV_PATH)
    load_uniprot_fasta(UNIPROT_FASTA_PATH)

    pdb_fasta = load_pdb_fasta(PDB_FASTA_PATH)

    for structure_data in combined_data:
        pdb_id = structure_data["pdb_id"].lower()
        for chain in structure_data["atom_data"]:
            chain_id = chain["chain_id"].upper()
            name, mol_type, uniprot_id = determine_molecule_info(pdb_id, chain_id, pdb_fasta)

            chain["molecule_name"] = name
            chain["molecule_type"] = mol_type
            chain["uniprot_id"] = uniprot_id

    # Optional: brief sample output for debugging/verification (limited to first ~20 structures)
    print("\nUniProt assignments (example for up to 20 input files):")
    for i, structure_data in enumerate(combined_data):
        if i >= 20:
            break
        pdb_id = structure_data["pdb_id"]
        for chain in structure_data["atom_data"]:
            print(
                f"  {pdb_id}_{chain['chain_id']}: {chain['molecule_name']} "
                f"({chain['molecule_type']}) UniProt-ID: {chain['uniprot_id']}"
            )
