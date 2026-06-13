"""Cached loaders for external reference data files."""

from __future__ import annotations

import csv
from functools import lru_cache
from typing import Dict, List, Optional

from Bio import SeqIO


@lru_cache(maxsize=None)
def load_valid_pdb_ids(pdb_fasta_path: str) -> set[str]:
    """Load all uppercase PDB IDs observed in a pdb_seqres-style FASTA."""
    valid_pdb_ids: set[str] = set()
    with open(pdb_fasta_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue

            token = line.split()[0][1:]
            parts = token.split("_")
            if parts and len(parts[0]) == 4:
                valid_pdb_ids.add(parts[0].upper())

    return valid_pdb_ids


@lru_cache(maxsize=None)
def load_pdb_fasta(pdb_fasta_path: str) -> Dict[str, Dict[str, str]]:
    """Parse pdb_seqres.txt-like FASTA and return per-chain info/sequences."""
    pdb_sequences: Dict[str, Dict[str, str]] = {}
    with open(pdb_fasta_path, "r", encoding="utf-8", errors="ignore") as handle:
        current_key: Optional[str] = None
        current_seq: List[str] = []
        for line in handle:
            if line.startswith(">"):
                if current_key and current_seq:
                    pdb_sequences[current_key]["sequence"] = "".join(current_seq)

                parts = line.split()
                fasta_header = parts[0][1:]
                if "_" in fasta_header:
                    pdb_id, chain_id = fasta_header.split("_", 1)
                    formatted_key = f"{pdb_id.lower()}_{chain_id.upper()}"
                    pdb_sequences[formatted_key] = {"info": " ".join(parts[1:]), "sequence": ""}
                    current_key = formatted_key
                    current_seq = []
                else:
                    current_key = None
                    current_seq = []
            elif current_key is not None:
                current_seq.append(line.strip())

        if current_key and current_seq:
            pdb_sequences[current_key]["sequence"] = "".join(current_seq)

    return pdb_sequences


@lru_cache(maxsize=None)
def load_sifts_mapping(tsv_path: str) -> Dict[str, str]:
    """Load SIFTS PDB-to-UniProt chain mapping from a TSV file."""
    pdb_to_uniprot: Dict[str, str] = {}
    with open(tsv_path, "r", encoding="utf-8", errors="ignore") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or len(row) < 3:
                continue
            pdb_id = row[0].strip().lower()
            chain = row[1].strip().upper()
            uniprot_id = row[2].strip()
            key = f"{pdb_id}_{chain}"
            pdb_to_uniprot[key] = uniprot_id

    return pdb_to_uniprot


@lru_cache(maxsize=None)
def load_uniprot_names(fasta_path: str) -> Dict[str, str]:
    """Load UniProt FASTA to map UniProt ID to protein display name."""
    uniprot_names: Dict[str, str] = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        parts = record.id.split("|")
        if len(parts) >= 2:
            uniprot_id = parts[1]
        else:
            uniprot_id = record.id

        protein_name = record.description.split(" ", 1)[1] if " " in record.description else record.description
        uniprot_names[uniprot_id] = protein_name

    return uniprot_names
