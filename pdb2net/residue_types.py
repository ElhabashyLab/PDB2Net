"""Shared residue constants and lightweight classification helpers."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Iterable


AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "SEC", "PYL",
}

DNA_RESIDUES = {"DA", "DT", "DG", "DC", "DI"}
RNA_RESIDUES = {"A", "U", "G", "C", "I"}
NUCLEIC_ACID_TYPES = {"Nucleic Acid", "DNA", "RNA", "DNA/RNA"}

MODIFIED_RESIDUES_3TO3 = {
    "MSE": "MET",
    "SEP": "SER",
    "TPO": "THR",
    "PTR": "TYR",
    "CSO": "CYS",
    "HYP": "PRO",
}


def normalize_residue_name(name: object) -> str:
    """Return an uppercase residue code with common modified residues normalized."""
    residue = str(name or "").strip().upper()
    return MODIFIED_RESIDUES_3TO3.get(residue, residue)


def count_polymer_lengths(residue_names: Iterable[object]) -> tuple[int, int]:
    """Count observed residue records by polymer role/type, or legacy names.

    Explicit nonpolymers and other polymer types do not contribute. Without
    usable type metadata, retain the existing name-based interpretation.
    Atom alternatives and occupancy do not affect the number of records.
    """
    aa_count = 0
    nt_count = 0
    for raw_name in residue_names:
        if isinstance(raw_name, Mapping):
            if raw_name.get("is_polymer") is False:
                continue
            kind = raw_name.get("polymer_type")
            if kind in {"PeptideL", "PeptideD", "CyclicPseudoPeptide"}:
                aa_count += 1
                continue
            if kind in {"Dna", "Rna", "DnaRnaHybrid"}:
                nt_count += 1
                continue
            if kind not in (None, "", "Unknown"):
                continue
            raw_name = raw_name.get("residue_name")
        name = normalize_residue_name(raw_name)
        if name in AMINO_ACIDS:
            aa_count += 1
        elif name in DNA_RESIDUES or name in RNA_RESIDUES:
            nt_count += 1
    return aa_count, nt_count
