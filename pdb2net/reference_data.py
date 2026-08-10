"""Cached loaders for external reference data files."""

from __future__ import annotations

import csv
from functools import lru_cache
import re
from typing import Dict, List, Optional

from Bio import SeqIO


def _split_structure_chain_token(token: str) -> tuple[str, str] | None:
    """Split FASTA tokens for legacy and extended IDs at the final underscore."""
    if "_" not in token:
        return None
    structure_id, chain_id = token.rsplit("_", 1)
    if not structure_id or not chain_id:
        return None
    return structure_id, chain_id


@lru_cache(maxsize=None)
def load_valid_pdb_ids(pdb_fasta_path: str) -> set[str]:
    """Load all uppercase PDB IDs observed in a pdb_seqres-style FASTA."""
    valid_pdb_ids: set[str] = set()
    with open(pdb_fasta_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue

            token = line.split()[0][1:]
            parsed = _split_structure_chain_token(token)
            if parsed:
                valid_pdb_ids.add(parsed[0].upper())

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
                parsed = _split_structure_chain_token(fasta_header)
                if parsed:
                    pdb_id, chain_id = parsed
                    formatted_key = f"{pdb_id.lower()}_{chain_id}"
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
def load_pdb_fasta_headers(
    pdb_fasta_path: str,
    pdb_ids: tuple[str, ...] | None = None,
) -> Dict[str, Dict[str, str]]:
    """Return per-chain header metadata without retaining reference sequences.

    Molecule classification only consumes the FASTA header text. Keeping a
    dedicated lightweight lookup prevents every web job from holding the full
    PDB SEQRES sequence corpus in memory unnecessarily.
    """
    pdb_headers: Dict[str, Dict[str, str]] = {}
    wanted_ids = {value.lower() for value in pdb_ids} if pdb_ids is not None else None
    with open(pdb_fasta_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            parts = line.split()
            fasta_header = parts[0][1:]
            parsed = _split_structure_chain_token(fasta_header)
            if not parsed:
                continue
            pdb_id, chain_id = parsed
            if wanted_ids is not None and pdb_id.lower() not in wanted_ids:
                continue
            formatted_key = f"{pdb_id.lower()}_{chain_id}"
            pdb_headers[formatted_key] = {"info": " ".join(parts[1:])}
    return pdb_headers


_UNIPROT_ACCESSION_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[0-9]+)?$"
)


def is_uniprot_accession(value: object) -> bool:
    """Return whether a complete UniProtKB accession (including isoforms) was supplied."""
    return _UNIPROT_ACCESSION_RE.fullmatch(str(value or "").strip().upper()) is not None


@lru_cache(maxsize=None)
def load_sifts_segments(tsv_path: str) -> Dict[str, tuple[Dict[str, str], ...]]:
    """Load every valid external SIFTS segment for each PDB/author-chain pair."""
    segments: Dict[str, List[Dict[str, str]]] = {}
    header: list[str] | None = None
    with open(tsv_path, "r", encoding="utf-8", errors="ignore") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or len(row) < 3:
                continue
            cleaned = [value.strip() for value in row]
            first = cleaned[0].lstrip("#").strip().upper()
            second = cleaned[1].strip().upper()
            if first == "PDB" and second == "CHAIN":
                header = [value.lstrip("#").strip().lower() for value in cleaned]
                continue
            if cleaned[0].startswith("#"):
                continue
            pdb_id = cleaned[0].lower()
            chain = cleaned[1]
            uniprot_id = cleaned[2].upper()
            if not pdb_id or not chain or not is_uniprot_accession(uniprot_id):
                continue
            key = f"{pdb_id}_{chain}"
            if header is None:
                names = [
                    "pdb",
                    "chain",
                    "sp_primary",
                    "res_beg",
                    "res_end",
                    "pdb_beg",
                    "pdb_end",
                    "sp_beg",
                    "sp_end",
                ]
            else:
                names = header
            segment = {
                names[index] if index < len(names) and names[index] else f"column_{index + 1}": value
                for index, value in enumerate(cleaned)
            }
            segment["pdb"] = pdb_id
            segment["chain"] = chain
            segment["accession"] = uniprot_id
            segments.setdefault(key, []).append(segment)

    result: Dict[str, tuple[Dict[str, str], ...]] = {}
    for key, values in segments.items():
        unique = {
            tuple(sorted(segment.items())): segment
            for segment in values
        }
        result[key] = tuple(unique[item] for item in sorted(unique))
    return result


@lru_cache(maxsize=None)
def load_sifts_mapping(tsv_path: str) -> Dict[str, str]:
    """Compatibility view containing only unambiguous chain mappings."""
    pdb_to_uniprot: Dict[str, str] = {}
    for key, segments in load_sifts_segments(tsv_path).items():
        accessions = {segment["accession"] for segment in segments}
        if len(accessions) == 1:
            pdb_to_uniprot[key] = next(iter(accessions))

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
