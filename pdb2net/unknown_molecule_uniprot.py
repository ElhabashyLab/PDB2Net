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

import re
from typing import Any, Dict, List, Optional, Tuple

from .config_loader import config
from .network_annotations import network_annotation_config
from .reference_data import (
    load_pdb_fasta_headers as _load_pdb_fasta_headers,
    load_sifts_mapping as _load_sifts_mapping,
    load_sifts_segments as _load_sifts_segments,
    load_uniprot_names,
)

# --- In-memory lookup tables ---
pdb_to_uniprot: Dict[str, str] = {}   # "pdbid_CHAIN" → UniProt ID
pdb_to_sifts_segments: Dict[str, tuple[Dict[str, str], ...]] = {}
uniprot_dict: Dict[str, str] = {}     # UniProt ID → Protein display name

# --- Path-keyed lazy-loading guards ---
_sifts_loaded_path: str | None = None
_uniprot_loaded_path: str | None = None


def load_sifts_mapping(tsv_path: str) -> None:
    """Load SIFTS PDB→UniProt chain mapping from a TSV file (once)."""
    global pdb_to_uniprot, pdb_to_sifts_segments, _sifts_loaded_path
    if _sifts_loaded_path == tsv_path:
        return

    # The reference loader already owns an immutable-by-convention cached
    # mapping. Reuse it instead of duplicating a potentially large hash table.
    pdb_to_uniprot = _load_sifts_mapping(tsv_path)
    pdb_to_sifts_segments = _load_sifts_segments(tsv_path)
    _sifts_loaded_path = tsv_path


def load_uniprot_fasta(fasta_path: str) -> None:
    """Load UniProt FASTA to map UniProt ID → protein name (once)."""
    global uniprot_dict, _uniprot_loaded_path
    if _uniprot_loaded_path == fasta_path:
        return

    # Reuse the cached mapping directly. A shallow ``dict(...)`` copy keeps all
    # strings alive and duplicates the large hash table for no functional gain.
    uniprot_dict = load_uniprot_names(fasta_path)
    _uniprot_loaded_path = fasta_path


def load_pdb_fasta(
    pdb_fasta_path: str,
    pdb_ids: tuple[str, ...] | None = None,
) -> Dict[str, Dict[str, str]]:
    """Parse PDB SEQRES headers needed for per-chain molecule classification.

    Returns
    -------
    dict
        {
          "pdbid_CHAIN": {"info": <header tail>},
          ...
        }
    """
    return _load_pdb_fasta_headers(pdb_fasta_path, pdb_ids)


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


def determine_molecule_info(
    pdb_id: str,
    chain_id: str,
    pdb_fasta: Dict[str, Dict[str, str]],
    *,
    pdb_aliases: tuple[str, ...] | None = None,
) -> Tuple[str, str, Optional[str]]:
    """Combine SIFTS and PDB FASTA to determine name/type/UniProt for a chain."""
    details = _determine_molecule_info_details(
        pdb_id,
        chain_id,
        pdb_fasta,
        pdb_aliases=pdb_aliases,
    )
    return details["name"], details["molecule_type"], details["uniprot_id"]


def _determine_molecule_info_details(
    pdb_id: str,
    chain_id: str,
    pdb_fasta: Dict[str, Dict[str, str]],
    *,
    pdb_aliases: tuple[str, ...] | None = None,
) -> Dict[str, Any]:
    """Return FASTA provenance and the full external-SIFTS mapping state."""
    aliases = pdb_aliases or (pdb_id,)
    search_keys = [f"{alias.lower()}_{chain_id}" for alias in aliases if alias]
    fasta_name, mol_type, _ = "Unknown", "Unknown", None
    for search_key in search_keys:
        fasta_name, mol_type, _ = determine_from_fasta(search_key, pdb_fasta)
        if search_key in pdb_fasta:
            break

    accessions: set[str] = set()
    for search_key in search_keys:
        for segment in pdb_to_sifts_segments.get(search_key, ()):
            accession = str(segment.get("accession") or "").strip()
            if accession:
                accessions.add(accession)
    # Compatibility for callers/tests that inject only the historical mapping.
    if not accessions:
        accessions.update(
            str(pdb_to_uniprot[key])
            for key in search_keys
            if key in pdb_to_uniprot and pdb_to_uniprot[key]
        )

    uniprot_id = next(iter(accessions)) if len(accessions) == 1 else None
    name = fasta_name
    status = "none"
    if uniprot_id:
        better_name = uniprot_dict.get(uniprot_id)
        if not better_name and "-" in uniprot_id:
            better_name = uniprot_dict.get(uniprot_id.split("-", 1)[0])
        if better_name and better_name != "Unknown Protein":
            name = better_name
        elif name == "Unknown":
            name = f"UniProt: {uniprot_id}"
        mol_type = "Protein"
        status = "unique_external_mapping"
    elif len(accessions) > 1:
        mol_type = "Protein"
        status = "ambiguous_external_mapping"
    return {
        "name": name,
        "fasta_name": fasta_name,
        "molecule_type": mol_type,
        "uniprot_id": uniprot_id,
        "external_status": status,
        "external_accessions": sorted(accessions),
    }


def process_molecule_info(
    combined_data: List[Dict[str, Any]],
    *,
    pdb_fasta_headers: Dict[str, Dict[str, str]] | None = None,
) -> None:
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
    pdb_fasta_path = str(config.get("pdb_fasta_path") or "")
    sifts_tsv_path = str(config.get("sifts_tsv_path") or "")
    uniprot_fasta_path = str(config.get("uniprot_fasta_path") or "")
    load_sifts_mapping(sifts_tsv_path)
    load_uniprot_fasta(uniprot_fasta_path)

    requested: set[str] = set()
    for structure in combined_data:
        requested.add(str(structure.get("pdb_id") or "").lower())
        identity = structure.get("structure_identity", {})
        if isinstance(identity, dict):
            for field in ("canonical_id", "legacy_id"):
                if identity.get(field):
                    requested.add(str(identity[field]).lower())
    requested_pdb_ids = tuple(sorted(requested))
    pdb_fasta = (
        pdb_fasta_headers
        if pdb_fasta_headers is not None
        else load_pdb_fasta(pdb_fasta_path, requested_pdb_ids)
    )

    annotation_cfg = network_annotation_config()
    for structure_data in combined_data:
        pdb_id = structure_data["pdb_id"].lower()
        identity = structure_data.get("structure_identity", {})
        aliases: list[str] = [pdb_id]
        if isinstance(identity, dict):
            for field in ("legacy_id", "canonical_id"):
                value = str(identity.get(field) or "").strip().lower()
                if value and value not in aliases:
                    aliases.append(value)
        for chain in structure_data["atom_data"]:
            chain_id = chain["chain_id"]
            details = _determine_molecule_info_details(
                pdb_id,
                chain_id,
                pdb_fasta,
                pdb_aliases=tuple(aliases),
            )

            name = str(details["name"])
            fasta_name = str(details["fasta_name"])
            mol_type = str(details["molecule_type"])
            uniprot_id = details["uniprot_id"]
            external_status = str(details["external_status"])
            external_accessions = list(details["external_accessions"])
            chain["external_sifts_status"] = external_status
            if external_accessions:
                chain["external_sifts_accessions"] = external_accessions

            embedded_status = str(chain.get("embedded_uniprot_status") or "")
            if annotation_cfg["use_embedded_sifts"] and embedded_status == "unique_best_mapping":
                embedded_accession = str(chain.get("uniprot_id") or "")
                external_mismatch = bool(
                    external_accessions and external_accessions != [embedded_accession]
                )
                embedded_name = uniprot_dict.get(embedded_accession)
                if not embedded_name and "-" in embedded_accession:
                    embedded_name = uniprot_dict.get(embedded_accession.split("-", 1)[0])
                chain["molecule_name"] = (
                    fasta_name
                    if external_mismatch and fasta_name != "Unknown"
                    else (
                        f"UniProt: {embedded_accession}"
                        if external_mismatch
                        else (
                            embedded_name
                            if embedded_name and embedded_name != "Unknown Protein"
                            else (
                                fasta_name
                                if fasta_name != "Unknown"
                                else f"UniProt: {embedded_accession}"
                            )
                        )
                    )
                )
                chain["molecule_type"] = "Protein"
                if external_mismatch:
                    chain.setdefault("annotation_warnings", []).append(
                        {
                            "code": "EMBEDDED_EXTERNAL_UNIPROT_MISMATCH",
                            "message": (
                                f"{chain.get('unique_chain_id')} maps to {embedded_accession} in embedded SIFTS "
                                f"but external SIFTS reports {', '.join(external_accessions)}; the embedded "
                                "accession was retained without borrowing an unrelated UniProt name."
                            ),
                        }
                    )
                continue
            if annotation_cfg["use_embedded_sifts"] and embedded_status == "ambiguous_multi_mapping":
                # An ambiguous embedded mapping blocks lower-priority identity
                # sources.  Do not borrow a UniProt display name discovered via
                # external SIFTS when no accession is actually assigned.
                chain["molecule_name"] = fasta_name
                chain["molecule_type"] = "Protein"
                chain["uniprot_id"] = None
                continue

            if external_status == "ambiguous_external_mapping":
                chain["molecule_name"] = fasta_name
                chain["molecule_type"] = "Protein"
                chain["uniprot_id"] = None
                chain["annotation_source"] = "external_sifts"
                chain["matched_database"] = "UniProtKB"
                chain["matched_id"] = ",".join(external_accessions)
                chain["annotation_confidence"] = "ambiguous"
                continue

            chain["molecule_name"] = name
            chain["molecule_type"] = mol_type
            chain["uniprot_id"] = uniprot_id
            if uniprot_id:
                chain["annotation_source"] = "sifts"
                chain["matched_database"] = "UniProtKB"
                chain["matched_id"] = uniprot_id
                chain["representative_accession"] = uniprot_id
                chain["annotation_confidence"] = "high"

    # # Optional: brief sample output for debugging/verification (limited to first ~20 structures)
    # print("\nUniProt assignments (example for up to 20 input files):")
    # for i, structure_data in enumerate(combined_data):
    #     if i >= 20:
    #         break
    #     pdb_id = structure_data["pdb_id"]
    #     for chain in structure_data["atom_data"]:
    #         print(
    #             f"  {pdb_id}_{chain['chain_id']}: {chain['molecule_name']} "
    #             f"({chain['molecule_type']}) UniProt-ID: {chain['uniprot_id']}"
    #         )
