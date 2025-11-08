"""UniProt matching via BLASTP.

This module:
- Builds (if needed) a local BLAST database from a UniProt FASTA.
- Extracts protein sequences from parsed chain data (3-letter → 1-letter).
- Runs BLASTP to find the best UniProt match for eligible chains.
- Updates chain dictionaries in-place with 'uniprot_id', 'molecule_name',
  and 'molecule_type' when a match is found.

Notes
-----
- Sequences are generated directly from residue names in the parsed structure.
- Non-canonical amino acids mapped to 'O' are replaced with 'X' before BLAST.
- DNA/RNA/DNA-RNA chains are skipped; existing specific NA labels are preserved.
"""

from __future__ import annotations

import os
import subprocess
import uuid
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Dict, List, Optional, Tuple

from Bio.Data import IUPACData
from config_loader import config

# --- Paths from configuration ---
BLAST_DB_PATH: str = config["blast_db_path"]
BLAST_EXECUTABLE: str = config["blastp_executable"]
UNIPROT_FASTA_PATH: str = config["uniprot_fasta_path"]

# 3-letter → 1-letter amino-acid mapping (Biopython)
three_to_one: Dict[str, str] = IUPACData.protein_letters_3to1

# Set of valid amino acid residue codes (3-letter)
AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "SEC", "PYL",
}


def create_blast_database() -> None:
    """Create the UniProt BLAST database if it does not already exist.

    The database is created at: {BLAST_DB_PATH}/uniprot_db.*
    """
    db_path = os.path.join(BLAST_DB_PATH, "uniprot_db")
    if not os.path.exists(db_path + ".pin"):
        result = subprocess.run(
            ["makeblastdb", "-in", UNIPROT_FASTA_PATH, "-dbtype", "prot", "-out", db_path],
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            print(f"Error while creating BLAST database:\n{result.stderr}")
        else:
            print("BLAST database successfully created.")
    else:
        print("BLAST database already exists.")


def extract_sequence_from_parsed_data(chain_data: Dict[str, Any]) -> Optional[str]:
    """Convert residue names (3-letter) to a 1-letter protein sequence.

    Parameters
    ----------
    chain_data : dict
        Chain dictionary with a 'residues' list; each residue dict
        should contain 'residue_name'.

    Returns
    -------
    str | None
        One-letter sequence string if any residues are present; otherwise None.
    """
    seq = "".join(three_to_one.get(res["residue_name"].capitalize(), "X") for res in chain_data["residues"])
    return seq or None


def run_blast_search(query_sequence: str) -> Optional[str]:
    """Run BLASTP for a given protein sequence and return the best UniProt match.

    Parameters
    ----------
    query_sequence : str
        One-letter amino-acid sequence to query.

    Returns
    -------
    str | None
        Matched UniProt accession (e.g., 'P12345') or None if no match is found.
    """
    unique_id = str(uuid.uuid4())[:8]
    query_file = f"query_{unique_id}.fasta"
    output_file = f"blast_results_{unique_id}.txt"

    try:
        with open(query_file, "w", encoding="utf-8") as f:
            f.write(f">query\n{query_sequence}\n")

        blast_cmd = [
            BLAST_EXECUTABLE,
            "-query", query_file,
            "-db", os.path.join(BLAST_DB_PATH, "uniprot_db"),
            "-out", output_file,
            "-evalue", "1e-5",
            "-max_target_seqs", "8",
            "-outfmt", "6",
        ]

        result = subprocess.run(blast_cmd, capture_output=True, text=True)
        if result.returncode != 0:
            return None

        if not os.path.exists(output_file):
            return None

        best_match: Optional[str] = None
        with open(output_file, "r", encoding="utf-8") as f:
            for line in f:
                # BLAST outfmt 6: first hit line; subject id usually like 'sp|P12345|...'
                subject = line.split("\t")[1]
                parts = subject.split("|")
                if len(parts) >= 2:
                    best_match = parts[1]
                break

        return best_match

    except Exception:
        return None

    finally:
        # Clean up temporary files
        if os.path.exists(query_file):
            try:
                os.remove(query_file)
            except Exception:
                pass
        if os.path.exists(output_file):
            try:
                os.remove(output_file)
            except Exception:
                pass


def classify_molecule_type(chain: Dict[str, Any]) -> str:
    """Classify the molecule type of a chain based on its residues.

    Respects existing DNA/RNA subcategories ('DNA', 'RNA', 'DNA/RNA').
    If those are present, they are returned unchanged.

    Parameters
    ----------
    chain : dict
        Chain dictionary with a 'residues' list and an optional 'molecule_type'.

    Returns
    -------
    str
        One of: "Protein", "Nucleic Acid", "DNA", "RNA", "DNA/RNA", or "Unknown".
    """
    # Keep specific NA types if already set
    if chain.get("molecule_type") in ["DNA", "RNA", "DNA/RNA"]:
        return chain["molecule_type"]

    residues = chain.get("residues", [])
    if not residues:
        return "Unknown"

    sample_names = {str(res.get("residue_name", "")).strip().upper() for res in residues}
    return "Protein" if all(res in AMINO_ACIDS for res in sample_names) else "Nucleic Acid"


def parallel_blast_search(parsed_data: List[Dict[str, Any]], max_workers: int = 16) -> None:
    """Run BLASTP in parallel for all eligible chains and update them in-place.

    Eligibility
    -----------
    - Skip chains that already have a 'uniprot_id' (not None/'Unknown').
    - Skip chains classified as nucleic acid ('Nucleic Acid', 'DNA', 'RNA', 'DNA/RNA').
    - Build 1-letter sequences from residues; skip if empty.

    Parameters
    ----------
    parsed_data : list[dict]
        List of structures; each structure contains 'atom_data' (chains).
    max_workers : int, default 16
        Maximum number of threads used for BLAST queries.
    """
    tasks: List[Tuple[int, str, Dict[str, Any]]] = []

    # Prepare candidate sequences
    for structure in parsed_data:
        for chain in structure["atom_data"]:
            if chain.get("uniprot_id") not in [None, "Unknown"]:
                continue
            if chain.get("molecule_type") == "Nucleic Acid":
                continue

            sequence = extract_sequence_from_parsed_data(chain)
            if not sequence:
                continue

            # Replace 'O' with 'X' (rare/non-standard residue mapping)
            sequence = sequence.replace("O", "X")
            tasks.append((len(sequence), sequence, chain))

    # Sort tasks by sequence length when many tasks exist (small → large)
    if len(tasks) > 40:
        tasks.sort(key=lambda t: t[0])

    # Run BLAST in parallel
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [(executor.submit(run_blast_search, seq), chain) for _, seq, chain in tasks]

        for future, chain in futures:
            uniprot_id = future.result()
            if uniprot_id:
                chain["uniprot_id"] = uniprot_id
                chain["molecule_name"] = f"Matched UniProt: {uniprot_id}"
                chain["molecule_type"] = "Protein"
            else:
                # Only override with a broad class if no specific DNA/RNA category exists
                if chain.get("molecule_type") not in ["DNA", "RNA", "DNA/RNA"]:
                    chain["molecule_type"] = classify_molecule_type(chain)
