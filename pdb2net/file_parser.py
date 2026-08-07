"""File parsing utilities for PDB/mmCIF inputs.

This module:
- Validates input paths/extensions.
- Derives PDB IDs from filenames or file contents, cross-checking against
  `pdb_seqres.txt` (configured via `config["pdb_fasta_path"]`).
- Parses structures with Gemmi and returns a list of dicts:
  {"file_path", "pdb_id", "structure"}.

Notes
-----
- If no valid 4-character PDB ID can be found (or it is not present in
  `pdb_seqres.txt`), the filename (without extension, uppercased) is used
  as a fallback PDB ID with a warning.
- The allow-list of valid PDB IDs is loaded once on import for efficiency.
"""

from __future__ import annotations

import csv
import gzip
import os
import re
from typing import IO, Any, Dict, List, Optional

import gemmi

from .config_loader import config
from .reference_data import load_valid_pdb_ids as _load_valid_pdb_ids

# --- Allowed file extensions for structural inputs ---
ALLOWED_EXTENSIONS = {".pdb", ".ent", ".cif", ".mmcif"}

# Path to `pdb_seqres.txt` (set in configuration)
PDB_FASTA_PATH = config["pdb_fasta_path"]


def load_valid_pdb_ids() -> set[str]:
    """Load all valid PDB IDs from `pdb_seqres.txt`.

    Returns
    -------
    set[str]
        Uppercase PDB IDs observed in the FASTA header lines.
    """
    try:
        return _load_valid_pdb_ids(PDB_FASTA_PATH)
    except Exception as e:
        print(f"Error loading pdb_seqres.txt: {e}")
        return set()


# Loaded only when raw/archive filenames need canonical-ID validation.  A pure
# precomputed cache-hit process must not scan the full PDB FASTA at import time.
# Keeping this public name preserves the existing test/legacy monkeypatch hook.
VALID_PDB_IDS: set[str] | None = None


def _valid_pdb_ids() -> set[str]:
    global VALID_PDB_IDS
    if VALID_PDB_IDS is None:
        VALID_PDB_IDS = load_valid_pdb_ids()
    return VALID_PDB_IDS


def _structure_extension(file_path: str) -> str:
    """Return the structure extension, ignoring one optional gzip suffix."""
    filename = os.path.basename(file_path)
    if filename.lower().endswith(".gz"):
        filename = filename[:-3]
    return os.path.splitext(filename)[1].lower()


def _open_structure_text(file_path: str) -> IO[str]:
    """Open plain or gzip-compressed structure content as text."""
    if file_path.lower().endswith(".gz"):
        return gzip.open(file_path, "rt", encoding="utf-8", errors="ignore")
    return open(file_path, "r", encoding="utf-8", errors="ignore")


def is_valid_file(file_path: str) -> bool:
    """Check whether the path has a supported structure-file extension.

    Parameters
    ----------
    file_path : str
        Path to the input file.

    Returns
    -------
    bool
        True if the extension is one of {'.pdb', '.ent', '.cif', '.mmcif'},
        optionally followed by '.gz'.
    """
    return _structure_extension(file_path) in ALLOWED_EXTENSIONS


def extract_pdb_id_from_filename(file_path: str) -> Optional[str]:
    """Try to extract a valid 4-character PDB ID from the filename.

    The filename (without extension) may be a bare PDB ID or one of the
    canonical archive forms ``pdb1abc`` and ``pdb_00001abc``. The resolved ID
    is still checked against ``VALID_PDB_IDS``.

    Parameters
    ----------
    file_path : str
        Path to the input file.

    Returns
    -------
    str | None
        Uppercased PDB ID if both pattern and membership in VALID_PDB_IDS hold;
        otherwise None. Emits a warning if the pattern matches but the ID is
        not present in `pdb_seqres.txt`.
    """
    filename = os.path.basename(file_path)
    stem = filename.split(".")[0]
    patterns = (
        r"^([0-9][A-Za-z0-9]{3})$",
        r"^pdb([0-9][A-Za-z0-9]{3})$",
        r"^pdb_0000([0-9][A-Za-z0-9]{3})$",
    )
    for pattern in patterns:
        match = re.match(pattern, stem, flags=re.IGNORECASE)
        if not match:
            continue
        pdb_id = match.group(1).upper()
        if pdb_id in _valid_pdb_ids():
            return pdb_id
        print(f"Warning: File {filename} contains a canonical-looking PDB ID not found in pdb_seqres.txt: {pdb_id}.")
        break
    return None


def extract_pdb_id_from_file(file_path: str) -> Optional[str]:
    """Extract a PDB ID by scanning file contents (HEADER/data_/loop constructs).

    For PDB (.pdb/.ent):
        - Uses the last token on 'HEADER' lines as PDB ID candidate.

    For CIF/mmCIF (.cif/.mmcif):
        - Uses the string after 'data_' in a data-block header, or
        - the value on lines containing '_entry.id'.

    The first 4-character candidate found is validated against VALID_PDB_IDS.

    Parameters
    ----------
    file_path : str
        Path to the input file.

    Returns
    -------
    str | None
        Uppercased PDB ID when 4 characters and present in VALID_PDB_IDS;
        otherwise None (a warning is printed if a 4-char candidate is invalid).
    """
    ext = _structure_extension(file_path)

    try:
        with _open_structure_text(file_path) as f:
            for line in f:
                pdb_id: Optional[str] = None

                if ext in {".pdb", ".ent"} and line.startswith("HEADER"):
                    # PDB convention: last token is typically the 4-char id
                    pdb_id = line.split()[-1].strip().upper()

                elif ext in {".cif", ".mmcif"}:
                    lower = line.lower()
                    if lower.startswith("data_"):
                        # Example: data_1abc → take '1abc'
                        pdb_id = line.split("_", 1)[1].strip().upper()
                    elif "_entry.id" in lower:
                        # Example: _entry.id 1abc
                        parts = line.split()
                        if parts:
                            pdb_id = parts[-1].strip().upper()

                if pdb_id:
                    if len(pdb_id) == 4 and pdb_id in _valid_pdb_ids():
                        return pdb_id
                    if len(pdb_id) == 4:
                        print(
                            "Warning: File "
                            f"{file_path} contains a canonical-looking PDB ID not found in pdb_seqres.txt: {pdb_id}."
                        )
                    return None

    except Exception as e:
        print(f"Error reading {file_path}: {e}")

    return None


def get_pdb_id(file_path: str) -> str:
    """Determine the PDB ID for a file, preferring the filename over contents.

    If neither yields a valid ID (or the ID is not part of VALID_PDB_IDS),
    the function falls back to using the uppercased filename stem and prints
    a warning.

    Parameters
    ----------
    file_path : str
        Path to the PDB/mmCIF file.

    Returns
    -------
    str
        Chosen PDB ID (may be the uppercased filename as a fallback).
    """
    pdb_id = extract_pdb_id_from_filename(file_path)
    if pdb_id:
        return pdb_id

    pdb_id = extract_pdb_id_from_file(file_path)
    if pdb_id:
        return pdb_id

    # Fallback: use filename stem (uppercased) if no valid PDB ID was found
    filename_stem = os.path.basename(file_path).split(".")[0].upper()
    print(f"Info: No canonical PDB ID found. Using filename as structure ID: {filename_stem}")
    return filename_stem


def parse_structure(file_path: str, pdb_id: str) -> Optional[gemmi.Structure]:
    """Load a plain or gzip-compressed PDB/mmCIF structure with Gemmi.

    Parameters
    ----------
    file_path : str
        Path to the structural file.
    pdb_id : str
        The (resolved) PDB identifier; currently not used by Gemmi here, but
        kept to preserve the function signature used upstream.

    Returns
    -------
    gemmi.Structure | None
        Parsed structure on success, otherwise None with an error message.
    """
    try:
        structure = gemmi.read_structure(file_path)
        structure.setup_entities()
        return structure
    except Exception as e:
        print(f"Error parsing {file_path}: {e}")
        return None


def read_files_from_csv(csv_path: str) -> List[Dict[str, Any]]:
    """Read input file paths from a CSV and parse valid structures.

    The CSV must contain a column named 'file_path'. For each entry, a structure
    is parsed when:
      - the path exists and has a valid extension, and
      - a PDB ID can be derived (filename → content → fallback to filename).

    Parameters
    ----------
    csv_path : str
        Path to the CSV file.

    Returns
    -------
    list[dict]
        Each dict has keys: {"file_path", "pdb_id", "structure"}.
    """
    file_paths: List[str] = []

    with open(csv_path, newline="", encoding="utf-8") as csvfile:
        reader = csv.DictReader(csvfile)
        if "file_path" not in reader.fieldnames:
            raise ValueError("CSV file must contain a 'file_path' column.")

        for row in reader:
            file_paths.append(row["file_path"])

    structures: List[Dict[str, Any]] = []
    for file_path in file_paths:
        if is_valid_file(file_path):
            pdb_id = get_pdb_id(file_path)
            if pdb_id:
                structure = parse_structure(file_path, pdb_id)
                if structure:
                    structures.append({"file_path": file_path, "pdb_id": pdb_id, "structure": structure})
            else:
                print(f"Skipping file (no valid PDB ID found): {file_path}")
        else:
            print(f"Skipping file (invalid extension): {file_path}")

    return structures


def read_files_from_folder(folder_path: str) -> List[Dict[str, Any]]:
    """Scan a folder for PDB/mmCIF files and parse valid structures.

    Parameters
    ----------
    folder_path : str
        Path to the input directory.

    Returns
    -------
    list[dict]
        Each dict has keys: {"file_path", "pdb_id", "structure"}.
    """
    structures: List[Dict[str, Any]] = []

    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)

        if os.path.isfile(file_path) and is_valid_file(file_path):
            pdb_id = get_pdb_id(file_path)
            if pdb_id:
                structure = parse_structure(file_path, pdb_id)
                if structure:
                    structures.append({"file_path": file_path, "pdb_id": pdb_id, "structure": structure})
            else:
                print(f"Skipping file (no valid PDB ID found): {file_path}")
        else:
            print(f"Skipping file (invalid extension or not a file): {file_path}")

    return structures
