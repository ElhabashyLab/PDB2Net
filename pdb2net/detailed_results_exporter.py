"""Bounded detailed atom-interaction CSV export utilities.

Detailed interaction tables can be substantially larger than the network
artifacts they accompany.  Rows are therefore serialized incrementally and
published atomically instead of being accumulated in a DataFrame.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
import io
import os
from pathlib import Path
import shutil
from typing import Any, BinaryIO, Dict, List, Sequence

import numpy as np
from scipy.spatial import cKDTree

from .artifact_names import portable_artifact_stem
from .input_contract import InputValidationError

DETAILED_INTERACTION_COLUMNS = (
    "PDB_ID",
    "Chain_A",
    "Residue_A",
    "Atom_A",
    "Chain_B",
    "Residue_B",
    "Atom_B",
    "Distance",
    "UniProt_A",
    "UniProt_B",
    "Interaction_Type",
)
DETAILED_INTERACTION_FILENAME_SUFFIX = "_detailed_interactions"


_WRITE_BUFFER_BYTES = 1024 * 1024


@dataclass
class DetailedInteractionBudget:
    """Cumulative detailed-export budget shared by every structure in a run."""

    max_rows: int | None = None
    max_bytes: int | None = None
    min_free_bytes: int | None = None
    rows: int = 0
    bytes: int = 0

    def ensure_rows_available(self, additional_rows: int) -> None:
        """Fail before materializing a neighbour list that cannot fit."""
        if additional_rows < 0:
            raise ValueError("additional_rows must not be negative")
        if self.max_rows is not None and self.rows + additional_rows > self.max_rows:
            raise InputValidationError(
                "DETAILED_INTERACTION_ROW_LIMIT_EXCEEDED",
                (
                    "Detailed interaction export exceeds the configured per-run "
                    f"limit of {self.max_rows} rows."
                ),
            )

    def reserve(self, serialized_bytes: int, *, rows: int = 0) -> None:
        """Reserve exact serialized bytes and data rows before buffering them."""
        if serialized_bytes < 0 or rows < 0:
            raise ValueError("Detailed interaction reservations must not be negative")
        self.ensure_rows_available(rows)
        if self.max_bytes is not None and self.bytes + serialized_bytes > self.max_bytes:
            raise InputValidationError(
                "DETAILED_INTERACTION_BYTE_LIMIT_EXCEEDED",
                (
                    "Detailed interaction export exceeds the configured per-run "
                    f"limit of {self.max_bytes} bytes."
                ),
            )
        self.rows += rows
        self.bytes += serialized_bytes

    def ensure_storage_available(self, path: Path, pending_bytes: int) -> None:
        """Keep the configured free-space reserve before each bounded write."""
        if self.min_free_bytes is None:
            return
        try:
            free_bytes = shutil.disk_usage(path).free
        except OSError as exc:
            raise InputValidationError(
                "DETAILED_INTERACTION_STORAGE_RESERVE_LOW",
                "Free space for detailed interaction export could not be determined.",
            ) from exc
        if free_bytes - pending_bytes < self.min_free_bytes:
            raise InputValidationError(
                "DETAILED_INTERACTION_STORAGE_RESERVE_LOW",
                (
                    "Detailed interaction export would cross the configured "
                    f"free-space reserve of {self.min_free_bytes} bytes."
                ),
            )

    def as_dict(self, *, enabled: bool = True) -> dict[str, int | bool | None]:
        return {
            "enabled": enabled,
            "rows": self.rows,
            "bytes": self.bytes,
            "max_rows": self.max_rows,
            "max_bytes": self.max_bytes,
            "min_free_bytes": self.min_free_bytes,
        }


def _csv_record(values: Sequence[object]) -> bytes:
    """Serialize one stable UTF-8 CSV record with Unix newlines."""
    buffer = io.StringIO(newline="")
    csv.writer(buffer, lineterminator="\n").writerow(values)
    return buffer.getvalue().encode("utf-8")


class _AtomicBoundedCsvWriter:
    def __init__(self, output_file: Path, budget: DetailedInteractionBudget) -> None:
        self.output_file = output_file
        self.budget = budget
        self.temporary_file = output_file.with_name(
            f".{output_file.name}.{os.getpid()}.tmp"
        )
        self._handle: BinaryIO | None = None
        self._buffer = bytearray()

    def __enter__(self) -> "_AtomicBoundedCsvWriter":
        if self.output_file.exists() or self.output_file.is_symlink():
            raise FileExistsError(
                f"Detailed interaction output already exists: {self.output_file}"
            )
        try:
            self._handle = self.temporary_file.open("xb")
            self.write(DETAILED_INTERACTION_COLUMNS, rows=0)
        except Exception:
            self._abort()
            raise
        return self

    def write(self, values: Sequence[object], *, rows: int = 1) -> None:
        payload = _csv_record(values)
        self.budget.reserve(len(payload), rows=rows)
        if self._buffer and len(self._buffer) + len(payload) > _WRITE_BUFFER_BYTES:
            self._flush()
        self._buffer.extend(payload)
        if len(self._buffer) >= _WRITE_BUFFER_BYTES:
            self._flush()

    def _flush(self) -> None:
        if not self._buffer:
            return
        if self._handle is None:
            raise RuntimeError("Detailed CSV writer is not open")
        self.budget.ensure_storage_available(self.output_file.parent, len(self._buffer))
        self._handle.write(self._buffer)
        self._buffer.clear()

    def __exit__(self, exc_type, exc, traceback) -> bool:
        if exc_type is not None:
            self._abort()
            return False
        try:
            self._flush()
            if self._handle is None:
                raise RuntimeError("Detailed CSV writer is not open")
            self._handle.flush()
            os.fsync(self._handle.fileno())
            self._handle.close()
            self._handle = None
            os.replace(self.temporary_file, self.output_file)
            try:
                directory = os.open(self.output_file.parent, os.O_RDONLY)
            except OSError:
                directory = None
            if directory is not None:
                try:
                    os.fsync(directory)
                finally:
                    os.close(directory)
        except Exception:
            self._abort()
            raise
        return False

    def _abort(self) -> None:
        if self._handle is not None:
            try:
                self._handle.close()
            except OSError:
                pass
            self._handle = None
        try:
            self.temporary_file.unlink(missing_ok=True)
        except OSError:
            pass


def export_detailed_interactions(
    structure_data: Dict[str, Any],
    interactions: List[Dict[str, Any]],
    run_output_path: str,
    *,
    budget: DetailedInteractionBudget | None = None,
) -> DetailedInteractionBudget:
    """Write one bounded per-structure detailed interaction table."""
    pdb_id = structure_data["pdb_id"]
    atom_data = structure_data["atom_data"]

    from .config_loader import config

    radius = float(config["distance_thresholds"]["all_atoms_radius"])
    include_models = str(config.get("structure_model_policy", "first")).strip().lower() == "all"
    active_budget = budget or DetailedInteractionBudget()

    residues_atoms_lookup: Dict[str, Dict[str, Any]] = {}
    chains_by_unique_id: Dict[str, Dict[str, Any]] = {}
    for chain in atom_data:
        unique_id = str(chain.get("unique_chain_id") or chain["chain_id"])
        chains_by_unique_id[unique_id] = chain
        atoms: List[Dict[str, str]] = []
        atom_coords: List[List[float]] = []

        for residue in chain.get("residues", []):
            res_id = residue.get("residue_name", "?")
            res_number = residue.get("residue_number", "?")
            for atom in residue.get("atoms", []):
                atoms.append(
                    {
                        "residue": f"{res_number}:{res_id}",
                        "atom_name": atom.get("atom_name", "?"),
                    }
                )
                atom_coords.append(atom.get("coordinates"))

        coords_arr = (
            np.array(atom_coords, dtype=float)
            if atom_coords
            else np.empty((0, 3), dtype=float)
        )
        residues_atoms_lookup[unique_id] = {"atoms": atoms, "coords": coords_arr}

    output_root = Path(run_output_path)
    output_root.mkdir(parents=True, exist_ok=True)
    output_stem = portable_artifact_stem(
        f"{pdb_id}{DETAILED_INTERACTION_FILENAME_SUFFIX}"
    )
    output_file = output_root / f"{output_stem}.csv"

    with _AtomicBoundedCsvWriter(output_file, active_budget) as writer:
        for interaction in interactions:
            chain_a_raw = str(interaction.get("chain_a", ""))
            chain_b_raw = str(interaction.get("chain_b", ""))
            chain_a_data = chains_by_unique_id.get(chain_a_raw)
            chain_b_data = chains_by_unique_id.get(chain_b_raw)
            if not chain_a_data or not chain_b_data:
                continue
            chain_a_id = chain_a_raw if include_models else str(chain_a_data.get("chain_id") or "")
            chain_b_id = chain_b_raw if include_models else str(chain_b_data.get("chain_id") or "")

            data_a = residues_atoms_lookup.get(
                chain_a_raw, {"atoms": [], "coords": np.empty((0, 3))}
            )
            data_b = residues_atoms_lookup.get(
                chain_b_raw, {"atoms": [], "coords": np.empty((0, 3))}
            )
            if data_a["coords"].size == 0 or data_b["coords"].size == 0:
                continue

            uniprot_a = chain_a_data.get("uniprot_id", "UNKNOWN")
            uniprot_b = chain_b_data.get("uniprot_id", "UNKNOWN")
            interaction_type = interaction.get("interaction_type", "AA")
            tree_a = cKDTree(data_a["coords"])

            # Query one atom at a time.  Count first so a row cap fails before
            # SciPy materializes a neighbour list that cannot be exported.
            for idx_b, coordinate_b in enumerate(data_b["coords"]):
                neighbour_count = int(
                    tree_a.query_ball_point(coordinate_b, r=radius, return_length=True)
                )
                active_budget.ensure_rows_available(neighbour_count)
                indices_a = tree_a.query_ball_point(coordinate_b, r=radius)
                for idx_a in indices_a:
                    atom_a = data_a["atoms"][idx_a]
                    atom_b = data_b["atoms"][idx_b]
                    distance = float(
                        np.linalg.norm(data_a["coords"][idx_a] - coordinate_b)
                    )
                    writer.write(
                        (
                            pdb_id,
                            chain_a_id,
                            atom_a["residue"],
                            atom_a["atom_name"],
                            chain_b_id,
                            atom_b["residue"],
                            atom_b["atom_name"],
                            round(distance, 2),
                            uniprot_a,
                            uniprot_b,
                            interaction_type,
                        )
                    )

    return active_budget
