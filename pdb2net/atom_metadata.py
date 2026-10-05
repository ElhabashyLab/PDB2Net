"""Retain occupancy information that Gemmi's numeric defaults would lose."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterator

import gemmi

from .input_contract import InputValidationError


AtomKey = tuple[int, str, int | None, str, str, str, str]


def atom_key(model_number: int, chain: str, residue, atom) -> AtomKey:
    return (model_number, chain.strip(), residue.seqid.num, residue.seqid.icode.strip(),
            residue.name, atom.name, atom.altloc.replace("\x00", "").strip())


def _pdb_number(value: str) -> int | None:
    """Read a PDB residue number, including the hybrid-36 form Gemmi accepts."""
    stripped = value.strip()
    if not stripped:
        return None
    if stripped[0] in "-+0123456789":
        return int(stripped)
    width = len(value)
    offset = -10 if stripped[0].isupper() else 16
    return int(stripped, 36) + offset * 36 ** (width - 1) + 10 ** width


def _override(token: str) -> tuple[bool, float | None]:
    number = gemmi.cif.as_number(token) if token else None
    needed = (number is None or not math.isfinite(number)
              or (number != 0 and (abs(number) < 1e-37 or abs(number) > 1e38)))
    return needed, number


def _pdb_records(text: str) -> Iterator[tuple[AtomKey, float | None]]:
    model = 1
    for line in text.splitlines():
        if line.startswith("MODEL "):
            model = int(line[10:14].strip() or "0")
        elif line.startswith(("ATOM  ", "HETATM")):
            needed, number = _override(line[54:60].strip())
            if not needed:
                continue
            # PDB permits Gemmi's two-character chain extension in columns 21–22.
            key = (model, line[20:22].strip(), _pdb_number(line[22:26]), line[26:27].strip(),
                   line[17:20].strip(), line[12:16].strip(), line[16:17].strip())
            yield key, number


def _cif_records(block: gemmi.cif.Block) -> Iterator[tuple[AtomKey, float | None]]:
    # Column views avoid copying the entire atom_site table into Python.
    names = ["occupancy", "Cartn_x", "pdbx_PDB_model_num", "pdbx_PDB_ins_code", "label_alt_id"]
    names += [prefix + field for prefix in ("auth_", "label_")
              for field in ("asym_id", "seq_id", "comp_id", "atom_id")]
    columns = {name: block.find_values("_atom_site." + name) for name in names}

    def value(name: str, index: int, fallback: str = "") -> str:
        column = columns[name]
        return gemmi.cif.as_string(column[index]) if len(column) else fallback

    for index in range(len(columns["Cartn_x"])):
        needed, occupancy = _override(value("occupancy", index))
        if not needed:
            continue
        def author(name: str) -> str:
            return value("auth_" + name, index) or value("label_" + name, index)

        number = author("seq_id") or value("label_seq_id", index)
        seqid = gemmi.SeqId(number) if number else None
        # SeqId's string constructor lowercases suffixes; the mmCIF reader
        # preserves their source case. Do not use the constructor's icode.
        suffix = number[-1] if seqid and seqid.icode.strip() else ""
        insertion = value("pdbx_PDB_ins_code", index).strip() or suffix
        key = (int(value("pdbx_PDB_model_num", index, "1") or "0"), author("asym_id").strip(),
               seqid.num if seqid else None, insertion, author("comp_id"), author("atom_id"),
               value("label_alt_id", index))
        yield key, occupancy


def occupancy_overrides(structure: gemmi.Structure, *, text: str, block, file_path: str) -> dict[AtomKey, float | None]:
    """Map only unknown/invalid or float32-extreme values by stable atom identity.

    Normal numeric occupancies stay on Gemmi atoms. Sparse overrides avoid a
    second complete Python atom table for ordinary deposited structures.
    NaN marks invalid source text; contact selection diagnoses it when relevant.
    """
    overrides: dict[AtomKey, float | None] = {}
    try:
        records = _cif_records(block) if block is not None else _pdb_records(text)
        for key, number in records:
            if key in overrides:
                raise ValueError("duplicate atom metadata")
            overrides[key] = number
        if overrides:
            seen: set[AtomKey] = set()
            for model in structure:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            key = atom_key(model.num, chain.name, residue, atom)
                            if key in overrides:
                                if key in seen:
                                    raise ValueError("ambiguous atom metadata")
                                seen.add(key)
            if seen != overrides.keys():
                raise ValueError("unmatched atom metadata")
    except (ValueError, RuntimeError, IndexError) as exc:
        raise InputValidationError(
            "AMBIGUOUS_ATOM_IDENTITY",
            f"Cannot associate atom metadata unambiguously: {Path(file_path).name}",
        ) from exc
    return overrides
