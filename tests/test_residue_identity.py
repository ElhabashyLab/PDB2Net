"""Preserve author residue identities from independent PDB/mmCIF inputs to CSV."""

import csv

import pytest

from pdb2net import distances, pipeline
from pdb2net.config_loader import config
from pdb2net.detailed_results_exporter import (
    DETAILED_INTERACTION_COLUMNS,
    export_detailed_interactions,
)


# chain, author number, insertion code, residue name, x coordinate.
RECORDS = [
    ("A", -1, "", "GLY", 0.0),
    ("A", 0, "", "ALA", 1.0),
    ("A", 10, "", "ALA", 2.0),
    ("A", 10, "A", "ALA", 3.0),
    ("B", 1, "", "ALA", 4.0),
]


def _pdb():
    lines = ["HEADER".ljust(62) + "1ABC\n"]
    for model in [1, 2]:
        lines.append(f"MODEL     {model:4d}\n")
        for serial, (chain, number, code, residue, x) in enumerate(RECORDS, 1):
            lines.append(
                f"ATOM  {serial:5d}  CA  {residue:3s} {chain}{number:4d}{code or ' '}   "
                f"{x:8.3f}{0.:8.3f}{0.:8.3f}  1.00 20.00           C  \n"
            )
        lines.append("ENDMDL\n")
    return "".join(lines) + "END\n"


def _mmcif():
    # This text is not converted from the PDB by Gemmi or production code.
    fields = [
        "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id",
        "label_comp_id", "label_asym_id", "label_entity_id", "label_seq_id",
        "pdbx_PDB_ins_code", "Cartn_x", "Cartn_y", "Cartn_z", "occupancy",
        "B_iso_or_equiv", "auth_seq_id", "auth_comp_id", "auth_asym_id",
        "auth_atom_id", "pdbx_PDB_model_num",
    ]
    lines = ["data_1ABC\n_entry.id 1ABC\nloop_\n"]
    lines.extend(f"_atom_site.{field}\n" for field in fields)
    for model in [1, 2]:
        for index, (chain, number, code, residue, x) in enumerate(RECORDS, 1):
            entity = 1 if chain == "A" else 2
            label_number = index if chain == "A" else 1
            lines.append(
                f"ATOM {(model - 1) * 5 + index} C CA . {residue} L{chain} "
                f"{entity} {label_number} {code or '?'} {x} 0 0 1.0 20.0 "
                f"{number} {residue} {chain} CA {model}\n"
            )
    return "".join(lines)


@pytest.mark.parametrize("extension", ["pdb", "cif"])
@pytest.mark.parametrize("policy", ["first", "all"])
def test_author_residue_identity_survives_parsing_and_csv(tmp_path, monkeypatch, extension, policy):
    source = tmp_path / f"1abc.{extension}"
    source.write_text(_pdb() if extension == "pdb" else _mmcif())
    monkeypatch.setitem(config, "structure_model_policy", policy)
    monkeypatch.setitem(config, "distance_thresholds", {"ca_radius": 12.0, "all_atoms_radius": 5.0})
    monkeypatch.setitem(config, "interaction_filters", {"protein_protein_min_ca_neighbors": 1})
    structure = pipeline.process_single_file(str(source))
    for chain in structure["atom_data"]:
        chain["molecule_type"] = "Protein"
        if chain["chain_id"] == "A":
            assert [
                (res["residue_number"], res.get("insertion_code", ""))
                for res in chain["residues"]
            ] == [(-1, ""), (0, ""), (10, ""), (10, "A")]
            assert all(
                "insertion_code" not in res
                for res in chain["residues"][:3]
            )
    distances.coords_cache.clear()
    distances.tree_cache.clear()
    try:
        contacts = distances.calculate_distances_with_ckdtree([structure])
        export_detailed_interactions(structure, contacts, str(tmp_path / "csv"))
    finally:
        distances.coords_cache.clear()
        distances.tree_cache.clear()

    with (tmp_path / "csv/1ABC_detailed_interactions.csv").open(newline="") as handle:
        reader = csv.DictReader(handle)
        assert tuple(reader.fieldnames) == DETAILED_INTERACTION_COLUMNS
        rows = list(reader)
    expected = set()
    for model in ([1] if policy == "first" else [1, 2]):
        prefix = "" if policy == "first" else f"1ABC:model{model}:"
        expected.update(
            (prefix + "A", residue, "CA", prefix + "B", "1:ALA", "CA", distance)
            for residue, distance in [
                ("-1:GLY", 4.0), ("0:ALA", 3.0), ("10:ALA", 2.0), ("10A:ALA", 1.0),
            ]
        )
    assert len(rows) == len(expected)
    assert {
        (row["Chain_A"], row["Residue_A"], row["Atom_A"],
         row["Chain_B"], row["Residue_B"], row["Atom_B"], float(row["Distance"]))
        for row in rows
    } == expected
