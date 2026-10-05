"""Original occupancy and altloc metadata must survive structure parsing."""

import math

import gemmi
import pytest

from pdb2net import data_processor, file_parser, pipeline
from pdb2net.input_contract import InputValidationError


def pdb_atom(serial, chain="A", number="   1", name="CA", alt=" ", occupancy="  1.00", residue="ALA", x=0):
    return (f"ATOM  {serial:5d} {name:^4s}{alt}{residue:3s}{chain:>2s}{number}    "
            f"{x:8.3f}{0.:8.3f}{0.:8.3f}{occupancy} 20.00           C  \n")


def cif_atoms(rows, *, occupancy_column=True):
    fields = ["group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id", "label_comp_id",
              "label_asym_id", "label_entity_id", "label_seq_id", "pdbx_PDB_ins_code",
              "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "auth_seq_id", "auth_comp_id",
              "auth_asym_id", "auth_atom_id", "pdbx_PDB_model_num"]
    lines = ["data_1abc\nloop_\n"]
    lines.extend(f"_atom_site.{field}\n" for field in fields if field != "occupancy" or occupancy_column)
    for index, (chain, number, atom, alt, occupancy, x) in enumerate(rows, 1):
        values = ["ATOM", f"atom{index}", "C", atom, alt, "ALA", chain, "1", str(number), "?",
                  str(x), "0", "0", occupancy, str(number), "ALA", chain, atom, "1"]
        if not occupancy_column:
            values.pop(13)
        lines.append(" ".join(values) + "\n")
    return "".join(lines)


def raw_metadata(source):
    parsed = file_parser.parse_structure_input(str(source))
    return data_processor.process_structure({**parsed, "file_path": str(source), "pdb_id": "1ABC"})


@pytest.mark.parametrize("extension", ["pdb", "cif"])
def test_missing_zero_and_positive_occupancy_are_distinct(tmp_path, extension):
    rows = [("A", 1, "CA", "A", "?", 0), ("A", 1, "CA", "B", "0.00", 1),
            ("A", 2, "CA", ".", "0.20", 2)]
    text = cif_atoms(rows) if extension == "cif" else "".join([
        pdb_atom(1, alt="A", occupancy="      "),
        pdb_atom(2, alt="B", occupancy="  0.00", x=1),
        pdb_atom(3, number="   2", occupancy="  0.20", x=2),
    ])
    source = tmp_path / f"1abc.{extension}"
    source.write_text(text)
    structure = pipeline.process_single_file(str(source))
    atoms = [a for c in structure["atom_data"] for r in c["residues"] for a in r["atoms"]]
    assert [a["altloc"] for a in atoms] == ["A", "B", ""]
    assert atoms[0]["occupancy"] is None
    assert atoms[1]["occupancy"] == 0
    assert atoms[2]["occupancy"] == pytest.approx(0.2)
    assert [a["coordinates"] for a in atoms] == [[0., 0., 0.], [1., 0., 0.], [2., 0., 0.]]


def test_absent_cif_occupancy_column_is_unknown(tmp_path):
    source = tmp_path / "1abc.cif"
    source.write_text(cif_atoms([("A", 1, "CA", ".", "1", 0)], occupancy_column=False))
    parsed = pipeline.process_single_file(str(source))
    assert parsed["atom_data"][0]["residues"][0]["atoms"][0]["occupancy"] is None


def test_occupancy_follows_identity_when_chain_parts_merge(tmp_path, monkeypatch):
    monkeypatch.setitem(pipeline.config, "structure_model_policy", "all")
    source = tmp_path / "1abc.pdb"
    source.write_text("MODEL        7\n" + pdb_atom(1, occupancy="      ")
                      + pdb_atom(2, chain="B", occupancy="  0.00")
                      + pdb_atom(3, number="   2", occupancy="  0.50")
                      + "ENDMDL\nMODEL       11\n" + pdb_atom(1, occupancy="  1.00") + "ENDMDL\n")
    parsed = raw_metadata(source)
    observed = {(c["model_index"], c["chain_id"], r["residue_number"]): a["occupancy"]
                for c in parsed["atom_data"] for r in c["residues"] for a in r["atoms"]}
    assert observed == {(1, "A", 1): None, (1, "A", 2): 0.5, (1, "B", 1): 0., (2, "A", 1): 1.}


def test_hybrid36_residue_and_two_character_chain_keep_unknown_occupancy(tmp_path):
    source = tmp_path / "1abc.pdb"
    source.write_text(pdb_atom(1, chain="XY", number="A000", occupancy="      "))
    parsed = pipeline.process_single_file(str(source))
    chain = parsed["atom_data"][0]
    assert chain["chain_id"] == "XY"
    assert chain["residues"][0]["residue_number"] == 10000
    assert chain["residues"][0]["atoms"][0]["occupancy"] is None


def test_invalid_occupancy_is_not_misread_as_zero_or_unknown(tmp_path):
    source = tmp_path / "1abc.pdb"
    source.write_text(pdb_atom(1, occupancy="   xyz"))
    parsed = raw_metadata(source)
    assert math.isnan(parsed["atom_data"][0]["residues"][0]["atoms"][0]["occupancy"])
    with pytest.raises(InputValidationError, match="INVALID_ATOM_OCCUPANCY"):
        pipeline.process_single_file(str(source))


def test_conflicting_duplicate_atom_metadata_is_diagnostic(tmp_path):
    source = tmp_path / "1abc.pdb"
    source.write_text(pdb_atom(1, occupancy="      ") + pdb_atom(2, occupancy="  0.00"))
    with pytest.raises(InputValidationError, match="AMBIGUOUS_ATOM_IDENTITY"):
        file_parser.parse_structure_input(str(source))


@pytest.mark.parametrize("column", ["auth_asym_id", "auth_comp_id", "auth_atom_id"])
@pytest.mark.parametrize("missing", ["?", "."])
def test_null_author_fields_use_label_identity(tmp_path, column, missing):
    block = gemmi.cif.read_string(cif_atoms([("A", 1, "CA", ".", "?", 0)]))[0]
    block.find_values("_atom_site." + column)[0] = missing
    source = tmp_path / "1abc.cif"
    source.write_text(block.as_string())
    parsed = pipeline.process_single_file(str(source))
    assert parsed["atom_data"][0]["residues"][0]["atoms"][0]["occupancy"] is None


@pytest.mark.parametrize(("number", "insertion", "expected"), [("1A", "?", "A"), ("1", "' '", "")])
def test_cif_insertion_identity_preserves_case_and_normalizes_blank(tmp_path, number, insertion, expected):
    block = gemmi.cif.read_string(cif_atoms([("A", 1, "CA", ".", "?", 0)]))[0]
    block.find_values("_atom_site.auth_seq_id")[0] = number
    block.find_values("_atom_site.pdbx_PDB_ins_code")[0] = insertion
    source = tmp_path / "1abc.cif"
    source.write_text(block.as_string())
    parsed = pipeline.process_single_file(str(source))
    residue = parsed["atom_data"][0]["residues"][0]
    assert residue.get("insertion_code", "") == expected
    assert residue["atoms"][0]["occupancy"] is None
