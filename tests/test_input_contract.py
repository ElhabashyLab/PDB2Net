import gzip
from pathlib import Path

import gemmi
import pytest

from pdb2net import file_parser
from pdb2net.file_parser import (
    extract_pdb_id_from_file,
    extract_pdb_id_from_filename,
    get_pdb_id,
    is_valid_file,
    parse_structure,
)
from pdb2net.pipeline import InputValidationError, discover_input_files


def _touch(path: Path) -> None:
    path.write_text("data_test\n", encoding="utf-8")


@pytest.mark.parametrize(
    "filename",
    [
        "1abc.pdb",
        "1abc.ent",
        "1abc.cif",
        "1abc.mmcif",
        "1abc.pdb.gz",
        "1abc.ent.gz",
        "1abc.cif.gz",
        "1abc.mmcif.gz",
        "1abc.PDB",
        "1abc.ENT",
        "1abc.CIF",
        "1abc.MMCIF",
        "1abc.PDB.GZ",
        "1abc.ENT.GZ",
        "1abc.CIF.GZ",
        "1abc.MMCIF.GZ",
    ],
)
def test_supported_structure_extensions_are_case_insensitive(filename: str) -> None:
    assert is_valid_file(filename)


def test_folder_input_discovers_only_supported_structure_files(tmp_path: Path) -> None:
    for name in [
        "a.pdb",
        "b.cif",
        "c.mmcif",
        "d.PDB",
        "e.ent",
        "f.cif.gz",
        "g.ent.GZ",
        "notes.txt",
        "image.png",
        "archive.gz",
    ]:
        _touch(tmp_path / name)
    nested = tmp_path / "nested"
    nested.mkdir()
    _touch(nested / "not_discovered.pdb")

    discovered = {Path(path).name for path in discover_input_files(str(tmp_path))}

    assert discovered == {"a.pdb", "b.cif", "c.mmcif", "d.PDB", "e.ent", "f.cif.gz", "g.ent.GZ"}


@pytest.mark.parametrize(
    "filename",
    [
        "1abc.cif.gz",
        "pdb1abc.ent.gz",
        "pdb_00001abc.cif.gz",
        "PDB1ABC.ENT.GZ",
        "PDB_00001ABC.MMCIF.GZ",
    ],
)
def test_archive_filename_patterns_normalize_to_valid_pdb_id(filename: str, monkeypatch) -> None:
    monkeypatch.setattr(file_parser, "VALID_PDB_IDS", {"1ABC"})

    assert extract_pdb_id_from_filename(filename) == "1ABC"


def test_archive_filename_candidate_does_not_require_legacy_fasta_membership(monkeypatch, capsys) -> None:
    monkeypatch.setattr(file_parser, "VALID_PDB_IDS", {"2XYZ"})

    assert extract_pdb_id_from_filename("pdb_00001abc.cif.gz") == "1ABC"
    assert "pdb_seqres" not in capsys.readouterr().out


@pytest.mark.parametrize(
    ("filename", "content"),
    [
        ("structure-copy.ent.gz", "HEADER    TEST                                    01-JAN-00   1ABC\nEND\n"),
        ("structure-copy.cif.gz", "data_1abc\n_entry.id 1abc\n"),
    ],
)
def test_pdb_id_is_read_from_gzip_content(
    tmp_path: Path,
    monkeypatch,
    filename: str,
    content: str,
) -> None:
    monkeypatch.setattr(file_parser, "VALID_PDB_IDS", {"1ABC"})
    archive = tmp_path / filename
    with gzip.open(archive, "wt", encoding="utf-8") as handle:
        handle.write(content)

    assert extract_pdb_id_from_file(str(archive)) == "1ABC"


@pytest.mark.parametrize("extension", ["pdb.gz", "ent.gz", "cif.gz", "mmcif.gz"])
def test_gemmi_parses_gzip_structure_files(tmp_path: Path, extension: str) -> None:
    pdb_text = (
        "HEADER    TEST                                    01-JAN-00   1ABC\n"
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n"
        "TER\n"
        "END\n"
    )
    if extension.startswith(("cif", "mmcif")):
        content = gemmi.read_pdb_string(pdb_text).make_mmcif_document().as_string().replace(
            "data_string", "data_1abc", 1
        )
    else:
        content = pdb_text

    archive = tmp_path / f"1abc.{extension}"
    with gzip.open(archive, "wt", encoding="utf-8") as handle:
        handle.write(content)

    structure = parse_structure(str(archive), "1ABC")

    assert structure is not None
    assert len(structure) == 1
    assert len(structure[0]) == 1


def test_empty_input_folder_raises_clear_error(tmp_path: Path) -> None:
    with pytest.raises(InputValidationError) as exc_info:
        discover_input_files(str(tmp_path))

    assert exc_info.value.code == "NO_VALID_INPUT_FILES"


def test_missing_input_folder_raises_clear_error(tmp_path: Path) -> None:
    missing = tmp_path / "missing"

    with pytest.raises(InputValidationError) as exc_info:
        discover_input_files(str(missing))

    assert exc_info.value.code == "INPUT_PATH_NOT_FOUND"


def test_file_input_path_raises_clear_error(tmp_path: Path) -> None:
    file_path = tmp_path / "single.pdb"
    _touch(file_path)

    with pytest.raises(InputValidationError) as exc_info:
        discover_input_files(str(file_path))

    assert exc_info.value.code == "INPUT_PATH_NOT_DIRECTORY"


def test_custom_cif_without_canonical_pdb_id_uses_filename_as_structure_id(tmp_path: Path, capsys) -> None:
    custom_cif = tmp_path / "AF_Q9BYF1_model_v4.cif"
    custom_cif.write_text("data_AF_Q9BYF1_model_v4\n", encoding="utf-8")

    structure_id = get_pdb_id(str(custom_cif))

    assert structure_id == "AF_Q9BYF1_MODEL_V4"
    captured = capsys.readouterr()
    assert "Using filename as structure ID" in captured.out
    assert "invalid PDB ID" not in captured.out
