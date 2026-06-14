from pathlib import Path

import pytest

from pdb2net.file_parser import get_pdb_id, is_valid_file
from pdb2net.pipeline import InputValidationError, discover_input_files


def _touch(path: Path) -> None:
    path.write_text("data_test\n", encoding="utf-8")


@pytest.mark.parametrize(
    "filename",
    [
        "1abc.pdb",
        "1abc.cif",
        "1abc.mmcif",
        "1abc.PDB",
        "1abc.CIF",
        "1abc.MMCIF",
    ],
)
def test_supported_structure_extensions_are_case_insensitive(filename: str) -> None:
    assert is_valid_file(filename)


def test_folder_input_discovers_only_supported_structure_files(tmp_path: Path) -> None:
    for name in ["a.pdb", "b.cif", "c.mmcif", "d.PDB", "notes.txt", "image.png"]:
        _touch(tmp_path / name)

    discovered = {Path(path).name for path in discover_input_files(str(tmp_path))}

    assert discovered == {"a.pdb", "b.cif", "c.mmcif", "d.PDB"}


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
