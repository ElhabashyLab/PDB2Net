from pathlib import Path

import pytest

from pdb2net import batching
from pdb2net.batching import batch_run
from pdb2net.input_contract import InputValidationError


def test_batch_run_rejects_missing_input_folder_before_processing(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setitem(batching.config, "output_path", str(tmp_path / "outputs"))
    missing = tmp_path / "missing"

    with pytest.raises(InputValidationError) as exc_info:
        batch_run(str(missing), lambda files: None)

    assert exc_info.value.code == "INPUT_PATH_NOT_FOUND"
    assert list((tmp_path / "outputs").glob("20*/manifest.json"))


def test_batch_run_rejects_file_input_before_processing(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setitem(batching.config, "output_path", str(tmp_path / "outputs"))
    file_path = tmp_path / "single.cif"
    file_path.write_text("data_single\n", encoding="utf-8")

    with pytest.raises(InputValidationError) as exc_info:
        batch_run(str(file_path), lambda files: None)

    assert exc_info.value.code == "INPUT_PATH_NOT_DIRECTORY"


def test_batch_run_rejects_folder_without_supported_structure_files(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setitem(batching.config, "output_path", str(tmp_path / "outputs"))
    (tmp_path / "notes.txt").write_text("not a structure\n", encoding="utf-8")

    with pytest.raises(InputValidationError) as exc_info:
        batch_run(str(tmp_path), lambda files: None)

    assert exc_info.value.code == "NO_VALID_INPUT_FILES"
