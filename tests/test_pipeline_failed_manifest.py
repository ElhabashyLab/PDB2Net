import json
from pathlib import Path

import pytest

from pdb2net import pipeline
from pdb2net.input_contract import InputValidationError


def test_run_pipeline_writes_failed_manifest_for_missing_input_folder(tmp_path: Path, monkeypatch) -> None:
    output_root = tmp_path / "outputs"
    monkeypatch.setitem(pipeline.config, "output_path", str(output_root))
    missing = tmp_path / "missing"

    with pytest.raises(InputValidationError):
        pipeline.run_pipeline(str(missing))

    manifests = list(output_root.glob("20*/manifest.json"))
    assert len(manifests) == 1
    manifest = json.loads(manifests[0].read_text(encoding="utf-8"))
    assert manifest["status"] == "failed"
    assert manifest["errors"][0]["code"] == "INPUT_PATH_NOT_FOUND"
    assert manifest["input_files"] == []
    assert manifest["input_path"] == str(missing)


def test_run_pipeline_writes_failed_manifest_for_empty_input_folder(tmp_path: Path, monkeypatch) -> None:
    output_root = tmp_path / "outputs"
    empty_input = tmp_path / "inputs"
    empty_input.mkdir()
    monkeypatch.setitem(pipeline.config, "output_path", str(output_root))

    with pytest.raises(InputValidationError):
        pipeline.run_pipeline(str(empty_input))

    manifests = list(output_root.glob("20*/manifest.json"))
    assert len(manifests) == 1
    manifest = json.loads(manifests[0].read_text(encoding="utf-8"))
    assert manifest["status"] == "failed"
    assert manifest["errors"][0]["code"] == "NO_VALID_INPUT_FILES"
