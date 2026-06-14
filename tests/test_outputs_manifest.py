import json
from pathlib import Path

from pdb2net.input_contract import InputValidationError
from pdb2net.outputs import create_run_output_paths, write_failed_run_manifest, write_run_manifest


def test_write_run_manifest_is_additive_and_machine_readable(tmp_path: Path) -> None:
    paths = create_run_output_paths(str(tmp_path), timestamp="2026-01-02_03-04-05")

    write_run_manifest(
        paths.manifest_file,
        input_files=["/inputs/a.pdb"],
        output_paths=paths,
        config_snapshot={"open_in_cytoscape": False},
        status="success",
        started_at="2026-01-02T03:04:05",
        finished_at="2026-01-02T03:04:06",
        total_time=1.25,
    )

    manifest = json.loads(Path(paths.manifest_file).read_text(encoding="utf-8"))
    assert manifest["status"] == "success"
    assert manifest["input_files"] == ["/inputs/a.pdb"]
    assert manifest["input_path"] is None
    assert manifest["outputs"]["runtime_analysis"] == paths.log_file
    assert manifest["config"]["open_in_cytoscape"] is False


def test_write_failed_run_manifest_records_error_code(tmp_path: Path) -> None:
    error = InputValidationError("NO_VALID_INPUT_FILES", "no valid files")

    paths = write_failed_run_manifest(
        str(tmp_path),
        input_path="/inputs",
        config_snapshot={"open_in_cytoscape": False},
        error=error,
        started_at="2026-01-02T03:04:05",
    )

    manifest = json.loads(Path(paths.manifest_file).read_text(encoding="utf-8"))
    assert manifest["status"] == "failed"
    assert manifest["input_path"] == "/inputs"
    assert manifest["errors"][0]["code"] == "NO_VALID_INPUT_FILES"
    assert manifest["errors"][0]["message"].startswith("NO_VALID_INPUT_FILES")
