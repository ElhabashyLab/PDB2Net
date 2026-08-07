import json
from pathlib import Path

from pdb2net import __version__
from pdb2net.input_contract import InputValidationError
from pdb2net.outputs import (
    OUTPUT_CONTRACT_VERSION,
    collect_web_outputs,
    create_run_output_paths,
    write_failed_run_manifest,
    write_run_manifest,
    write_run_summary,
)


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
        annotations={"chains_total": 1, "by_source": {"sifts": 1}},
        references={"manifest_id": "mini-ref-v1"},
        resources={"input": {"total_bytes": 123}},
        skipped_outputs=[{"name": "Combined_Network_X", "reason": "combined_graph_limit_exceeded"}],
    )

    manifest = json.loads(Path(paths.manifest_file).read_text(encoding="utf-8"))
    assert manifest["output_contract_version"] == OUTPUT_CONTRACT_VERSION
    assert manifest["pdb2net_version"] == __version__
    assert manifest["status"] == "success"
    assert manifest["started_at"] == "2026-01-02T03:04:05"
    assert manifest["finished_at"] == "2026-01-02T03:04:06"
    assert manifest["input_files"] == ["/inputs/a.pdb"]
    assert manifest["input_path"] is None
    assert manifest["outputs"]["runtime_analysis"] == paths.log_file
    assert manifest["config"]["open_in_cytoscape"] is False
    assert manifest["counts"]["input_files"] == 1
    assert manifest["annotations"]["by_source"] == {"sifts": 1}
    assert manifest["references"]["manifest_id"] == "mini-ref-v1"
    assert manifest["resources"]["input"]["total_bytes"] == 123
    assert manifest["skipped_outputs"][0]["name"] == "Combined_Network_X"


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
    summary = json.loads(Path(paths.summary_file).read_text(encoding="utf-8"))
    assert manifest["status"] == "failed"
    assert summary["status"] == "failed"
    assert summary["output_contract_version"] == OUTPUT_CONTRACT_VERSION
    assert summary["pdb2net_version"] == __version__
    assert manifest["input_path"] == "/inputs"
    assert manifest["errors"][0]["code"] == "NO_VALID_INPUT_FILES"
    assert manifest["errors"][0]["message"].startswith("NO_VALID_INPUT_FILES")
    assert Path(paths.summary_file).exists()


def test_collect_web_outputs_creates_stable_summary_networks_and_interactions(tmp_path: Path) -> None:
    paths = create_run_output_paths(str(tmp_path / "internal"), timestamp="2026-01-02_03-04-05")
    Path(paths.combined_dir, "A.cx2").write_text("combined", encoding="utf-8")
    Path(paths.protein_dir, "B.cx2").write_text("protein", encoding="utf-8")
    Path(paths.chain_dir, "C.cx2").write_text("chain", encoding="utf-8")
    Path(paths.distances_dir, "D.csv").write_text("interactions", encoding="utf-8")
    Path(paths.log_file).write_text("runtime", encoding="utf-8")

    write_run_manifest(
        paths.manifest_file,
        input_files=["/inputs/a.pdb"],
        output_paths=paths,
        config_snapshot={"open_in_cytoscape": False},
        status="success",
        started_at="2026-01-02T03:04:05",
        finished_at="2026-01-02T03:04:06",
        total_time=1.25,
        annotations={"chains_total": 1},
        references={"manifest_id": "mini-ref-v1"},
        resources={"input": {"total_bytes": 123}},
        skipped_outputs=[{"name": "Combined_Network_X", "reason": "combined_graph_limit_exceeded"}],
    )
    write_run_summary(paths)

    collect_web_outputs(paths, str(tmp_path / "outputs"))

    web_root = tmp_path / "outputs"
    assert (web_root / "summary.json").exists()
    assert (web_root / "networks" / "A.cx2").read_text(encoding="utf-8") == "combined"
    assert (web_root / "networks" / "B.cx2").read_text(encoding="utf-8") == "protein"
    assert (web_root / "networks" / "C.cx2").read_text(encoding="utf-8") == "chain"
    assert (web_root / "interactions" / "D.csv").read_text(encoding="utf-8") == "interactions"

    summary = json.loads((web_root / "summary.json").read_text(encoding="utf-8"))
    assert summary["output_contract_version"] == OUTPUT_CONTRACT_VERSION
    assert summary["pdb2net_version"] == __version__
    assert summary["status"] == "success"
    assert summary["started_at"] == "2026-01-02T03:04:05"
    assert summary["finished_at"] == "2026-01-02T03:04:06"
    assert summary["input_files"] == ["/inputs/a.pdb"]
    assert summary["input_path"] is None
    assert summary["networks"]
    assert summary["interactions"]
    assert summary["counts"] == {"networks": 3, "interactions": 1, "skipped_outputs": 1}
    assert summary["config"]["open_in_cytoscape"] is False
    assert summary["annotations"] == {"chains_total": 1}
    assert summary["references"]["manifest_id"] == "mini-ref-v1"
    assert summary["resources"]["input"]["total_bytes"] == 123
    assert summary["skipped_outputs"][0]["name"] == "Combined_Network_X"
