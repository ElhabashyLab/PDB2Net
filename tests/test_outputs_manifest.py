import json
import hashlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from pdb2net import __version__
from pdb2net import outputs
from pdb2net.input_contract import InputValidationError
from pdb2net.outputs import (
    OUTPUT_CONTRACT_VERSION,
    collect_web_outputs,
    create_run_output_paths,
    reserve_web_output_directory,
    write_failed_run_manifest,
    write_run_manifest,
    write_run_summary,
)
from pdb2net.structure_identity import identity_from_official_id


def _literal_cx2() -> str:
    return json.dumps(
        [
            {"CXVersion": "2.0", "hasFragments": False},
            {
                "metaData": [
                    {"name": "attributeDeclarations", "elementCount": 1},
                    {"name": "networkAttributes", "elementCount": 1},
                    {"name": "nodes", "elementCount": 2},
                    {"name": "edges", "elementCount": 1},
                    {"name": "visualProperties", "elementCount": 1},
                ]
            },
            {
                "attributeDeclarations": [
                    {
                        "networkAttributes": {"name": {"d": "string"}},
                        "nodes": {"name": {"d": "string"}},
                        "edges": {"interaction": {"d": "string"}},
                    }
                ]
            },
            {"networkAttributes": [{"name": "fixture"}]},
            {
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.0, "v": {"name": "A"}},
                    {"id": 1, "x": 1.0, "y": 1.0, "v": {"name": "B"}},
                ]
            },
            {"edges": [{"id": 0, "s": 0, "t": 1, "v": {"interaction": "x"}}]},
            {"visualProperties": [{}]},
            {"status": [{"error": "", "success": True}]},
        ]
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


def test_create_run_output_paths_reserves_a_unique_directory(tmp_path: Path) -> None:
    first = create_run_output_paths(str(tmp_path), timestamp="2026-01-02_03-04-05")
    Path(first.chain_dir, "first.cx2").write_text("first", encoding="utf-8")

    second = create_run_output_paths(str(tmp_path), timestamp="2026-01-02_03-04-05")

    assert Path(first.run_output_path).name == "2026-01-02_03-04-05"
    assert Path(second.run_output_path).name == "2026-01-02_03-04-05_2"
    assert first.run_output_path != second.run_output_path
    assert not Path(second.chain_dir, "first.cx2").exists()


def test_create_run_output_paths_is_safe_under_concurrency(tmp_path: Path) -> None:
    def create_path(_index: int) -> str:
        return create_run_output_paths(
            str(tmp_path), timestamp="2026-01-02_03-04-05"
        ).run_output_path

    with ThreadPoolExecutor(max_workers=8) as executor:
        paths = list(executor.map(create_path, range(8)))

    assert len(set(paths)) == 8


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
    Path(paths.combined_dir, "A.cx2").write_text(_literal_cx2(), encoding="utf-8")
    Path(paths.protein_dir, "B.cx2").write_text(_literal_cx2(), encoding="utf-8")
    Path(paths.chain_dir, "C.cx2").write_text(_literal_cx2(), encoding="utf-8")
    Path(paths.distances_dir, "D.csv").write_text("source,target\nA,B\n", encoding="utf-8")
    Path(paths.log_file).write_text("runtime", encoding="utf-8")
    identity = identity_from_official_id("1abc").as_dict()

    write_run_manifest(
        paths.manifest_file,
        input_files=["/inputs/a.pdb"],
        output_paths=paths,
        config_snapshot={
            "networks": {"chain_per_pdb": True},
            "structure_model_policy": "first",
            "export_detailed_interactions": True,
        },
        status="success",
        started_at="2026-01-02T03:04:05",
        finished_at="2026-01-02T03:04:06",
        total_time=1.25,
        annotations={"chains_total": 1},
        references={"manifest_id": "mini-ref-v1"},
        resources={"input": {"total_bytes": 123}},
        skipped_outputs=[{"name": "Combined_Network_X", "reason": "combined_graph_limit_exceeded"}],
        extra_counts={"structures": 1, "chains": 2},
        identities=[identity],
        structure_inputs=[
            {
                "file": "/inputs/a.pdb",
                    "identity": identity,
                    "format": "pdb",
                    "kind": "pdb",
                    "embedded_annotation_counts": {
                        "uniprot": 0,
                        "pfam": 0,
                        "cath": 0,
                        "scop2": 0,
                    },
            }
        ],
    )
    write_run_summary(paths)

    (tmp_path / "outputs").mkdir()
    collect_web_outputs(paths, str(tmp_path / "outputs"))

    web_root = tmp_path / "outputs"
    assert (web_root / "summary.json").exists()
    assert (web_root / "networks" / "A.cx2").read_text(encoding="utf-8") == _literal_cx2()
    assert (web_root / "networks" / "B.cx2").read_text(encoding="utf-8") == _literal_cx2()
    assert (web_root / "networks" / "C.cx2").read_text(encoding="utf-8") == _literal_cx2()
    assert (web_root / "interactions" / "D.csv").read_text(encoding="utf-8") == "source,target\nA,B\n"

    summary = json.loads((web_root / "summary.json").read_text(encoding="utf-8"))
    assert summary["output_contract_version"] == OUTPUT_CONTRACT_VERSION
    assert summary["pdb2net_version"] == __version__
    assert summary["status"] == "success"
    assert summary["started_at"] == "2026-01-02T03:04:05"
    assert summary["finished_at"] == "2026-01-02T03:04:06"
    assert summary["input_files"] == ["a.pdb"]
    assert summary["identities"] == [identity]
    assert summary["structure_inputs"][0]["file"] == "a.pdb"
    assert summary["networks"] == ["networks/A.cx2", "networks/B.cx2", "networks/C.cx2"]
    assert summary["interactions"] == ["interactions/D.csv"]
    assert summary["counts"] == {
        "networks": 3,
        "interactions": 1,
        "structures": 1,
        "chains": 2,
        "skipped_outputs": 1,
    }
    assert summary["config"]["structure_model_policy"] == "first"
    assert "open_in_cytoscape" not in summary["config"]
    assert summary["annotations"] == {"chains_total": 1}
    assert summary["references"]["manifest_id"] == "mini-ref-v1"
    assert summary["resources"]["input"]["total_bytes"] == 123
    assert summary["skipped_outputs"][0]["name"] == "Combined_Network_X"
    assert "run_summary" not in summary
    assert "internal_run_output_path" not in summary
    network_record = summary["artifacts"]["networks"][0]
    assert network_record["nodes"] == 2
    assert network_record["edges"] == 1
    assert network_record["sha256"] == hashlib.sha256(_literal_cx2().encode()).hexdigest()
    csv_record = summary["artifacts"]["interactions"][0]
    assert csv_record["rows"] == 1
    assert csv_record["columns"] == ["source", "target"]


def test_collect_web_outputs_rejects_a_nonempty_target(tmp_path: Path) -> None:
    paths = create_run_output_paths(str(tmp_path / "internal"), timestamp="2026-01-02_03-04-05")
    write_run_manifest(
        paths.manifest_file,
        input_files=[],
        output_paths=paths,
        config_snapshot={},
        status="failed",
        started_at="2026-01-02T03:04:05",
        finished_at="2026-01-02T03:04:06",
        total_time=1.0,
        errors=[{"code": "CORE_FAILED", "message": "failed"}],
    )
    write_run_summary(paths)
    web_root = tmp_path / "outputs"
    web_root.mkdir()
    sentinel = web_root / "existing.txt"
    sentinel.write_text("keep", encoding="utf-8")

    with pytest.raises(InputValidationError) as error:
        collect_web_outputs(paths, str(web_root))

    assert error.value.code == "WEB_OUTPUT_DIR_NOT_EMPTY"
    assert sentinel.read_text(encoding="utf-8") == "keep"
    assert set(web_root.iterdir()) == {sentinel}


def test_web_output_reservation_is_exclusive(tmp_path: Path) -> None:
    web_root = tmp_path / "outputs"
    web_root.mkdir()

    reserve_web_output_directory(str(web_root))

    with pytest.raises(InputValidationError) as error:
        reserve_web_output_directory(str(web_root))
    assert error.value.code == "WEB_OUTPUT_DIR_NOT_EMPTY"


def test_failed_web_summary_is_diagnostic_only_and_never_copies_partial_artifacts(
    tmp_path: Path,
) -> None:
    paths = create_run_output_paths(str(tmp_path / "internal"), timestamp="2026-01-02_03-04-05")
    Path(paths.chain_dir, "partial.cx2").write_text(_literal_cx2(), encoding="utf-8")
    write_run_manifest(
        paths.manifest_file,
        input_files=["/private/input.cif"],
        output_paths=paths,
        config_snapshot={},
        status="failed",
        started_at="2026-01-02T03:04:05",
        finished_at="2026-01-02T03:04:06",
        total_time=1.0,
        errors=[{"code": "CORE_FAILED", "message": "failed at /private/reference/db"}],
    )
    write_run_summary(paths)

    summary = collect_web_outputs(paths, str(tmp_path / "outputs"))

    assert summary["status"] == "failed"
    assert summary["networks"] == []
    assert summary["interactions"] == []
    assert summary["artifacts"] == {"networks": [], "interactions": []}
    assert list((tmp_path / "outputs" / "networks").iterdir()) == []
    assert "/private" not in json.dumps(summary)


@pytest.mark.parametrize("failure", ["second_copy", "partial_copy", "summary"])
def test_publication_failure_removes_only_its_own_artifacts(tmp_path, monkeypatch, failure) -> None:
    paths = create_run_output_paths(str(tmp_path / "internal"))
    for name in ("A.cx2", "B.cx2"):
        Path(paths.chain_dir, name).write_text(_literal_cx2(), encoding="utf-8")
    Path(paths.distances_dir, "contacts.csv").write_text("a,b\n1,2\n", encoding="utf-8")
    Path(paths.log_file).write_text("runtime", encoding="utf-8")
    identity = identity_from_official_id("1abc").as_dict()
    write_run_manifest(
        paths.manifest_file, input_files=["1abc.pdb"], output_paths=paths,
        config_snapshot={}, status="success", started_at="2026-01-01T00:00:00",
        finished_at="2026-01-01T00:00:01", total_time=1,
        references={"manifest_id": "review-fixture"}, identities=[identity],
        structure_inputs=[{"file": "1abc.pdb", "identity": identity, "format": "pdb", "kind": "pdb",
                           "embedded_annotation_counts": {key: 0 for key in ("uniprot", "pfam", "cath", "scop2")}}],
    )
    write_run_summary(paths)
    web = tmp_path / "web"
    reserve_web_output_directory(str(web))
    sentinel = web / "networks" / "keep.txt"
    sentinel.write_text("unrelated", encoding="utf-8")
    original_copy = outputs.shutil.copy2
    original_summary = outputs._write_public_summary
    calls = 0

    def copy(source, destination, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2 and failure != "summary":
            if failure == "partial_copy":
                Path(destination).write_bytes(b"incomplete")
            raise OSError("injected publication failure")
        return original_copy(source, destination, **kwargs)

    def summary(path, value):
        original_summary(path, value)
        if value["status"] == "success":
            raise OSError("injected publication failure")

    with monkeypatch.context() as patch:
        patch.setattr(outputs.shutil, "copy2", copy)
        if failure == "summary":
            patch.setattr(outputs, "_write_public_summary", summary)
        with pytest.raises(OSError, match="injected publication failure"):
            collect_web_outputs(paths, str(web), web_output_prepared=True)

    assert list((web / "networks").iterdir()) == [sentinel]
    assert list((web / "interactions").iterdir()) == []
    assert not (web / "runtime_analysis.txt").exists()
    assert not (web / "summary.json").exists()
    assert sentinel.read_text(encoding="utf-8") == "unrelated"
    assert Path(paths.chain_dir, "A.cx2").read_text(encoding="utf-8") == _literal_cx2()
