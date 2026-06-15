import json
from pathlib import Path

import pytest

from pdb2net import __version__
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
    run_summary = json.loads(manifests[0].with_name("run_summary.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "failed"
    assert run_summary["status"] == "failed"
    assert manifest["output_contract_version"] == "1.0"
    assert manifest["pdb2net_version"] == __version__
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


def test_run_pipeline_writes_failed_web_summary_for_empty_input_folder(tmp_path: Path, monkeypatch) -> None:
    output_root = tmp_path / "outputs"
    web_root = tmp_path / "web_outputs"
    empty_input = tmp_path / "inputs"
    empty_input.mkdir()
    monkeypatch.setitem(pipeline.config, "output_path", str(output_root))

    with pytest.raises(InputValidationError):
        pipeline.run_pipeline(str(empty_input), web_output_dir=str(web_root))

    summary = json.loads((web_root / "summary.json").read_text(encoding="utf-8"))
    assert summary["status"] == "failed"
    assert summary["output_contract_version"] == "1.0"
    assert summary["pdb2net_version"] == __version__
    assert summary["input_files"] == []
    assert summary["input_path"] == str(empty_input)
    assert summary["errors"][0]["code"] == "NO_VALID_INPUT_FILES"
    assert summary["networks"] == []
    assert summary["interactions"] == []


def test_reference_preflight_reports_missing_files(tmp_path: Path, monkeypatch) -> None:
    output_root = tmp_path / "outputs"
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    (input_dir / "1abc.cif").write_text("data_1abc\n", encoding="utf-8")
    monkeypatch.setitem(pipeline.config, "output_path", str(output_root))
    monkeypatch.setitem(pipeline.config, "pdb_fasta_path", str(tmp_path / "missing_pdb_seqres.txt"))
    monkeypatch.setitem(pipeline.config, "sifts_tsv_path", str(tmp_path / "missing_sifts.tsv"))
    monkeypatch.setitem(pipeline.config, "uniprot_fasta_path", str(tmp_path / "missing_uniprot.fasta"))

    with pytest.raises(InputValidationError) as exc_info:
        pipeline.run_pipeline(str(input_dir))

    assert exc_info.value.code == "REFERENCE_FILE_MISSING"
    message = str(exc_info.value)
    assert "pdb_seqres.txt" in message
    assert "pdb_chain_uniprot.tsv" in message
    assert "uniprot_sprot.fasta" in message


def test_output_preflight_rejects_missing_output_path(monkeypatch) -> None:
    monkeypatch.setitem(pipeline.config, "output_path", "")

    with pytest.raises(InputValidationError) as exc_info:
        pipeline._validate_output_root()

    assert exc_info.value.code == "OUTPUT_PATH_MISSING"


def test_run_pipeline_fails_when_no_structure_parses(tmp_path: Path, monkeypatch) -> None:
    output_root = tmp_path / "outputs"
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    (input_dir / "bad.cif").write_text("data_bad\n", encoding="utf-8")
    monkeypatch.setitem(pipeline.config, "output_path", str(output_root))
    monkeypatch.setattr(pipeline, "_parse_input_files", lambda _files: [])

    with pytest.raises(InputValidationError) as exc_info:
        pipeline.run_pipeline(str(input_dir))

    assert exc_info.value.code == "NO_PARSEABLE_STRUCTURES"
    manifest = json.loads(next(output_root.glob("20*/manifest.json")).read_text(encoding="utf-8"))
    assert manifest["status"] == "failed"
    assert manifest["errors"][0]["code"] == "NO_PARSEABLE_STRUCTURES"


def test_blast_preflight_fails_when_fallback_needs_missing_database(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setitem(pipeline.config, "blastp_executable", "blastp")
    monkeypatch.setitem(pipeline.config, "blast_db_path", str(tmp_path / "missing_db"))
    monkeypatch.setattr(pipeline.shutil, "which", lambda _name: "/usr/bin/blastp")
    combined_data = [
        {
            "pdb_id": "TST1",
            "atom_data": [
                {
                    "chain_id": "A",
                    "unique_chain_id": "TST1:A",
                    "molecule_type": "Protein",
                    "uniprot_id": None,
                    "residues": [{"residue_name": "ALA"}],
                }
            ],
        }
    ]

    with pytest.raises(InputValidationError) as exc_info:
        pipeline._validate_blast_ready_if_needed(combined_data)

    assert exc_info.value.code == "BLAST_DATABASE_MISSING"


def test_blast_preflight_skips_when_no_fallback_needed(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setitem(pipeline.config, "blast_db_path", str(tmp_path / "missing_db"))
    combined_data = [
        {
            "pdb_id": "TST1",
            "atom_data": [
                {
                    "chain_id": "D",
                    "unique_chain_id": "TST1:D",
                    "molecule_type": "DNA",
                    "uniprot_id": None,
                    "residues": [{"residue_name": "DA"}],
                }
            ],
        }
    ]

    pipeline._validate_blast_ready_if_needed(combined_data)


def test_blast_preflight_skips_single_chain_direct_uniprot_filename(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setitem(pipeline.config, "blast_db_path", str(tmp_path / "missing_db"))
    combined_data = [
        {
            "pdb_id": "AF_Q9BYF1_F1_MODEL_V4",
            "file_path": "/inputs/AF-Q9BYF1-F1-model_v4.cif",
            "atom_data": [
                {
                    "chain_id": "A",
                    "unique_chain_id": "AF_Q9BYF1_F1_MODEL_V4:A",
                    "molecule_type": "Unknown",
                    "uniprot_id": None,
                    "residues": [{"residue_name": "ALA"}],
                }
            ],
        }
    ]

    pipeline._validate_blast_ready_if_needed(combined_data)
