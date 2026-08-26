import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from pdb2net import __version__
from pdb2net.__main__ import main


def test_cli_help_includes_run_command(capsys) -> None:
    with pytest.raises(SystemExit) as exc_info:
        main(["--help"])

    assert exc_info.value.code == 0
    captured = capsys.readouterr()
    assert "run" in captured.out
    assert "capabilities" in captured.out


def test_cli_version_is_available_without_a_subcommand(capsys) -> None:
    with pytest.raises(SystemExit) as exc_info:
        main(["--version"])

    assert exc_info.value.code == 0
    assert capsys.readouterr().out.strip() == f"pdb2net {__version__}"


def test_capabilities_json_reports_stable_server_contracts(capsys) -> None:
    assert main(["capabilities", "--json"]) == 0

    document = json.loads(capsys.readouterr().out)
    assert document == {
        "capabilities_schema_version": "3",
        "pdb2net_version": __version__,
        "cli_contract_version": "2",
        "output_contract_version": "2.0",
        "precomputed_artifact_schema_version": "3",
        "commands": ["run", "precompute", "assemble", "capabilities"],
        "input_formats": [
            ".pdb",
            ".cif",
            ".mmcif",
            ".pdb.gz",
            ".cif.gz",
            ".mmcif.gz",
        ],
        "network_outputs": [
            "chain_per_pdb",
            "protein_per_pdb",
            "combined_chain_network",
            "combined_protein_network",
        ],
        "structure_model_policies": ["first", "all"],
        "features": [
            "headless_cx2",
            "detailed_interactions",
            "embedded_sifts",
            "diamond_uniref90",
            "precomputed_store",
        ],
    }


def test_capabilities_are_configuration_free(tmp_path: Path) -> None:
    environment = dict(os.environ)
    environment["PDB2NET_CONFIG_FILE"] = str(tmp_path / "does-not-exist.json")
    result = subprocess.run(
        [sys.executable, "-m", "pdb2net", "capabilities", "--json"],
        cwd=Path(__file__).resolve().parents[1],
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["capabilities_schema_version"] == "3"


def test_run_help_includes_backend_output_option(capsys) -> None:
    with pytest.raises(SystemExit) as exc_info:
        main(["run", "--help"])

    assert exc_info.value.code == 0
    captured = capsys.readouterr()
    assert "--web-output-dir" in captured.out


def test_run_help_includes_optional_diamond_flags(capsys) -> None:
    with pytest.raises(SystemExit) as exc_info:
        main(["run", "--help"])

    assert exc_info.value.code == 0
    captured = capsys.readouterr()
    assert "--diamond-uniref90" in captured.out
    assert "--diamond-uniref90-db" in captured.out
    assert "--diamond-threads" in captured.out
    assert "--diamond-iterate" in captured.out
    assert "--diamond-sensitivity" in captured.out


def test_run_returns_nonzero_for_missing_input_dir(tmp_path, capsys) -> None:
    exit_code = main([
        "run",
        "--input-dir",
        str(tmp_path / "missing"),
        "--output-dir",
        str(tmp_path / "outputs"),
        "--headless",
    ])

    assert exit_code == 1
    captured = capsys.readouterr()
    assert "PDB2Net run failed" in captured.err


def test_run_returns_nonzero_and_web_summary_for_empty_input_dir(tmp_path) -> None:
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    web_root = tmp_path / "web"

    exit_code = main([
        "run",
        "--input-dir",
        str(input_dir),
        "--output-dir",
        str(tmp_path / "outputs"),
        "--web-output-dir",
        str(web_root),
        "--headless",
    ])

    assert exit_code == 1
    assert (web_root / "summary.json").exists()
