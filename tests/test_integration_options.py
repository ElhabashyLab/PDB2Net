"""Default real tests are reproducible; overrides cannot change their gold data."""

import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
DATASET = ROOT / "tests/fixtures/core_reference"
ONE_CASE = "tests/integration/test_real_inputs.py::test_real_file_through_public_cli[1a8o.pdb]"


def invoke(*args, env=None):
    return subprocess.run([sys.executable, "-m", "pytest", *args], cwd=ROOT, env=env,
                          capture_output=True, text=True, timeout=60)


@pytest.mark.parametrize("options", [[], ["--run-integration"]])
def test_default_real_suite_collects_without_external_options(options):
    result = invoke("tests/integration", "--collect-only", "-q", *options)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "test_real_file_through_public_cli[1ehz.cif]" in result.stdout
    assert "test_modified_nucleotide_lengths" in result.stdout


@pytest.mark.parametrize(("options", "message"), [
    (["--integration-data", "unused"], "require --run-integration"),
    (["--run-integration", "--integration-data", "unused"], "require both"),
])
def test_incomplete_external_options_fail_clearly(options, message):
    result = invoke("tests/integration", "--collect-only", "-q", *options)
    assert result.returncode == pytest.ExitCode.USAGE_ERROR
    assert message in result.stderr


@pytest.mark.parametrize("damage", ["missing", "empty", "structure", "reference"])
@pytest.mark.integration
def test_external_snapshots_must_match_frozen_payloads(tmp_path, damage):
    data, refs = DATASET / "structures", DATASET / "references"
    if damage in {"missing", "empty"}:
        data = tmp_path / "structures"
        if damage == "empty":
            data.mkdir()
    elif damage == "structure":
        data = tmp_path / "structures"
        shutil.copytree(DATASET / "structures", data)
        with (data / "4hhb.pdb").open("ab") as stream:
            stream.write(b"REMARK changed fixture\n")
    else:
        refs = tmp_path / "references"
        shutil.copytree(DATASET / "references", refs)
        with (refs / "sifts.tsv").open("ab") as stream:
            stream.write(b"# changed reference\n")
    result = invoke(ONE_CASE, "-q", "--tb=short", "--run-integration", "--integration-data", str(data),
                    "--integration-references", str(refs))
    assert result.returncode != 0
    message = {"missing": "must be existing directories", "empty": "Integration fixture unavailable",
               "structure": "Frozen fixture checksum mismatch", "reference": "Frozen fixture checksum mismatch"}[damage]
    assert message in result.stdout + result.stderr


@pytest.mark.parametrize(("expression", "has_real_cases"), [
    ("(integration)", True), ("integration or not integration", True),
    ("not integration", False), ("not (integration)", False),
])
def test_marker_selection_preserves_explicit_fast_run(expression, has_real_cases):
    result = invoke("tests", "--collect-only", "-q", "-m", expression)
    assert result.returncode == 0, result.stdout + result.stderr
    assert ("test_real_file_through_public_cli[1ehz.cif]" in result.stdout) == has_real_cases


@pytest.mark.integration
def test_tool_setup_is_lazy_but_missing_blast_does_not_skip_real_tests(tmp_path):
    env = {k: v for k, v in os.environ.items() if not k.startswith("PDB2NET_")}
    env["PATH"] = str(tmp_path)
    quick = invoke("tests/test_data_processor.py", "-q", "-m", "not integration", env=env)
    assert quick.returncode == 0, quick.stdout + quick.stderr
    real = invoke(ONE_CASE, "-q", "--tb=short", env=env)
    assert real.returncode != 0
    assert "Integration tools missing from PATH" in real.stdout + real.stderr


@pytest.mark.integration
def test_identical_external_snapshot_is_supported():
    if not shutil.which("blastp") or not shutil.which("makeblastdb"):
        pytest.skip("explicit external real run needs BLAST tools")
    result = invoke(ONE_CASE, "-q", "--run-integration", "--integration-data", str(DATASET / "structures"),
                    "--integration-references", str(DATASET / "references"))
    assert result.returncode == 0, result.stdout + result.stderr
