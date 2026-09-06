"""Explicit integration requests must not silently succeed without fixtures."""

from pathlib import Path
import subprocess
import sys

import pytest


@pytest.mark.parametrize(
    ("options", "message"),
    [
        ([], "require --run-integration"),
        (
            ["--run-integration"],
            "requires --integration-data and --integration-references",
        ),
        (["--integration-data", "{data}"], "require --run-integration"),
        (
            [
                "--run-integration",
                "--integration-data",
                "{data}/missing",
                "--integration-references",
                "{references}",
            ],
            "must be existing directories",
        ),
        (
            [
                "--run-integration",
                "--integration-data",
                "{data}",
                "--integration-references",
                "{references}",
            ],
            "Integration fixture unavailable",
        ),
    ],
)
def test_integration_preflight_fails_clearly(options, message, tmp_path):
    data = tmp_path / "structures"
    references = tmp_path / "references"
    data.mkdir()
    references.mkdir()
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "tests/integration",
            "--collect-only",
            "-q",
            *(option.format(data=data, references=references) for option in options),
        ],
        cwd=Path(__file__).resolve().parents[1],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == pytest.ExitCode.USAGE_ERROR
    assert message in result.stderr


@pytest.mark.parametrize(
    ("expression", "expected"),
    [
        ("(integration)", pytest.ExitCode.USAGE_ERROR),
        ("integration or not integration", pytest.ExitCode.USAGE_ERROR),
        ("not integration", pytest.ExitCode.OK),
        ("not (integration)", pytest.ExitCode.OK),
    ],
)
def test_marker_selection_respects_integration_opt_in(expression, expected):
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "tests",
            "--collect-only",
            "-q",
            "-m",
            expression,
        ],
        cwd=Path(__file__).resolve().parents[1],
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert result.returncode == expected, result.stdout + result.stderr
    if expected == pytest.ExitCode.USAGE_ERROR:
        assert "require --run-integration" in result.stderr
