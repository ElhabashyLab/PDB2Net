import pytest

from pdb2net.__main__ import main


def test_cli_help_includes_run_command(capsys) -> None:
    with pytest.raises(SystemExit) as exc_info:
        main(["--help"])

    assert exc_info.value.code == 0
    captured = capsys.readouterr()
    assert "run" in captured.out


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
