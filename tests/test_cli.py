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
