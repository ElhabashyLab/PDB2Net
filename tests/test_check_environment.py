from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from scripts import check_environment


def _configure_ready_environment(tmp_path: Path, monkeypatch) -> dict[str, Path]:
    for name in list(os.environ):
        if name.startswith("PDB2NET_"):
            monkeypatch.delenv(name, raising=False)

    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    (input_dir / "mini.cif").write_text("data_mini\n", encoding="utf-8")
    output_dir = tmp_path / "outputs"
    output_dir.mkdir()
    reference_dir = tmp_path / "reference"
    reference_dir.mkdir()
    paths = {
        "pdb_fasta_path": reference_dir / "pdb_seqres.txt",
        "uniprot_fasta_path": reference_dir / "uniprot_sprot.fasta",
        "sifts_tsv_path": reference_dir / "pdb_chain_uniprot.tsv",
    }
    for path in paths.values():
        path.write_text("fixture\n", encoding="utf-8")

    blast_db = reference_dir / "blast_db"
    blast_db.mkdir()
    for suffix in (".pin", ".phr", ".psq"):
        Path(str(blast_db / "uniprot_db") + suffix).write_bytes(b"fixture")

    config_dir = tmp_path / "configs"
    config_dir.mkdir()
    (config_dir / "config.base.json").write_text(
        json.dumps(
            {
                "input_folder_path": str(input_dir),
                "output_path": str(output_dir),
                **{key: str(path) for key, path in paths.items()},
                "blast_db_path": str(blast_db),
                "blastp_executable": "blastp",
                "open_in_cytoscape": False,
                "diamond": {"enabled": False},
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setenv("PDB2NET_CONFIG_DIR", str(config_dir))
    monkeypatch.setattr(
        check_environment,
        "_import_available",
        lambda _name: (True, None),
    )
    monkeypatch.setattr(
        check_environment,
        "_find_executable",
        lambda value: f"/tools/{Path(str(value)).name}",
    )
    return {**paths, "blast_db_path": blast_db}


def test_ready_environment_succeeds_without_obsolete_java_check(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    _configure_ready_environment(tmp_path, monkeypatch)

    assert check_environment.main() == 0

    output = capsys.readouterr().out
    assert "java:" not in output
    assert "blast_db_path:" in output


def test_missing_required_reference_fails(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    paths = _configure_ready_environment(tmp_path, monkeypatch)
    paths["sifts_tsv_path"].unlink()

    assert check_environment.main() == 1

    output = capsys.readouterr().out
    assert "sifts_tsv_path:" in output
    assert "missing; SIFTS TSV" in output


def test_incomplete_blast_database_is_a_nonfatal_fallback_warning(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    paths = _configure_ready_environment(tmp_path, monkeypatch)
    Path(str(paths["blast_db_path"] / "uniprot_db") + ".psq").unlink()

    assert check_environment.main() == 0

    assert "BLAST fallback is not ready" in capsys.readouterr().out


def test_missing_blast_commands_are_nonfatal_fallback_warnings(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    _configure_ready_environment(tmp_path, monkeypatch)
    monkeypatch.setattr(
        check_environment,
        "_find_executable",
        lambda value: None if str(value) == "blastp" else "/tools/makeblastdb",
    )

    assert check_environment.main() == 0

    assert "BLAST fallback is not ready" in capsys.readouterr().out


def test_invalid_environment_override_uses_authoritative_config_loader(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    _configure_ready_environment(tmp_path, monkeypatch)
    monkeypatch.setenv("PDB2NET_OPEN_IN_CYTOSCAPE", "maybe")

    assert check_environment.main() == 1

    assert "Configuration error:" in capsys.readouterr().out


def test_unsupported_python_version_fails(tmp_path: Path, monkeypatch) -> None:
    _configure_ready_environment(tmp_path, monkeypatch)
    monkeypatch.setattr(check_environment.sys, "version_info", (3, 10, 0))

    assert check_environment.main() == 1


def test_missing_networkx_fails(tmp_path: Path, monkeypatch) -> None:
    _configure_ready_environment(tmp_path, monkeypatch)
    monkeypatch.setattr(
        check_environment,
        "_import_available",
        lambda name: (False, "broken import") if name == "networkx" else (True, None),
    )

    assert check_environment.main() == 1


@pytest.mark.skipif(os.name == "nt", reason="POSIX executable bits are not portable to Windows")
def test_non_executable_file_is_not_an_executable(tmp_path: Path) -> None:
    command = tmp_path / "blastp"
    command.write_text("not executable\n", encoding="utf-8")
    command.chmod(0o644)

    assert check_environment._find_executable(str(command)) is None


def test_non_command_file_is_not_an_executable_on_any_platform(tmp_path: Path) -> None:
    command = tmp_path / "blastp.txt"
    command.write_text("not a command\n", encoding="utf-8")

    assert check_environment._find_executable(str(command)) is None


def test_output_below_a_regular_file_fails(tmp_path: Path, monkeypatch) -> None:
    _configure_ready_environment(tmp_path, monkeypatch)
    blocking_parent = tmp_path / "not-a-directory"
    blocking_parent.write_text("file\n", encoding="utf-8")
    monkeypatch.setenv("PDB2NET_OUTPUT", str(blocking_parent / "outputs"))

    assert check_environment.main() == 1


def test_absent_cytoscape_setting_uses_runtime_live_default() -> None:
    assert check_environment._cytoscape_enabled({}) is True
    assert check_environment._cytoscape_enabled({"open_in_cytoscape": None}) is False


def test_enabled_cytoscape_requires_optional_python_extra(
    tmp_path: Path, monkeypatch
) -> None:
    _configure_ready_environment(tmp_path, monkeypatch)
    monkeypatch.setenv("PDB2NET_OPEN_IN_CYTOSCAPE", "true")
    monkeypatch.setattr(
        check_environment,
        "_import_available",
        lambda name: (False, "missing") if name == "py4cytoscape" else (True, None),
    )

    assert check_environment.main() == 1


def test_unready_enabled_diamond_is_a_nonfatal_fallback_warning(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    _configure_ready_environment(tmp_path, monkeypatch)
    monkeypatch.setenv("PDB2NET_DIAMOND_ENABLED", "true")
    monkeypatch.setenv("PDB2NET_DIAMOND_UNIREF90_DB", str(tmp_path / "missing-diamond-db"))
    monkeypatch.setattr(
        check_environment,
        "_find_executable",
        lambda value: None if str(value) == "diamond" else f"/tools/{value}",
    )

    assert check_environment.main() == 0

    assert "DIAMOND fallback is enabled but not ready" in capsys.readouterr().out


def test_existing_unwritable_output_directory_fails(
    tmp_path: Path, monkeypatch
) -> None:
    output_dir = tmp_path / "outputs"
    output_dir.mkdir()
    monkeypatch.setattr(check_environment.os, "access", lambda _path, _mode: False)

    assert check_environment._output_path_status(output_dir)[0] is False


def test_dangling_output_symlink_fails(tmp_path: Path) -> None:
    output_path = tmp_path / "outputs"
    output_path.symlink_to(tmp_path / "missing-target", target_is_directory=True)

    assert check_environment._output_path_status(output_path)[0] is False


def test_enabled_cytoscape_requires_running_or_executable_application(
    tmp_path: Path, monkeypatch
) -> None:
    _configure_ready_environment(tmp_path, monkeypatch)
    monkeypatch.setenv("PDB2NET_OPEN_IN_CYTOSCAPE", "true")
    monkeypatch.setattr(
        check_environment,
        "_find_executable",
        lambda value: None if "cytoscape" in str(value).lower() else f"/tools/{value}",
    )
    monkeypatch.setattr(check_environment, "_cytoscape_running", lambda: False)

    assert check_environment.main() == 1


def test_missing_makeblastdb_does_not_mark_fallback_unready(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    _configure_ready_environment(tmp_path, monkeypatch)
    monkeypatch.setattr(
        check_environment,
        "_find_executable",
        lambda value: None if str(value) == "makeblastdb" else f"/tools/{value}",
    )

    assert check_environment.main() == 0

    output = capsys.readouterr().out
    assert "BLAST fallback is not ready" not in output
    assert "BLAST database rebuild is unavailable" in output


def test_broken_module_import_is_unavailable(monkeypatch) -> None:
    monkeypatch.setattr(
        check_environment.importlib,
        "import_module",
        lambda _name: (_ for _ in ()).throw(ImportError("broken binary module")),
    )

    assert check_environment._import_available("networkx") == (
        False,
        "broken binary module",
    )


def test_unreadable_reference_file_is_not_ready(tmp_path: Path, monkeypatch) -> None:
    reference = tmp_path / "reference.tsv"
    reference.write_text("fixture\n", encoding="utf-8")
    path_open = Path.open

    def deny_reference(path, *args, **kwargs):
        if path == reference:
            raise PermissionError("denied")
        return path_open(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", deny_reference)

    assert check_environment._readable_file(reference) is False


def test_unreadable_input_directory_is_not_ready(tmp_path: Path, monkeypatch) -> None:
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    monkeypatch.setattr(
        check_environment.os,
        "scandir",
        lambda _path: (_ for _ in ()).throw(PermissionError("denied")),
    )

    assert check_environment._input_directory_ready(input_dir) is False


def test_input_directory_symlink_is_not_ready(tmp_path: Path) -> None:
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    input_link = tmp_path / "input-link"
    input_link.symlink_to(input_dir, target_is_directory=True)

    assert check_environment._input_directory_ready(input_link) is False


def test_supported_input_symlink_is_not_ready(tmp_path: Path) -> None:
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    target = tmp_path / "target.cif"
    target.write_text("data_target\n", encoding="utf-8")
    (input_dir / "link.cif").symlink_to(target)

    assert check_environment._input_directory_ready(input_dir) is False


def test_output_below_dangling_parent_symlink_fails(tmp_path: Path) -> None:
    parent = tmp_path / "dangling-parent"
    parent.symlink_to(tmp_path / "missing-parent", target_is_directory=True)

    assert check_environment._output_path_status(parent / "outputs")[0] is False


def test_enabled_cytoscape_accepts_reachable_server_without_path(
    tmp_path: Path, monkeypatch
) -> None:
    _configure_ready_environment(tmp_path, monkeypatch)
    monkeypatch.setenv("PDB2NET_OPEN_IN_CYTOSCAPE", "true")
    monkeypatch.setattr(check_environment, "_cytoscape_running", lambda: True)

    assert check_environment.main() == 0


def test_enabled_cytoscape_accepts_startable_configured_path_without_ping(
    tmp_path: Path, monkeypatch
) -> None:
    _configure_ready_environment(tmp_path, monkeypatch)
    monkeypatch.setenv("PDB2NET_OPEN_IN_CYTOSCAPE", "true")
    monkeypatch.setenv("PDB2NET_CYTO_PATH", "/apps/cytoscape")
    monkeypatch.setattr(
        check_environment,
        "_cytoscape_running",
        lambda: (_ for _ in ()).throw(AssertionError("ping must not run")),
    )

    assert check_environment.main() == 0
