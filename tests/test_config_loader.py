from concurrent.futures import ProcessPoolExecutor
import json
import multiprocessing
import os
from pathlib import Path

import pytest

from pdb2net import config_loader


def _spawned_config_values() -> tuple[str, bool]:
    from pdb2net.config_loader import config

    return (
        str(config["structure_model_policy"]),
        bool(config["network_annotations"]["use_embedded_sifts"]),
    )


def test_load_config_does_not_create_output_path(tmp_path: Path, monkeypatch) -> None:
    cfg_dir = tmp_path / "configs"
    cfg_dir.mkdir()
    output_path = tmp_path / "not_created_by_import"
    (cfg_dir / "config.base.json").write_text(
        json.dumps(
            {
                "output_path": str(output_path),
                "input_folder_path": str(tmp_path / "inputs"),
                "open_in_cytoscape": False,
            }
        ),
        encoding="utf-8",
    )

    monkeypatch.setenv("PDB2NET_CONFIG_DIR", str(cfg_dir))

    cfg = config_loader.load_config()

    assert cfg["output_path"] == str(output_path)
    assert not output_path.exists()


def test_diamond_batch_defaults_are_safe_and_bounded(tmp_path: Path, monkeypatch) -> None:
    cfg_dir = tmp_path / "configs"
    cfg_dir.mkdir()
    (cfg_dir / "config.base.json").write_text("{}", encoding="utf-8")
    monkeypatch.setenv("PDB2NET_CONFIG_DIR", str(cfg_dir))

    diamond = config_loader.load_config()["diamond"]

    assert diamond["enabled"] is False
    assert diamond["assign_uniprot_id"] == "never"
    assert diamond["threads"] == 6
    assert diamond["iterate"] is True
    assert diamond["sensitivity"] == "sensitive"
    assert diamond["block_size"] == 1.0
    assert diamond["index_chunks"] == 4
    assert diamond["max_target_seqs"] == 50
    assert diamond["batch_max_sequences"] == 5000
    assert diamond["batch_max_fasta_bytes"] == 50 * 1024 * 1024


def test_detailed_interaction_limits_default_to_unbounded_for_standalone_use(
    tmp_path: Path, monkeypatch
) -> None:
    cfg_dir = tmp_path / "configs"
    cfg_dir.mkdir()
    (cfg_dir / "config.base.json").write_text("{}", encoding="utf-8")
    monkeypatch.setenv("PDB2NET_CONFIG_DIR", str(cfg_dir))

    limits = config_loader.load_config()["resource_limits"]

    assert limits["max_detailed_interaction_rows"] is None
    assert limits["max_detailed_interaction_bytes"] is None
    assert limits["min_output_free_bytes"] is None


def test_detailed_interaction_limit_environment_overrides_are_numeric(
    monkeypatch,
) -> None:
    monkeypatch.setenv("PDB2NET_MAX_DETAILED_INTERACTION_ROWS", "123")
    monkeypatch.setenv("PDB2NET_MAX_DETAILED_INTERACTION_BYTES", "456")
    monkeypatch.setenv("PDB2NET_MIN_OUTPUT_FREE_BYTES", "789")
    loaded: dict[str, object] = {}

    config_loader._apply_env_overrides(loaded)

    assert loaded["resource_limits"] == {
        "max_detailed_interaction_rows": 123,
        "max_detailed_interaction_bytes": 456,
        "min_output_free_bytes": 789,
    }


def test_load_config_returns_a_fresh_snapshot(monkeypatch) -> None:
    first = config_loader.load_config()
    first["distance_thresholds"]["ca_radius"] = 7.5

    second = config_loader.load_config()

    assert second is not first
    assert second["distance_thresholds"]["ca_radius"] == 12.0


def test_environment_paths_are_normalized_after_overrides(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("PDB2NET_PDB_FASTA", "~/references/pdb_seqres.txt")

    loaded = config_loader.load_config()

    assert loaded["pdb_fasta_path"] == str(tmp_path / "references" / "pdb_seqres.txt")


@pytest.mark.parametrize(
    ("name", "value"),
    [
        ("PDB2NET_DIAMOND_THREADS", "not-an-integer"),
        ("PDB2NET_OPEN_IN_CYTOSCAPE", "maybe"),
        ("PDB2NET_CA_RADIUS", "not-a-number"),
    ],
)
def test_invalid_environment_values_fail_as_configuration_errors(
    name: str, value: str, monkeypatch
) -> None:
    monkeypatch.setenv(name, value)

    with pytest.raises(config_loader.ConfigError, match=name):
        config_loader.load_config()


def test_explicit_config_is_direct_and_does_not_mutate_environment(
    tmp_path: Path, monkeypatch
) -> None:
    explicit = tmp_path / "job.json"
    explicit.write_text(
        json.dumps({"distance_thresholds": {"ca_radius": 8.0}}),
        encoding="utf-8",
    )
    monkeypatch.delenv("PDB2NET_CONFIG_FILE", raising=False)

    loaded = config_loader.load_config(explicit_file=explicit)

    assert loaded["distance_thresholds"]["ca_radius"] == 8.0
    assert "PDB2NET_CONFIG_FILE" not in os.environ


def test_explicit_config_must_exist(tmp_path: Path) -> None:
    with pytest.raises(config_loader.ConfigError, match="does not exist"):
        config_loader.load_config(explicit_file=tmp_path / "missing.json")


def test_removed_layout_backend_fails_instead_of_silent_fallback(
    tmp_path: Path,
) -> None:
    explicit = tmp_path / "old-layout.json"
    explicit.write_text(json.dumps({"layout_mode": "prefuse_headless"}), encoding="utf-8")

    with pytest.raises(config_loader.ConfigError, match="layout_mode"):
        config_loader.load_config(explicit_file=explicit)


def test_higher_priority_explicit_config_can_replace_an_old_layout_value(
    tmp_path: Path, monkeypatch
) -> None:
    cfg_dir = tmp_path / "configs"
    cfg_dir.mkdir()
    (cfg_dir / "config.base.json").write_text(
        json.dumps({"layout_mode": "prefuse_headless"}),
        encoding="utf-8",
    )
    explicit = tmp_path / "current.json"
    explicit.write_text(json.dumps({"layout_mode": "python_fast"}), encoding="utf-8")
    monkeypatch.setenv("PDB2NET_CONFIG_DIR", str(cfg_dir))

    loaded = config_loader.load_config(explicit_file=explicit)

    assert loaded["layout_mode"] == "python_fast"


def test_spawned_worker_activates_the_parent_config_snapshot() -> None:
    from pdb2net import pipeline

    snapshot = config_loader.load_config()
    snapshot["structure_model_policy"] = "all"
    snapshot["network_annotations"]["use_embedded_sifts"] = False
    context = multiprocessing.get_context("spawn")

    with ProcessPoolExecutor(
        max_workers=1,
        mp_context=context,
        initializer=pipeline._activate_parsing_worker,
        initargs=(snapshot,),
    ) as executor:
        observed = executor.submit(_spawned_config_values).result(timeout=30)

    assert observed == ("all", False)
