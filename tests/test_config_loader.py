import json
from pathlib import Path

from pdb2net import config_loader


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
    monkeypatch.setattr(config_loader, "_config_cache", None)

    cfg = config_loader.load_config()

    assert cfg["output_path"] == str(output_path)
    assert not output_path.exists()


def test_diamond_batch_defaults_are_safe_and_bounded(tmp_path: Path, monkeypatch) -> None:
    cfg_dir = tmp_path / "configs"
    cfg_dir.mkdir()
    (cfg_dir / "config.base.json").write_text("{}", encoding="utf-8")
    monkeypatch.setenv("PDB2NET_CONFIG_DIR", str(cfg_dir))
    monkeypatch.setattr(config_loader, "_config_cache", None)

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
    monkeypatch.setattr(config_loader, "_config_cache", None)

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
