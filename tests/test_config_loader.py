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
