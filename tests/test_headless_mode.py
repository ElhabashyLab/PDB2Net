from pathlib import Path

import pandas as pd
import pytest

from pdb2net import cytoscape, layout_engine


def test_headless_export_does_not_import_or_require_live_cytoscape(monkeypatch, tmp_path: Path) -> None:
    def fail_if_called():
        raise AssertionError("py4cytoscape should not be requested in headless mode")

    monkeypatch.setattr(cytoscape, "_get_py4cytoscape", fail_if_called)
    monkeypatch.setitem(cytoscape.config, "open_in_cytoscape", False)
    monkeypatch.setitem(cytoscape.config, "layout_mode", "python_fast")

    cytoscape.create_cytoscape_network(
        [
            {
                "chain_a": "TST1:A",
                "chain_b": "TST1:B",
                "all_atoms_count": 2,
                "interaction_type": "Protein-Protein",
            }
        ],
        "Headless_Test_Network",
        str(tmp_path),
        nodes_data=[
            {"id": "TST1:A", "name": "TST1:A", "tooltip": "A", "color_group": "Protein"},
            {"id": "TST1:B", "name": "TST1:B", "tooltip": "B", "color_group": "Protein"},
        ],
    )

    assert (tmp_path / "Headless_Test_Network.cx2").exists()


def test_live_cytoscape_is_only_requested_when_enabled(monkeypatch, tmp_path: Path) -> None:
    calls = {"requested": False}

    class FakeP4C:
        @staticmethod
        def get_network_list():
            return []

        @staticmethod
        def create_network_from_data_frames(nodes, edges, title, collection):
            return 1

        @staticmethod
        def get_visual_style_names():
            return []

        @staticmethod
        def create_visual_style(style_name, mappings, defaults):
            return None

        @staticmethod
        def set_visual_style(style_name):
            return None

    def fake_get_py4cytoscape():
        calls["requested"] = True
        return FakeP4C

    monkeypatch.setattr(cytoscape, "_get_py4cytoscape", fake_get_py4cytoscape)
    monkeypatch.setitem(cytoscape.config, "open_in_cytoscape", True)
    monkeypatch.setitem(cytoscape.config, "keep_last_n_networks", 46)
    monkeypatch.setitem(cytoscape.config, "layout_mode", "python_fast")

    nodes_df = pd.DataFrame([{"id": "A", "name": "A", "tooltip": "A", "color_group": "Protein"}])
    assert not nodes_df.empty

    cytoscape.create_cytoscape_network([], "Live_Test_Network", str(tmp_path), nodes_data=nodes_df.to_dict("records"))

    assert calls["requested"]
    assert (tmp_path / "Live_Test_Network.cx2").exists()


def test_headless_export_failure_propagates(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setitem(cytoscape.config, "open_in_cytoscape", False)
    monkeypatch.setattr(
        cytoscape,
        "export_cx2_headless",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(OSError("disk full")),
    )

    with pytest.raises(OSError, match="disk full"):
        cytoscape.create_cytoscape_network(
            [],
            "Broken_Export",
            str(tmp_path),
            nodes_data=[{"id": "A", "name": "A", "color_group": "Protein"}],
        )


@pytest.mark.parametrize(
    "layout_mode", [None, "prefuse_headless", "cytoscape_live", "", "unknown"]
)
def test_removed_or_unknown_layout_mode_does_not_silently_fallback(
    layout_mode: str | None,
) -> None:
    nodes = pd.DataFrame([{"id": "A"}])

    with pytest.raises(ValueError, match="Unsupported layout_mode"):
        layout_engine.calculate_positions(
            nodes,
            pd.DataFrame(),
            "Invalid_Layout",
            layout_mode,
        )


def test_live_creation_failure_still_writes_portable_cx2(monkeypatch, tmp_path: Path) -> None:
    class BrokenLiveCytoscape:
        @staticmethod
        def get_network_list():
            return []

        @staticmethod
        def create_network_from_data_frames(**_kwargs):
            raise RuntimeError("live Cytoscape unavailable")

    monkeypatch.setitem(cytoscape.config, "open_in_cytoscape", True)
    monkeypatch.setitem(cytoscape.config, "keep_last_n_networks", 46)
    monkeypatch.setattr(cytoscape, "_get_py4cytoscape", lambda: BrokenLiveCytoscape)

    cytoscape.create_cytoscape_network(
        [],
        "Live_Fallback",
        str(tmp_path),
        nodes_data=[{"id": "A", "name": "A", "color_group": "Protein"}],
    )

    assert (tmp_path / "Live_Fallback.cx2").is_file()
