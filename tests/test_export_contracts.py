import json
import math
from pathlib import Path

import pandas as pd
import pytest

from pdb2net.cytoscape import generate_nodes_from_atom_data
from pdb2net.cx2_export import export_cx2_headless
from pdb2net.detailed_results_exporter import DETAILED_INTERACTION_COLUMNS, export_detailed_interactions
from pdb2net.visual_style import get_network_visual_profile


def _aspect_map(cx: list[dict]) -> dict[str, list]:
    aspects: dict[str, list] = {}
    for block in cx:
        if not isinstance(block, dict):
            continue
        for key, value in block.items():
            if key in {"CXVersion", "hasFragments"}:
                continue
            aspects[key] = value if isinstance(value, list) else [value]
    return aspects


def test_detailed_interactions_csv_columns_are_stable(tmp_path: Path) -> None:
    structure_data = {
        "pdb_id": "TST1",
        "atom_data": [
            {
                "chain_id": "A",
                "uniprot_id": "PAAAAA",
                "residues": [
                    {
                        "residue_name": "ALA",
                        "residue_number": 1,
                        "atoms": [{"atom_name": "CA", "coordinates": [0.0, 0.0, 0.0]}],
                    }
                ],
            },
            {
                "chain_id": "B",
                "uniprot_id": "PBBBBB",
                "residues": [
                    {
                        "residue_name": "GLY",
                        "residue_number": 2,
                        "atoms": [{"atom_name": "CA", "coordinates": [1.0, 0.0, 0.0]}],
                    }
                ],
            },
        ],
    }
    interactions = [
        {
            "chain_a": "TST1:A",
            "chain_b": "TST1:B",
            "all_atoms_count": 1,
            "interaction_type": "Protein-Protein",
        }
    ]

    export_detailed_interactions(structure_data, interactions, str(tmp_path))

    out_file = tmp_path / "TST1_detailed_interactions.csv"
    df = pd.read_csv(out_file)
    assert list(df.columns) == DETAILED_INTERACTION_COLUMNS


def test_headless_cx2_uses_only_native_inline_attributes_and_layout(tmp_path: Path) -> None:
    nodes_df = pd.DataFrame(
        [
            {
                "id": "A",
                "name": "A",
                "tooltip": "node A",
                "color_group": "Protein",
                "pdb_count": "",
                "source_chains": "",
                "node_kind": "",
            },
            {
                "id": "B",
                "name": "B",
                "tooltip": "node B",
                "color_group": "Protein",
                "pdb_count": 2,
                "source_chains": "PDB1:B, PDB2:B",
                "node_kind": "protein",
            },
        ]
    )
    edges_df = pd.DataFrame(
        [
            {
                "source": "A",
                "target": "B",
                "interaction": "Protein-Protein",
                "all_atoms_count": 3,
            }
        ]
    )
    positions = {"A": {"x": 1.0, "y": 2.0}, "B": {"x": 3.0, "y": 4.0}}

    export_cx2_headless(
        "Mini_Network",
        str(tmp_path),
        nodes_df,
        edges_df,
        {"Protein": "#1f77b4"},
        positions,
    )

    cx = json.loads((tmp_path / "Mini_Network.cx2").read_text(encoding="utf-8"))
    aspects = _aspect_map(cx)

    for aspect in ["nodes", "edges", "visualProperties"]:
        assert aspect in aspects
    assert {"nodeAttributes", "edgeAttributes", "cartesianLayout"}.isdisjoint(aspects)

    assert len(aspects["nodes"]) == 2
    assert len(aspects["edges"]) == 1
    assert [node["id"] for node in aspects["nodes"]] == [0, 1]
    assert aspects["nodes"][0]["v"]["name"] == "A"
    assert aspects["nodes"][0]["x"] == 1.0
    assert aspects["nodes"][0]["y"] == 2.0
    assert aspects["nodes"][1]["v"]["pdb_count"] == 2
    assert aspects["nodes"][1]["v"]["source_chains"] == "PDB1:B, PDB2:B"
    assert aspects["nodes"][1]["v"]["node_kind"] == "protein"
    assert aspects["edges"] == [
        {
            "id": 0,
            "s": 0,
            "t": 1,
            "v": {"interaction": "Protein-Protein", "all_atoms_count": 3},
        }
    ]
    metadata = {entry["name"]: entry["elementCount"] for entry in aspects["metaData"]}
    assert metadata["nodes"] == 2
    assert metadata["edges"] == 1
    assert aspects["status"] == [{"error": "", "success": True}]


def test_headless_cx2_is_deterministic_and_sorts_nodes_and_edges(tmp_path: Path) -> None:
    nodes = pd.DataFrame(
        [
            {"id": "B", "name": "B", "tooltip": "β", "color_group": "Protein"},
            {"id": "A", "name": "A", "tooltip": "α", "color_group": "Protein"},
        ]
    )
    edges = pd.DataFrame(
        [{"source": "B", "target": "A", "interaction": "interacts_with", "all_atoms_count": 1}]
    )
    positions = {"A": {"x": 0, "y": 1}, "B": {"x": 2, "y": 3}}
    export_cx2_headless("stable", str(tmp_path), nodes, edges, {}, positions)
    first = (tmp_path / "stable.cx2").read_bytes()

    export_cx2_headless(
        "stable",
        str(tmp_path),
        nodes.iloc[::-1],
        edges.iloc[::-1],
        {},
        positions,
    )

    assert (tmp_path / "stable.cx2").read_bytes() == first


@pytest.mark.parametrize("bad_value", [math.nan, math.inf, -math.inf])
def test_headless_cx2_rejects_nonfinite_coordinates_and_attributes(
    tmp_path: Path, bad_value: float
) -> None:
    nodes = pd.DataFrame(
        [{"id": "A", "name": "A", "tooltip": "", "color_group": "Protein"}]
    )
    with pytest.raises(ValueError, match="finite"):
        export_cx2_headless(
            "bad-coordinate",
            str(tmp_path),
            nodes,
            pd.DataFrame(columns=["source", "target", "interaction", "all_atoms_count"]),
            {},
            {"A": {"x": bad_value, "y": 0.0}},
        )

    nodes["score"] = [bad_value]
    with pytest.raises(ValueError, match="finite"):
        export_cx2_headless(
            "bad-attribute",
            str(tmp_path),
            nodes,
            pd.DataFrame(columns=["source", "target", "interaction", "all_atoms_count"]),
            {},
            {"A": {"x": 0.0, "y": 0.0}},
        )


def test_headless_cx2_rejects_duplicate_nodes_and_edges(tmp_path: Path) -> None:
    duplicate_nodes = pd.DataFrame(
        [
            {"id": "A", "name": "A", "tooltip": "", "color_group": "Protein"},
            {"id": "A", "name": "A2", "tooltip": "", "color_group": "Protein"},
        ]
    )
    with pytest.raises(ValueError, match="node IDs"):
        export_cx2_headless("duplicate-nodes", str(tmp_path), duplicate_nodes, pd.DataFrame(), {}, {})

    nodes = pd.DataFrame(
        [
            {"id": "A", "name": "A", "tooltip": "", "color_group": "Protein"},
            {"id": "B", "name": "B", "tooltip": "", "color_group": "Protein"},
        ]
    )
    duplicate_edges = pd.DataFrame(
        [
            {"source": "A", "target": "B", "interaction": "x", "all_atoms_count": 1},
            {"source": "B", "target": "A", "interaction": "x", "all_atoms_count": 2},
        ]
    )
    with pytest.raises(ValueError, match="duplicate"):
        export_cx2_headless("duplicate-edges", str(tmp_path), nodes, duplicate_edges, {}, {})


def test_headless_cx2_omits_missing_optional_attributes_and_rejects_mixed_types(
    tmp_path: Path,
) -> None:
    nodes = pd.DataFrame(
        [
            {"id": "A", "name": "A", "tooltip": "", "color_group": "Protein", "score": 1},
            {"id": "B", "name": "B", "tooltip": "", "color_group": "Protein", "score": ""},
        ]
    )
    export_cx2_headless(
        "optional",
        str(tmp_path),
        nodes,
        pd.DataFrame(columns=["source", "target", "interaction", "all_atoms_count"]),
        {},
        {"A": {"x": 0.0, "y": 0.0}, "B": {"x": 1.0, "y": 1.0}},
    )
    aspects = _aspect_map(json.loads((tmp_path / "optional.cx2").read_text(encoding="utf-8")))
    assert aspects["nodes"][0]["v"]["score"] == 1
    assert "score" not in aspects["nodes"][1]["v"]

    nodes.loc[1, "score"] = "high"
    with pytest.raises(ValueError, match="incompatible scalar types"):
        export_cx2_headless(
            "mixed",
            str(tmp_path),
            nodes,
            pd.DataFrame(columns=["source", "target", "interaction", "all_atoms_count"]),
            {},
            {"A": {"x": 0.0, "y": 0.0}, "B": {"x": 1.0, "y": 1.0}},
        )


def test_component_combined_protein_titles_use_combined_protein_profile() -> None:
    profile = get_network_visual_profile("Combined_Protein_Network_P0DTC2_Q9BYF1")

    assert profile["is_combined_protein_network"] is True
    assert profile["is_combined_network"] is True


def test_chain_nodes_include_filterable_source_metadata() -> None:
    nodes = generate_nodes_from_atom_data(
        [
            {
                "unique_chain_id": "6M17:A",
                "chain_id": "A",
                "_parent_pdb_id": "6M17",
                "_parent_file_label": "6m17.cif",
                "molecule_type": "Protein",
                "molecule_name": "Spike protein",
                "uniprot_id": "P0DTC2",
                "residues": [{"residue_name": "ALA"}, {"residue_name": "DA"}],
            }
        ],
        pdb_id="6M17",
    )

    node = nodes[0]
    assert node["node_kind"] == "chain"
    assert node["pdb_id"] == "6M17"
    assert node["chain_id"] == "A"
    assert node["source_file"] == "6m17.cif"
    assert node["molecule_type"] == "Protein"
    assert node["uniprot_id"] == "P0DTC2"
    assert node["aa_len"] == 1
    assert node["nt_len"] == 1


def test_chain_nodes_include_diamond_annotation_metadata() -> None:
    nodes = generate_nodes_from_atom_data(
        [
            {
                "unique_chain_id": "TST1:A",
                "chain_id": "A",
                "_parent_pdb_id": "TST1",
                "_parent_file_label": "model.cif",
                "molecule_type": "Protein",
                "molecule_name": "Example oxidase",
                "uniprot_id": None,
                "annotation_source": "diamond_uniref90",
                "matched_database": "UniRef90",
                "matched_id": "UniRef90_P12345",
                "representative_accession": "P12345",
                "annotation_confidence": "high",
                "residues": [{"residue_name": "ALA"}],
            }
        ],
        pdb_id="TST1",
    )

    node = nodes[0]
    assert node["annotation_source"] == "diamond_uniref90"
    assert node["matched_database"] == "UniRef90"
    assert node["matched_id"] == "UniRef90_P12345"
    assert node["representative_accession"] == "P12345"
    assert node["annotation_confidence"] == "high"
    assert "Annotation: diamond_uniref90 / UniRef90 (high)" in node["tooltip"]
    assert "Matched ID: UniRef90_P12345" in node["tooltip"]
    assert "Representative accession: P12345" in node["tooltip"]
