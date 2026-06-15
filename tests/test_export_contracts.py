import json
from pathlib import Path

import pandas as pd

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


def test_headless_cx2_contains_core_aspects_and_attributes(tmp_path: Path) -> None:
    nodes_df = pd.DataFrame(
        [
            {"id": "A", "name": "A", "tooltip": "node A", "color_group": "Protein"},
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

    for aspect in [
        "nodes",
        "edges",
        "nodeAttributes",
        "edgeAttributes",
        "cartesianLayout",
        "visualProperties",
    ]:
        assert aspect in aspects

    assert len(aspects["nodes"]) == 2
    assert len(aspects["edges"]) == 1
    node_attr_names = {record["n"] for record in aspects["nodeAttributes"]}
    edge_attr_names = {record["n"] for record in aspects["edgeAttributes"]}
    assert {"name", "tooltip", "color_group", "pdb_count", "source_chains", "node_kind"}.issubset(node_attr_names)
    assert {"interaction", "all_atoms_count"}.issubset(edge_attr_names)


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
