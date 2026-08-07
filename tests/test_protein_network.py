from pdb2net import protein_network
from pdb2net.protein_network import create_protein_network


def _protein_chain(chain_id: str, uniprot_id: str) -> dict:
    return {
        "chain_id": chain_id,
        "molecule_type": "Protein",
        "molecule_name": f"Protein {uniprot_id}",
        "uniprot_id": uniprot_id,
        "residues": [{"residue_name": "ALA"}],
    }


def _dna_chain(chain_id: str) -> dict:
    return {
        "chain_id": chain_id,
        "molecule_type": "DNA",
        "molecule_name": "DNA chain",
        "uniprot_id": None,
        "residues": [{"residue_name": "DA"}, {"residue_name": "DT"}],
    }


def test_combined_protein_networks_follow_linked_chain_components(monkeypatch) -> None:
    exported = []

    def fake_create_cytoscape_network(edges, title, output_path, nodes_data=None):
        exported.append({
            "title": title,
            "edges": edges,
            "nodes": nodes_data or [],
            "output_path": output_path,
        })

    monkeypatch.setattr(protein_network, "create_cytoscape_network", fake_create_cytoscape_network)

    combined_data = [
        {"pdb_id": "PDB1", "file_path": "/inputs/pdb1.cif", "atom_data": [_protein_chain("A", "U1"), _protein_chain("B", "U2")]},
        {"pdb_id": "PDB2", "file_path": "/inputs/pdb2.cif", "atom_data": [_protein_chain("A", "U1"), _protein_chain("B", "U3")]},
        {"pdb_id": "PDB3", "file_path": "/inputs/pdb3.cif", "atom_data": [_protein_chain("A", "U4")]},
        {"pdb_id": "PDB4", "file_path": "/inputs/pdb4.cif", "atom_data": [_protein_chain("A", "U4")]},
    ]
    results = [
        {"chain_a": "PDB1:A", "chain_b": "PDB1:B", "all_atoms_count": 2},
        {"chain_a": "PDB2:A", "chain_b": "PDB2:B", "all_atoms_count": 3},
    ]

    create_protein_network(
        results,
        combined_data,
        "/tmp/unused",
        {"protein_per_pdb": False, "combined_protein_network": True},
        combined_output_path="/tmp/combined",
    )

    assert len(exported) == 2
    assert {entry["title"] for entry in exported} == {
        "Combined_Protein_Network_U1_U2_U3",
        "Combined_Protein_Network_U4",
    }

    linked = next(entry for entry in exported if entry["title"] == "Combined_Protein_Network_U1_U2_U3")
    linked_nodes = {node["id"]: node for node in linked["nodes"]}
    assert set(linked_nodes) == {"U1", "U2", "U3"}
    assert linked_nodes["U1"]["name"] == "U1\n2 PDBs"
    assert linked_nodes["U1"]["color_group"] == "Multi"
    assert linked_nodes["U1"]["pdb_count"] == 2
    assert linked_nodes["U1"]["pdb_ids"] == "PDB1, PDB2"
    assert linked_nodes["U1"]["source_chains"] == "PDB1:A, PDB2:A"
    assert linked_nodes["U1"]["source_files"] == "pdb1.cif, pdb2.cif"
    assert "Source chains: PDB1:A, PDB2:A" in linked_nodes["U1"]["tooltip"]
    assert {
        tuple(sorted((edge["chain_a"], edge["chain_b"]))) + (edge["all_atoms_count"],)
        for edge in linked["edges"]
    } == {("U1", "U2", 2), ("U1", "U3", 3)}

    nodes_only = next(entry for entry in exported if entry["title"] == "Combined_Protein_Network_U4")
    assert [node["id"] for node in nodes_only["nodes"]] == ["U4"]
    assert nodes_only["edges"] == []


def test_protein_na_nodes_include_filterable_source_metadata(monkeypatch) -> None:
    exported = []

    def fake_create_cytoscape_network(edges, title, output_path, nodes_data=None):
        exported.append({"title": title, "edges": edges, "nodes": nodes_data or []})

    monkeypatch.setattr(protein_network, "create_cytoscape_network", fake_create_cytoscape_network)

    combined_data = [
        {
            "pdb_id": "PDB1",
            "file_path": "/inputs/pdb1.cif",
            "atom_data": [_protein_chain("A", "U1"), _dna_chain("D")],
        }
    ]
    results = [{"chain_a": "PDB1:A", "chain_b": "PDB1:D", "all_atoms_count": 5}]

    create_protein_network(
        results,
        combined_data,
        "/tmp/unused",
        {"protein_per_pdb": True, "combined_protein_network": False},
        per_pdb_output_path="/tmp/protein",
    )

    network = exported[0]
    nodes = {node["id"]: node for node in network["nodes"]}

    assert nodes["U1"]["node_kind"] == "protein"
    assert nodes["PDB1:D"]["name"] == "PDB1:D"
    assert nodes["PDB1:D"]["node_kind"] == "nucleic_acid"
    assert nodes["PDB1:D"]["pdb_id"] == "PDB1"
    assert nodes["PDB1:D"]["chain_id"] == "D"
    assert nodes["PDB1:D"]["molecule_type"] == "DNA"
    assert nodes["PDB1:D"]["source_file"] == "pdb1.cif"
    assert nodes["PDB1:D"]["nt_len"] == 2


def test_combined_protein_network_over_limit_is_reported_and_not_exported(monkeypatch) -> None:
    exported = []
    monkeypatch.setattr(
        protein_network,
        "create_cytoscape_network",
        lambda edges, title, output_path, nodes_data=None: exported.append(title),
    )
    combined_data = [
        {"pdb_id": "PDB1", "atom_data": [_protein_chain("A", "U1"), _protein_chain("B", "U2")]},
        {"pdb_id": "PDB2", "atom_data": [_protein_chain("A", "U1"), _protein_chain("B", "U3")]},
    ]
    results = [
        {"chain_a": "PDB1:A", "chain_b": "PDB1:B", "all_atoms_count": 2},
        {"chain_a": "PDB2:A", "chain_b": "PDB2:B", "all_atoms_count": 3},
    ]

    skipped = create_protein_network(
        results,
        combined_data,
        "/tmp/unused",
        {"protein_per_pdb": False, "combined_protein_network": True},
        combined_graph_limits={"max_nodes": 2, "max_edges": 10},
    )

    assert exported == []
    assert skipped == [
        {
            "output_type": "network",
            "network_kind": "combined_protein_network",
            "name": "Combined_Protein_Network_U1_U2_U3",
            "reason": "combined_graph_limit_exceeded",
            "exceeded_limits": ["max_nodes"],
            "counts": {"nodes": 3, "edges": 2},
            "limits": {"max_nodes": 2, "max_edges": 10},
        }
    ]
