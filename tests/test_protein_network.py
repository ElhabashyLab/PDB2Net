from pdb2net import protein_network
from pdb2net.protein_network import create_protein_network
from pdb2net import pipeline
from pdb2net.config_loader import config
from pdb2net.outputs import _cx2_counts
from pdb2net.structure_identity import identity_from_official_id

import json
import ntpath
import posixpath

import pytest


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


@pytest.mark.parametrize("reverse", [False, True])
def test_four_export_views_preserve_their_declared_scientific_scope(tmp_path, monkeypatch, reverse):
    """Hand-counted aggregation after annotation, through actual validated CX2."""
    monkeypatch.setitem(config, "open_in_cytoscape", False)
    monkeypatch.setitem(config, "layout_mode", "python_fast")
    monkeypatch.setitem(config, "structure_model_policy", "first")
    rna = {**_dna_chain("E"), "molecule_type": "RNA"}
    structures = []
    for pdb_id, chains in [
        ("1ABC", [_protein_chain("A", "P12345"), _protein_chain("B", "P12345"),
                  _protein_chain("C", "Q99999"), _dna_chain("D"), rna, _protein_chain("U", None)]),
        ("2XYZ", [_protein_chain("A", "P12345"), _protein_chain("B", "Q99999")]),
        ("3ABC", [_protein_chain("A", "O12345")]),
    ]:
        for c in chains:
            c["unique_chain_id"] = f"{pdb_id}:{c['chain_id']}"
        structures.append({"pdb_id": pdb_id, "file_path": f"{pdb_id}.pdb", "atom_data": chains,
                           "structure_identity": identity_from_official_id(pdb_id).as_dict()})
    edges = [{"chain_a": a, "chain_b": b, "all_atoms_count": count, "interaction_type": kind}
             for a, b, count, kind in [
                 ("1ABC:A", "1ABC:B", 20, "Protein-Protein"),
                 ("1ABC:A", "1ABC:C", 5, "Protein-Protein"),
                 ("1ABC:B", "1ABC:C", 7, "Protein-Protein"),
                 ("1ABC:B", "1ABC:D", 4, "Protein-DNA"),
                 ("1ABC:B", "1ABC:U", 2, "Protein-Protein"),
                 ("1ABC:D", "1ABC:E", 3, "DNA-RNA"),
                 ("2XYZ:A", "2XYZ:B", 9, "Protein-Protein"),
             ]]
    if reverse:
        structures.reverse()
        edges.reverse()
    pipeline._export_chain_networks(structures, edges, str(tmp_path / "chain"))
    pipeline._create_linked_identity_network(edges, structures, str(tmp_path / "linked"))
    create_protein_network(edges, structures, str(tmp_path / "protein"),
                           {"protein_per_pdb": True, "combined_protein_network": True})
    exported = {}
    for path in tmp_path.rglob("*.cx2"):
        _cx2_counts(path)
        aspects = {key: value for block in json.loads(path.read_text()) for key, value in block.items()}
        names = {n["id"]: n["v"]["name"].split("\n")[0] for n in aspects["nodes"]}
        pairs = {(*sorted((names[e["s"]], names[e["t"]])), e["v"]["interaction"]): e["v"]["all_atoms_count"]
                 for e in aspects["edges"]}
        exported[path.stem] = set(names.values()), pairs
    chain_nodes, chain_edges = exported["Chain_Interaction_Network_1ABC"]
    assert chain_nodes == {f"1ABC:{c}" for c in "ABCDEU"}
    assert len(chain_edges) == 6
    assert chain_edges["1ABC:A", "1ABC:B", "Protein-Protein"] == 20
    assert chain_edges["1ABC:B", "1ABC:U", "Protein-Protein"] == 2
    assert chain_edges["1ABC:D", "1ABC:E", "DNA-RNA"] == 3
    assert exported["Protein_Network_1ABC"] == (
        {"P12345", "Q99999", "1ABC:D"},
        {("P12345", "Q99999", "interacts_with"): 12, ("1ABC:D", "P12345", "interacts_with"): 4},
    )
    assert exported["Combined_Protein_Network_P12345_Q99999"] == (
        {"P12345", "Q99999"}, {("P12345", "Q99999", "interacts_with"): 21},
    )
    linked_nodes, linked_edges = exported["Combined_Network_P12345_Q99999"]
    assert linked_nodes == chain_nodes | {"2XYZ:A", "2XYZ:B"}
    assert {key: value for key, value in linked_edges.items() if key[2] != "identity"} == (
        chain_edges | {("2XYZ:A", "2XYZ:B", "Protein-Protein"): 9}
    )
    assert {key: value for key, value in linked_edges.items() if key[2] == "identity"} == {
        ("1ABC:A", "2XYZ:A", "identity"): 0,
        ("1ABC:B", "2XYZ:A", "identity"): 0,
        ("1ABC:C", "2XYZ:B", "identity"): 0,
    }
    assert exported["Chain_Interaction_Network_3ABC"] == ({"3ABC:A"}, {})


@pytest.mark.parametrize("reverse", [False, True])
def test_per_pdb_protein_metadata_does_not_leak_other_inputs(monkeypatch, reverse) -> None:
    exported = {}
    monkeypatch.setattr(
        protein_network,
        "create_cytoscape_network",
        lambda edges, title, output_path, nodes_data=None: exported.update({title: nodes_data}),
    )
    structures = [
        {"pdb_id": pdb_id, "file_path": f"/inputs/{pdb_id}.cif", "atom_data": [
            {**_protein_chain("A", "P12345"), "molecule_name": f"Name from {pdb_id}"}
        ]}
        for pdb_id in ("ONE1", "TWO2")
    ]
    create_protein_network(
        [], structures[::-1] if reverse else structures, "/tmp/unused",
        {"protein_per_pdb": True, "combined_protein_network": True},
    )
    for pdb_id in ("ONE1", "TWO2"):
        node, = exported[f"Protein_Network_{pdb_id}"]
        assert node["pdb_count"] == 1
        assert node["pdb_ids"] == pdb_id
        assert node["source_files"] == f"{pdb_id}.cif"
        assert node["source_chains"] == f"{pdb_id}:A"
        assert node["molecule_name"] == f"Name from {pdb_id}"
        assert node["tooltip"].splitlines()[0] == f"Name from {pdb_id}"
    combined, = exported["Combined_Protein_Network_P12345"]
    assert combined["pdb_count"] == 2
    assert combined["pdb_ids"] == "ONE1, TWO2"
    assert combined["source_files"] == "ONE1.cif, TWO2.cif"
    assert combined["source_chains"] == "ONE1:A, TWO2:A"


@pytest.mark.parametrize(("path_module", "paths", "labels"), [
    (ntpath, (r"C:\inputs\one.cif", r"D:\inputs\two.cif"), ("one.cif", "two.cif")),
    (ntpath, (r"\\server\share\one.cif", r"\\server\share\two.cif"), ("one.cif", "two.cif")),
    (ntpath, (r"C:one.cif", r"inputs\two.cif"), ("one.cif", "two.cif")),
    (posixpath, ("/inputs/one.cif", "inputs/two.cif"), ("one.cif", "two.cif")),
    (posixpath, (r"/inputs/part\one.cif", r"/inputs/part\two.cif"),
     (r"part\one.cif", r"part\two.cif")),
    (ntpath, ("pdb:pdb_00001abc", "pdb:pdb_00002xyz"), ("pdb:pdb_00001abc", "pdb:pdb_00002xyz")),
    (posixpath, ("", None), ("1ABC", "2XYZ")),
])
def test_protein_cx2_source_labels_follow_native_path_rules(tmp_path, monkeypatch, path_module, paths, labels):
    # Exercise the actual stdlib Windows/POSIX implementations on any test host,
    # without changing os.path globally or requiring a Windows filesystem.
    monkeypatch.setattr(protein_network, "basename", path_module.basename, raising=False)
    monkeypatch.setitem(config, "open_in_cytoscape", False)
    monkeypatch.setitem(config, "layout_mode", "python_fast")
    structures = [
        {"pdb_id": "1ABC", "file_path": paths[0], "atom_data": [_protein_chain("A", "P12345"), _dna_chain("D")]},
        {"pdb_id": "2XYZ", "file_path": paths[1], "atom_data": [_protein_chain("A", "P12345")]},
    ]
    create_protein_network(
        [{"chain_a": "1ABC:A", "chain_b": "1ABC:D", "all_atoms_count": 5}],
        structures, str(tmp_path), {"protein_per_pdb": True, "combined_protein_network": True},
    )
    exported = {}
    for path in tmp_path.glob("*.cx2"):
        _cx2_counts(path)
        aspects = {key: value for block in json.loads(path.read_text()) for key, value in block.items()}
        exported[path.stem] = aspects
    assert len(exported) == 3
    for pdb_id, label in zip(("1ABC", "2XYZ"), labels):
        nodes = {n["v"]["name"]: n["v"] for n in exported[f"Protein_Network_{pdb_id}"]["nodes"]}
        assert nodes["P12345"]["source_files"] == label
        if pdb_id == "1ABC":
            assert nodes["1ABC:D"]["source_file"] == label
            edge, = exported[f"Protein_Network_{pdb_id}"]["edges"]
            assert edge["v"]["all_atoms_count"] == 5
    node, = exported["Combined_Protein_Network_P12345"]["nodes"]
    expected = ", ".join(sorted(labels))
    assert node["v"]["source_files"] == expected
    assert f"Source files: {expected}" in node["v"]["tooltip"]


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
