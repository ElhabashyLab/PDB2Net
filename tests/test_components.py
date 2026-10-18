from pdb2net.components import build_identity_edges, find_linked_components

import json

from pdb2net import pipeline
from pdb2net.config_loader import config
from pdb2net.structure_identity import identity_from_official_id
from pdb2net.protein_network import create_protein_network

import pytest


def test_identity_bridges_include_every_chain_without_same_pdb_edges() -> None:
    chains = {
        "P12345": ["1AAA:A", "1AAA:B", "2BBB:C", "2BBB:D", "3CCC:E"],
    }
    chain_to_pdb = {
        "1AAA:A": "1AAA",
        "1AAA:B": "1AAA",
        "2BBB:C": "2BBB",
        "2BBB:D": "2BBB",
        "3CCC:E": "3CCC",
    }

    edges = build_identity_edges(chains, chain_to_pdb)

    assert len(edges) == len(chain_to_pdb) - 1
    assert all(edge["all_atoms_count"] == 0 and edge["interaction_type"] == "identity" for edge in edges)
    assert all(
        chain_to_pdb[edge["chain_a"]] != chain_to_pdb[edge["chain_b"]]
        for edge in edges
    )
    components = find_linked_components([], edges, valid_nodes=chain_to_pdb)
    assert components == [set(chain_to_pdb)]


def test_identity_bridges_require_at_least_two_distinct_pdbs() -> None:
    edges = build_identity_edges(
        {"P12345": ["1AAA:A", "1AAA:B"]},
        {"1AAA:A": "1AAA", "1AAA:B": "1AAA"},
    )

    assert edges == []


def test_reference_oracle_accepts_a_spanning_tree_and_rejects_broken_bridges(tmp_path, monkeypatch):
    from integration.test_real_inputs import assert_combined_exports

    monkeypatch.setitem(config, "open_in_cytoscape", False)
    monkeypatch.setitem(config, "layout_mode", "python_fast")
    structures = []
    expected = {"structures": {}, "networks": {}}
    for pdb_id in ("1ABC", "2XYZ"):
        filename = f"{pdb_id}.cif"
        chains = [{"chain_id": c, "unique_chain_id": f"{pdb_id}:{c}",
                   "molecule_type": "Protein", "molecule_name": "Shared protein",
                   "uniprot_id": "P12345", "aa_len": 1} for c in "AB"]
        structures.append({"pdb_id": pdb_id, "file_path": filename, "atom_data": chains})
        expected["structures"][filename] = {
            "display_id": pdb_id, "chains": {c: {"uniprot_id": "P12345"} for c in "AB"},
        }
        expected["networks"][filename] = []
    destination = tmp_path / "web/networks"
    pipeline._create_linked_identity_network([], structures, str(destination))
    create_protein_network([], structures, str(destination),
                           {"protein_per_pdb": False, "combined_protein_network": True})
    assert_combined_exports(tmp_path, list(expected["structures"]), expected)
    path, = destination.glob("Combined_Network_*.cx2")
    original = path.read_text()
    for damage in ("missing", "extra", "nonzero"):
        document = json.loads(original)
        aspects = {key: value for block in document for key, value in block.items()}
        edges = aspects["edges"]
        assert len(edges) == 3
        if damage == "missing":
            edges.pop()
        elif damage == "nonzero":
            edges[0]["v"]["all_atoms_count"] = 1
        else:
            names = {n["v"]["name"].split("\n")[0]: n["id"] for n in aspects["nodes"]}
            existing = {frozenset((e["s"], e["t"])) for e in edges}
            pair, = [(names[f"1ABC:{a}"], names[f"2XYZ:{b}"]) for a in "AB" for b in "AB"
                     if frozenset((names[f"1ABC:{a}"], names[f"2XYZ:{b}"])) not in existing]
            edges.append({**edges[0], "id": max(e["id"] for e in edges) + 1, "s": pair[0], "t": pair[1]})
        path.write_text(json.dumps(document))
        with pytest.raises(AssertionError):
            assert_combined_exports(tmp_path, list(expected["structures"]), expected)
    path.write_text(original)


def test_combined_chain_color_groups_do_not_depend_on_input_paths(tmp_path, monkeypatch) -> None:
    monkeypatch.setitem(config, "open_in_cytoscape", False)
    monkeypatch.setitem(config, "layout_mode", "python_fast")
    groups_by_run = []
    for run, filenames in enumerate((("a.cif", "z.cif"), ("z.cif", "a.cif"))):
        structures = []
        for pdb_id, filename in zip(("1ABC", "2XYZ"), filenames):
            identity = identity_from_official_id(pdb_id)
            structures.append({
                "pdb_id": pdb_id,
                "structure_identity": identity.as_dict(),
                "file_path": f"/private/jobs/{filename}",
                "atom_data": [{
                    "chain_id": "A", "unique_chain_id": f"{pdb_id}:A",
                    "molecule_type": "Protein", "molecule_name": "Shared protein",
                    "uniprot_id": "P12345", "aa_len": 1,
                }],
            })
        destination = tmp_path / str(run)
        pipeline._create_linked_identity_network([], structures, str(destination))
        path, = destination.glob("*.cx2")
        text = path.read_text(encoding="utf-8")
        document = {key: value for block in json.loads(text) for key, value in block.items()}
        groups_by_run.append({node["v"]["name"]: node["v"]["color_group"] for node in document["nodes"]})
        assert "/private/jobs" not in text
    assert groups_by_run[0] == groups_by_run[1] == {
        "1ABC:A\nP12345": "pdb_00001abc",
        "2XYZ:A\nP12345": "pdb_00002xyz",
    }
