from pdb2net.components import build_identity_edges, find_linked_components

import json

from pdb2net import pipeline
from pdb2net.config_loader import config
from pdb2net.structure_identity import identity_from_official_id


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
