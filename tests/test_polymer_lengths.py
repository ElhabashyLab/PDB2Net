"""Observed polymer lengths use source roles, independently of contact selection."""

import copy

import pytest

from pdb2net import cytoscape, pipeline, protein_network
from pdb2net.residue_types import count_polymer_lengths


def residue(name, kind, *, polymer=True, number=1, insertion="", atoms=None):
    return {"residue_name": name, "polymer_type": kind, "is_polymer": polymer,
            "residue_number": number, "insertion_code": insertion, "atoms": atoms or []}


@pytest.mark.parametrize(("kind", "name", "expected"), [
    ("PeptideL", "MSE", (1, 0)), ("PeptideD", "UNK", (1, 0)),
    ("CyclicPseudoPeptide", "HYP", (1, 0)), ("Dna", "5MC", (0, 1)),
    ("Rna", "PSU", (0, 1)), ("DnaRnaHybrid", "CAR", (0, 1)),
])
def test_polymer_type_counts_modified_and_unknown_residues(kind, name, expected):
    assert count_polymer_lengths([residue(name, kind)]) == expected


def test_nonpolymer_names_and_other_polymer_types_are_not_lengths():
    records = [residue("LEU", "PeptideL", polymer=False),
               residue("DA", "Dna", polymer=False), residue("HOH", "Unknown", polymer=False),
               residue("ALA", "SaccharideD"), residue("DA", "Other")]
    assert count_polymer_lengths(records) == (0, 0)


def test_legacy_names_and_records_without_usable_types_keep_name_counting():
    names = ["ala", "MSE", "DG", "PSU", "UNK", None]
    assert count_polymer_lengths(residue_names=iter(names)) == (2, 1)
    records = [{"residue_name": name} for name in names]
    assert count_polymer_lengths(records) == (2, 1)
    for unknown in (None, "", "Unknown"):
        records = [residue(name, unknown) for name in names]
        assert count_polymer_lengths(records) == (2, 1)


def test_lengths_ignore_atom_alternatives_and_occupancy_without_mutating_input():
    records = [residue("PSU", "Rna", number=10,
                       atoms=[{"altloc": "A", "occupancy": 0.5}, {"altloc": "B", "occupancy": 0.5}]),
               residue("U", "Rna", number=10, insertion="A", atoms=[{"occupancy": 0}])]
    before = copy.deepcopy(records)
    assert count_polymer_lengths(records) == (0, 2)
    assert records == before


@pytest.mark.parametrize("compact", [False, True])
def test_chain_and_protein_exports_share_raw_and_compact_lengths(tmp_path, monkeypatch, compact):
    chains = [
        {"chain_id": "A", "unique_chain_id": "1ABC:A", "molecule_type": "Protein",
         "uniprot_id": "P12345", "residues": [residue("MSE", "PeptideL"), residue("LEU", "Unknown", polymer=False)]},
        {"chain_id": "B", "unique_chain_id": "1ABC:B", "molecule_type": "RNA", "uniprot_id": None,
         "residues": [residue("PSU", "Rna"), residue("U", "Rna", number=2), residue("HOH", "Unknown", polymer=False)]},
    ]
    structures = [{"pdb_id": "1ABC", "file_path": "1abc.cif", "atom_data": chains}]
    if compact:
        structures = pipeline._compact_structure_summaries(structures)
    nodes = cytoscape.generate_nodes_from_atom_data(structures[0]["atom_data"], "1ABC")
    assert {n["chain_id"]: (n["aa_len"], n["nt_len"]) for n in nodes} == {"A": (1, 0), "B": (0, 2)}
    exported = []
    monkeypatch.setattr(protein_network, "create_cytoscape_network",
                        lambda *_args, nodes_data=None, **_kwargs: exported.extend(nodes_data))
    protein_network.create_protein_network(
        [{"chain_a": "1ABC:A", "chain_b": "1ABC:B", "all_atoms_count": 3}], structures, str(tmp_path),
        {"protein_per_pdb": True, "combined_protein_network": False})
    by_name = {n["name"]: n for n in exported}
    assert "Length: 1 aa" in by_name["P12345"]["tooltip"]
    assert by_name["1ABC:B"]["nt_len"] == 2
