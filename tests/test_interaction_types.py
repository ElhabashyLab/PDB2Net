import numpy as np
import pytest

from pdb2net import distances
from pdb2net.distances import calculate_distances_with_ckdtree, determine_interaction_type
from pdb2net.uniprot_matcher import classify_molecule_type


@pytest.mark.parametrize(
    ("left", "right", "expected"),
    [
        ("Protein", "Protein", "Protein-Protein"),
        ("Protein", "RNA", "Protein-RNA"),
        ("RNA", "Protein", "Protein-RNA"),
        ("Protein", "DNA", "Protein-DNA"),
        ("RNA", "RNA", "RNA-RNA"),
        ("DNA", "DNA", "DNA-DNA"),
        ("Unknown", "Protein", None),
        ("Protein", "Unknown", None),
    ],
)
def test_determine_interaction_type_contract(left: str, right: str, expected: str | None) -> None:
    assert determine_interaction_type(left, right) == expected


def test_dense_kdtree_contact_count_does_not_change_scientific_result() -> None:
    left = distances.cKDTree(np.zeros((1_000, 3), dtype=float))
    right = distances.cKDTree(np.zeros((1_000, 3), dtype=float))

    assert distances.count_nearby_atoms(left, right, 1.0) == 1_000_000


@pytest.mark.parametrize(
    ("residues", "expected"),
    [
        (["ALA", "GLY", "LYS"], "Protein"),
        (["A", "U", "G"], "RNA"),
        (["DA", "DT", "DG"], "DNA"),
        (["DA", "U"], "DNA/RNA"),
        (["UNK", "LIG"], "Unknown"),
    ],
)
def test_classify_molecule_type_from_residue_composition(residues: list[str], expected: str) -> None:
    chain = {"residues": [{"residue_name": name} for name in residues]}
    assert classify_molecule_type(chain) == expected


def _protein_chain(unique_id: str, offset: float = 0.0) -> dict:
    residues = []
    for index in range(10):
        residues.append({
            "residue_name": "ALA",
            "atoms": [{"atom_name": "CA", "coordinates": [float(index * 10) + offset, 0.0, 0.0]}],
        })
    return {
        "chain_id": unique_id.split(":", 1)[1],
        "unique_chain_id": unique_id,
        "molecule_type": "Protein",
        "residues": residues,
    }


def _dna_chain(unique_id: str, offset: float = 0.0) -> dict:
    return {
        "chain_id": unique_id.split(":", 1)[1],
        "unique_chain_id": unique_id,
        "molecule_type": "DNA",
        "residues": [
            {"residue_name": "DA", "atoms": [{"atom_name": "P", "coordinates": [offset, 0.0, 0.0]}]},
            {"residue_name": "DT", "atoms": [{"atom_name": "P", "coordinates": [offset + 100.0, 0.0, 0.0]}]},
        ],
    }


def test_protein_contact_filters_are_configurable(monkeypatch) -> None:
    distances.coords_cache.clear()
    distances.tree_cache.clear()
    monkeypatch.setitem(distances.config, "distance_thresholds", {"ca_radius": 1.1, "all_atoms_radius": 1.1})
    monkeypatch.setitem(
        distances.config,
        "interaction_filters",
        {
            "protein_protein_min_ca_neighbors": 11,
            "protein_protein_min_all_atom_contacts": 1,
            "protein_nucleic_acid_min_all_atom_contacts": 1,
            "nucleic_acid_min_all_atom_contacts": 1,
        },
    )

    no_edges = calculate_distances_with_ckdtree([
        {"file_path": "/tmp/tst.cif", "atom_data": [_protein_chain("TST:A"), _protein_chain("TST:B")]}
    ])

    assert no_edges == []

    monkeypatch.setitem(
        distances.config,
        "interaction_filters",
        {
            "protein_protein_min_ca_neighbors": 10,
            "protein_protein_min_all_atom_contacts": 1,
            "protein_nucleic_acid_min_all_atom_contacts": 1,
            "nucleic_acid_min_all_atom_contacts": 1,
        },
    )
    distances.coords_cache.clear()
    distances.tree_cache.clear()

    edges = calculate_distances_with_ckdtree([
        {"file_path": "/tmp/tst.cif", "atom_data": [_protein_chain("TST:A"), _protein_chain("TST:B")]}
    ])

    assert len(edges) == 1
    assert edges[0]["interaction_type"] == "Protein-Protein"


def test_nucleic_acid_contact_filters_are_configurable(monkeypatch) -> None:
    distances.coords_cache.clear()
    distances.tree_cache.clear()
    monkeypatch.setitem(distances.config, "distance_thresholds", {"ca_radius": 1.1, "all_atoms_radius": 1.1})
    monkeypatch.setitem(
        distances.config,
        "interaction_filters",
        {
            "protein_protein_min_ca_neighbors": 10,
            "protein_protein_min_all_atom_contacts": 1,
            "protein_nucleic_acid_min_all_atom_contacts": 1,
            "nucleic_acid_min_all_atom_contacts": 3,
        },
    )

    edges = calculate_distances_with_ckdtree([
        {"file_path": "/tmp/tst.cif", "atom_data": [_dna_chain("TST:A"), _dna_chain("TST:B")]}
    ])

    assert edges == []


def test_kdtree_cache_does_not_reuse_same_chain_ids_across_structures(monkeypatch) -> None:
    distances.coords_cache.clear()
    distances.tree_cache.clear()
    monkeypatch.setitem(distances.config, "distance_thresholds", {"ca_radius": 1.1, "all_atoms_radius": 1.1})
    monkeypatch.setitem(
        distances.config,
        "interaction_filters",
        {
            "protein_protein_min_ca_neighbors": 10,
            "protein_protein_min_all_atom_contacts": 1,
            "protein_nucleic_acid_min_all_atom_contacts": 1,
            "nucleic_acid_min_all_atom_contacts": 1,
        },
    )

    edges = calculate_distances_with_ckdtree([
        {
            "file_path": "/tmp/first.cif",
            "atom_data": [_protein_chain("DUP:A"), _protein_chain("DUP:B")],
        },
        {
            "file_path": "/tmp/second.cif",
            "atom_data": [_protein_chain("DUP:A", offset=1000.0), _protein_chain("DUP:B", offset=2000.0)],
        },
    ])

    assert len(edges) == 1
    assert edges[0]["file_path"] == "/tmp/first.cif"
