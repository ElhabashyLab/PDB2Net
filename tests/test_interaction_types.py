import pytest

from pdb2net.distances import determine_interaction_type
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
