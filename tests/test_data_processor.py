import gemmi
import pytest

from pdb2net import data_processor
from pdb2net.data_processor import process_structure
from pdb2net.input_contract import InputValidationError


class _Element:
    def __init__(self, name: str) -> None:
        self.name = name


class _Atom:
    def __init__(self, name: str, element: str, pos: tuple[float, float, float]) -> None:
        self.name = name
        self.element = _Element(element)
        self.pos = pos


class _SeqId:
    def __init__(self, num: int) -> None:
        self.num = num


class _Residue:
    def __init__(self, name: str, num: int, atoms: list[_Atom]) -> None:
        self.name = name
        self.seqid = _SeqId(num)
        self.entity_type = gemmi.EntityType.Polymer
        self._atoms = atoms

    def __iter__(self):
        return iter(self._atoms)


class _Chain:
    def __init__(self, name: str, residues: list[_Residue]) -> None:
        self.name = name
        self._residues = residues

    def __iter__(self):
        return iter(self._residues)


class _Model:
    def __init__(self, chains: list[_Chain]) -> None:
        self._chains = chains

    def __iter__(self):
        return iter(self._chains)


class _Structure:
    def __init__(self, models: list[_Model]) -> None:
        self._models = models

    def __iter__(self):
        return iter(self._models)


def _protein_chain(chain_id: str, residue_name: str = "ALA") -> _Chain:
    return _Chain(chain_id, [_Residue(residue_name, 1, [_Atom("CA", "C", (1.0, 2.0, 3.0))])])


def test_modified_residue_keeps_chain_and_records_original_name(monkeypatch) -> None:
    monkeypatch.setitem(data_processor.config, "structure_model_policy", "first")
    structure = _Structure([_Model([_protein_chain("A", "MSE")])])

    processed = process_structure({"file_path": "/tmp/test.cif", "pdb_id": "TST1", "structure": structure})

    chain = processed["atom_data"][0]
    residue = chain["residues"][0]
    assert residue["residue_name"] == "MET"
    assert residue["original_residue_name"] == "MSE"
    assert chain["unique_chain_id"] == "TST1:A"


def test_default_model_policy_uses_first_model_only(monkeypatch) -> None:
    monkeypatch.setitem(data_processor.config, "structure_model_policy", "first")
    structure = _Structure([
        _Model([_protein_chain("A")]),
        _Model([_protein_chain("B")]),
    ])

    processed = process_structure({"file_path": "/tmp/test.cif", "pdb_id": "TST1", "structure": structure})

    assert [chain["chain_id"] for chain in processed["atom_data"]] == ["A"]


def test_all_model_policy_keeps_model_in_unique_chain_id(monkeypatch) -> None:
    monkeypatch.setitem(data_processor.config, "structure_model_policy", "all")
    structure = _Structure([
        _Model([_protein_chain("A")]),
        _Model([_protein_chain("A")]),
    ])

    processed = process_structure({"file_path": "/tmp/test.cif", "pdb_id": "TST1", "structure": structure})

    assert [chain["unique_chain_id"] for chain in processed["atom_data"]] == [
        "TST1:model1:A",
        "TST1:model2:A",
    ]


@pytest.mark.parametrize("chain_id", ["", "  "])
@pytest.mark.parametrize("residue_name", ["ALA", "HOH"])
def test_missing_chain_id_is_a_file_error_before_residue_selection(monkeypatch, chain_id, residue_name):
    monkeypatch.setitem(data_processor.config, "structure_model_policy", "first")
    structure = _Structure([_Model([_protein_chain("A"), _protein_chain(chain_id, residue_name)])])

    with pytest.raises(InputValidationError) as error:
        process_structure({"file_path": "/tmp/test.cif", "pdb_id": "TST1", "structure": structure})

    assert error.value.code == "MISSING_CHAIN_ID"
    assert "model 1" in error.value.message


@pytest.mark.parametrize("extension", ["pdb", "cif"])
def test_missing_chain_id_in_later_model_obeys_model_selection(tmp_path, monkeypatch, extension):
    from pdb2net import pipeline

    text = (
        "MODEL        1\n"
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n"
        "ENDMDL\nMODEL        2\n"
        "ATOM      1  CA  ALA     1       3.000   0.000   0.000  1.00 20.00           C  \n"
        "ENDMDL\nEND\n"
    )
    if extension == "cif":
        text = gemmi.read_pdb_string(text).make_mmcif_document().as_string()
    source = tmp_path / f"models.{extension}"
    source.write_text(text)
    monkeypatch.setitem(data_processor.config, "structure_model_policy", "first")
    processed = pipeline.process_single_file(str(source))
    assert [chain["chain_id"] for chain in processed["atom_data"]] == ["A"]

    monkeypatch.setitem(data_processor.config, "structure_model_policy", "all")
    with pytest.raises(InputValidationError) as error:
        pipeline.process_single_file(str(source))
    assert error.value.code == "MISSING_CHAIN_ID"
    assert "model 2" in error.value.message


@pytest.mark.parametrize("extension", ["pdb", "cif"])
def test_assembly_metadata_does_not_generate_coordinate_chains(tmp_path, monkeypatch, extension):
    from pdb2net import file_parser, pipeline

    monkeypatch.setitem(pipeline.config, "structure_model_policy", "first")
    text = "HEADER".ljust(62) + "1ABC\n"
    text += "REMARK 350 BIOMOLECULE: 1\nREMARK 350 APPLY THE FOLLOWING TO CHAINS: A\n"
    for operation in (1, 2):
        for row in range(3):
            text += (f"REMARK 350   BIOMT{row+1} {operation:3d}  "
                     + " ".join(f"{float(i == row):9.6f}" for i in range(3))
                     + f"     {3. if operation == 2 and row == 1 else 0.:10.5f}\n")
    text += "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \nEND\n"
    if extension == "cif":
        document = gemmi.read_pdb_string(text).make_mmcif_document()
        document[0].name = "1abc"
        text = document.as_string()
    source = tmp_path / f"1abc.{extension}"
    source.write_text(text)
    parsed = file_parser.parse_structure_input(str(source))["structure"]
    assert len(parsed.assemblies) == 1
    assert len(parsed.assemblies[0].generators[0].operators) == 2
    processed = pipeline.process_single_file(str(source))
    assert [c["chain_id"] for c in processed["atom_data"]] == ["A"]
