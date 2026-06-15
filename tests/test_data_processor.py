from pdb2net import data_processor
from pdb2net.data_processor import process_structure


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
