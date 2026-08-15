from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from pdb2net import detailed_results_exporter
from pdb2net.detailed_results_exporter import (
    DetailedInteractionBudget,
    export_detailed_interactions,
)
from pdb2net.input_contract import InputValidationError
from pdb2net.server_interface import DETAILED_INTERACTION_COLUMNS


def _structure(pdb_id: str, *, atoms_a: int = 1) -> dict:
    return {
        "pdb_id": pdb_id,
        "atom_data": [
            {
                "chain_id": "A",
                "unique_chain_id": f"{pdb_id}:A",
                "uniprot_id": "PAAAAA",
                "residues": [
                    {
                        "residue_name": "ALA",
                        "residue_number": 1,
                        "atoms": [
                            {
                                "atom_name": f"C{index}",
                                "coordinates": [float(index), 0.0, 0.0],
                            }
                            for index in range(atoms_a)
                        ],
                    }
                ],
            },
            {
                "chain_id": "B",
                "unique_chain_id": f"{pdb_id}:B",
                "uniprot_id": "PBBBBB",
                "residues": [
                    {
                        "residue_name": "GLY",
                        "residue_number": 2,
                        "atoms": [
                            {"atom_name": "CA", "coordinates": [0.5, 0.0, 0.0]}
                        ],
                    }
                ],
            },
        ],
    }


def _interactions(pdb_id: str) -> list[dict]:
    return [
        {
            "chain_a": f"{pdb_id}:A",
            "chain_b": f"{pdb_id}:B",
            "interaction_type": "Protein-Protein",
        }
    ]


def test_detailed_export_streams_stable_csv_and_records_exact_budget(tmp_path: Path) -> None:
    budget = DetailedInteractionBudget(max_rows=2, max_bytes=10_000)

    returned = export_detailed_interactions(
        _structure("TST1"), _interactions("TST1"), str(tmp_path), budget=budget
    )

    output = tmp_path / "TST1_detailed_interactions.csv"
    frame = pd.read_csv(output)
    assert returned is budget
    assert list(frame.columns) == list(DETAILED_INTERACTION_COLUMNS)
    assert len(frame) == budget.rows == 1
    assert budget.bytes == output.stat().st_size
    assert not list(tmp_path.glob(".*.tmp"))


def test_row_limit_fails_before_publication_and_is_cumulative(tmp_path: Path) -> None:
    budget = DetailedInteractionBudget(max_rows=1, max_bytes=10_000)
    export_detailed_interactions(
        _structure("ONE1"), _interactions("ONE1"), str(tmp_path), budget=budget
    )

    with pytest.raises(InputValidationError) as error:
        export_detailed_interactions(
            _structure("TWO2"), _interactions("TWO2"), str(tmp_path), budget=budget
        )

    assert error.value.code == "DETAILED_INTERACTION_ROW_LIMIT_EXCEEDED"
    assert (tmp_path / "ONE1_detailed_interactions.csv").is_file()
    assert not (tmp_path / "TWO2_detailed_interactions.csv").exists()
    assert not list(tmp_path.glob(".*.tmp"))


def test_row_limit_is_checked_before_large_neighbour_list_is_exported(tmp_path: Path) -> None:
    budget = DetailedInteractionBudget(max_rows=1, max_bytes=10_000)

    with pytest.raises(InputValidationError) as error:
        export_detailed_interactions(
            _structure("MANY", atoms_a=2),
            _interactions("MANY"),
            str(tmp_path),
            budget=budget,
        )

    assert error.value.code == "DETAILED_INTERACTION_ROW_LIMIT_EXCEEDED"
    assert budget.rows == 0
    assert not (tmp_path / "MANY_detailed_interactions.csv").exists()


def test_byte_limit_includes_header_and_removes_temporary_file(tmp_path: Path) -> None:
    budget = DetailedInteractionBudget(max_rows=10, max_bytes=1)

    with pytest.raises(InputValidationError) as error:
        export_detailed_interactions(
            _structure("BYTE"), _interactions("BYTE"), str(tmp_path), budget=budget
        )

    assert error.value.code == "DETAILED_INTERACTION_BYTE_LIMIT_EXCEEDED"
    assert not any(tmp_path.iterdir())


def test_free_space_reserve_is_checked_before_buffer_flush(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        detailed_results_exporter.shutil,
        "disk_usage",
        lambda _path: SimpleNamespace(free=100),
    )
    budget = DetailedInteractionBudget(
        max_rows=10,
        max_bytes=10_000,
        min_free_bytes=100,
    )

    with pytest.raises(InputValidationError) as error:
        export_detailed_interactions(
            _structure("DISK"), _interactions("DISK"), str(tmp_path), budget=budget
        )

    assert error.value.code == "DETAILED_INTERACTION_STORAGE_RESERVE_LOW"
    assert not (tmp_path / "DISK_detailed_interactions.csv").exists()
    assert not list(tmp_path.glob(".*.tmp"))
