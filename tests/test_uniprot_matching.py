import subprocess

import pytest

from pdb2net import uniprot_matcher
from pdb2net.uniprot_matcher import extract_direct_uniprot_accession, parallel_blast_search


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("AF-Q9BYF1-F1-model_v4.cif", "Q9BYF1"),
        ("AF-F4HVG8-2-F1-model_v4.cif", "F4HVG8-2"),
        ("AF_Q9BYF1_model_v4.cif", "Q9BYF1"),
        ("Q8WZ42-F1.cif", "Q8WZ42"),
        ("123456.cif", None),
        ("my_model.cif", None),
    ],
)
def test_extract_direct_uniprot_accession_from_alphafold_style_names(name: str, expected: str | None) -> None:
    assert extract_direct_uniprot_accession(name) == expected


def test_alphafold_single_chain_direct_uniprot_skips_blast(monkeypatch: pytest.MonkeyPatch) -> None:
    def fail_blast(*_args, **_kwargs):
        raise AssertionError("BLAST should not run when a direct AlphaFold UniProt ID is available")

    monkeypatch.setattr(uniprot_matcher, "run_blast_search", fail_blast)
    monkeypatch.setattr(uniprot_matcher, "_lookup_direct_uniprot_name", lambda _accession: "Test protein")

    parsed_data = [
        {
            "file_path": "/tmp/AF-Q9BYF1-F1-model_v4.cif",
            "pdb_id": "AF_Q9BYF1_F1_MODEL_V4",
            "atom_data": [
                {
                    "chain_id": "A",
                    "unique_chain_id": "AF_Q9BYF1_F1_MODEL_V4:A",
                    "molecule_type": "Unknown",
                    "uniprot_id": None,
                    "residues": [{"residue_name": "ALA"}, {"residue_name": "GLY"}],
                }
            ],
        }
    ]

    parallel_blast_search(parsed_data, max_workers=1)

    chain = parsed_data[0]["atom_data"][0]
    assert chain["uniprot_id"] == "Q9BYF1"
    assert chain["molecule_type"] == "Protein"
    assert chain["molecule_name"] == "Test protein"


def test_short_blast_query_uses_more_targets_one_hsp_and_blastp_short(
    tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    commands: list[list[str]] = []

    def fake_run(cmd, capture_output, text):
        commands.append(cmd)
        output_path = cmd[cmd.index("-out") + 1]
        with open(output_path, "w", encoding="utf-8") as handle:
            handle.write("")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(uniprot_matcher.subprocess, "run", fake_run)
    monkeypatch.setattr(uniprot_matcher, "BLAST_DB_PATH", str(tmp_path))
    monkeypatch.setattr(uniprot_matcher, "BLAST_EXECUTABLE", "blastp")

    hit = uniprot_matcher.run_blast_search("ACDEFGHIKLMNPQRSTVWY")

    assert hit is None
    assert commands
    cmd = commands[0]
    assert cmd[cmd.index("-max_target_seqs") + 1] == "100"
    assert cmd[cmd.index("-max_hsps") + 1] == "1"
    assert cmd[cmd.index("-task") + 1] == "blastp-short"
