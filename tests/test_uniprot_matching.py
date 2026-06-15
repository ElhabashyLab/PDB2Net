import subprocess

import pytest

from pdb2net import uniprot_matcher
from pdb2net.uniprot_matcher import BlastHit, extract_direct_uniprot_accession, parallel_blast_search


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


def test_swissprot_hit_skips_diamond_fallback(monkeypatch: pytest.MonkeyPatch) -> None:
    swissprot_hit = BlastHit(
        accession="P12345",
        reviewed=True,
        bitscore=200.0,
        evalue=1e-50,
        qcov=0.99,
        pident=98.0,
        title="sp|P12345| Example protein",
    )

    def fail_diamond(*_args, **_kwargs):
        raise AssertionError("DIAMOND should not run when Swiss-Prot BLAST accepted a hit")

    monkeypatch.setattr(uniprot_matcher, "_run_blastp_swissprot_search", lambda *_args, **_kwargs: swissprot_hit)
    monkeypatch.setattr(uniprot_matcher, "_run_diamond_uniref90_search", fail_diamond)

    assert uniprot_matcher.run_blast_search("A" * 100) == swissprot_hit


def test_diamond_uniref90_fallback_parses_tabular_output(tmp_path, monkeypatch: pytest.MonkeyPatch) -> None:
    commands: list[list[str]] = []

    def fake_run(cmd, capture_output, text):
        commands.append(cmd)
        if cmd[0] == "blastp":
            output_path = cmd[cmd.index("-out") + 1]
            with open(output_path, "w", encoding="utf-8") as handle:
                handle.write("")
            return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

        output_path = cmd[cmd.index("--out") + 1]
        with open(output_path, "w", encoding="utf-8") as handle:
            handle.write(
                "query\tUniRef90_P12345\t98.5\t100\t1\t100\t100\t1e-80\t250\t"
                "UniRef90_P12345 Cluster: Example oxidase n=12 Tax=Testus exampleus\n"
            )
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(uniprot_matcher.subprocess, "run", fake_run)
    monkeypatch.setattr(uniprot_matcher, "_CACHE_ENABLED", False)
    monkeypatch.setattr(uniprot_matcher, "BLAST_DB_PATH", str(tmp_path))
    monkeypatch.setattr(uniprot_matcher, "BLAST_EXECUTABLE", "blastp")
    monkeypatch.setitem(
        uniprot_matcher.config,
        "diamond",
        {
            "enabled": True,
            "executable": "diamond",
            "uniref90_db_path": str(tmp_path / "uniref90.dmnd"),
            "max_target_seqs": 50,
            "assign_uniprot_id": "never",
        },
    )

    hit = uniprot_matcher.run_blast_search("A" * 100)

    assert hit is not None
    assert hit.source == "diamond_uniref90"
    assert hit.database == "UniRef90"
    assert hit.matched_id == "UniRef90_P12345"
    assert hit.representative_accession == "P12345"
    assert hit.accession == "P12345"
    assert hit.confidence == "high"
    assert any(cmd[0] == "diamond" for cmd in commands)


def test_diamond_uniref90_annotation_keeps_uniprot_id_conservative(monkeypatch: pytest.MonkeyPatch) -> None:
    hit = BlastHit(
        accession="P12345",
        reviewed=False,
        bitscore=250.0,
        evalue=1e-80,
        qcov=1.0,
        pident=98.5,
        title="UniRef90_P12345 Cluster: Example oxidase n=12",
        source="diamond_uniref90",
        database="UniRef90",
        matched_id="UniRef90_P12345",
        representative_accession="P12345",
        confidence="high",
    )
    monkeypatch.setattr(uniprot_matcher, "run_blast_search", lambda *_args, **_kwargs: hit)
    monkeypatch.setattr(uniprot_matcher, "_CACHE_ENABLED", False)
    monkeypatch.setitem(
        uniprot_matcher.config,
        "diamond",
        {
            "enabled": True,
            "executable": "diamond",
            "uniref90_db_path": "/tmp/uniref90.dmnd",
            "max_target_seqs": 50,
            "assign_uniprot_id": "never",
        },
    )
    parsed_data = [
        {
            "file_path": "/tmp/model.cif",
            "pdb_id": "TST1",
            "atom_data": [
                {
                    "chain_id": "A",
                    "unique_chain_id": "TST1:A",
                    "molecule_type": "Unknown",
                    "uniprot_id": None,
                    "residues": [{"residue_name": "ALA"}, {"residue_name": "GLY"}],
                }
            ],
        }
    ]

    parallel_blast_search(parsed_data, max_workers=1)

    chain = parsed_data[0]["atom_data"][0]
    assert chain.get("uniprot_id") is None
    assert chain["molecule_type"] == "Protein"
    assert chain["annotation_source"] == "diamond_uniref90"
    assert chain["matched_database"] == "UniRef90"
    assert chain["matched_id"] == "UniRef90_P12345"
    assert chain["representative_accession"] == "P12345"
    assert chain["annotation_confidence"] == "high"
