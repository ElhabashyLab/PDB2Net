import json
import itertools
import shutil
import subprocess
from pathlib import Path

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
    monkeypatch.setitem(uniprot_matcher.config, "blast_db_path", str(tmp_path))
    monkeypatch.setitem(uniprot_matcher.config, "blastp_executable", "blastp")

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

    monkeypatch.setattr(
        uniprot_matcher,
        "_run_blastp_swissprot_batch",
        lambda queries, **_kwargs: {
            key: uniprot_matcher._hit_outcome(swissprot_hit)
            for key, _sequence, _label in queries
        },
    )
    monkeypatch.setattr(uniprot_matcher, "_run_diamond_uniref90_batch", fail_diamond)

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
    monkeypatch.setitem(uniprot_matcher.config, "blast_db_path", str(tmp_path))
    monkeypatch.setitem(uniprot_matcher.config, "blastp_executable", "blastp")
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
    monkeypatch.setattr(
        uniprot_matcher,
        "_run_blastp_swissprot_batch",
        lambda queries, **_kwargs: {key: None for key, _sequence, _label in queries},
    )
    monkeypatch.setattr(
        uniprot_matcher,
        "_run_diamond_uniref90_batch",
        lambda queries, **_kwargs: {key: hit for key, _sequence, _label in queries},
    )
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


def test_diamond_batch_uses_one_multifasta_process_and_configured_limits(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    commands: list[list[str]] = []
    query_fastas: list[str] = []

    def fake_run(cmd, capture_output, text):
        commands.append(cmd)
        query_path = Path(cmd[cmd.index("--query") + 1])
        query_fastas.append(query_path.read_text(encoding="utf-8"))
        output_path = Path(cmd[cmd.index("--out") + 1])
        output_path.write_text(
            "q000000\tUniRef90_P12345\t98.5\t100\t1\t100\t100\t1e-80\t250\t"
            "UniRef90_P12345 Synthetic alpha enzyme n=1 Tax=Testus exampleus\n",
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(uniprot_matcher.subprocess, "run", fake_run)
    monkeypatch.setitem(
        uniprot_matcher.config,
        "diamond",
        {
            "enabled": True,
            "executable": "diamond",
            "uniref90_db_path": str(tmp_path / "mini_uniref90"),
            "threads": 6,
            "iterate": True,
            "sensitivity": "sensitive",
            "block_size": 1.0,
            "index_chunks": 4,
            "max_target_seqs": 50,
            "batch_max_sequences": 5000,
            "batch_max_fasta_bytes": 50 * 1024 * 1024,
            "assign_uniprot_id": "high_confidence",
        },
    )

    results = uniprot_matcher._run_diamond_uniref90_batch(
        [
            ("first", "A" * 100, "first chain"),
            ("second", "C" * 100, "second chain"),
        ]
    )

    assert len(commands) == 1
    assert query_fastas == [f">q000000\n{'A' * 100}\n>q000001\n{'C' * 100}\n"]
    cmd = commands[0]
    assert cmd[cmd.index("--threads") + 1] == "6"
    assert cmd[cmd.index("--block-size") + 1] == "1.0"
    assert cmd[cmd.index("--index-chunks") + 1] == "4"
    assert cmd[cmd.index("--max-target-seqs") + 1] == "50"
    assert "--iterate" in cmd
    assert "--sensitive" in cmd
    assert results["first"].hit is not None
    assert results["first"].hit.matched_id == "UniRef90_P12345"
    assert results["second"].status == "no_hit"


def test_swissprot_batch_uses_one_multifasta_process_for_same_task_group(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    commands: list[list[str]] = []
    query_fastas: list[str] = []

    def fake_run(cmd, capture_output, text):
        commands.append(cmd)
        query_path = Path(cmd[cmd.index("-query") + 1])
        query_fastas.append(query_path.read_text(encoding="utf-8"))
        output_path = Path(cmd[cmd.index("-out") + 1])
        output_path.write_text(
            "q000000\tsp|P11111|SYNTH_TEST\t99.0\t100\t1\t100\t100\t100\t1e-70\t230\t"
            "sp|P11111|SYNTH_TEST Synthetic Swiss-Prot enzyme OS=Testus exampleus\n",
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(uniprot_matcher.subprocess, "run", fake_run)
    monkeypatch.setitem(uniprot_matcher.config, "blast_db_path", str(tmp_path))
    monkeypatch.setitem(uniprot_matcher.config, "blastp_executable", "blastp")
    monkeypatch.setitem(
        uniprot_matcher.config,
        "blast",
        {"batch_max_sequences": 5000, "batch_max_fasta_bytes": 50 * 1024 * 1024},
    )

    results = uniprot_matcher._run_blastp_swissprot_batch(
        [
            ("first", "A" * 100, "first chain"),
            ("second", "C" * 100, "second chain"),
        ],
        max_workers=2,
    )

    assert len(commands) == 1
    assert query_fastas == [f">q000000\n{'A' * 100}\n>q000001\n{'C' * 100}\n"]
    assert results["first"].hit is not None
    assert results["first"].hit.accession == "P11111"
    assert results["second"].status == "no_hit"


def test_diamond_batch_chunks_at_configured_sequence_limit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    commands: list[list[str]] = []

    def fake_run(cmd, capture_output, text):
        commands.append(cmd)
        Path(cmd[cmd.index("--out") + 1]).write_text("", encoding="utf-8")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(uniprot_matcher.subprocess, "run", fake_run)
    monkeypatch.setitem(
        uniprot_matcher.config,
        "diamond",
        {
            "enabled": True,
            "uniref90_db_path": str(tmp_path / "mini_uniref90"),
            "batch_max_sequences": 2,
            "batch_max_fasta_bytes": 50 * 1024 * 1024,
        },
    )

    results = uniprot_matcher._run_diamond_uniref90_batch(
        [(str(index), "A" * 100, str(index)) for index in range(5)]
    )

    assert len(commands) == 3
    assert {key: outcome.status for key, outcome in results.items()} == {
        str(index): "no_hit" for index in range(5)
    }


def test_parallel_search_deduplicates_and_sends_only_swissprot_misses_to_diamond(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    swiss_hit = BlastHit(
        accession="P11111",
        reviewed=True,
        bitscore=200.0,
        evalue=1e-50,
        qcov=1.0,
        pident=100.0,
        title="sp|P11111| Swiss hit",
    )
    diamond_hit = BlastHit(
        accession="P12345",
        reviewed=False,
        bitscore=220.0,
        evalue=1e-60,
        qcov=1.0,
        pident=98.0,
        title="UniRef90_P12345 Synthetic alpha enzyme",
        source="diamond_uniref90",
        database="UniRef90",
        matched_id="UniRef90_P12345",
        representative_accession="P12345",
        confidence="high",
    )
    swiss_queries: list[str] = []
    diamond_batches: list[list[tuple[str, str, str]]] = []

    def fake_swiss(queries, **_kwargs):
        swiss_queries.extend(sequence for _key, sequence, _label in queries)
        return {
            key: swiss_hit if sequence == "AG" else None
            for key, sequence, _label in queries
        }

    def fake_diamond(queries, **_kwargs):
        diamond_batches.append(queries)
        return {key: diamond_hit for key, _sequence, _label in queries}

    monkeypatch.setattr(uniprot_matcher, "_CACHE_ENABLED", False)
    monkeypatch.setattr(uniprot_matcher, "_run_blastp_swissprot_batch", fake_swiss)
    monkeypatch.setattr(uniprot_matcher, "_run_diamond_uniref90_batch", fake_diamond)
    monkeypatch.setitem(
        uniprot_matcher.config,
        "diamond",
        {"enabled": True, "assign_uniprot_id": "high_confidence"},
    )

    def chain(chain_id: str, residues: list[str]) -> dict:
        return {
            "chain_id": chain_id,
            "unique_chain_id": f"TST1:{chain_id}",
            "molecule_type": "Unknown",
            "uniprot_id": None,
            "residues": [{"residue_name": residue} for residue in residues],
        }

    parsed_data = [
        {
            "file_path": "/tmp/model.cif",
            "pdb_id": "TST1",
            "atom_data": [
                chain("A", ["ALA", "GLY"]),
                chain("B", ["VAL", "LEU"]),
                chain("C", ["VAL", "LEU"]),
            ],
        }
    ]

    parallel_blast_search(parsed_data, max_workers=2)

    assert sorted(swiss_queries) == ["AG", "VL"]
    assert len(diamond_batches) == 1
    assert [sequence for _key, sequence, _label in diamond_batches[0]] == ["VL"]
    assert parsed_data[0]["atom_data"][0]["uniprot_id"] == "P11111"
    assert parsed_data[0]["atom_data"][1]["uniprot_id"] == "P12345"
    assert parsed_data[0]["atom_data"][2]["uniprot_id"] == "P12345"


def test_real_diamond_search_against_original_format_mini_uniref90(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    diamond = shutil.which("diamond")
    if diamond is None:
        pytest.skip("DIAMOND is not installed")

    fixture_dir = Path(__file__).parent / "fixtures" / "uniref90"
    fasta_path = fixture_dir / "mini_uniref90.fasta"
    profile = json.loads((fixture_dir / "config.test.json").read_text(encoding="utf-8"))
    db_prefix = tmp_path / "mini_uniref90"
    subprocess.run(
        [diamond, "makedb", "--in", str(fasta_path), "--db", str(db_prefix)],
        check=True,
        capture_output=True,
        text=True,
    )
    diamond_config = profile["diamond"]
    diamond_config["executable"] = diamond
    diamond_config["uniref90_db_path"] = str(db_prefix)
    monkeypatch.setitem(uniprot_matcher.config, "diamond", diamond_config)

    first_sequence = "".join(
        line.strip()
        for line in fasta_path.read_text(encoding="utf-8").splitlines()[1:2]
    )
    outcome = uniprot_matcher._run_diamond_uniref90_batch(
        [("fixture-query", first_sequence, "fixture query")]
    )["fixture-query"]
    result = outcome.hit

    assert result is not None
    assert result.matched_id == "UniRef90_P12345"
    assert result.representative_accession == "P12345"
    assert result.confidence == "high"


@pytest.mark.parametrize("batch_kind", ["swissprot", "diamond"])
@pytest.mark.parametrize("failure_mode", ["nonzero", "missing_output", "exception"])
def test_batch_failures_are_distinct_from_no_hit(
    batch_kind: str,
    failure_mode: str,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_run(cmd, capture_output, text):
        if failure_mode == "exception":
            raise OSError("search process could not start")
        if failure_mode == "nonzero":
            return subprocess.CompletedProcess(cmd, 2, stdout="", stderr="fatal search error")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(uniprot_matcher.subprocess, "run", fake_run)
    query = [("query-key", "ACDEFGHIKLMNPQRSTVWY" * 5, "test query")]
    if batch_kind == "swissprot":
        monkeypatch.setitem(uniprot_matcher.config, "blast_db_path", str(tmp_path))
        outcome = uniprot_matcher._run_blastp_swissprot_batch(query)["query-key"]
    else:
        monkeypatch.setitem(
            uniprot_matcher.config,
            "diamond",
            {
                "enabled": True,
                "uniref90_db_path": str(tmp_path / "mini_uniref90"),
            },
        )
        outcome = uniprot_matcher._run_diamond_uniref90_batch(query)["query-key"]

    assert outcome.status == "error"
    assert outcome.hit is None
    assert outcome.error


@pytest.mark.parametrize("batch_kind", ["swissprot", "diamond"])
def test_successful_empty_batch_output_is_a_real_no_hit(
    batch_kind: str,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_run(cmd, capture_output, text):
        output_flag = "-out" if batch_kind == "swissprot" else "--out"
        Path(cmd[cmd.index(output_flag) + 1]).write_text("", encoding="utf-8")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    monkeypatch.setattr(uniprot_matcher.subprocess, "run", fake_run)
    query = [("query-key", "ACDEFGHIKLMNPQRSTVWY" * 5, "test query")]
    if batch_kind == "swissprot":
        monkeypatch.setitem(uniprot_matcher.config, "blast_db_path", str(tmp_path))
        outcome = uniprot_matcher._run_blastp_swissprot_batch(query)["query-key"]
    else:
        monkeypatch.setitem(
            uniprot_matcher.config,
            "diamond",
            {
                "enabled": True,
                "uniref90_db_path": str(tmp_path / "mini_uniref90"),
            },
        )
        outcome = uniprot_matcher._run_diamond_uniref90_batch(query)["query-key"]

    assert outcome.status == "no_hit"
    assert outcome.hit is None
    assert not outcome.error


def test_search_error_is_not_cached_or_forwarded_to_diamond(monkeypatch: pytest.MonkeyPatch) -> None:
    def fail_diamond(*_args, **_kwargs):
        raise AssertionError("DIAMOND must not run after a Swiss-Prot process error")

    monkeypatch.setattr(uniprot_matcher, "_cache_get", lambda _key: (False, None))
    monkeypatch.setattr(
        uniprot_matcher,
        "_cache_put",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("search errors must not be cached")),
    )
    monkeypatch.setattr(
        uniprot_matcher,
        "_run_blastp_swissprot_batch",
        lambda queries, **_kwargs: {
            key: uniprot_matcher._error_outcome("blast failed")
            for key, _sequence, _label in queries
        },
    )
    monkeypatch.setattr(uniprot_matcher, "_run_diamond_uniref90_batch", fail_diamond)
    monkeypatch.setitem(uniprot_matcher.config, "diamond", {"enabled": True})
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

    with pytest.raises(uniprot_matcher.SequenceSearchError, match="blast failed"):
        parallel_blast_search(parsed_data, max_workers=1)

    assert parsed_data[0]["atom_data"][0].get("uniprot_id") is None


def test_diamond_error_after_real_swissprot_miss_is_not_cached(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(uniprot_matcher, "_cache_get", lambda _key: (False, None))
    monkeypatch.setattr(
        uniprot_matcher,
        "_cache_put",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("search errors must not be cached")),
    )
    monkeypatch.setattr(
        uniprot_matcher,
        "_run_blastp_swissprot_batch",
        lambda queries, **_kwargs: {
            key: uniprot_matcher._no_hit_outcome()
            for key, _sequence, _label in queries
        },
    )
    monkeypatch.setattr(
        uniprot_matcher,
        "_run_diamond_uniref90_batch",
        lambda queries, **_kwargs: {
            key: uniprot_matcher._error_outcome("diamond failed")
            for key, _sequence, _label in queries
        },
    )
    monkeypatch.setitem(uniprot_matcher.config, "diamond", {"enabled": True})
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

    with pytest.raises(uniprot_matcher.SequenceSearchError, match="diamond failed"):
        parallel_blast_search(parsed_data, max_workers=1)


def test_legitimate_no_hit_is_persistently_cacheable(monkeypatch: pytest.MonkeyPatch) -> None:
    cache_writes: list[tuple[str, BlastHit | None]] = []
    monkeypatch.setattr(uniprot_matcher, "_cache_get", lambda _key: (False, None))
    monkeypatch.setattr(uniprot_matcher, "_cache_put", lambda key, hit: cache_writes.append((key, hit)))
    monkeypatch.setattr(
        uniprot_matcher,
        "_run_blastp_swissprot_batch",
        lambda queries, **_kwargs: {
            key: uniprot_matcher._no_hit_outcome()
            for key, _sequence, _label in queries
        },
    )
    monkeypatch.setitem(uniprot_matcher.config, "diamond", {"enabled": False})
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

    assert len(cache_writes) == 1
    assert cache_writes[0][1] is None


def test_cache_signature_tracks_result_policies_but_not_chunk_sizes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    swissprot = tmp_path / "uniprot_sprot.fasta"
    swissprot.write_text(">sp|P12345|TEST\nACDEFGHIK\n", encoding="utf-8")
    diamond_db = tmp_path / "mini_uniref90.dmnd"
    diamond_db.write_bytes(b"test database signature")
    monkeypatch.setitem(uniprot_matcher.config, "uniprot_fasta_path", str(swissprot))
    baseline_config = {
        "enabled": True,
        "uniref90_db_path": str(diamond_db),
        "sensitivity": "sensitive",
        "iterate": True,
        "max_target_seqs": 50,
        "assign_uniprot_id": "never",
        "batch_max_sequences": 5000,
        "batch_max_fasta_bytes": 50 * 1024 * 1024,
    }

    def signature(**overrides):
        monkeypatch.setitem(
            uniprot_matcher.config,
            "diamond",
            {**baseline_config, **overrides},
        )
        return uniprot_matcher._db_signature()

    baseline = signature()
    assert signature(sensitivity="very-sensitive") != baseline
    assert signature(iterate=False) != baseline
    assert signature(block_size=2.0) != baseline
    assert signature(max_target_seqs=25) != baseline
    assert signature(assign_uniprot_id="high_confidence") != baseline
    assert signature(batch_max_sequences=1, batch_max_fasta_bytes=1024) == baseline

    monkeypatch.setitem(uniprot_matcher.config, "diamond", baseline_config)
    monkeypatch.setitem(
        uniprot_matcher.config,
        "blast",
        {**uniprot_matcher.config["blast"], "max_target_seqs_default": 51},
    )
    assert uniprot_matcher._db_signature() != baseline

    payload = uniprot_matcher._search_policy_payload()
    assert "thresholds" in payload["swissprot"]
    assert "thresholds" in payload["diamond"]
    assert "batch_max_sequences" not in payload["swissprot"]
    assert "batch_max_sequences" not in payload["diamond"]


def test_cache_signature_tracks_every_blast_component_and_reference_manifest(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fasta = tmp_path / "uniprot_sprot.fasta"
    fasta.write_text(">sp|P11111|TEST\nAAAA\n", encoding="utf-8")
    blast_dir = tmp_path / "blast"
    blast_dir.mkdir()
    for suffix in ("pin", "phr", "psq"):
        (blast_dir / f"uniprot_db.{suffix}").write_bytes(suffix.encode("ascii"))
    monkeypatch.setitem(uniprot_matcher.config, "uniprot_fasta_path", str(fasta))
    monkeypatch.setitem(uniprot_matcher.config, "blast_db_path", str(blast_dir))
    monkeypatch.setitem(uniprot_matcher.config, "diamond", {"enabled": False})
    monkeypatch.setitem(uniprot_matcher.config, "reference_manifest_id", "release-a")
    baseline = uniprot_matcher._db_signature()

    (blast_dir / "uniprot_db.psq").write_bytes(b"changed component")
    component_changed = uniprot_matcher._db_signature()
    monkeypatch.setitem(uniprot_matcher.config, "reference_manifest_id", "release-b")

    assert component_changed != baseline
    assert uniprot_matcher._db_signature() != component_changed
    assert uniprot_matcher._search_policy_payload()["swissprot"]["selection_policy"]["version"] == 2


def _swissprot_row(accession: str, *, pident: float = 100.0, bitscore: float = 250.0) -> str:
    return (
        f"query\tsp|{accession}|TEST\t{pident}\t100\t1\t100\t100\t100\t"
        f"1e-80\t{bitscore}\tsp|{accession}|TEST Synthetic protein"
    )


def test_perfect_swissprot_ties_are_deterministic_across_all_row_orders() -> None:
    rows = [_swissprot_row(accession) for accession in ("P33333", "P22222", "P11111")]
    selected = {
        uniprot_matcher._select_blastp_swissprot_hit(
            "A" * 100,
            list(permutation),
            use_qcovs=True,
        ).accession
        for permutation in itertools.permutations(rows)
    }

    assert selected == {"P11111"}


def test_perfect_tie_escape_requires_every_contender_to_be_near_perfect() -> None:
    hit = uniprot_matcher._select_blastp_swissprot_hit(
        "A" * 100,
        [
            _swissprot_row("P11111", pident=100.0, bitscore=250.0),
            _swissprot_row("P22222", pident=99.5, bitscore=249.0),
            _swissprot_row("P33333", pident=98.0, bitscore=248.0),
        ],
        use_qcovs=True,
    )

    assert hit is None


def test_high_confidence_uniref_upi_does_not_assign_canonical_uniprot_id(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    upi_hit = BlastHit(
        accession="UPI000012345678",
        reviewed=False,
        bitscore=240.0,
        evalue=1e-70,
        qcov=1.0,
        pident=99.0,
        title="UniRef90_UPI000012345678 Synthetic cluster",
        source="diamond_uniref90",
        database="UniRef90",
        matched_id="UniRef90_UPI000012345678",
        representative_accession="UPI000012345678",
        confidence="high",
    )
    monkeypatch.setattr(uniprot_matcher, "_CACHE_ENABLED", False)
    monkeypatch.setattr(
        uniprot_matcher,
        "_run_blastp_swissprot_batch",
        lambda queries, **_kwargs: {
            key: uniprot_matcher._no_hit_outcome()
            for key, _sequence, _label in queries
        },
    )
    monkeypatch.setattr(
        uniprot_matcher,
        "_run_diamond_uniref90_batch",
        lambda queries, **_kwargs: {
            key: uniprot_matcher._hit_outcome(upi_hit)
            for key, _sequence, _label in queries
        },
    )
    monkeypatch.setitem(
        uniprot_matcher.config,
        "diamond",
        {"enabled": True, "assign_uniprot_id": "high_confidence"},
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
    assert chain["annotation_source"] == "diamond_uniref90"
    assert chain["matched_id"] == "UniRef90_UPI000012345678"
    assert chain["representative_accession"] == "UPI000012345678"
    assert chain["annotation_confidence"] == "high"


@pytest.mark.parametrize(
    ("accession", "expected"),
    [
        ("P12345", True),
        ("Q9BYF1", True),
        ("A0A1B2", True),
        ("P12345-2", True),
        ("UPI000012345678", False),
        ("UniRef90_P12345", False),
        ("prefix-P12345-suffix", False),
    ],
)
def test_uniprotkb_accession_validation_requires_a_complete_accession(
    accession: str,
    expected: bool,
) -> None:
    assert uniprot_matcher._is_uniprotkb_accession(accession) is expected


def test_fasta_byte_cap_uses_exact_record_boundary() -> None:
    queries = [("a", "AAAA", "first"), ("b", "CCCC", "second")]
    exact_bytes = sum(len(f">{key}\n{sequence}\n".encode("utf-8")) for key, sequence, _ in queries)

    exact_chunks = uniprot_matcher._chunk_sequence_queries(
        queries,
        max_sequences=5000,
        max_fasta_bytes=exact_bytes,
    )
    below_boundary_chunks = uniprot_matcher._chunk_sequence_queries(
        queries,
        max_sequences=5000,
        max_fasta_bytes=exact_bytes - 1,
    )

    assert exact_chunks == [queries]
    assert below_boundary_chunks == [[queries[0]], [queries[1]]]
