from __future__ import annotations

import copy
import gzip
import json
import threading
import time
from pathlib import Path
from types import SimpleNamespace

import pytest

from pdb2net import precomputed_store as store
from pdb2net.__main__ import main
from pdb2net.config_loader import config
from pdb2net.input_contract import InputValidationError


@pytest.fixture(autouse=True)
def _stable_precomputed_config(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    pdb_fasta = tmp_path / "pdb_seqres.txt"
    pdb_fasta.write_text(">1abc_A mol:protein length:1\nA\n", encoding="utf-8")
    monkeypatch.setitem(config, "pdb_fasta_path", str(pdb_fasta))
    monkeypatch.setitem(config, "reference_manifest_id", "test-references-v1")
    monkeypatch.setitem(config, "structure_model_policy", "first")
    monkeypatch.setitem(config, "export_detailed_interactions", False)
    monkeypatch.setitem(config, "open_in_cytoscape", False)
    monkeypatch.setitem(config, "output_path", str(tmp_path / "outputs"))


def _raw_structure(source: Path, pdb_id: str = "1ABC") -> dict:
    return {
        "file_path": str(source),
        "pdb_id": pdb_id,
        "atom_data": [
            {
                "chain_id": "A",
                "unique_chain_id": f"{pdb_id}:A",
                "model_index": 1,
                "molecule_name": "Protein Alpha",
                "molecule_type": "Protein",
                "sequence": "A",
                "is_hetatm": False,
                "uniprot_id": "P12345",
                "annotation_source": "sifts",
                "matched_database": "UniProtKB",
                "annotation_confidence": "high",
                "residues": [
                    {
                        "residue_name": "ALA",
                        "residue_number": 1,
                        "atoms": [{"atom_name": "CA", "coordinates": [0.0, 0.0, 0.0]}],
                    }
                ],
            },
            {
                "chain_id": "B",
                "unique_chain_id": f"{pdb_id}:B",
                "model_index": 1,
                "molecule_name": "Protein Beta",
                "molecule_type": "Protein",
                "sequence": "G",
                "is_hetatm": False,
                "uniprot_id": "Q99999",
                "annotation_source": "sifts",
                "matched_database": "UniProtKB",
                "annotation_confidence": "high",
                "residues": [
                    {
                        "residue_name": "GLY",
                        "residue_number": 1,
                        "atoms": [{"atom_name": "CA", "coordinates": [1.0, 0.0, 0.0]}],
                    }
                ],
            },
        ],
    }


def _raw_edges(source: Path, pdb_id: str = "1ABC") -> list[dict]:
    return [
        {
            "file_path": str(source),
            "chain_a": f"{pdb_id}:A",
            "chain_b": f"{pdb_id}:B",
            "ca_neighbors": 11,
            "all_atoms_count": 17,
            "interaction_type": "Protein-Protein",
        }
    ]


def _write_fixture_entry(tmp_path: Path) -> tuple[Path, Path, dict]:
    source = tmp_path / "1abc.pdb"
    source.write_text("HEADER TEST 1ABC\n", encoding="utf-8")
    store_root = tmp_path / "store"
    path, changed = store.write_entry(
        store_root,
        source,
        _raw_structure(source),
        _raw_edges(source),
    )
    assert changed
    return store_root, source, store.load_entry(store_root, "1abc")


def _rewrite_payload(store_root: Path, pdb_id: str, payload: dict) -> None:
    store.entry_path(store_root, pdb_id).write_bytes(
        gzip.compress(json.dumps(payload).encode("utf-8"), mtime=0)
    )


def _scientific_view(structure: dict, interactions: list[dict]) -> tuple[list[dict], list[dict]]:
    chains = [
        {
            key: value
            for key, value in chain.items()
            if key not in {"_parent_file_path", "_parent_file_label", "_parent_pdb_id"}
        }
        for chain in structure["atom_data"]
    ]
    edges = [
        {key: value for key, value in edge.items() if key != "file_path"}
        for edge in interactions
    ]
    return chains, edges


def test_raw_and_precomputed_scientific_chain_and_edge_data_are_equal(tmp_path: Path) -> None:
    from pdb2net.pipeline import _compact_structure_summaries

    store_root, source, cached = _write_fixture_entry(tmp_path)
    raw_compact = _compact_structure_summaries([_raw_structure(source)])[0]

    assert _scientific_view(raw_compact, _raw_edges(source)) == _scientific_view(
        cached["structure"], cached["interactions"]
    )
    assert "residues" not in cached["structure"]["atom_data"][0]
    assert str(source.parent) not in gzip.decompress(
        store.entry_path(store_root, "1abc").read_bytes()
    ).decode("utf-8")


def test_entry_is_profile_namespaced_and_same_source_is_a_cache_hit(tmp_path: Path) -> None:
    store_root, source, cached = _write_fixture_entry(tmp_path)
    expected = (
        store_root
        / "profiles"
        / store.profile_id()
        / "entries"
        / "ab"
        / "1abc.json.gz"
    )
    assert expected.is_file()
    assert cached["source"]["scope"] == "asymmetric_unit"

    path, changed = store.write_entry(
        store_root,
        source,
        _raw_structure(source),
        _raw_edges(source),
    )
    assert path == expected
    assert not changed


def test_precompute_uses_normal_scientific_stages_once_then_hits_cache(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    from pdb2net import pipeline
    from pdb2net import unknown_molecule_uniprot

    source_dir = tmp_path / "inputs"
    source_dir.mkdir()
    source = source_dir / "1abc.pdb"
    source.write_text("HEADER TEST 1ABC\n", encoding="utf-8")
    store_root = tmp_path / "store"
    calls = {"parse": 0, "annotate": 0, "search": 0, "distance": 0}

    monkeypatch.setattr(pipeline, "_validate_required_reference_files", lambda: None)

    def fake_parse(paths):
        assert paths == [str(source.resolve())]
        calls["parse"] += 1
        return [_raw_structure(source.resolve())]

    def fake_annotate(structures, **_kwargs):
        assert structures[0]["atom_data"][0]["uniprot_id"] == "P12345"
        calls["annotate"] += 1

    def fake_search(structures):
        assert structures
        calls["search"] += 1

    def fake_distances(structures):
        calls["distance"] += 1
        return _raw_edges(source.resolve())

    monkeypatch.setattr(pipeline, "_parse_input_files", fake_parse)
    monkeypatch.setattr(unknown_molecule_uniprot, "process_molecule_info", fake_annotate)
    monkeypatch.setattr(pipeline, "_run_blast_annotation", fake_search)
    monkeypatch.setattr(store, "calculate_distances_with_ckdtree", fake_distances)

    first = store.precompute_directory(store_root, source_dir, recursive=False)
    assert first["written"] == 1
    assert first["failed"] == 0
    assert calls == {"parse": 1, "annotate": 1, "search": 1, "distance": 1}

    def forbidden(*_args, **_kwargs):
        raise AssertionError("an unchanged cache entry must skip every scientific stage")

    monkeypatch.setattr(pipeline, "_parse_input_files", forbidden)
    monkeypatch.setattr(unknown_molecule_uniprot, "process_molecule_info", forbidden)
    monkeypatch.setattr(pipeline, "_run_blast_annotation", forbidden)
    monkeypatch.setattr(store, "calculate_distances_with_ckdtree", forbidden)
    second = store.precompute_directory(store_root, source_dir, recursive=False)
    assert second["cache_hits"] == 1
    assert second["written"] == 0


def test_precompute_rejects_source_changed_during_processing(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    from pdb2net import pipeline, unknown_molecule_uniprot

    source_dir = tmp_path / "inputs"
    source_dir.mkdir()
    source = source_dir / "1abc.pdb"
    source.write_text("HEADER ORIGINAL 1ABC\n", encoding="utf-8")
    monkeypatch.setattr(pipeline, "_validate_required_reference_files", lambda: None)

    def fake_parse(_paths):
        parsed = _raw_structure(source.resolve())
        source.write_text("HEADER REVISED 1ABC\n", encoding="utf-8")
        return [parsed]

    monkeypatch.setattr(pipeline, "_parse_input_files", fake_parse)
    monkeypatch.setattr(
        unknown_molecule_uniprot,
        "process_molecule_info",
        lambda _structures, **_kwargs: None,
    )
    monkeypatch.setattr(pipeline, "_run_blast_annotation", lambda _structures: None)
    monkeypatch.setattr(
        store,
        "calculate_distances_with_ckdtree",
        lambda _structures: _raw_edges(source.resolve()),
    )

    report = store.precompute_directory(tmp_path / "store", source_dir, recursive=False)

    assert report["failed"] == 1
    assert report["errors"][0]["code"] == "PRECOMPUTE_SOURCE_CHANGED"
    assert not store.entry_path(tmp_path / "store", "1abc").exists()


def test_changed_scientific_profile_uses_a_new_namespace(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    store_root, _, _ = _write_fixture_entry(tmp_path)
    old_path = store.entry_path(store_root, "1abc")
    old_profile = old_path.parents[2].name

    monkeypatch.setitem(config["distance_thresholds"], "ca_radius", 13.0)
    assert store.profile_id() != old_profile
    with pytest.raises(InputValidationError, match="PRECOMPUTED_ENTRY_MISSING"):
        store.load_entry(store_root, "1abc")
    assert old_path.is_file()


def test_profile_tracks_result_policy_but_not_batch_chunk_sizes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    original = store.profile_id()
    monkeypatch.setitem(config["diamond"], "block_size", 2.0)
    assert store.profile_id() != original

    monkeypatch.setitem(config["diamond"], "block_size", 1.0)
    restored = store.profile_id()
    monkeypatch.setitem(config["diamond"], "batch_max_sequences", 17)
    monkeypatch.setitem(config["diamond"], "batch_max_fasta_bytes", 12345)
    assert store.profile_id() == restored


def test_cache_validator_rejects_reference_chain_and_contact_corruption(tmp_path: Path) -> None:
    store_root, _, payload = _write_fixture_entry(tmp_path)

    bad_reference = copy.deepcopy(payload)
    bad_reference["references"]["manifest_id"] = "other-release"
    _rewrite_payload(store_root, "1abc", bad_reference)
    with pytest.raises(InputValidationError) as exc_info:
        store.load_entry(store_root, "1abc")
    assert exc_info.value.code == "CACHE_REFERENCE_MISMATCH"

    _rewrite_payload(store_root, "1abc", payload)
    bad_chain = copy.deepcopy(payload)
    bad_chain["structure"]["atom_data"][0]["uniprot_id"] = "Unknown"
    _rewrite_payload(store_root, "1abc", bad_chain)
    with pytest.raises(InputValidationError, match="UniProt accession"):
        store.load_entry(store_root, "1abc")

    bad_contacts = copy.deepcopy(payload)
    del bad_contacts["interactions"][0]["ca_neighbors"]
    _rewrite_payload(store_root, "1abc", bad_contacts)
    with pytest.raises(InputValidationError, match="protein contact counts"):
        store.load_entry(store_root, "1abc")


def test_cache_write_rejects_parent_symlink(tmp_path: Path) -> None:
    store_root = tmp_path / "store"
    _selected, profile_directory = store.ensure_profile_manifest(store_root)
    outside = tmp_path / "outside"
    outside.mkdir()
    (profile_directory / "entries").symlink_to(outside, target_is_directory=True)
    source = tmp_path / "1abc.pdb"
    source.write_text("HEADER TEST 1ABC\n", encoding="utf-8")

    with pytest.raises(InputValidationError) as exc_info:
        store.write_entry(store_root, source, _raw_structure(source), _raw_edges(source))

    assert exc_info.value.code == "UNSAFE_CACHE_PATH"
    assert list(outside.iterdir()) == []


def test_cache_hit_assemble_never_parses_or_searches(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    store_root, _, cached = _write_fixture_entry(tmp_path)
    from pdb2net import file_parser, pipeline

    def forbidden(*_args, **_kwargs):
        raise AssertionError("cache hit must not parse, annotate, search, or calculate distances")

    monkeypatch.setattr(pipeline, "_parse_input_files", forbidden)
    monkeypatch.setattr(pipeline, "_run_blast_annotation", forbidden)
    monkeypatch.setattr(store, "calculate_distances_with_ckdtree", forbidden)
    monkeypatch.setattr(file_parser, "VALID_PDB_IDS", None)
    monkeypatch.setattr(file_parser, "load_valid_pdb_ids", forbidden)
    captured: dict = {}

    def capture_exports(structures, interactions, *_args, **_kwargs):
        captured["structures"] = copy.deepcopy(structures)
        captured["interactions"] = copy.deepcopy(interactions)
        return []

    monkeypatch.setattr(pipeline, "_export_network_outputs", capture_exports)
    result = store.run_assemble_pipeline(store_root, ["1ABC"])

    assert Path(result.summary_file).is_file()
    assert captured["structures"] == [cached["structure"]]
    assert captured["interactions"] == cached["interactions"]
    summary = json.loads(Path(result.summary_file).read_text(encoding="utf-8"))
    assert summary["status"] == "success"
    assert summary["resources"]["processing"]["cache_hits"] == 1
    assert summary["input_files"] == ["pdb:1ABC"]


def test_missing_entry_can_be_populated_from_exact_sharded_gzip_source(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    archive = tmp_path / "archive"
    source = archive / "ab" / "pdb_00001abc.cif.gz"
    source.parent.mkdir(parents=True)
    source.write_bytes(gzip.compress(b"data_1abc\n"))
    store_root = tmp_path / "store"

    def fake_precompute(target_store, paths):
        assert paths == [source.resolve()]
        store.write_entry(
            target_store,
            source,
            _raw_structure(source),
            _raw_edges(source),
        )
        return {"failed": 0, "errors": []}

    monkeypatch.setattr(store, "precompute_sources", fake_precompute)
    payload, was_hit = store._load_or_populate_entry(
        store_root,
        "1abc",
        source_dir=archive,
        populate_missing=True,
    )
    assert not was_hit
    assert payload["pdb_id"] == "1abc"


def test_multiple_lazy_misses_are_precomputed_in_one_batch(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    archive = tmp_path / "archive"
    archive.mkdir()
    first = archive / "1abc.cif"
    second = archive / "2xyz.cif"
    first.write_text("data_1abc\n", encoding="utf-8")
    second.write_text("data_2xyz\n", encoding="utf-8")
    calls: list[list[Path]] = []

    def fake_precompute(target_store, paths):
        calls.append(paths)
        for source, pdb_id in ((first, "1ABC"), (second, "2XYZ")):
            store.write_entry(
                target_store,
                source,
                _raw_structure(source, pdb_id),
                _raw_edges(source, pdb_id),
            )
        return {"failed": 0, "errors": []}

    monkeypatch.setattr(store, "precompute_sources", fake_precompute)
    loaded = store._load_or_populate_entries(
        tmp_path / "store",
        ["1abc", "2xyz"],
        source_dir=archive,
        populate_missing=True,
    )

    assert calls == [[first.resolve(), second.resolve()]]
    assert [payload["pdb_id"] for payload, _hit in loaded] == ["1abc", "2xyz"]
    assert [hit for _payload, hit in loaded] == [False, False]


def test_population_lock_serializes_same_profile_entry(tmp_path: Path) -> None:
    store_root = tmp_path / "store"
    waiting = threading.Event()
    acquired = threading.Event()

    def contender() -> None:
        waiting.set()
        with store._entry_population_lock(store_root, "1abc"):
            acquired.set()

    with store._entry_population_lock(store_root, "1abc"):
        thread = threading.Thread(target=contender)
        thread.start()
        assert waiting.wait(timeout=1)
        time.sleep(0.1)
        assert not acquired.is_set()
    thread.join(timeout=2)
    assert acquired.is_set()


def test_aggregate_chain_limit_is_enforced_while_loading(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    store_root = tmp_path / "store"
    for pdb_id in ("1ABC", "2XYZ"):
        source = tmp_path / f"{pdb_id.lower()}.pdb"
        source.write_text(f"HEADER TEST {pdb_id}\n", encoding="utf-8")
        store.write_entry(
            store_root,
            source,
            _raw_structure(source, pdb_id),
            _raw_edges(source, pdb_id),
        )
    monkeypatch.setattr(store, "MAX_REQUEST_CHAINS", 3)

    with pytest.raises(InputValidationError) as exc_info:
        store._load_or_populate_entries(
            store_root,
            ["1abc", "2xyz"],
            source_dir=None,
            populate_missing=False,
        )

    assert exc_info.value.code == "PRECOMPUTED_REQUEST_LIMIT_EXCEEDED"


def test_ambiguous_archive_sources_are_rejected(tmp_path: Path) -> None:
    archive = tmp_path / "archive"
    (archive / "ab").mkdir(parents=True)
    (archive / "1abc.cif").write_text("data_1abc\n", encoding="utf-8")
    (archive / "ab" / "1abc.cif.gz").write_bytes(gzip.compress(b"data_1abc\n"))

    with pytest.raises(InputValidationError, match="PDB_ARCHIVE_ENTRY_AMBIGUOUS"):
        store.resolve_archive_source(archive, "1ABC")


def test_corrupt_cache_edge_is_rejected(tmp_path: Path) -> None:
    store_root, _, payload = _write_fixture_entry(tmp_path)
    path = store.entry_path(store_root, "1abc")
    payload["interactions"][0]["chain_b"] = "1ABC:DOES_NOT_EXIST"
    path.write_bytes(gzip.compress(json.dumps(payload).encode("utf-8"), mtime=0))

    with pytest.raises(InputValidationError, match="endpoints are invalid"):
        store.load_entry(store_root, "1abc")


def test_missing_assemble_writes_failed_contract_summary(tmp_path: Path) -> None:
    web_output = tmp_path / "web"
    with pytest.raises(InputValidationError, match="PRECOMPUTED_ENTRY_MISSING"):
        store.run_assemble_pipeline(
            tmp_path / "store",
            ["1abc"],
            web_output_dir=str(web_output),
        )

    summary = json.loads((web_output / "summary.json").read_text(encoding="utf-8"))
    assert summary["output_contract_version"] == "1.1"
    assert summary["status"] == "failed"
    assert summary["errors"][0]["code"] == "PRECOMPUTED_ENTRY_MISSING"


def test_assemble_rejects_detailed_interaction_export(tmp_path: Path) -> None:
    config["export_detailed_interactions"] = True
    with pytest.raises(InputValidationError, match="DETAILED_INTERACTIONS_UNAVAILABLE"):
        store.run_assemble_pipeline(tmp_path / "store", ["1abc"])


def test_cli_help_exposes_precompute_and_assemble_contract(capsys) -> None:
    with pytest.raises(SystemExit) as exc_info:
        main(["precompute", "--help"])
    assert exc_info.value.code == 0
    precompute_help = capsys.readouterr().out
    assert "--input-dir" in precompute_help
    assert "--store" in precompute_help
    assert "--recursive" in precompute_help

    with pytest.raises(SystemExit) as exc_info:
        main(["assemble", "--help"])
    assert exc_info.value.code == 0
    assemble_help = capsys.readouterr().out
    assert "--pdb-id" in assemble_help
    assert "--web-output-dir" in assemble_help
    assert "--source-dir" in assemble_help
    assert "--populate-missing" in assemble_help


def test_cli_precompute_dispatches_without_changing_run_flags(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    called = {}

    def fake_precompute(store_path, input_dir, *, recursive):
        called.update(store=store_path, input=input_dir, recursive=recursive)
        return {
            "written": 1,
            "cache_hits": 0,
            "failed": 0,
            "profile_id": "a" * 64,
            "errors": [],
        }

    monkeypatch.setattr(store, "precompute_directory", fake_precompute)
    exit_code = main(
        [
            "precompute",
            "--input-dir",
            str(tmp_path / "input"),
            "--store",
            str(tmp_path / "store"),
            "--recursive",
            "--headless",
        ]
    )
    assert exit_code == 0
    assert called["recursive"] is True

    # The established local command keeps the same required input/output flags.
    with pytest.raises(SystemExit) as exc_info:
        main(["run", "--help"])
    assert exc_info.value.code == 0


def test_cli_assemble_passes_repeated_ids_and_lazy_options(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    called = {}
    summary = tmp_path / "summary.json"

    def fake_assemble(store_path, pdb_ids, **kwargs):
        called.update(store=store_path, ids=pdb_ids, kwargs=kwargs)
        return SimpleNamespace(run_output_path=str(tmp_path / "run"), summary_file=str(summary))

    monkeypatch.setattr(store, "run_assemble_pipeline", fake_assemble)
    exit_code = main(
        [
            "assemble",
            "--store",
            str(tmp_path / "store"),
            "--pdb-id",
            "1abc",
            "--pdb-id",
            "2XYZ",
            "--output-dir",
            str(tmp_path / "output"),
            "--source-dir",
            str(tmp_path / "archive"),
            "--populate-missing",
            "--headless",
        ]
    )
    assert exit_code == 0
    assert called["ids"] == ["1abc", "2XYZ"]
    assert called["kwargs"]["populate_missing"] is True
