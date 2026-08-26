from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from pdb2net import precomputed
from pdb2net.__main__ import main
from pdb2net.config_loader import config
from pdb2net.input_contract import InputValidationError
from pdb2net.precomputed import assemble, build, io, schema
from pdb2net.structure_identity import ChainIdentity, identity_from_official_id


@pytest.fixture(autouse=True)
def _stable_precomputed_config(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    pdb_fasta = tmp_path / "pdb_seqres.txt"
    pdb_fasta.write_text(
        ">1abc_A mol:protein length:1\nA\n"
        ">1abc_B mol:protein length:1\nG\n"
        ">2xyz_A mol:protein length:1\nA\n"
        ">2xyz_B mol:protein length:1\nG\n",
        encoding="utf-8",
    )
    monkeypatch.setitem(config, "pdb_fasta_path", str(pdb_fasta))
    monkeypatch.setitem(config, "reference_manifest_id", "test-references-v1")
    monkeypatch.setitem(config, "structure_model_policy", "first")
    monkeypatch.setitem(config, "export_detailed_interactions", False)
    monkeypatch.setitem(config, "open_in_cytoscape", False)
    monkeypatch.setitem(config, "output_path", str(tmp_path / "outputs"))
    monkeypatch.setitem(
        config,
        "distance_thresholds",
        {"ca_radius": 12.0, "all_atoms_radius": 5.0},
    )
    monkeypatch.setitem(
        config,
        "interaction_filters",
        {
            "protein_protein_min_ca_neighbors": 10,
            "protein_protein_min_all_atom_contacts": 1,
            "protein_nucleic_acid_min_all_atom_contacts": 1,
            "nucleic_acid_min_all_atom_contacts": 1,
        },
    )
    monkeypatch.setitem(
        config,
        "network_annotations",
        {
            "use_embedded_sifts": True,
            "tooltip_fields": ["uniprot"],
            "max_tooltip_segments_per_database": 20,
        },
    )


def _raw_structure(source: Path, pdb_id: str = "1ABC") -> dict:
    identity = identity_from_official_id(pdb_id)
    embedded_counts = {"uniprot": 0, "pfam": 0, "cath": 0, "scop2": 0}
    return {
        "file_path": str(source),
        "pdb_id": identity.display_id,
        "structure_identity": identity.as_dict(),
        "input_format": "pdb",
        "input_kind": "pdb",
        "embedded_annotation_counts": embedded_counts,
        "atom_data": [
            {
                "chain_id": "A",
                "unique_chain_id": f"{identity.display_id}:A",
                "chain_identity": ChainIdentity(
                    identity.canonical_id, identity.display_id, "A", 1
                ).as_dict(),
                "structure_key": identity.canonical_id,
                "structure_display_id": identity.display_id,
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
                        "atoms": [
                            {"atom_name": "CA", "coordinates": [0.0, 0.0, 0.0]}
                        ],
                    }
                ],
            },
            {
                "chain_id": "B",
                "unique_chain_id": f"{identity.display_id}:B",
                "chain_identity": ChainIdentity(
                    identity.canonical_id, identity.display_id, "B", 1
                ).as_dict(),
                "structure_key": identity.canonical_id,
                "structure_display_id": identity.display_id,
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
                        "atoms": [
                            {"atom_name": "CA", "coordinates": [1.0, 0.0, 0.0]}
                        ],
                    }
                ],
            },
        ],
    }


def _raw_edges(source: Path, pdb_id: str = "1ABC") -> list[dict]:
    identity = identity_from_official_id(pdb_id)
    return [
        {
            "file_path": str(source),
            "structure_key": identity.canonical_id,
            "model_index": 1,
            "chain_a": f"{identity.display_id}:A",
            "chain_b": f"{identity.display_id}:B",
            "ca_neighbors": 11,
            "all_atoms_count": 17,
            "interaction_type": "Protein-Protein",
        }
    ]


def _source(tmp_path: Path, pdb_id: str) -> Path:
    path = tmp_path / f"{pdb_id.lower()}.pdb"
    path.write_text(f"HEADER TEST {pdb_id.upper()}\n", encoding="utf-8")
    return path


def _install_fake_science(
    monkeypatch: pytest.MonkeyPatch,
    *,
    omit: set[str] | None = None,
) -> None:
    from pdb2net import pipeline, unknown_molecule_uniprot

    omitted = omit or set()
    monkeypatch.setattr(pipeline, "_validate_required_reference_files", lambda: None)

    def fake_parse(paths, **_kwargs):
        structures = []
        for raw_path in paths:
            path = Path(raw_path)
            pdb_id = path.stem.upper()
            canonical = identity_from_official_id(pdb_id).canonical_id
            if canonical not in omitted:
                structures.append(_raw_structure(path, pdb_id))
        return structures

    monkeypatch.setattr(pipeline, "_parse_input_files", fake_parse)
    monkeypatch.setattr(
        unknown_molecule_uniprot, "process_molecule_info", lambda *_args, **_kwargs: None
    )
    monkeypatch.setattr(pipeline, "_run_blast_annotation", lambda _structures: None)
    monkeypatch.setattr(
        build,
        "calculate_distances_with_ckdtree",
        lambda structures: [
            edge
            for structure in structures
            for edge in _raw_edges(
                Path(structure["file_path"]), structure["pdb_id"]
            )
        ],
    )


def _build_store(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    *pdb_ids: str,
) -> tuple[Path, list[Path], dict]:
    _install_fake_science(monkeypatch)
    sources = [_source(tmp_path, pdb_id) for pdb_id in (pdb_ids or ("1abc",))]
    root = tmp_path / "store"
    report = precomputed.precompute_sources(root, sources)
    assert report["failed"] == 0
    return root, sources, report


def _read_raw_entry(root: Path, pdb_id: str = "1abc") -> dict:
    return json.loads(gzip.decompress(precomputed.entry_path(root, pdb_id).read_bytes()))


def _rewrite_entry(root: Path, payload: dict, pdb_id: str = "1abc") -> None:
    precomputed.entry_path(root, pdb_id).write_bytes(
        gzip.compress(schema.canonical_json(payload), mtime=0)
    )


def _store_digest(root: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        digest.update(path.relative_to(root).as_posix().encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def test_profile_is_exact_fixed_shape_and_tooltips_do_not_affect_hash(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    profile = precomputed.scientific_profile()
    assert set(profile) == {
        "artifact_schema_version",
        "scientific_pipeline_version",
        "interaction_pipeline_version",
        "annotation_pipeline_version",
        "source_scope",
        "structure_model_policy",
        "parser",
        "distance_thresholds",
        "interaction_filters",
        "annotations",
    }
    assert profile["artifact_schema_version"] == "3"
    assert profile["distance_thresholds"] == {
        "ca_radius": 12.0,
        "all_atoms_radius": 5.0,
    }
    before = precomputed.profile_id()
    monkeypatch.setitem(
        config["network_annotations"], "tooltip_fields", ["uniprot", "pfam"]
    )
    assert precomputed.profile_id() == before
    monkeypatch.setitem(config, "reference_manifest_id", "test-references-v2")
    assert precomputed.profile_id() != before


@pytest.mark.parametrize(
    ("section", "key", "value"),
    [
        ("distance_thresholds", "ca_radius", 11.0),
        ("distance_thresholds", "ca_radius", "12.0"),
        ("interaction_filters", "protein_protein_min_ca_neighbors", 9),
        ("interaction_filters", "protein_protein_min_ca_neighbors", 10.9),
    ],
)
def test_profile_rejects_nonstandard_scientific_values(
    monkeypatch: pytest.MonkeyPatch, section: str, key: str, value: object
) -> None:
    monkeypatch.setitem(config[section], key, value)
    with pytest.raises(InputValidationError) as error:
        precomputed.scientific_profile()
    assert error.value.code == "PRECOMPUTED_PROFILE_UNSUPPORTED"


def test_successful_build_uses_one_profile_layout_and_manifest_is_last(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    root, sources, report = _build_store(monkeypatch, tmp_path, "1abc")
    assert report["written"] == 1
    assert report["cache_hits"] == 0
    assert precomputed.entry_path(root, "1abc") == (
        root / "entries" / "ab" / "pdb_00001abc.json.gz"
    )
    assert not (root / "profiles").exists()
    manifest = precomputed.load_manifest(root)
    assert set(manifest) == {
        "artifact_schema_version",
        "created_at",
        "producer",
        "profile_id",
        "profile",
        "entry_count",
    }
    assert manifest["entry_count"] == 1
    assert manifest["producer"]["name"] == "pdb2net"

    entry = precomputed.load_entry(root, "1ABC")
    assert set(entry) == {
        "artifact_schema_version",
        "created_at",
        "producer",
        "profile_id",
        "pdb_id",
        "structure_identity",
        "source",
        "geometry",
        "annotations",
        "counts",
    }
    assert set(entry["source"]) == {"sha256", "size_bytes", "scope"}
    serialized = gzip.decompress(precomputed.entry_path(root, "1abc").read_bytes())
    assert str(sources[0].parent).encode() not in serialized
    assert b'"residues"' not in serialized
    assert b'"coordinates"' not in serialized


def test_failed_build_leaves_no_manifest_and_retry_reuses_valid_entry(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    first = _source(tmp_path, "1abc")
    second = _source(tmp_path, "2xyz")
    _install_fake_science(
        monkeypatch, omit={identity_from_official_id("2xyz").canonical_id}
    )
    root = tmp_path / "store"
    failed = precomputed.precompute_sources(root, [first, second])
    assert failed["written"] == 1
    assert failed["failed"] == 1
    assert not (root / "manifest.json").exists()
    first_bytes = precomputed.entry_path(root, "1abc").read_bytes()
    with pytest.raises(InputValidationError) as error:
        precomputed.load_entry(root, "1abc")
    assert error.value.code == "PRECOMPUTED_STORE_UNPUBLISHED"

    _install_fake_science(monkeypatch)
    retried = precomputed.precompute_sources(root, [first, second])
    assert retried["cache_hits"] == 1
    assert retried["written"] == 1
    assert retried["failed"] == 0
    assert precomputed.entry_path(root, "1abc").read_bytes() == first_bytes
    assert precomputed.load_manifest(root)["entry_count"] == 2


def test_unpublished_retry_recomputes_an_entry_when_its_source_hash_changes(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    first = _source(tmp_path, "1abc")
    second = _source(tmp_path, "2xyz")
    _install_fake_science(
        monkeypatch, omit={identity_from_official_id("2xyz").canonical_id}
    )
    root = tmp_path / "store"
    assert precomputed.precompute_sources(root, [first, second])["failed"] == 1
    original = precomputed.entry_path(root, "1abc").read_bytes()
    first.write_text("HEADER CHANGED 1ABC\n", encoding="utf-8")

    _install_fake_science(monkeypatch)
    retried = precomputed.precompute_sources(root, [first, second])
    assert retried["cache_hits"] == 0
    assert retried["written"] == 2
    assert precomputed.entry_path(root, "1abc").read_bytes() != original


def test_published_store_is_immutable(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    root, sources, _report = _build_store(monkeypatch, tmp_path, "1abc")
    with pytest.raises(InputValidationError) as error:
        precomputed.precompute_sources(root, sources)
    assert error.value.code == "PRECOMPUTED_STORE_PUBLISHED"


@pytest.mark.parametrize(
    "relative_path",
    (
        "entries/zz/pdb_00001abc.json.gz",
        "entries/ab/1abc.json.gz",
    ),
)
def test_build_rejects_noncanonical_or_duplicate_entry_paths(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    relative_path: str,
) -> None:
    _install_fake_science(monkeypatch)
    source = _source(tmp_path, "1abc")
    root = tmp_path / "store"
    unexpected = root / relative_path
    unexpected.parent.mkdir(parents=True)
    unexpected.write_bytes(b"unexpected")

    report = precomputed.precompute_sources(root, [source])

    assert report["failed"] == 1
    assert report["errors"][-1]["code"] == "PRECOMPUTED_ENTRY_LAYOUT_INVALID"
    assert not (root / "manifest.json").exists()


def test_raw_and_precomputed_scientific_data_are_equivalent(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    from pdb2net.pipeline import _compact_structure_summaries

    root, sources, _report = _build_store(monkeypatch, tmp_path, "1abc")
    payload = precomputed.load_entry(root, "1abc")
    cached_structure, cached_edges, _references = schema.materialize_entry(payload)
    raw_compact = _compact_structure_summaries(
        [_raw_structure(sources[0], "1ABC")]
    )[0]

    for chain in raw_compact["atom_data"]:
        chain.pop("residues", None)
    assert cached_structure["atom_data"] == raw_compact["atom_data"]
    assert [
        {key: value for key, value in edge.items() if key != "file_path"}
        for edge in cached_edges
    ] == [
        {key: value for key, value in edge.items() if key != "file_path"}
        for edge in _raw_edges(sources[0])
    ]


def test_assemble_is_read_only_and_writes_contract_two_outputs(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    root, _sources, _report = _build_store(monkeypatch, tmp_path, "1abc")
    before = _store_digest(root)
    monkeypatch.setitem(
        config,
        "networks",
        {
            "chain_per_pdb": True,
            "protein_per_pdb": False,
            "combined_chain_network": False,
            "combined_protein_network": False,
        },
    )
    web_output = tmp_path / "web"
    precomputed.run_assemble_pipeline(root, ["1ABC"], web_output_dir=str(web_output))
    assert _store_digest(root) == before
    summary = json.loads((web_output / "summary.json").read_text(encoding="utf-8"))
    assert summary["output_contract_version"] == "2.0"
    assert summary["status"] == "success"
    assert summary["references"]["precomputed_store"] == {
        "artifact_schema_version": "3",
        "profile_id": precomputed.profile_id(),
        "source_scope": "asymmetric_unit",
    }
    assert summary["interactions"] == []


def test_normalization_rejects_duplicate_aliases() -> None:
    assert precomputed.normalize_pdb_id("1ABC") == "pdb_00001abc"
    with pytest.raises(InputValidationError) as error:
        precomputed.normalize_requested_ids(["1ABC", "pdb_00001abc"])
    assert error.value.code == "DUPLICATE_STRUCTURE_IDENTITY"


def test_missing_entry_fails_without_modifying_store(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    root, _sources, _report = _build_store(monkeypatch, tmp_path, "1abc")
    before = _store_digest(root)
    with pytest.raises(InputValidationError) as error:
        precomputed.load_entry(root, "2xyz")
    assert error.value.code == "PRECOMPUTED_ENTRY_MISSING"
    assert _store_digest(root) == before


@pytest.mark.parametrize("content", [b"not-gzip", gzip.compress(b"not-json", mtime=0)])
def test_corrupt_entry_fails_closed(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, content: bytes
) -> None:
    root, _sources, _report = _build_store(monkeypatch, tmp_path, "1abc")
    precomputed.entry_path(root, "1abc").write_bytes(content)
    with pytest.raises(InputValidationError) as error:
        precomputed.load_entry(root, "1abc")
    assert error.value.code == "CORRUPT_PRECOMPUTED_ENTRY"


@pytest.mark.parametrize(
    ("mutation", "code"),
    [
        (lambda value: value.update(artifact_schema_version="2"), "PRECOMPUTED_SCHEMA_MISMATCH"),
        (lambda value: value.update(profile_id="0" * 64), "PRECOMPUTED_PROFILE_MISMATCH"),
        (lambda value: value.update(pdb_id="pdb_00002xyz"), "PRECOMPUTED_PDB_ID_MISMATCH"),
    ],
)
def test_wrong_schema_profile_or_identity_fails_closed(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    mutation,
    code: str,
) -> None:
    root, _sources, _report = _build_store(monkeypatch, tmp_path, "1abc")
    payload = _read_raw_entry(root)
    mutation(payload)
    _rewrite_entry(root, payload)
    with pytest.raises(InputValidationError) as error:
        precomputed.load_entry(root, "1abc")
    assert error.value.code == code


def test_entry_size_and_count_ceilings_fail_closed(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    root, _sources, _report = _build_store(monkeypatch, tmp_path, "1abc")
    entry_size = precomputed.entry_path(root, "1abc").stat().st_size
    monkeypatch.setattr(io, "MAX_COMPRESSED_ENTRY_BYTES", entry_size - 1)
    with pytest.raises(InputValidationError) as size_error:
        precomputed.load_entry(root, "1abc")
    assert size_error.value.code == "PRECOMPUTED_ENTRY_TOO_LARGE"

    monkeypatch.setattr(io, "MAX_COMPRESSED_ENTRY_BYTES", entry_size + 1)
    monkeypatch.setattr(schema, "MAX_CHAINS_PER_ENTRY", 1)
    with pytest.raises(InputValidationError) as count_error:
        precomputed.load_entry(root, "1abc")
    assert count_error.value.code == "CORRUPT_PRECOMPUTED_ENTRY"

    monkeypatch.setattr(schema, "MAX_CHAINS_PER_ENTRY", 50_000)
    monkeypatch.setattr(schema, "MAX_EDGES_PER_ENTRY", 0)
    with pytest.raises(InputValidationError) as edge_error:
        precomputed.load_entry(root, "1abc")
    assert edge_error.value.code == "CORRUPT_PRECOMPUTED_ENTRY"

    monkeypatch.setattr(schema, "MAX_EDGES_PER_ENTRY", 500_000)
    monkeypatch.setattr(io, "MAX_DECOMPRESSED_ENTRY_BYTES", 10)
    with pytest.raises(InputValidationError) as expanded_error:
        precomputed.load_entry(root, "1abc")
    assert expanded_error.value.code == "PRECOMPUTED_ENTRY_TOO_LARGE"


def test_request_aggregate_ceiling_fails_closed(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    root, _sources, _report = _build_store(monkeypatch, tmp_path, "1abc")
    monkeypatch.setattr(assemble, "MAX_REQUEST_CHAINS", 1)
    with pytest.raises(InputValidationError) as error:
        assemble._load_requested_entries(
            root,
            ["pdb_00001abc"],
            manifest=precomputed.load_manifest(root),
            resource_limits={},
        )
    assert error.value.code == "PRECOMPUTED_REQUEST_LIMIT_EXCEEDED"


def test_symlinked_store_or_entry_is_rejected(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    root, _sources, _report = _build_store(monkeypatch, tmp_path, "1abc")
    link = tmp_path / "store-link"
    link.symlink_to(root, target_is_directory=True)
    with pytest.raises(InputValidationError) as root_error:
        precomputed.load_manifest(link)
    assert root_error.value.code == "UNSAFE_PRECOMPUTED_PATH"

    entry = precomputed.entry_path(root, "1abc")
    target = tmp_path / "outside.json.gz"
    target.write_bytes(entry.read_bytes())
    entry.unlink()
    entry.symlink_to(target)
    with pytest.raises(InputValidationError) as entry_error:
        precomputed.load_entry(root, "1abc")
    assert entry_error.value.code in {
        "UNSAFE_PRECOMPUTED_PATH",
        "UNSAFE_PRECOMPUTED_ENTRY",
    }


def test_manifest_profile_mismatch_fails_closed(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    root, _sources, _report = _build_store(monkeypatch, tmp_path, "1abc")
    manifest_path = root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["profile"]["distance_thresholds"]["ca_radius"] = 11.0
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(InputValidationError) as error:
        precomputed.load_manifest(root)
    assert error.value.code == "PRECOMPUTED_PROFILE_MISMATCH"


def test_cli_help_and_dispatch_have_no_lazy_surface(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, capsys
) -> None:
    with pytest.raises(SystemExit) as help_exit:
        main(["assemble", "--help"])
    assert help_exit.value.code == 0
    help_text = capsys.readouterr().out
    assert "--pdb-id" in help_text
    assert "--source-dir" not in help_text
    assert "--populate-missing" not in help_text

    called = {}

    def fake_assemble(store_path, pdb_ids, **kwargs):
        called.update(store=store_path, ids=pdb_ids, kwargs=kwargs)
        return SimpleNamespace(
            run_output_path=str(tmp_path / "run"),
            summary_file=str(tmp_path / "summary.json"),
        )

    monkeypatch.setattr(precomputed, "run_assemble_pipeline", fake_assemble)
    assert main(
        [
            "assemble",
            "--store",
            str(tmp_path / "store"),
            "--pdb-id",
            "1abc",
            "--pdb-id",
            "2xyz",
            "--output-dir",
            str(tmp_path / "output"),
            "--headless",
        ]
    ) == 0
    assert called["ids"] == ["1abc", "2xyz"]
    assert called["kwargs"] == {"web_output_dir": None}
