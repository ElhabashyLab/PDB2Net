from __future__ import annotations

import gzip
import os
from pathlib import Path

import pytest

from pdb2net import distances, file_parser, pipeline, precomputed_store
from pdb2net.config_loader import config
from pdb2net.input_contract import InputValidationError
from pdb2net.structure_identity import (
    StructureIdentity,
    canonical_pdb_id,
    legacy_pdb_id,
    pdb_archive_shard,
)


def _assert_code(error: pytest.ExceptionInfo[InputValidationError], code: str) -> None:
    assert error.value.code == code


def _minimal_mmcif(identity: str) -> bytes:
    return f"data_{identity}\n_entry.id {identity}\n".encode()


def test_official_id_grammar_normalization_and_archive_shard_are_exact() -> None:
    assert canonical_pdb_id("1AbC") == "pdb_00001abc"
    assert canonical_pdb_id("PDB_Ab12Cd34") == "pdb_ab12cd34"
    assert legacy_pdb_id("pdb_00001abc") == "1abc"
    assert legacy_pdb_id("pdb_ab12cd34") is None
    assert pdb_archive_shard("pdb_ab12cd34") == "d3"

    for invalid in ("AABC", "123", "12345", "pdb_1234567", "pdb_123456789", "pdb_1234-678"):
        with pytest.raises(ValueError):
            canonical_pdb_id(invalid)

    for unsafe in ("/private/input", r"C:\\private\\input", "bad\x00id"):
        with pytest.raises(ValueError):
            StructureIdentity("local", unsafe)


def test_precomputed_requests_reject_legacy_extended_identity_duplicates() -> None:
    with pytest.raises(InputValidationError) as error:
        precomputed_store.normalize_requested_ids(["1ABC", "pdb_00001abc"])
    _assert_code(error, "DUPLICATE_STRUCTURE_IDENTITY")


def test_equivalent_content_identity_claims_coalesce(tmp_path: Path) -> None:
    source = tmp_path / "pdb_00001abc.cif"
    source.write_text(
        "data_1abc\n"
        "_entry.id pdb_00001ABC\n"
        "loop_\n"
        "_database_2.database_id\n"
        "_database_2.database_code\n"
        "PDB 1AbC\n",
        encoding="utf-8",
    )

    identity, warnings, claims = file_parser.resolve_structure_identity(str(source))

    assert identity.canonical_id == "pdb_00001abc"
    assert warnings == []
    assert {claim["source"] for claim in claims} == {
        "data_block",
        "_entry.id",
        "_database_2[PDB].database_code",
        "filename",
    }


def test_conflicting_content_claims_fail_but_filename_mismatch_is_content_first(
    tmp_path: Path,
) -> None:
    conflict = tmp_path / "conflict.cif"
    conflict.write_text("data_1abc\n_entry.id 2xyz\n", encoding="utf-8")
    with pytest.raises(InputValidationError) as error:
        file_parser.resolve_structure_identity(str(conflict))
    _assert_code(error, "CONFLICTING_STRUCTURE_IDENTITY")

    mismatch = tmp_path / "2xyz.cif"
    mismatch.write_bytes(_minimal_mmcif("1abc"))
    identity, warnings, _claims = file_parser.resolve_structure_identity(str(mismatch))
    assert identity.canonical_id == "pdb_00001abc"
    assert [warning["code"] for warning in warnings] == ["STRUCTURE_FILENAME_ID_MISMATCH"]


@pytest.mark.parametrize(
    ("left_name", "left_content", "right_name", "right_content"),
    [
        ("first.cif", _minimal_mmcif("1abc"), "second.cif", _minimal_mmcif("pdb_00001abc")),
        ("local-one.cif", _minimal_mmcif("local_sample"), "local-two.cif", _minimal_mmcif("LOCAL_SAMPLE")),
        (
            "af-one.cif",
            _minimal_mmcif("AF_Q9BYF1_F1_MODEL_V4"),
            "af-two.cif",
            _minimal_mmcif("af_q9byf1_f1_model_v4"),
        ),
    ],
)
def test_duplicate_official_local_and_alphafold_identities_are_rejected(
    tmp_path: Path,
    left_name: str,
    left_content: bytes,
    right_name: str,
    right_content: bytes,
) -> None:
    left = tmp_path / left_name
    right = tmp_path / right_name
    left.write_bytes(left_content)
    right.write_bytes(right_content)
    paths = [str(left), str(right)]
    inventory = pipeline.inspect_input_files(paths)

    with pytest.raises(InputValidationError) as error:
        pipeline._preflight_structure_identities(paths, inventory)

    _assert_code(error, "DUPLICATE_STRUCTURE_IDENTITY")


def test_core_and_directory_modes_reject_structure_symlinks(tmp_path: Path) -> None:
    target = tmp_path / "target.cif"
    target.write_bytes(_minimal_mmcif("1abc"))
    link = tmp_path / "linked.cif"
    link.symlink_to(target)

    with pytest.raises(InputValidationError) as error:
        pipeline.inspect_input_files([str(link)])
    _assert_code(error, "SYMLINK_INPUT_NOT_ALLOWED")

    with pytest.raises(InputValidationError) as error:
        pipeline.discover_input_files(str(tmp_path))
    _assert_code(error, "SYMLINK_INPUT_NOT_ALLOWED")


def test_plain_and_gzip_magic_suffix_must_match(tmp_path: Path) -> None:
    fake_gzip = tmp_path / "plain.cif.gz"
    fake_gzip.write_bytes(_minimal_mmcif("1abc"))
    disguised = tmp_path / "compressed.cif"
    disguised.write_bytes(gzip.compress(_minimal_mmcif("1abc"), mtime=0))

    for source in (fake_gzip, disguised):
        with pytest.raises(InputValidationError) as error:
            file_parser.read_validated_structure_bytes(str(source))
        _assert_code(error, "INPUT_COMPRESSION_MISMATCH")


def test_empty_plain_and_empty_gzip_payloads_have_stable_errors(tmp_path: Path) -> None:
    plain = tmp_path / "empty.cif"
    plain.write_bytes(b"")
    with pytest.raises(InputValidationError) as error:
        file_parser.read_validated_structure_bytes(str(plain))
    _assert_code(error, "EMPTY_INPUT_FILE")

    compressed = tmp_path / "empty.cif.gz"
    compressed.write_bytes(gzip.compress(b"", mtime=0))
    with pytest.raises(InputValidationError) as error:
        file_parser.read_validated_structure_bytes(str(compressed))
    _assert_code(error, "INVALID_GZIP_INPUT")


def test_gzip_crc_nested_and_trailing_data_are_rejected(tmp_path: Path) -> None:
    payload = _minimal_mmcif("1abc")
    corrupt = bytearray(gzip.compress(payload, mtime=0))
    corrupt[-8] ^= 0x80
    cases = {
        "crc.cif.gz": bytes(corrupt),
        "nested.cif.gz": gzip.compress(gzip.compress(payload, mtime=0), mtime=0),
        "trailing.cif.gz": gzip.compress(payload, mtime=0) + b"not-a-member",
    }
    expected = {
        "crc.cif.gz": "INVALID_GZIP_INPUT",
        "nested.cif.gz": "NESTED_GZIP_INPUT",
        "trailing.cif.gz": "INVALID_GZIP_INPUT",
    }
    for name, raw in cases.items():
        source = tmp_path / name
        source.write_bytes(raw)
        with pytest.raises(InputValidationError) as error:
            file_parser.read_validated_structure_bytes(str(source))
        _assert_code(error, expected[name])


@pytest.mark.parametrize("removed", [1, 5, 8])
def test_every_truncated_gzip_trailer_is_rejected(tmp_path: Path, removed: int) -> None:
    source = tmp_path / f"truncated-{removed}.cif.gz"
    source.write_bytes(gzip.compress(_minimal_mmcif("1abc"), mtime=0)[:-removed])

    with pytest.raises(InputValidationError) as error:
        file_parser.read_validated_structure_bytes(str(source))

    _assert_code(error, "INVALID_GZIP_INPUT")


def test_concatenated_gzip_members_share_one_exact_expansion_limit(tmp_path: Path) -> None:
    payload = _minimal_mmcif("1abc")
    split = 7
    source = tmp_path / "members.cif.gz"
    source.write_bytes(
        gzip.compress(payload[:split], mtime=0) + gzip.compress(payload[split:], mtime=0)
    )

    expanded, size = file_parser.read_validated_structure_bytes(
        str(source), maximum_expanded_bytes=len(payload)
    )
    assert expanded == payload
    assert size == len(payload)
    with pytest.raises(InputValidationError) as error:
        file_parser.read_validated_structure_bytes(
            str(source), maximum_expanded_bytes=len(payload) - 1
        )
    _assert_code(error, "INPUT_FILE_EXPANDED_BYTES_LIMIT_EXCEEDED")


def test_concatenated_transport_cannot_smuggle_multiple_mmcif_blocks(tmp_path: Path) -> None:
    first = b"data_1abc\n_entry.id 1abc\n"
    second = b"data_2xyz\n_entry.id 2xyz\n"
    source = tmp_path / "two-blocks.cif.gz"
    source.write_bytes(gzip.compress(first, mtime=0) + gzip.compress(second, mtime=0))

    with pytest.raises(InputValidationError) as error:
        file_parser.parse_structure_input(str(source))

    _assert_code(error, "INVALID_MMCIF_DATA_BLOCK_COUNT")


def test_input_signature_change_is_detected_after_read(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    source = tmp_path / "stable.cif"
    source.write_bytes(_minimal_mmcif("1abc"))
    real_signature = file_parser.input_file_signature
    calls = 0

    def changing_signature(path: str | Path):
        nonlocal calls
        calls += 1
        signature = real_signature(path)
        if calls >= 2:
            return signature[:-1] + (signature[-1] + 1,)
        return signature

    monkeypatch.setattr(file_parser, "input_file_signature", changing_signature)
    with pytest.raises(InputValidationError) as error:
        file_parser.read_validated_structure_bytes(str(source))
    _assert_code(error, "INPUT_CHANGED_DURING_PROCESSING")


def test_same_length_input_rewrite_with_restored_mtime_changes_digest(tmp_path: Path) -> None:
    source = tmp_path / "stable.cif"
    original = _minimal_mmcif("1abc")
    replacement = _minimal_mmcif("2xyz")
    assert len(original) == len(replacement)
    source.write_bytes(original)
    before_stat = source.stat()
    expected_signature = file_parser.input_file_signature(source)
    expected_sha256 = file_parser.input_file_sha256(
        source,
        expected_signature=expected_signature,
    )

    source.write_bytes(replacement)
    os.utime(source, ns=(before_stat.st_atime_ns, before_stat.st_mtime_ns))

    assert file_parser.input_file_sha256(source) != expected_sha256
    with pytest.raises(InputValidationError) as error:
        file_parser.read_validated_structure_bytes(
            str(source),
            expected_signature=expected_signature,
            expected_sha256=expected_sha256,
        )
    _assert_code(error, "INPUT_CHANGED_DURING_PROCESSING")


def test_parse_worker_enforces_forwarded_gzip_expansion_limit(tmp_path: Path) -> None:
    source = tmp_path / "bounded.cif.gz"
    source.write_bytes(gzip.compress(_minimal_mmcif("1abc"), mtime=0))
    signature = file_parser.input_file_signature(source)

    with pytest.raises(InputValidationError) as error:
        pipeline.process_single_file(
            str(source),
            expected_signature=signature,
            maximum_expanded_bytes=8,
        )
    _assert_code(error, "INPUT_FILE_EXPANDED_BYTES_LIMIT_EXCEEDED")


def test_duplicate_public_chain_node_ids_are_rejected_across_batches() -> None:
    seen: dict[str, str] = {}
    pipeline._register_unique_chain_ids(
        [
            {
                "file_path": "/inputs/first.cif",
                "atom_data": [{"unique_chain_id": "FOO:BAR:X"}],
            }
        ],
        seen,
    )

    with pytest.raises(InputValidationError) as error:
        pipeline._register_unique_chain_ids(
            [
                {
                    "file_path": "/inputs/second.cif",
                    "atom_data": [{"unique_chain_id": "FOO:BAR:X"}],
                }
            ],
            seen,
        )
    _assert_code(error, "DUPLICATE_CHAIN_NODE_ID")


def _protein_chain(unique_id: str, chain_id: str, model_index: int, offset: float) -> dict:
    return {
        "chain_id": chain_id,
        "unique_chain_id": unique_id,
        "structure_key": "local:TST",
        "model_index": model_index,
        "molecule_type": "Protein",
        "molecule_name": f"Protein {chain_id}",
        "uniprot_id": f"P0000{model_index}{chain_id}",
        "residues": [
            {
                "residue_name": "ALA",
                "residue_number": 1,
                "atoms": [{"atom_name": "CA", "coordinates": [offset, 0.0, 0.0]}],
            }
        ],
    }


def test_all_model_chain_interactions_are_kept_and_exported_separately(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setitem(
        distances.config,
        "distance_thresholds",
        {"ca_radius": 2.0, "all_atoms_radius": 2.0},
    )
    monkeypatch.setitem(
        distances.config,
        "interaction_filters",
        {
            "protein_protein_min_ca_neighbors": 1,
            "protein_protein_min_all_atom_contacts": 1,
            "protein_nucleic_acid_min_all_atom_contacts": 1,
            "nucleic_acid_min_all_atom_contacts": 1,
        },
    )
    distances.coords_cache.clear()
    distances.tree_cache.clear()
    chains = [
        _protein_chain("TST:model1:A", "A", 1, 0.0),
        _protein_chain("TST:model1:B", "B", 1, 1.0),
        _protein_chain("TST:model2:A", "A", 2, 10.0),
        _protein_chain("TST:model2:B", "B", 2, 11.0),
    ]
    structure = {
        "file_path": "/inputs/tst.cif",
        "pdb_id": "TST",
        "structure_identity": {
            "source": "local",
            "canonical_id": "TST",
            "legacy_id": None,
            "key": "local:TST",
            "display_id": "TST",
        },
        "atom_data": chains,
    }

    edges = distances.calculate_distances_with_ckdtree([structure])

    assert [edge["model_index"] for edge in edges] == [1, 2]
    assert [(edge["chain_a"], edge["chain_b"]) for edge in edges] == [
        ("TST:model1:A", "TST:model1:B"),
        ("TST:model2:A", "TST:model2:B"),
    ]

    exports: list[tuple[str, list[dict]]] = []
    monkeypatch.setitem(config, "structure_model_policy", "all")
    monkeypatch.setattr(pipeline, "generate_nodes_from_atom_data", lambda *_args, **_kwargs: [])
    monkeypatch.setattr(
        pipeline,
        "create_cytoscape_network",
        lambda model_edges, title, _path, nodes_data=None: exports.append((title, model_edges)),
    )
    pipeline._export_chain_networks([structure], edges, "/tmp/unused")

    assert [title for title, _edges in exports] == [
        "Chain_Interaction_Network_TST_model1",
        "Chain_Interaction_Network_TST_model2",
    ]
    assert [len(model_edges) for _title, model_edges in exports] == [1, 1]


def test_all_model_policy_rejects_every_protein_network_before_analysis(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setitem(config, "structure_model_policy", "all")
    monkeypatch.setitem(
        config,
        "networks",
        {
            "chain_per_pdb": True,
            "protein_per_pdb": True,
            "combined_chain_network": False,
            "combined_protein_network": False,
        },
    )

    with pytest.raises(InputValidationError) as error:
        pipeline._validate_analysis_config()

    _assert_code(error, "PROTEIN_NETWORKS_UNSUPPORTED_FOR_ALL_MODELS")
