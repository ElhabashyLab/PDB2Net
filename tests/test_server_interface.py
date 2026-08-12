from __future__ import annotations

import argparse
import gzip
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

import pdb2net
from pdb2net import (
    artifact_names,
    components,
    config_loader,
    cx2_export,
    data_processor,
    detailed_results_exporter,
    distances,
    file_parser,
    network_annotations,
    pipeline,
    precomputed_store,
    protein_network,
    uniprot_matcher,
)
from pdb2net.__main__ import build_parser
from pdb2net.capabilities import capability_document
from pdb2net.config_loader import config
from pdb2net.contracts import (
    CAPABILITIES_SCHEMA_VERSION,
    OUTPUT_CONTRACT_VERSION,
    PRECOMPUTED_ARTIFACT_SCHEMA_VERSION,
)
from pdb2net.input_contract import InputValidationError
from pdb2net.outputs import collect_web_outputs, create_run_output_paths, write_run_manifest
from pdb2net.server_interface import (
    ALLOWED_INTERACTION_TYPES,
    ANNOTATION_DATABASES,
    ANNOTATION_PIPELINE_VERSION,
    COMBINED_COMPONENT_SEMANTICS,
    COMMAND_SEMANTICS_IDS,
    CX2_DECLARATION_SCOPES,
    CX2_HEADER,
    CX2_REQUIRED_ASPECT_ORDER,
    CX2_SUCCESS_STATUS,
    DETAILED_INTERACTION_COLUMNS,
    DETAILED_INTERACTION_FILENAME_SUFFIX,
    DISTANCE_THRESHOLD_RULES,
    INTERACTION_PIPELINE_VERSION,
    INTERACTION_FILTER_RULES,
    MAX_ARTIFACT_STEM_BYTES,
    NETWORK_OUTPUT_FIELDS,
    NETWORK_TITLE_BASES,
    PARSER_SEMANTICS,
    PORTABLE_ARTIFACT_STEM_SEMANTICS_ID,
    PRECOMPUTED_ANNOTATION_PIPELINE_VERSION,
    PRECOMPUTED_PARSER_SEMANTICS,
    PRECOMPUTED_SCIENTIFIC_PIPELINE_VERSION,
    PRECOMPUTED_SOURCE_SCOPE,
    PRECOMPUTED_ENTRY_FIELDS,
    RESOURCE_LIMIT_FIELDS,
    SERVER_ENVIRONMENT,
    UPLOAD_SOURCE_SCOPE,
    WEB_OUTPUT_SUMMARY_FIELDS,
    WEB_OUTPUT_VALIDATION_SEMANTICS,
    WEB_SERVER_INPUT_SUFFIXES,
)
from pdb2net.structure_identity import ChainIdentity, identity_from_official_id


ROOT = Path(__file__).resolve().parents[1]


def _subcommand_parsers() -> dict[str, argparse.ArgumentParser]:
    parser = build_parser()
    action = next(
        item
        for item in parser._actions
        if isinstance(item, argparse._SubParsersAction)
    )
    return dict(action.choices)


def _option_action(parser: argparse.ArgumentParser, option: str) -> argparse.Action:
    return next(action for action in parser._actions if option in action.option_strings)


def test_capability_v2_server_interface_is_json_native_and_configuration_free(
    tmp_path: Path,
) -> None:
    document = capability_document()

    assert document["capabilities_schema_version"] == CAPABILITIES_SCHEMA_VERSION == "2"
    assert set(document["server_interface"]) == {
        "commands",
        "contracts",
        "scientific_profiles",
    }
    assert json.loads(json.dumps(document["server_interface"])) == document["server_interface"]

    environment = dict(os.environ)
    environment["PDB2NET_CONFIG_FILE"] = str(tmp_path / "does-not-exist.json")
    result = subprocess.run(
        [sys.executable, "-m", "pdb2net", "capabilities", "--json"],
        cwd=ROOT,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["server_interface"] == document["server_interface"]


def test_declared_server_command_options_match_the_real_argparse_parser() -> None:
    commands = capability_document()["server_interface"]["commands"]
    parsers = _subcommand_parsers()

    assert set(commands) == {"capabilities", "run", "precompute", "assemble"}
    assert {name: value["semantics_id"] for name, value in commands.items()} == (
        COMMAND_SEMANTICS_IDS
    )
    assert len(set(COMMAND_SEMANTICS_IDS.values())) == len(COMMAND_SEMANTICS_IDS)
    for command, descriptor in commands.items():
        assert descriptor["success_exit_code"] == 0
        parser = parsers[command]
        for option, expected in descriptor["options"].items():
            action = _option_action(parser, option)
            actual_type = (
                "boolean"
                if isinstance(action, argparse._StoreTrueAction)
                else "string"
            )
            cardinality = (
                "one_or_more"
                if isinstance(action, argparse._AppendAction)
                else "one"
            )
            choices = list(action.choices) if action.choices is not None else None
            assert actual_type == expected["type"], (command, option)
            assert bool(action.required) is expected["required"], (command, option)
            assert cardinality == expected["cardinality"], (command, option)
            assert choices == expected["choices"], (command, option)


def test_server_interface_constants_are_used_by_runtime_validators() -> None:
    assert pipeline.NETWORK_OUTPUT_FIELDS is NETWORK_OUTPUT_FIELDS
    assert pipeline.DISTANCE_THRESHOLD_RULES is DISTANCE_THRESHOLD_RULES
    assert pipeline.INTERACTION_FILTER_RULES is INTERACTION_FILTER_RULES
    assert pipeline.RESOURCE_LIMIT_KEYS is RESOURCE_LIMIT_FIELDS
    assert network_annotations.ANNOTATION_DATABASES is ANNOTATION_DATABASES
    assert precomputed_store.SCIENTIFIC_PIPELINE_VERSION == PRECOMPUTED_SCIENTIFIC_PIPELINE_VERSION
    assert precomputed_store.ANNOTATION_PIPELINE_VERSION == PRECOMPUTED_ANNOTATION_PIPELINE_VERSION
    assert precomputed_store.SOURCE_SCOPE == PRECOMPUTED_SOURCE_SCOPE
    assert config_loader.SERVER_ENVIRONMENT is SERVER_ENVIRONMENT
    assert detailed_results_exporter.DETAILED_INTERACTION_COLUMNS is DETAILED_INTERACTION_COLUMNS
    assert precomputed_store.ALLOWED_INTERACTION_TYPES is ALLOWED_INTERACTION_TYPES
    assert artifact_names.MAX_ARTIFACT_STEM_BYTES == MAX_ARTIFACT_STEM_BYTES
    assert (
        artifact_names.PORTABLE_ARTIFACT_STEM_SEMANTICS_ID
        == PORTABLE_ARTIFACT_STEM_SEMANTICS_ID
    )
    assert cx2_export.CX2_HEADER is CX2_HEADER
    assert cx2_export.CX2_REQUIRED_ASPECT_ORDER is CX2_REQUIRED_ASPECT_ORDER
    assert cx2_export.CX2_DECLARATION_SCOPES is CX2_DECLARATION_SCOPES
    assert cx2_export.CX2_SUCCESS_STATUS is CX2_SUCCESS_STATUS
    assert pipeline.NETWORK_TITLE_BASES is NETWORK_TITLE_BASES
    assert protein_network.NETWORK_TITLE_BASES is NETWORK_TITLE_BASES
    assert data_processor.PARSER_SEMANTICS == PARSER_SEMANTICS
    assert distances.INTERACTION_PIPELINE_VERSION == INTERACTION_PIPELINE_VERSION
    assert network_annotations.ANNOTATION_PIPELINE_VERSION == ANNOTATION_PIPELINE_VERSION
    assert uniprot_matcher.ANNOTATION_PIPELINE_VERSION == ANNOTATION_PIPELINE_VERSION
    assert components.COMBINED_COMPONENT_SEMANTICS == COMBINED_COMPONENT_SEMANTICS


@pytest.mark.parametrize("model_indices", ([2], [1, 3]))
def test_all_model_validation_semantics_match_noncontiguous_runtime_exports(
    monkeypatch: pytest.MonkeyPatch, model_indices: list[int]
) -> None:
    selected = WEB_OUTPUT_VALIDATION_SEMANTICS["selected_outputs"]
    assert selected == {
        "first_model_per_identity_networks_cover_each_identity_exactly_once": True,
        "all_model_chain_networks_cover_each_identity_at_least_once": True,
        "all_model_chain_network_identity_model_pairs_are_unique": True,
        "all_model_chain_network_model_indices_need_not_be_contiguous": True,
        "interaction_csvs_cover_each_identity_exactly_once": True,
        "chain_network_node_sum_equals_summary_chains": True,
    }

    titles: list[str] = []
    monkeypatch.setitem(config, "structure_model_policy", "all")
    monkeypatch.setattr(
        pipeline, "generate_nodes_from_atom_data", lambda *_args, **_kwargs: []
    )
    monkeypatch.setattr(
        pipeline,
        "create_cytoscape_network",
        lambda _edges, title, _path, nodes_data=None: titles.append(title),
    )
    pipeline._export_chain_networks(
        [
            {
                "pdb_id": "TST1",
                "structure_identity": {"key": "local:TST1"},
                "atom_data": [
                    {
                        "model_index": model_index,
                        "unique_chain_id": f"TST1:model{model_index}:A",
                    }
                    for model_index in model_indices
                ],
            }
        ],
        [],
        "/tmp/unused",
    )

    assert titles == [
        f"Chain_Interaction_Network_TST1_model{model_index}"
        for model_index in model_indices
    ]


def test_declared_server_environment_is_the_real_config_override_map(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    declared = capability_document()["server_interface"]["contracts"]["environment"]

    assert declared == SERVER_ENVIRONMENT
    values: dict[str, str] = {}
    for index, (name, descriptor) in enumerate(SERVER_ENVIRONMENT.items(), start=1):
        value = f"/contract/value-{index}"
        monkeypatch.setenv(name, value)
        values[descriptor["config_path"]] = value
        assert descriptor["type"] in {"path", "executable"}
        assert descriptor["supply_condition"] == "real_worker"

    loaded: dict[str, object] = {}
    config_loader._apply_env_overrides(loaded)
    for path, value in values.items():
        assert loaded[path] == value


def test_declared_web_input_suffixes_and_single_gzip_layer_match_real_parser(
    tmp_path: Path,
) -> None:
    document = capability_document()
    declared = document["server_interface"]["contracts"]["input"]

    assert declared == {
        "supported_suffixes": list(WEB_SERVER_INPUT_SUFFIXES),
        "suffix_matching": "case_insensitive",
        "gzip": {"optional": True, "maximum_layers": 1},
    }
    assert document["input_formats"] == list(WEB_SERVER_INPUT_SUFFIXES)
    for suffix in WEB_SERVER_INPUT_SUFFIXES:
        assert file_parser.is_valid_file(f"tiny{suffix}")
        assert file_parser.is_valid_file(f"tiny{suffix.upper()}")
    for suffix in WEB_SERVER_INPUT_SUFFIXES[3:]:
        assert not file_parser.is_valid_file(f"tiny{suffix}.gz")

    nested = tmp_path / "nested.cif.gz"
    nested.write_bytes(gzip.compress(gzip.compress(b"data_tiny\n", mtime=0), mtime=0))
    with pytest.raises(InputValidationError) as error:
        file_parser.read_validated_structure_bytes(str(nested))
    assert error.value.code == "NESTED_GZIP_INPUT"


def test_distributable_scientific_defaults_match_the_declared_upload_profile() -> None:
    base = json.loads(
        (Path(pdb2net.__file__).parent / "configs" / "config.base.json").read_text(
            encoding="utf-8"
        )
    )
    profile = capability_document()["server_interface"]["scientific_profiles"]["upload"]

    assert profile["source_scope"] == UPLOAD_SOURCE_SCOPE
    assert profile["parser_semantics"] == PARSER_SEMANTICS
    assert profile["interaction_pipeline_version"] == INTERACTION_PIPELINE_VERSION
    assert profile["annotation_pipeline_version"] == ANNOTATION_PIPELINE_VERSION
    assert profile["combined_component_semantics"] == COMBINED_COMPONENT_SEMANTICS

    for name, rule in profile["distance_thresholds"].items():
        assert base["distance_thresholds"][name] == rule["default"]
    for name, rule in profile["interaction_filters"].items():
        assert base["interaction_filters"][name] == rule["default"]
    assert base["structure_model_policy"] == profile["structure_model_policy"]["default"]
    assert base["export_detailed_interactions"] is profile["detailed_interactions"]["default"]
    assert base["network_annotations"]["tooltip_fields"] == profile["network_annotations"][
        "tooltip_fields"
    ]["default"]


def test_precomputed_profile_is_content_addressed_not_fixed_to_default_thresholds(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setitem(
        config,
        "distance_thresholds",
        {"ca_radius": 9.5, "all_atoms_radius": 4.25},
    )
    monkeypatch.setitem(
        config,
        "interaction_filters",
        {
            name: int(rule["default"]) + 1
            for name, rule in INTERACTION_FILTER_RULES.items()
        },
    )
    manifest_hash, manifest_root = precomputed_store.ensure_profile_manifest(tmp_path)
    manifest = json.loads((manifest_root / "manifest.json").read_text(encoding="utf-8"))
    declared = capability_document()["server_interface"]["scientific_profiles"][
        "precomputed"
    ]

    assert manifest["profile_id"] == manifest_hash
    assert manifest["artifact_schema_version"] == PRECOMPUTED_ARTIFACT_SCHEMA_VERSION
    assert manifest["profile"]["distance_thresholds"] == {
        "ca_radius": 9.5,
        "all_atoms_radius": 4.25,
    }
    assert manifest["profile"]["parser"]["semantics"] == PRECOMPUTED_PARSER_SEMANTICS
    assert manifest["profile"]["interaction_pipeline_version"] == INTERACTION_PIPELINE_VERSION
    assert declared["interaction_pipeline_version"] == INTERACTION_PIPELINE_VERSION
    assert declared["distance_thresholds"]["content_addressed"] is True
    assert all(
        "default" not in rule
        for rule in declared["distance_thresholds"]["fields"].values()
    )
    assert set(declared["network_outputs"]) == set(NETWORK_OUTPUT_FIELDS)
    assert all(
        output["models"] == ["first"]
        for output in declared["network_outputs"].values()
    )
    assert declared["combined_component_semantics"] == COMBINED_COMPONENT_SEMANTICS


def test_precomputed_annotation_output_profile_matches_runtime_validation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    declared = capability_document()["server_interface"]["scientific_profiles"][
        "precomputed"
    ]["network_annotations"]
    monkeypatch.setitem(config, "reference_manifest_id", "refs-capability-test")
    monkeypatch.setitem(
        config,
        "network_annotations",
        {
            "use_embedded_sifts": declared["use_embedded_sifts"],
            "tooltip_fields": declared["tooltip_fields"],
            "max_tooltip_segments_per_database": declared[
                "max_tooltip_segments_per_database"
            ],
        },
    )

    assert network_annotations.network_annotation_config() == {
        "use_embedded_sifts": True,
        "tooltip_fields": [],
        "max_tooltip_segments_per_database": 20,
    }
    assert precomputed_store.annotation_profile()["use_embedded_sifts"] is True
    assert declared["cache_validation"] == "current_annotation_profile_required"


def test_emitted_summary_and_artifact_manifest_use_contract_version_constants(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    identity = identity_from_official_id("1abc")
    paths = create_run_output_paths(str(tmp_path / "internal"))
    write_run_manifest(
        paths.manifest_file,
        input_files=["1abc.pdb"],
        output_paths=paths,
        config_snapshot={},
        status="success",
        started_at="2026-08-11T00:00:00",
        finished_at="2026-08-11T00:00:01",
        total_time=1.0,
        extra_counts={"structures": 1, "chains": 1, "skipped_outputs": 0},
        references={"manifest_id": "tiny-reference"},
        identities=[identity.as_dict()],
        structure_inputs=[
            {
                "file": "1abc.pdb",
                "identity": identity.as_dict(),
                "format": "pdb",
                "kind": "pdb",
                "embedded_annotation_counts": {
                    "uniprot": 0,
                    "pfam": 0,
                    "cath": 0,
                    "scop2": 0,
                },
            }
        ],
    )
    public = collect_web_outputs(paths, str(tmp_path / "web"))
    monkeypatch.setitem(config, "reference_manifest_id", "tiny-reference")
    source = tmp_path / "1abc.pdb"
    source.write_text("HEADER TINY 1ABC\n", encoding="utf-8")
    structure = {
        "file_path": str(source),
        "pdb_id": "1ABC",
        "structure_identity": identity.as_dict(),
        "atom_data": [
            {
                "chain_id": "A",
                "chain_identity": ChainIdentity(
                    identity.canonical_id, identity.display_id, "A", 1
                ).as_dict(),
                "structure_key": identity.canonical_id,
                "structure_display_id": identity.display_id,
                "model_index": 1,
                "molecule_name": "Tiny protein",
                "molecule_type": "Protein",
                "sequence": "A",
                "is_hetatm": False,
                "residues": [{"residue_name": "ALA", "residue_number": 1, "atoms": []}],
            }
        ],
    }
    entry_path, changed = precomputed_store.write_entry(
        tmp_path / "store", source, structure, []
    )
    assert changed is True
    raw_entry = json.loads(gzip.decompress(entry_path.read_bytes()))
    profile_root = entry_path.parents[2]
    artifact_manifest = json.loads(
        (profile_root / "manifest.json").read_text(encoding="utf-8")
    )

    assert set(public) == set(WEB_OUTPUT_SUMMARY_FIELDS)
    assert public["output_contract_version"] == OUTPUT_CONTRACT_VERSION
    validation = capability_document()["server_interface"]["contracts"]["web_output"][
        "validation_semantics"
    ]
    assert validation == WEB_OUTPUT_VALIDATION_SEMANTICS
    assert set(public["structure_inputs"][0]) == set(
        validation["structure_inputs"]["fields"]
    )
    assert public["structure_inputs"][0]["format"] in validation[
        "structure_inputs"
    ]["formats"]
    assert public["structure_inputs"][0]["kind"] in validation[
        "structure_inputs"
    ]["kinds"]
    assert set(public["structure_inputs"][0]["embedded_annotation_counts"]) == set(
        validation["structure_inputs"]["embedded_annotation_count_fields"]
    )
    assert set(raw_entry) == set(PRECOMPUTED_ENTRY_FIELDS)
    assert artifact_manifest["artifact_schema_version"] == PRECOMPUTED_ARTIFACT_SCHEMA_VERSION


def test_populate_missing_requires_source_directory_for_a_cache_miss(
    tmp_path: Path,
) -> None:
    with pytest.raises(InputValidationError, match="PDB_ARCHIVE_REQUIRED"):
        precomputed_store._load_or_populate_entries(
            tmp_path / "store",
            ["pdb_00001abc"],
            source_dir=None,
            populate_missing=True,
            selected_profile="a" * 64,
            profile_directory=tmp_path / "profile",
        )


def test_populate_missing_allows_cache_hit_without_source_directory(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    selected_profile = "a" * 64
    store_root = tmp_path / "store"
    cached = precomputed_store.entry_path(
        store_root, "pdb_00001abc", selected_profile
    )
    cached.parent.mkdir(parents=True)
    cached.write_bytes(b"cached-entry")
    payload = {"counts": {"chains": 1, "interactions": 0}}
    monkeypatch.setattr(
        precomputed_store,
        "_read_entry_for_profile",
        lambda *_args, **_kwargs: (payload, 12),
    )

    loaded = precomputed_store._load_or_populate_entries(
        store_root,
        ["pdb_00001abc"],
        source_dir=None,
        populate_missing=True,
        selected_profile=selected_profile,
        profile_directory=cached.parents[2],
    )

    assert loaded == [(payload, True)]
