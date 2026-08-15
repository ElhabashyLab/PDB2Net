"""Stable, configuration-free semantics consumed by server integrations.

This module deliberately contains only JSON-native constants.  Importing it
must not load machine-local configuration, Gemmi, reference data, or the
precomputed store.
"""

from __future__ import annotations

import copy
from typing import Any


NETWORK_OUTPUT_FIELDS = (
    "chain_per_pdb",
    "protein_per_pdb",
    "combined_chain_network",
    "combined_protein_network",
)
CORE_STRUCTURE_EXTENSIONS = (".pdb", ".ent", ".cif", ".mmcif")
WEB_SERVER_STRUCTURE_EXTENSIONS = (".pdb", ".cif", ".mmcif")
WEB_SERVER_INPUT_SUFFIXES = (
    *WEB_SERVER_STRUCTURE_EXTENSIONS,
    *(f"{suffix}.gz" for suffix in WEB_SERVER_STRUCTURE_EXTENSIONS),
)
STRUCTURE_MODEL_POLICIES = ("first", "all")
ALL_MODEL_FORBIDDEN_NETWORKS = (
    "protein_per_pdb",
    "combined_protein_network",
)
RESOURCE_LIMIT_FIELDS = (
    "max_input_files",
    "max_total_input_bytes",
    "max_single_input_bytes",
    "max_processing_batch_bytes",
    "max_total_input_expanded_bytes",
    "max_single_input_expanded_bytes",
    "max_detailed_interaction_rows",
    "max_detailed_interaction_bytes",
    "min_output_free_bytes",
)
ANNOTATION_DATABASES = ("uniprot", "pfam", "cath", "scop2")
SERVER_ENVIRONMENT = {
    "PDB2NET_PDB_FASTA": {
        "config_path": "pdb_fasta_path",
        "type": "path",
        "supply_condition": "real_worker",
    },
    "PDB2NET_SIFTS_TSV": {
        "config_path": "sifts_tsv_path",
        "type": "path",
        "supply_condition": "real_worker",
    },
    "PDB2NET_UNIPROT_FASTA": {
        "config_path": "uniprot_fasta_path",
        "type": "path",
        "supply_condition": "real_worker",
    },
    "PDB2NET_BLAST_DB": {
        "config_path": "blast_db_path",
        "type": "path",
        "supply_condition": "real_worker",
    },
    "PDB2NET_BLASTP": {
        "config_path": "blastp_executable",
        "type": "executable",
        "supply_condition": "real_worker",
    },
    "PDB2NET_BLAST_CACHE_PATH": {
        "config_path": "blast_cache_path",
        "type": "path",
        "supply_condition": "real_worker",
    },
}

DETAILED_INTERACTION_COLUMNS = (
    "PDB_ID",
    "Chain_A",
    "Residue_A",
    "Atom_A",
    "Chain_B",
    "Residue_B",
    "Atom_B",
    "Distance",
    "UniProt_A",
    "UniProt_B",
    "Interaction_Type",
)
ALLOWED_INTERACTION_TYPES = frozenset(
    {
        "Protein-Protein",
        "Protein-DNA",
        "Protein-RNA",
        "Protein-DNA/RNA",
        "Protein-Nucleic Acid",
        "DNA-DNA",
        "DNA-RNA",
        "DNA-DNA/RNA",
        "RNA-RNA",
        "RNA-DNA/RNA",
        "DNA/RNA-DNA/RNA",
        "Nucleic Acid-Nucleic Acid",
    }
)
DETAILED_INTERACTION_FILENAME_SUFFIX = "_detailed_interactions"

CX2_HEADER = {"CXVersion": "2.0", "hasFragments": False}
CX2_REQUIRED_ASPECT_ORDER = (
    "metaData",
    "attributeDeclarations",
    "networkAttributes",
    "nodes",
    "edges",
    "visualProperties",
    "status",
)
CX2_FORBIDDEN_ASPECTS = ("cartesianLayout", "edgeAttributes", "nodeAttributes")
CX2_DECLARATION_SCOPES = ("networkAttributes", "nodes", "edges")
CX2_SUPPORTED_ATTRIBUTE_TYPES = ("boolean", "double", "integer", "long", "string")
CX2_SUCCESS_STATUS = ({"error": "", "success": True},)
NETWORK_TITLE_BASES = {
    "chain_per_pdb": "Chain_Interaction_Network",
    "protein_per_pdb": "Protein_Network",
    "combined_chain": "Combined_Network",
    "combined_protein": "Combined_Protein_Network",
}
ALL_MODEL_TITLE_SUFFIX_PATTERN = r"_model([1-9][0-9]*)"
PORTABLE_ARTIFACT_STEM_SEMANTICS_ID = "pdb2net-portable-artifact-stem-v1"
MAX_ARTIFACT_STEM_BYTES = 180

PARSER_SEMANTICS = (
    "validated-single-document-gemmi-heavy-atoms-no-hydrogen-or-deuterium-v3"
)
MMCIF_IDENTITY_POLICY = {
    "data_block_count": {"minimum": 1, "maximum": 1},
    "content_claims": [
        "data_block",
        "_entry.id",
        "_database_2[PDB].database_code",
    ],
    "content_claim_conflict": "reject_with_CONFLICTING_STRUCTURE_IDENTITY",
    "filename_mismatch": "content_wins_with_STRUCTURE_FILENAME_ID_MISMATCH_warning",
}
INTERACTION_PIPELINE_VERSION = "pdb2net-distance-classification-v1"
ANNOTATION_PIPELINE_VERSION = "pdb2net-sifts-fasta-search-annotations-v3"
COMBINED_COMPONENT_SEMANTICS = "pdb2net-cross-pdb-uniprot-linked-components-v1"
UPLOAD_SOURCE_SCOPE = "input_document"
COMMAND_SEMANTICS_IDS = {
    "capabilities": "pdb2net-capabilities-json-v2",
    "run": "pdb2net-run-structure-pipeline-v1",
    "precompute": "pdb2net-precompute-artifact-store-v2",
    "assemble": "pdb2net-assemble-precomputed-artifacts-v2",
}

DISTANCE_THRESHOLD_RULES = {
    "ca_radius": {"type": "number", "default": 12.0, "minimum": 2.0, "maximum": 30.0},
    "all_atoms_radius": {"type": "number", "default": 5.0, "minimum": 1.0, "maximum": 15.0},
}
INTERACTION_FILTER_RULES = {
    "protein_protein_min_ca_neighbors": {
        "type": "integer",
        "default": 10,
        "minimum": 1,
        "maximum": 100_000,
    },
    "protein_protein_min_all_atom_contacts": {
        "type": "integer",
        "default": 1,
        "minimum": 1,
        "maximum": 100_000,
    },
    "protein_nucleic_acid_min_all_atom_contacts": {
        "type": "integer",
        "default": 1,
        "minimum": 1,
        "maximum": 100_000,
    },
    "nucleic_acid_min_all_atom_contacts": {
        "type": "integer",
        "default": 1,
        "minimum": 1,
        "maximum": 100_000,
    },
}
NETWORK_ANNOTATION_RULES = {
    "use_embedded_sifts": {"type": "boolean", "default": True},
    "tooltip_fields": {
        "type": "array",
        "default": ["uniprot"],
        "items": {"type": "string", "choices": list(ANNOTATION_DATABASES)},
        "unique_items": True,
    },
    "max_tooltip_segments_per_database": {
        "type": "integer",
        "default": 20,
        "minimum": 1,
        "maximum": 1_000,
    },
}

PRECOMPUTED_SCIENTIFIC_PIPELINE_VERSION = "pdb2net-asu-first-standard-interactions-v3"
PRECOMPUTED_ANNOTATION_PIPELINE_VERSION = ANNOTATION_PIPELINE_VERSION
PRECOMPUTED_SOURCE_SCOPE = "asymmetric_unit"
PRECOMPUTED_PARSER_SEMANTICS = PARSER_SEMANTICS

WEB_OUTPUT_SUMMARY_FIELDS = (
    "output_contract_version",
    "pdb2net_version",
    "status",
    "started_at",
    "finished_at",
    "input_files",
    "identities",
    "structure_inputs",
    "networks",
    "interactions",
    "artifacts",
    "runtime_analysis",
    "counts",
    "annotations",
    "references",
    "resources",
    "skipped_outputs",
    "config",
    "errors",
    "warnings",
)
WEB_OUTPUT_COUNT_FIELDS = (
    "networks",
    "interactions",
    "structures",
    "chains",
    "skipped_outputs",
)
PRECOMPUTED_ENTRY_FIELDS = (
    "artifact_schema_version",
    "scientific_pipeline_version",
    "pdb2net_version",
    "profile_id",
    "pdb_id",
    "structure_identity",
    "created_at",
    "source",
    "geometry",
    "annotations",
    "counts",
)
STRUCTURE_IDENTITY_FIELDS = (
    "source",
    "canonical_id",
    "legacy_id",
    "key",
    "display_id",
)
WEB_OUTPUT_VALIDATION_SEMANTICS = {
    "semantics_id": "pdb2net-web-output-validation-v1",
    "structure_inputs": {
        "fields": [
            "file",
            "identity",
            "format",
            "kind",
            "embedded_annotation_counts",
        ],
        "formats": ["pdb", "mmcif"],
        "kinds": ["pdb", "mmcif", "nextgen_enriched_mmcif"],
        "embedded_annotation_count_fields": list(ANNOTATION_DATABASES),
    },
    "artifact_records": {
        "network_fields": ["path", "size_bytes", "sha256", "nodes", "edges"],
        "interaction_fields": [
            "path",
            "size_bytes",
            "sha256",
            "rows",
            "columns",
        ],
        "artifact_lists_match_summary_path_order": True,
    },
    "selected_outputs": {
        "first_model_per_identity_networks_cover_each_identity_exactly_once": True,
        "all_model_chain_networks_cover_each_identity_at_least_once": True,
        "all_model_chain_network_identity_model_pairs_are_unique": True,
        "all_model_chain_network_model_indices_need_not_be_contiguous": True,
        "interaction_csvs_cover_each_identity_exactly_once": True,
        "chain_network_node_sum_equals_summary_chains": True,
    },
}


def _command_option(
    value_type: str,
    *,
    required: bool = False,
    cardinality: str = "one",
    value_format: str | None = None,
    choices: list[str] | None = None,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "type": value_type,
        "required": required,
        "cardinality": cardinality,
        "choices": choices,
    }
    if value_format is not None:
        result["format"] = value_format
    return result


COMMANDS = {
    "capabilities": {
        "semantics_id": COMMAND_SEMANTICS_IDS["capabilities"],
        "success_exit_code": 0,
        "options": {
            "--json": _command_option("boolean"),
        },
    },
    "run": {
        "semantics_id": COMMAND_SEMANTICS_IDS["run"],
        "success_exit_code": 0,
        "options": {
            "--input-dir": _command_option("string", required=True, value_format="path"),
            "--output-dir": _command_option("string", required=True, value_format="path"),
            "--web-output-dir": _command_option("string", value_format="path"),
            "--config": _command_option("string", value_format="path"),
            "--headless": _command_option("boolean"),
        },
    },
    "precompute": {
        "semantics_id": COMMAND_SEMANTICS_IDS["precompute"],
        "success_exit_code": 0,
        "options": {
            "--input-dir": _command_option("string", required=True, value_format="path"),
            "--store": _command_option("string", required=True, value_format="path"),
            "--config": _command_option("string", value_format="path"),
            "--recursive": _command_option("boolean"),
            "--headless": _command_option("boolean"),
        },
    },
    "assemble": {
        "semantics_id": COMMAND_SEMANTICS_IDS["assemble"],
        "success_exit_code": 0,
        "options": {
            "--store": _command_option("string", required=True, value_format="path"),
            "--pdb-id": _command_option(
                "string", required=True, cardinality="one_or_more", value_format="pdb_id"
            ),
            "--output-dir": _command_option("string", required=True, value_format="path"),
            "--web-output-dir": _command_option("string", value_format="path"),
            "--config": _command_option("string", value_format="path"),
            "--source-dir": _command_option("string", value_format="path"),
            "--populate-missing": _command_option("boolean"),
            "--headless": _command_option("boolean"),
        },
        "cross_field_rules": [
            {
                "if_present": "--populate-missing",
                "when": "requested_entry_missing",
                "requires": ["--source-dir"],
            },
        ],
    },
}


def _web_config_fields() -> dict[str, dict[str, Any]]:
    fields: dict[str, dict[str, Any]] = {}
    for name, rule in DISTANCE_THRESHOLD_RULES.items():
        fields[f"distance_thresholds.{name}"] = copy.deepcopy(rule)
    for name, rule in INTERACTION_FILTER_RULES.items():
        fields[f"interaction_filters.{name}"] = copy.deepcopy(rule)
    for name in NETWORK_OUTPUT_FIELDS:
        fields[f"networks.{name}"] = {"type": "boolean"}
    fields["structure_model_policy"] = {
        "type": "string",
        "choices": list(STRUCTURE_MODEL_POLICIES),
    }
    for name in ("parsing", "blast_threads"):
        fields[f"workers.{name}"] = {"type": "integer", "minimum": 1}
    for name in RESOURCE_LIMIT_FIELDS:
        fields[f"resource_limits.{name}"] = {"type": "integer", "minimum": 1}
    for name, rule in NETWORK_ANNOTATION_RULES.items():
        fields[f"network_annotations.{name}"] = copy.deepcopy(rule)
    for name in ("max_nodes", "max_edges"):
        fields[f"combined_graph_limits.{name}"] = {"type": "integer", "minimum": 1}
    fields["reference_manifest_id"] = {"type": "string", "minimum_length": 1}
    fields.update(
        {
            "diamond.enabled": {"type": "boolean"},
            "diamond.executable": {"type": "string", "minimum_length": 1},
            "diamond.uniref90_db_path": {"type": "string", "minimum_length": 1},
            "diamond.temp_dir": {"type": "string", "minimum_length": 1},
            "diamond.threads": {"type": "integer", "minimum": 1},
            "diamond.iterate": {"type": "boolean"},
            "diamond.sensitivity": {
                "type": "string",
                "choices": [
                    "default",
                    "fast",
                    "mid-sensitive",
                    "sensitive",
                    "more-sensitive",
                    "very-sensitive",
                    "ultra-sensitive",
                ],
            },
            "diamond.block_size": {"type": "number", "exclusive_minimum": 0},
            "diamond.index_chunks": {"type": "integer", "minimum": 1},
            "diamond.max_target_seqs": {"type": "integer", "minimum": 1},
            "diamond.batch_max_sequences": {"type": "integer", "minimum": 1},
            "diamond.batch_max_fasta_bytes": {"type": "integer", "minimum": 1},
            "diamond.assign_uniprot_id": {
                "type": "string",
                "choices": ["never", "high_confidence", "always"],
            },
            "export_detailed_interactions": {"type": "boolean"},
            "open_in_cytoscape": {"type": "boolean", "const": False},
        }
    )
    return fields


def _constraints_without_defaults(
    rules: dict[str, dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    constraints = copy.deepcopy(rules)
    for rule in constraints.values():
        rule.pop("default", None)
    return constraints


CONTRACTS = {
    "input": {
        "supported_suffixes": list(WEB_SERVER_INPUT_SUFFIXES),
        "suffix_matching": "case_insensitive",
        "gzip": {
            "optional": True,
            "maximum_layers": 1,
        },
        "mmcif_identity_policy": copy.deepcopy(MMCIF_IDENTITY_POLICY),
    },
    "environment": copy.deepcopy(SERVER_ENVIRONMENT),
    "web_config": {
        "fields": _web_config_fields(),
        "cross_field_rules": [
            {
                "if": {"structure_model_policy": "all"},
                "requires": {
                    "networks.protein_per_pdb": False,
                    "networks.combined_protein_network": False,
                },
            },
            {
                "when_web_output_requested": True,
                "requires_nonempty": ["reference_manifest_id"],
            },
        ],
    },
    "web_output": {
        "layout": {
            "summary": "summary.json",
            "network_directory": "networks",
            "network_suffix": ".cx2",
            "interaction_directory": "interactions",
            "interaction_suffix": ".csv",
            "runtime_analysis": "runtime_analysis.txt",
        },
        "summary": {
            "fields": list(WEB_OUTPUT_SUMMARY_FIELDS),
            "success_status": "success",
            "failure_status": "failed",
            "success_errors": [],
            "counts_fields": list(WEB_OUTPUT_COUNT_FIELDS),
            "artifact_groups": ["networks", "interactions"],
            "identity_fields": list(STRUCTURE_IDENTITY_FIELDS),
            "successful_inputs_have_exactly_one_identity": True,
            "reference_manifest_id_required_on_success": True,
            "precomputed_reference_fields": [
                "artifact_schema_version",
                "profile_id",
                "annotation_profile_id",
                "source_scope",
            ],
        },
        "web_acceptance": {
            "scope": "successful_real_jobs_validated_with_selected_outputs",
            "selected_outputs_must_not_be_skipped": True,
        },
        "validation_semantics": copy.deepcopy(WEB_OUTPUT_VALIDATION_SEMANTICS),
        "artifacts": {
            "paths_are_relative": True,
            "directory_contents_are_exact": True,
            "sha256_and_size_required": True,
            "network_counts": ["nodes", "edges"],
            "interaction_counts": ["rows", "columns"],
            "portable_stem": {
                "semantics_id": PORTABLE_ARTIFACT_STEM_SEMANTICS_ID,
                "maximum_utf8_bytes": MAX_ARTIFACT_STEM_BYTES,
            },
        },
        "interaction_csv": {
            "columns": list(DETAILED_INTERACTION_COLUMNS),
            "allowed_interaction_types": sorted(ALLOWED_INTERACTION_TYPES),
            "filename": {
                "scope": "per_identity",
                "identity_field": "display_id",
                "identity_suffix": DETAILED_INTERACTION_FILENAME_SUFFIX,
                "extension": ".csv",
                "derived_via_portable_stem": True,
            },
        },
        "cx2": {
            "header": copy.deepcopy(CX2_HEADER),
            "required_aspect_order": list(CX2_REQUIRED_ASPECT_ORDER),
            "forbidden_aspects": list(CX2_FORBIDDEN_ASPECTS),
            "attribute_declarations": {
                "scopes": list(CX2_DECLARATION_SCOPES),
                "network_name": {"field": "name", "type": "string"},
                "supported_types": list(CX2_SUPPORTED_ATTRIBUTE_TYPES),
            },
            "success_status": [dict(item) for item in CX2_SUCCESS_STATUS],
            "network_titles": {
                "chain_per_pdb": {
                    "prefix": f"{NETWORK_TITLE_BASES['chain_per_pdb']}_",
                    "scope": "per_identity",
                    "identity_field": "display_id",
                    "first_model_suffix": "",
                    "all_model_suffix_pattern": ALL_MODEL_TITLE_SUFFIX_PATTERN,
                },
                "protein_per_pdb": {
                    "prefix": f"{NETWORK_TITLE_BASES['protein_per_pdb']}_",
                    "scope": "per_identity",
                    "identity_field": "display_id",
                    "model_policy": "first",
                },
                "combined_chain": {
                    "prefix": f"{NETWORK_TITLE_BASES['combined_chain']}_",
                    "scope": "inter_pdb_identity_component",
                    "suffix": "nonempty",
                },
                "combined_protein": {
                    "prefix": f"{NETWORK_TITLE_BASES['combined_protein']}_",
                    "scope": "inter_pdb_identity_component",
                    "suffix": "nonempty",
                },
            },
            "filename": {
                "source": "networkAttributes.name",
                "extension": ".cx2",
                "derived_via_portable_stem": True,
            },
        },
    },
    "precomputed_artifact": {
        "entry_fields": list(PRECOMPUTED_ENTRY_FIELDS),
        "identity_fields": list(STRUCTURE_IDENTITY_FIELDS),
        "source_fields": ["name", "sha256", "size_bytes", "scope"],
        "source_scope": PRECOMPUTED_SOURCE_SCOPE,
        "profile_id_format": "sha256",
        "annotation_profile_id_format": "sha256",
        "canonical_pdb_id_format": "pdb_XXXXXXXX_lowercase",
        "duplicate_canonical_ids_rejected": True,
    },
}


SCIENTIFIC_PROFILES = {
    "upload": {
        "source_scope": UPLOAD_SOURCE_SCOPE,
        "parser_semantics": PARSER_SEMANTICS,
        "interaction_pipeline_version": INTERACTION_PIPELINE_VERSION,
        "annotation_pipeline_version": ANNOTATION_PIPELINE_VERSION,
        "combined_component_semantics": COMBINED_COMPONENT_SEMANTICS,
        "distance_thresholds": copy.deepcopy(DISTANCE_THRESHOLD_RULES),
        "interaction_filters": copy.deepcopy(INTERACTION_FILTER_RULES),
        "structure_model_policy": {
            "default": "first",
            "choices": list(STRUCTURE_MODEL_POLICIES),
            "all_forbids_networks": list(ALL_MODEL_FORBIDDEN_NETWORKS),
        },
        "network_outputs": {
            "chain_per_pdb": {
                "scope": "per_identity",
                "models": ["first", "all"],
                "selected_emission": "required_per_identity",
            },
            "protein_per_pdb": {
                "scope": "per_identity",
                "models": ["first"],
                "selected_emission": "annotation_conditional",
                "requires": "at_least_one_uniprot_annotated_protein",
            },
            "combined_chain_network": {
                "scope": "inter_pdb_identity_component",
                "models": ["first", "all"],
                "selected_emission": "component_conditional",
                "may_be_skipped_by": ["combined_graph_limits"],
            },
            "combined_protein_network": {
                "scope": "inter_pdb_identity_component",
                "models": ["first"],
                "selected_emission": "component_conditional",
                "may_be_skipped_by": ["combined_graph_limits"],
            },
        },
        "detailed_interactions": {"supported": True, "default": False},
        "network_annotations": {
            "databases": list(ANNOTATION_DATABASES),
            **copy.deepcopy(NETWORK_ANNOTATION_RULES),
        },
    },
    "precomputed": {
        "scientific_pipeline_version": PRECOMPUTED_SCIENTIFIC_PIPELINE_VERSION,
        "annotation_pipeline_version": PRECOMPUTED_ANNOTATION_PIPELINE_VERSION,
        "interaction_pipeline_version": INTERACTION_PIPELINE_VERSION,
        "source_scope": PRECOMPUTED_SOURCE_SCOPE,
        "structure_model_policy": "first",
        "parser_semantics": PRECOMPUTED_PARSER_SEMANTICS,
        "parser_runtime_version_is_content_addressed": True,
        "combined_component_semantics": COMBINED_COMPONENT_SEMANTICS,
        "network_outputs": {
            "chain_per_pdb": {
                "scope": "per_identity",
                "models": ["first"],
                "selected_emission": "required_per_identity",
            },
            "protein_per_pdb": {
                "scope": "per_identity",
                "models": ["first"],
                "selected_emission": "annotation_conditional",
                "requires": "at_least_one_uniprot_annotated_protein",
            },
            "combined_chain_network": {
                "scope": "inter_pdb_identity_component",
                "models": ["first"],
                "selected_emission": "component_conditional",
                "may_be_skipped_by": ["combined_graph_limits"],
            },
            "combined_protein_network": {
                "scope": "inter_pdb_identity_component",
                "models": ["first"],
                "selected_emission": "component_conditional",
                "may_be_skipped_by": ["combined_graph_limits"],
            },
        },
        "network_annotations": {
            "use_embedded_sifts": True,
            "tooltip_fields": [],
            "max_tooltip_segments_per_database": NETWORK_ANNOTATION_RULES[
                "max_tooltip_segments_per_database"
            ]["default"],
            "cache_validation": "current_annotation_profile_required",
        },
        "distance_thresholds": {
            "content_addressed": True,
            "fields": _constraints_without_defaults(DISTANCE_THRESHOLD_RULES),
        },
        "interaction_filters": {
            "content_addressed": True,
            "fields": _constraints_without_defaults(INTERACTION_FILTER_RULES),
        },
        "network_selection_affects_geometry_profile": False,
        "detailed_interactions": False,
        "reference_manifest_id_required": True,
    },
}


SERVER_INTERFACE = {
    "commands": COMMANDS,
    "contracts": CONTRACTS,
    "scientific_profiles": SCIENTIFIC_PROFILES,
}


def server_interface_document() -> dict[str, Any]:
    """Return a mutable JSON-native copy of the stable server interface."""
    return copy.deepcopy(SERVER_INTERFACE)


__all__ = [
    "ALLOWED_INTERACTION_TYPES",
    "ALL_MODEL_FORBIDDEN_NETWORKS",
    "ALL_MODEL_TITLE_SUFFIX_PATTERN",
    "ANNOTATION_DATABASES",
    "ANNOTATION_PIPELINE_VERSION",
    "COMBINED_COMPONENT_SEMANTICS",
    "COMMAND_SEMANTICS_IDS",
    "CORE_STRUCTURE_EXTENSIONS",
    "CX2_DECLARATION_SCOPES",
    "CX2_FORBIDDEN_ASPECTS",
    "CX2_HEADER",
    "CX2_REQUIRED_ASPECT_ORDER",
    "CX2_SUCCESS_STATUS",
    "CX2_SUPPORTED_ATTRIBUTE_TYPES",
    "DETAILED_INTERACTION_COLUMNS",
    "DETAILED_INTERACTION_FILENAME_SUFFIX",
    "DISTANCE_THRESHOLD_RULES",
    "INTERACTION_PIPELINE_VERSION",
    "INTERACTION_FILTER_RULES",
    "MAX_ARTIFACT_STEM_BYTES",
    "MMCIF_IDENTITY_POLICY",
    "NETWORK_ANNOTATION_RULES",
    "NETWORK_OUTPUT_FIELDS",
    "NETWORK_TITLE_BASES",
    "PARSER_SEMANTICS",
    "PORTABLE_ARTIFACT_STEM_SEMANTICS_ID",
    "PRECOMPUTED_ANNOTATION_PIPELINE_VERSION",
    "PRECOMPUTED_PARSER_SEMANTICS",
    "PRECOMPUTED_SCIENTIFIC_PIPELINE_VERSION",
    "PRECOMPUTED_SOURCE_SCOPE",
    "RESOURCE_LIMIT_FIELDS",
    "SERVER_ENVIRONMENT",
    "SERVER_INTERFACE",
    "STRUCTURE_MODEL_POLICIES",
    "UPLOAD_SOURCE_SCOPE",
    "WEB_OUTPUT_COUNT_FIELDS",
    "WEB_OUTPUT_SUMMARY_FIELDS",
    "WEB_OUTPUT_VALIDATION_SEMANTICS",
    "WEB_SERVER_INPUT_SUFFIXES",
    "WEB_SERVER_STRUCTURE_EXTENSIONS",
    "server_interface_document",
]
