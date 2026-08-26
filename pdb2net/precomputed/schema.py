"""Schema-3 profile and entry validation for the precomputed store."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import math
import re
from typing import Any, Iterable, Mapping

from .. import __version__
from ..config_loader import config
from ..contracts import PRECOMPUTED_ARTIFACT_SCHEMA_VERSION
from ..data_processor import PARSER_SEMANTICS
from ..distances import (
    ALLOWED_INTERACTION_TYPES,
    DISTANCE_THRESHOLD_RULES,
    INTERACTION_FILTER_RULES,
    INTERACTION_PIPELINE_VERSION,
    determine_interaction_type,
)
from ..input_contract import InputValidationError
from ..network_annotations import (
    ANNOTATION_PIPELINE_VERSION,
    network_annotation_config,
)
from ..structure_identity import ChainIdentity, StructureIdentity, canonical_pdb_id


ARTIFACT_SCHEMA_VERSION = PRECOMPUTED_ARTIFACT_SCHEMA_VERSION
SCIENTIFIC_PIPELINE_VERSION = "pdb2net-asu-first-standard-interactions-v3"
SOURCE_SCOPE = "asymmetric_unit"

MAX_COMPRESSED_ENTRY_BYTES = 64 * 1024 * 1024
MAX_DECOMPRESSED_ENTRY_BYTES = 128 * 1024 * 1024
MAX_CHAINS_PER_ENTRY = 50_000
MAX_EDGES_PER_ENTRY = 500_000
MAX_REQUEST_EXPANDED_BYTES = 256 * 1024 * 1024
MAX_REQUEST_CHAINS = 200_000
MAX_REQUEST_EDGES = 2_000_000
MAX_EMBEDDED_SEGMENTS_PER_CHAIN = 100_000

EMBEDDED_SEGMENT_FIELDS = {
    "database",
    "source_database",
    "accession",
    "domain_name",
    "chain_id",
    "label_asym_id",
    "entity_id",
    "segment_id",
    "instance_id",
    "pdb_start",
    "pdb_end",
    "database_start",
    "database_end",
    "best_mapping",
    "identity",
}
ANNOTATION_CHAIN_FIELDS = {
    "molecule_name",
    "molecule_type",
    "uniprot_id",
    "annotation_source",
    "matched_database",
    "matched_id",
    "representative_accession",
    "annotation_confidence",
    "external_sifts_status",
    "external_sifts_accessions",
    "annotation_warnings",
}
GEOMETRY_CHAIN_FIELDS = {
    "chain_id",
    "unique_chain_id",
    "chain_identity",
    "structure_key",
    "structure_display_id",
    "model_index",
    "sequence",
    "is_hetatm",
    "aa_len",
    "nt_len",
    "embedded_annotations",
    "embedded_annotation_source",
    "embedded_uniprot_status",
    "embedded_uniprot_accessions",
}
ALLOWED_MOLECULE_TYPES = {
    "Protein",
    "DNA",
    "RNA",
    "DNA/RNA",
    "Nucleic Acid",
    "Unknown",
    "UNKNOWN",
}


def canonical_json(value: object) -> bytes:
    """Serialize one canonical, finite JSON value."""
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")


def utc_timestamp() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def producer() -> dict[str, str]:
    return {"name": "pdb2net", "version": __version__}


def normalize_pdb_id(value: object) -> str:
    """Return a canonical lowercase wwPDB Extended identifier."""
    raw = str(value or "").strip()
    try:
        return canonical_pdb_id(raw)
    except ValueError as exc:
        raise InputValidationError(
            "INVALID_PDB_ID",
            "PDB IDs must be a legacy four-character ID or an extended "
            f"pdb_XXXXXXXX ID: {raw!r}",
        ) from exc


def normalize_requested_ids(values: Iterable[object]) -> list[str]:
    """Normalize requested IDs and reject canonical duplicates."""
    normalized: list[str] = []
    seen: set[str] = set()
    for value in values:
        pdb_id = normalize_pdb_id(value)
        if pdb_id in seen:
            raise InputValidationError(
                "DUPLICATE_STRUCTURE_IDENTITY",
                f"More than one requested input resolves to structure identity {pdb_id!r}.",
            )
        seen.add(pdb_id)
        normalized.append(pdb_id)
    if not normalized:
        raise InputValidationError("NO_PDB_IDS", "At least one PDB ID is required.")
    return normalized


def _reference_manifest_id() -> str:
    value = str(config.get("reference_manifest_id") or "").strip()
    if not value:
        raise InputValidationError(
            "REFERENCE_MANIFEST_ID_REQUIRED",
            "Precomputed runs require a non-empty reference_manifest_id.",
        )
    return value


def _fixed_numeric_values(
    section: str,
    rules: Mapping[str, Mapping[str, Any]],
    *,
    integer: bool,
) -> dict[str, int | float]:
    raw = config.get(section, {})
    if not isinstance(raw, Mapping):
        raw = {}
    result: dict[str, int | float] = {}
    for name, rule in rules.items():
        expected = rule["default"]
        value = raw.get(name, expected)
        valid_type = (
            isinstance(value, int) and not isinstance(value, bool)
            if integer
            else isinstance(value, (int, float))
            and not isinstance(value, bool)
            and math.isfinite(float(value))
        )
        if not valid_type:
            raise InputValidationError(
                "PRECOMPUTED_PROFILE_UNSUPPORTED",
                f"Precomputed schema 3 requires the standard {section}.{name} value.",
            )
        normalized = value if integer else float(value)
        if normalized != expected:
            raise InputValidationError(
                "PRECOMPUTED_PROFILE_UNSUPPORTED",
                f"Precomputed schema 3 requires {section}.{name}={expected}.",
            )
        result[name] = normalized
    return result


def scientific_profile() -> dict[str, Any]:
    """Return the one complete fixed scientific/annotation store profile."""
    import gemmi

    policy = str(config.get("structure_model_policy") or "first").strip().lower()
    if policy != "first":
        raise InputValidationError(
            "PRECOMPUTED_MODEL_POLICY_UNSUPPORTED",
            "Precomputed schema 3 requires structure_model_policy='first'.",
        )
    annotation_config = network_annotation_config()
    if annotation_config["use_embedded_sifts"] is not True:
        raise InputValidationError(
            "PRECOMPUTED_PROFILE_UNSUPPORTED",
            "Precomputed schema 3 requires embedded SIFTS annotations.",
        )
    from ..uniprot_matcher import _search_policy_payload

    profile = {
        "artifact_schema_version": ARTIFACT_SCHEMA_VERSION,
        "scientific_pipeline_version": SCIENTIFIC_PIPELINE_VERSION,
        "interaction_pipeline_version": INTERACTION_PIPELINE_VERSION,
        "annotation_pipeline_version": ANNOTATION_PIPELINE_VERSION,
        "source_scope": SOURCE_SCOPE,
        "structure_model_policy": "first",
        "parser": {
            "semantics": PARSER_SEMANTICS,
            "gemmi_version": str(getattr(gemmi, "__version__", "unknown")),
        },
        "distance_thresholds": _fixed_numeric_values(
            "distance_thresholds", DISTANCE_THRESHOLD_RULES, integer=False
        ),
        "interaction_filters": _fixed_numeric_values(
            "interaction_filters", INTERACTION_FILTER_RULES, integer=True
        ),
        "annotations": {
            "reference_manifest_id": _reference_manifest_id(),
            "use_embedded_sifts": True,
            "search_policy": _search_policy_payload(),
        },
    }
    return json.loads(canonical_json(profile))


def profile_id(profile: Mapping[str, Any] | None = None) -> str:
    selected = dict(profile or scientific_profile())
    return hashlib.sha256(canonical_json(selected)).hexdigest()


def manifest_document(profile: Mapping[str, Any], entry_count: int) -> dict[str, Any]:
    selected_profile = dict(profile)
    return {
        "artifact_schema_version": ARTIFACT_SCHEMA_VERSION,
        "created_at": utc_timestamp(),
        "producer": producer(),
        "profile_id": profile_id(selected_profile),
        "profile": selected_profile,
        "entry_count": entry_count,
    }


def _valid_timestamp(value: object) -> bool:
    if not isinstance(value, str) or not value or len(value) > 100:
        return False
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return False
    return parsed.tzinfo is not None and parsed.utcoffset() == timezone.utc.utcoffset(parsed)


def _valid_producer(value: object) -> bool:
    return (
        isinstance(value, dict)
        and set(value) == {"name", "version"}
        and value.get("name") == "pdb2net"
        and isinstance(value.get("version"), str)
        and bool(value["version"])
        and len(value["version"]) <= 64
    )


def validate_manifest(
    payload: object,
    *,
    expected_profile: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Validate one published schema-3 manifest."""
    if not isinstance(payload, dict) or set(payload) != {
        "artifact_schema_version",
        "created_at",
        "producer",
        "profile_id",
        "profile",
        "entry_count",
    }:
        raise InputValidationError(
            "CORRUPT_PRECOMPUTED_MANIFEST",
            "Precomputed manifest fields are invalid.",
        )
    if payload.get("artifact_schema_version") != ARTIFACT_SCHEMA_VERSION:
        raise InputValidationError(
            "PRECOMPUTED_SCHEMA_MISMATCH",
            "Precomputed manifest schema is unsupported.",
        )
    if not _valid_timestamp(payload.get("created_at")) or not _valid_producer(
        payload.get("producer")
    ):
        raise InputValidationError(
            "CORRUPT_PRECOMPUTED_MANIFEST",
            "Precomputed manifest metadata is invalid.",
        )
    profile = payload.get("profile")
    if not isinstance(profile, dict):
        raise InputValidationError(
            "CORRUPT_PRECOMPUTED_MANIFEST",
            "Precomputed manifest profile is invalid.",
        )
    calculated_id = profile_id(profile)
    if payload.get("profile_id") != calculated_id:
        raise InputValidationError(
            "PRECOMPUTED_PROFILE_MISMATCH",
            "Precomputed manifest profile hash is invalid.",
        )
    if expected_profile is not None and profile != dict(expected_profile):
        raise InputValidationError(
            "PRECOMPUTED_PROFILE_MISMATCH",
            "Precomputed manifest does not match the active fixed profile.",
        )
    entry_count = payload.get("entry_count")
    if isinstance(entry_count, bool) or not isinstance(entry_count, int) or entry_count < 0:
        raise InputValidationError(
            "CORRUPT_PRECOMPUTED_MANIFEST",
            "Precomputed manifest entry count is invalid.",
        )
    return dict(payload)


def _valid_nonnegative_int(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value >= 0


def _bounded_optional_text(value: object, *, maximum: int) -> bool:
    return value is None or (isinstance(value, str) and len(value) <= maximum)


def _validate_embedded_annotations(chain: Mapping[str, Any], chain_id: str) -> None:
    embedded = chain.get("embedded_annotations", {})
    if not embedded:
        return
    if not isinstance(embedded, dict) or any(
        key not in {"uniprot", "pfam", "cath", "scop2"} for key in embedded
    ):
        raise InputValidationError(
            "CORRUPT_PRECOMPUTED_ENTRY", "Embedded annotations are invalid."
        )
    for database, segments in embedded.items():
        if not isinstance(segments, list) or len(segments) > MAX_EMBEDDED_SEGMENTS_PER_CHAIN:
            raise InputValidationError(
                "CORRUPT_PRECOMPUTED_ENTRY", "Embedded annotation segments are invalid."
            )
        for segment in segments:
            if not isinstance(segment, dict) or not set(segment).issubset(
                EMBEDDED_SEGMENT_FIELDS
            ):
                raise InputValidationError(
                    "CORRUPT_PRECOMPUTED_ENTRY", "Embedded annotation segment is invalid."
                )
            for value in segment.values():
                if (
                    not isinstance(value, (str, int, float, bool))
                    or isinstance(value, str)
                    and len(value) > 2_000
                    or isinstance(value, float)
                    and not math.isfinite(value)
                ):
                    raise InputValidationError(
                        "CORRUPT_PRECOMPUTED_ENTRY", "Embedded annotation value is invalid."
                    )
            if (
                segment.get("database") != database
                or segment.get("chain_id") != chain_id
                or not isinstance(segment.get("accession"), str)
            ):
                raise InputValidationError(
                    "CORRUPT_PRECOMPUTED_ENTRY", "Embedded annotation chain is invalid."
                )
            pdb_start, pdb_end = segment.get("pdb_start"), segment.get("pdb_end")
            if (
                isinstance(pdb_start, bool)
                or not isinstance(pdb_start, int)
                or isinstance(pdb_end, bool)
                or not isinstance(pdb_end, int)
                or pdb_start < 1
                or pdb_end < pdb_start
            ):
                raise InputValidationError(
                    "CORRUPT_PRECOMPUTED_ENTRY", "Embedded annotation range is invalid."
                )
            if database == "uniprot":
                from ..reference_data import is_uniprot_accession

                database_start = segment.get("database_start")
                database_end = segment.get("database_end")
                identity_value = segment.get("identity")
                if (
                    not is_uniprot_accession(segment.get("accession"))
                    or not isinstance(segment.get("best_mapping"), bool)
                    or isinstance(identity_value, bool)
                    or not isinstance(identity_value, (int, float))
                    or not math.isfinite(float(identity_value))
                    or not 0 <= float(identity_value) <= 1
                    or isinstance(database_start, bool)
                    or not isinstance(database_start, int)
                    or isinstance(database_end, bool)
                    or not isinstance(database_end, int)
                    or database_start < 1
                    or database_end < database_start
                ):
                    raise InputValidationError(
                        "CORRUPT_PRECOMPUTED_ENTRY",
                        "Embedded UniProt segment is invalid.",
                    )


def _materialized_chains(payload: Mapping[str, Any]) -> list[dict[str, Any]]:
    geometry = payload["geometry"]
    annotations = payload["annotations"]
    annotation_map = {
        item["unique_chain_id"]: item
        for item in annotations["chains"]
    }
    chains: list[dict[str, Any]] = []
    for geometry_chain in geometry["structure"]["atom_data"]:
        chain = dict(geometry_chain)
        annotation = dict(annotation_map[chain["unique_chain_id"]])
        annotation.pop("unique_chain_id", None)
        chain.update(annotation)
        chains.append(chain)
    return chains


def validate_entry(
    payload: object,
    *,
    expected_pdb_id: str,
    expected_profile_id: str,
    expected_profile: Mapping[str, Any],
) -> dict[str, Any]:
    """Validate one untrusted schema-3 entry without changing its shape."""
    if not isinstance(payload, dict) or set(payload) != {
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
    }:
        raise InputValidationError(
            "CORRUPT_PRECOMPUTED_ENTRY", "Precomputed entry fields are invalid."
        )
    if payload.get("artifact_schema_version") != ARTIFACT_SCHEMA_VERSION:
        raise InputValidationError(
            "PRECOMPUTED_SCHEMA_MISMATCH", "Precomputed entry schema is unsupported."
        )
    if not _valid_timestamp(payload.get("created_at")) or not _valid_producer(
        payload.get("producer")
    ):
        raise InputValidationError(
            "CORRUPT_PRECOMPUTED_ENTRY", "Precomputed entry metadata is invalid."
        )
    if payload.get("profile_id") != expected_profile_id:
        raise InputValidationError(
            "PRECOMPUTED_PROFILE_MISMATCH", "Precomputed entry profile is stale."
        )
    if payload.get("pdb_id") != expected_pdb_id:
        raise InputValidationError(
            "PRECOMPUTED_PDB_ID_MISMATCH",
            "Precomputed entry ID does not match its requested path.",
        )

    identity_raw = payload.get("structure_identity")
    if not isinstance(identity_raw, dict):
        raise InputValidationError(
            "CORRUPT_PRECOMPUTED_ENTRY", "Precomputed structure identity is missing."
        )
    try:
        identity = StructureIdentity.from_mapping(identity_raw)
    except (TypeError, ValueError) as exc:
        raise InputValidationError(
            "CORRUPT_PRECOMPUTED_ENTRY", "Precomputed structure identity is invalid."
        ) from exc
    if (
        identity.source != "pdb"
        or identity.canonical_id != expected_pdb_id
        or identity_raw != identity.as_dict()
    ):
        raise InputValidationError(
            "PRECOMPUTED_PDB_ID_MISMATCH",
            "Precomputed identity does not match its requested path.",
        )

    source = payload.get("source")
    if (
        not isinstance(source, dict)
        or set(source) != {"sha256", "size_bytes", "scope"}
        or not re.fullmatch(r"[a-f0-9]{64}", str(source.get("sha256") or ""))
        or not _valid_nonnegative_int(source.get("size_bytes"))
        or source.get("scope") != SOURCE_SCOPE
    ):
        raise InputValidationError(
            "CORRUPT_PRECOMPUTED_ENTRY", "Precomputed source metadata is invalid."
        )

    geometry = payload.get("geometry")
    if not isinstance(geometry, dict) or set(geometry) != {"structure", "interactions"}:
        raise InputValidationError(
            "CORRUPT_PRECOMPUTED_ENTRY", "Precomputed geometry is invalid."
        )
    structure = geometry.get("structure")
    allowed_structure_fields = {
        "file_path",
        "pdb_id",
        "atom_data",
        "structure_identity",
        "input_format",
        "input_kind",
        "embedded_annotation_counts",
    }
    if not isinstance(structure, dict) or not set(structure).issubset(
        allowed_structure_fields
    ):
        raise InputValidationError(
            "CORRUPT_PRECOMPUTED_ENTRY", "Precomputed structure metadata is invalid."
        )
    if (
        structure.get("pdb_id") != identity.display_id
        or structure.get("file_path") != f"pdb:{expected_pdb_id}"
        or structure.get("structure_identity") != identity.as_dict()
    ):
        raise InputValidationError(
            "CORRUPT_PRECOMPUTED_ENTRY", "Precomputed structure labels are invalid."
        )

    annotations = payload.get("annotations")
    if not isinstance(annotations, dict) or set(annotations) != {"references", "chains"}:
        raise InputValidationError(
            "CORRUPT_PRECOMPUTED_ENTRY", "Precomputed annotations are invalid."
        )
    references = annotations.get("references")
    expected_manifest_id = expected_profile["annotations"]["reference_manifest_id"]
    if (
        not isinstance(references, dict)
        or set(references) != {"manifest_id"}
        or references.get("manifest_id") != expected_manifest_id
    ):
        raise InputValidationError(
            "PRECOMPUTED_PROFILE_MISMATCH", "Precomputed references are stale."
        )

    geometry_chains = structure.get("atom_data")
    annotation_chains = annotations.get("chains")
    if (
        not isinstance(geometry_chains, list)
        or not geometry_chains
        or len(geometry_chains) > MAX_CHAINS_PER_ENTRY
        or not isinstance(annotation_chains, list)
        or len(annotation_chains) != len(geometry_chains)
    ):
        raise InputValidationError(
            "CORRUPT_PRECOMPUTED_ENTRY", "Precomputed chain collection is invalid."
        )
    annotation_map: dict[str, dict[str, Any]] = {}
    for annotation in annotation_chains:
        if (
            not isinstance(annotation, dict)
            or not set(annotation).issubset(ANNOTATION_CHAIN_FIELDS | {"unique_chain_id"})
        ):
            raise InputValidationError(
                "CORRUPT_PRECOMPUTED_ENTRY", "Precomputed chain annotation is invalid."
            )
        unique_id = annotation.get("unique_chain_id")
        if not isinstance(unique_id, str) or unique_id in annotation_map:
            raise InputValidationError(
                "CORRUPT_PRECOMPUTED_ENTRY", "Precomputed annotation IDs are invalid."
            )
        annotation_map[unique_id] = annotation

    chains: list[dict[str, Any]] = []
    chain_ids: set[str] = set()
    for raw_chain in geometry_chains:
        if not isinstance(raw_chain, dict) or not set(raw_chain).issubset(
            GEOMETRY_CHAIN_FIELDS
        ):
            raise InputValidationError(
                "CORRUPT_PRECOMPUTED_ENTRY", "Precomputed geometry chain is invalid."
            )
        chain_id = raw_chain.get("chain_id")
        unique_id = raw_chain.get("unique_chain_id")
        if (
            not isinstance(chain_id, str)
            or len(chain_id) > 256
            or unique_id != f"{identity.display_id}:{chain_id}"
            or unique_id in chain_ids
            or unique_id not in annotation_map
            or raw_chain.get("model_index") != 1
        ):
            raise InputValidationError(
                "CORRUPT_PRECOMPUTED_ENTRY", "Precomputed chain ID is invalid."
            )
        try:
            structured_chain = ChainIdentity.from_mapping(raw_chain.get("chain_identity", {}))
        except (TypeError, ValueError) as exc:
            raise InputValidationError(
                "CORRUPT_PRECOMPUTED_ENTRY", "Structured chain identity is invalid."
            ) from exc
        if (
            structured_chain.structure_key != expected_pdb_id
            or structured_chain.structure_display_id != identity.display_id
            or structured_chain.chain_id != chain_id
            or structured_chain.model_index != 1
            or raw_chain.get("chain_identity") != structured_chain.as_dict()
            or raw_chain.get("structure_key") != expected_pdb_id
            or raw_chain.get("structure_display_id") != identity.display_id
        ):
            raise InputValidationError(
                "CORRUPT_PRECOMPUTED_ENTRY", "Structured chain identity is inconsistent."
            )
        if (
            not isinstance(raw_chain.get("sequence"), str)
            or len(raw_chain["sequence"]) > 1_000_000
            or not isinstance(raw_chain.get("is_hetatm"), bool)
            or not _valid_nonnegative_int(raw_chain.get("aa_len"))
            or not _valid_nonnegative_int(raw_chain.get("nt_len"))
        ):
            raise InputValidationError(
                "CORRUPT_PRECOMPUTED_ENTRY", "Precomputed chain data is invalid."
            )
        _validate_embedded_annotations(raw_chain, chain_id)
        accessions = raw_chain.get("embedded_uniprot_accessions", [])
        if accessions and (
            not isinstance(accessions, list)
            or len(accessions) > 1_000
            or any(not isinstance(value, str) or len(value) > 100 for value in accessions)
        ):
            raise InputValidationError(
                "CORRUPT_PRECOMPUTED_ENTRY", "Embedded UniProt accessions are invalid."
            )

        chain = dict(raw_chain)
        annotation = dict(annotation_map[unique_id])
        annotation.pop("unique_chain_id", None)
        chain.update(annotation)
        if chain.get("molecule_type") not in ALLOWED_MOLECULE_TYPES:
            raise InputValidationError(
                "CORRUPT_PRECOMPUTED_ENTRY", "Precomputed molecule type is invalid."
            )
        if not _bounded_optional_text(chain.get("molecule_name"), maximum=10_000):
            raise InputValidationError(
                "CORRUPT_PRECOMPUTED_ENTRY", "Precomputed molecule name is invalid."
            )
        uniprot_id = chain.get("uniprot_id")
        if uniprot_id not in (None, ""):
            from ..uniprot_matcher import _is_uniprotkb_accession

            if not isinstance(uniprot_id, str) or not _is_uniprotkb_accession(uniprot_id):
                raise InputValidationError(
                    "CORRUPT_PRECOMPUTED_ENTRY", "Precomputed UniProt accession is invalid."
                )
        external_accessions = chain.get("external_sifts_accessions", [])
        if external_accessions and (
            not isinstance(external_accessions, list)
            or len(external_accessions) > 1_000
            or any(
                not isinstance(value, str) or len(value) > 100
                for value in external_accessions
            )
        ):
            raise InputValidationError(
                "CORRUPT_PRECOMPUTED_ENTRY", "External SIFTS accessions are invalid."
            )
        warnings = chain.get("annotation_warnings", [])
        if warnings and (
            not isinstance(warnings, list)
            or len(warnings) > 100
            or any(
                not isinstance(value, dict)
                or set(value) != {"code", "message"}
                or not isinstance(value.get("code"), str)
                or len(value["code"]) > 100
                or not isinstance(value.get("message"), str)
                or len(value["message"]) > 2_000
                for value in warnings
            )
        ):
            raise InputValidationError(
                "CORRUPT_PRECOMPUTED_ENTRY", "Annotation warnings are invalid."
            )
        for field in ANNOTATION_CHAIN_FIELDS - {
            "molecule_name",
            "molecule_type",
            "uniprot_id",
            "external_sifts_accessions",
            "annotation_warnings",
        }:
            if not _bounded_optional_text(chain.get(field), maximum=1_000):
                raise InputValidationError(
                    "CORRUPT_PRECOMPUTED_ENTRY", f"Precomputed chain field {field} is invalid."
                )
        chains.append(chain)
        chain_ids.add(unique_id)

    interactions = geometry.get("interactions")
    if not isinstance(interactions, list) or len(interactions) > MAX_EDGES_PER_ENTRY:
        raise InputValidationError(
            "CORRUPT_PRECOMPUTED_ENTRY", "Precomputed interactions are invalid."
        )
    chains_by_id = {chain["unique_chain_id"]: chain for chain in chains}
    seen_edges: set[tuple[str, str]] = set()
    filters = expected_profile["interaction_filters"]
    for edge in interactions:
        if not isinstance(edge, dict) or not set(edge).issubset(
            {
                "file_path",
                "structure_key",
                "model_index",
                "chain_a",
                "chain_b",
                "all_atoms_count",
                "ca_neighbors",
                "interaction_type",
            }
        ):
            raise InputValidationError(
                "CORRUPT_PRECOMPUTED_ENTRY", "Precomputed interaction is invalid."
            )
        chain_a, chain_b = edge.get("chain_a"), edge.get("chain_b")
        unordered = tuple(sorted((str(chain_a), str(chain_b))))
        if (
            chain_a not in chain_ids
            or chain_b not in chain_ids
            or chain_a == chain_b
            or unordered in seen_edges
            or edge.get("interaction_type") not in ALLOWED_INTERACTION_TYPES
            or not _valid_nonnegative_int(edge.get("all_atoms_count"))
        ):
            raise InputValidationError(
                "CORRUPT_PRECOMPUTED_ENTRY", "Precomputed interaction content is invalid."
            )
        seen_edges.add(unordered)
        expected_type = determine_interaction_type(
            chains_by_id[chain_a].get("molecule_type"),
            chains_by_id[chain_b].get("molecule_type"),
        )
        if edge.get("interaction_type") != expected_type:
            raise InputValidationError(
                "CORRUPT_PRECOMPUTED_ENTRY", "Interaction type disagrees with its chains."
            )
        if expected_type == "Protein-Protein":
            if (
                not _valid_nonnegative_int(edge.get("ca_neighbors"))
                or edge["ca_neighbors"] < filters["protein_protein_min_ca_neighbors"]
                or edge["all_atoms_count"]
                < filters["protein_protein_min_all_atom_contacts"]
            ):
                raise InputValidationError(
                    "CORRUPT_PRECOMPUTED_ENTRY", "Protein contact counts are below the profile."
                )
        else:
            minimum = (
                filters["protein_nucleic_acid_min_all_atom_contacts"]
                if str(expected_type).startswith("Protein-")
                else filters["nucleic_acid_min_all_atom_contacts"]
            )
            if "ca_neighbors" in edge or edge["all_atoms_count"] < minimum:
                raise InputValidationError(
                    "CORRUPT_PRECOMPUTED_ENTRY", "Contact counts disagree with the profile."
                )
        if (
            edge.get("file_path") != f"pdb:{expected_pdb_id}"
            or edge.get("structure_key") not in (None, expected_pdb_id)
            or edge.get("model_index", 1) != 1
        ):
            raise InputValidationError(
                "CORRUPT_PRECOMPUTED_ENTRY", "Precomputed interaction labels are invalid."
            )

    counts = payload.get("counts")
    if (
        not isinstance(counts, dict)
        or set(counts) != {"chains", "interactions"}
        or counts.get("chains") != len(chains)
        or counts.get("interactions") != len(interactions)
    ):
        raise InputValidationError(
            "CORRUPT_PRECOMPUTED_ENTRY", "Precomputed counts do not match content."
        )
    return dict(payload)


def materialize_entry(
    payload: Mapping[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, Any]]:
    """Merge validated annotation chains into the compact runtime structure."""
    structure = dict(payload["geometry"]["structure"])
    structure["atom_data"] = _materialized_chains(payload)
    interactions = [dict(edge) for edge in payload["geometry"]["interactions"]]
    references = dict(payload["annotations"]["references"])
    return structure, interactions, references


__all__ = [
    "ARTIFACT_SCHEMA_VERSION",
    "ANNOTATION_CHAIN_FIELDS",
    "GEOMETRY_CHAIN_FIELDS",
    "MAX_CHAINS_PER_ENTRY",
    "MAX_COMPRESSED_ENTRY_BYTES",
    "MAX_DECOMPRESSED_ENTRY_BYTES",
    "MAX_EDGES_PER_ENTRY",
    "MAX_REQUEST_CHAINS",
    "MAX_REQUEST_EDGES",
    "MAX_REQUEST_EXPANDED_BYTES",
    "SCIENTIFIC_PIPELINE_VERSION",
    "SOURCE_SCOPE",
    "canonical_json",
    "manifest_document",
    "materialize_entry",
    "normalize_pdb_id",
    "normalize_requested_ids",
    "producer",
    "profile_id",
    "scientific_profile",
    "utc_timestamp",
    "validate_entry",
    "validate_manifest",
]
