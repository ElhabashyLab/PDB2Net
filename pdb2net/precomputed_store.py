"""Versioned, portable per-PDB graph cache and assembly pipeline.

The store is intentionally independent of the webserver.  It persists the
small scientific primitives needed to recreate per-PDB and combined networks.
Schema 2 separates reference-independent compact geometry/interactions from
refreshable external annotations and retains compact embedded SIFTS segments.
Raw atom coordinates and detailed atom-pair tables are deliberately not cached.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import os
import re
import stat
import tempfile
import time
from collections import Counter
from contextlib import ExitStack, contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from . import __version__
from .config_loader import config
from .distances import calculate_distances_with_ckdtree, coords_cache, tree_cache
from .file_parser import get_structure_identity, is_valid_file
from .input_contract import InputValidationError
from .network_annotations import network_annotation_config
from .structure_identity import (
    ChainIdentity,
    StructureIdentity,
    canonical_pdb_id,
    identity_from_official_id,
    pdb_archive_shard,
)


ARTIFACT_SCHEMA_VERSION = "2"
SCIENTIFIC_PIPELINE_VERSION = "pdb2net-asu-first-standard-interactions-v3"
ANNOTATION_PIPELINE_VERSION = "pdb2net-sifts-fasta-search-annotations-v3"
SOURCE_SCOPE = "asymmetric_unit"
MAX_COMPRESSED_ENTRY_BYTES = 64 * 1024 * 1024
MAX_DECOMPRESSED_ENTRY_BYTES = 128 * 1024 * 1024
MAX_CHAINS_PER_ENTRY = 50_000
MAX_EDGES_PER_ENTRY = 500_000
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
MAX_REQUEST_EXPANDED_BYTES = 256 * 1024 * 1024
MAX_REQUEST_CHAINS = 200_000
MAX_REQUEST_EDGES = 2_000_000
POPULATION_LOCK_TIMEOUT_SECONDS = 300.0
MAX_SIMULTANEOUS_ENTRY_LOCKS = 64

ALLOWED_INTERACTION_TYPES = {
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
ALLOWED_MOLECULE_TYPES = {
    "Protein",
    "DNA",
    "RNA",
    "DNA/RNA",
    "Nucleic Acid",
    "Unknown",
    "UNKNOWN",
}
CHAIN_FIELDS = {
    "chain_id",
    "unique_chain_id",
    "chain_identity",
    "structure_key",
    "structure_display_id",
    "model_index",
    "molecule_name",
    "molecule_type",
    "sequence",
    "is_hetatm",
    "uniprot_id",
    "annotation_source",
    "matched_database",
    "matched_id",
    "representative_accession",
    "annotation_confidence",
    "aa_len",
    "nt_len",
    "embedded_annotations",
    "embedded_annotation_source",
    "embedded_uniprot_status",
    "embedded_uniprot_accessions",
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


def normalize_pdb_id(value: object) -> str:
    """Return a canonical lowercase wwPDB extended identifier."""
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
    """Normalize requested PDB IDs and reject canonical identity collisions."""
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


def _canonical_json(value: object) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")


def _require_first_model() -> None:
    policy = str(config.get("structure_model_policy") or "first").strip().lower()
    if policy != "first":
        raise InputValidationError(
            "PRECOMPUTED_MODEL_POLICY_UNSUPPORTED",
            "Precomputed store schema 2 supports structure_model_policy='first' only.",
        )


def _reference_manifest_id() -> str:
    value = str(config.get("reference_manifest_id") or "").strip()
    if not value:
        raise InputValidationError(
            "REFERENCE_MANIFEST_ID_REQUIRED",
            "Precomputed runs require a non-empty reference_manifest_id so stale annotations cannot be reused.",
        )
    return value


def _profile_interaction_filters() -> dict[str, int]:
    filters = config.get("interaction_filters", {})
    if not isinstance(filters, Mapping):
        filters = {}
    return {
        "protein_protein_min_ca_neighbors": int(
            filters.get("protein_protein_min_ca_neighbors", 10)
        ),
        "protein_protein_min_all_atom_contacts": int(
            filters.get("protein_protein_min_all_atom_contacts", 1)
        ),
        "protein_nucleic_acid_min_all_atom_contacts": int(
            filters.get("protein_nucleic_acid_min_all_atom_contacts", 1)
        ),
        "nucleic_acid_min_all_atom_contacts": int(
            filters.get("nucleic_acid_min_all_atom_contacts", 1)
        ),
    }


def scientific_profile() -> dict[str, Any]:
    """Return the geometry profile, independent of weekly references/tooltips."""
    import gemmi

    _require_first_model()
    thresholds = config.get("distance_thresholds", {})
    if not isinstance(thresholds, Mapping):
        thresholds = {}
    profile = {
        "artifact_schema_version": ARTIFACT_SCHEMA_VERSION,
        "scientific_pipeline_version": SCIENTIFIC_PIPELINE_VERSION,
        "source_scope": SOURCE_SCOPE,
        "structure_model_policy": "first",
        "parser": {
            "semantics": "validated-single-document-gemmi-heavy-atoms-no-hydrogen-or-deuterium-v3",
            "gemmi_version": str(getattr(gemmi, "__version__", "unknown")),
        },
        "distance_thresholds": {
            "ca_radius": float(thresholds.get("ca_radius", 12.0)),
            "all_atoms_radius": float(thresholds.get("all_atoms_radius", 5.0)),
        },
        "interaction_filters": _profile_interaction_filters(),
    }
    # Search thresholds are naturally represented as tuples in Python.  Store
    # one canonical JSON-native shape so an in-memory profile compares exactly
    # with the same profile read back from manifest.json.
    return json.loads(_canonical_json(profile))


def annotation_profile() -> dict[str, Any]:
    """Return settings that invalidate annotations but not cached geometry."""
    from .uniprot_matcher import _search_policy_payload

    annotation_cfg = network_annotation_config()
    profile = {
        "annotation_pipeline_version": ANNOTATION_PIPELINE_VERSION,
        "reference_manifest_id": _reference_manifest_id(),
        "use_embedded_sifts": annotation_cfg["use_embedded_sifts"],
        "search_policy": _search_policy_payload(),
    }
    return json.loads(_canonical_json(profile))


def annotation_profile_id(profile: Mapping[str, Any] | None = None) -> str:
    return hashlib.sha256(_canonical_json(dict(profile or annotation_profile()))).hexdigest()


def profile_id(profile: Mapping[str, Any] | None = None) -> str:
    return hashlib.sha256(_canonical_json(dict(profile or scientific_profile()))).hexdigest()


def profile_root(store: Path | str, profile_hash: str | None = None) -> Path:
    selected = profile_hash or profile_id()
    if not re.fullmatch(r"[a-f0-9]{64}", selected):
        raise InputValidationError("INVALID_CACHE_PROFILE", "Cache profile ID is malformed.")
    raw_store = Path(store).expanduser()
    if raw_store.is_symlink():
        raise InputValidationError("UNSAFE_CACHE_PATH", "Precomputed store root must not be a symlink.")
    store_root = raw_store.resolve()
    profiles = store_root / "profiles"
    selected_root = profiles / selected
    for candidate in (profiles, selected_root):
        if candidate.is_symlink() or (candidate.exists() and not candidate.is_dir()):
            raise InputValidationError("UNSAFE_CACHE_PATH", "Cache profile path is unsafe.")
    return selected_root


def entry_path(store: Path | str, pdb_id: object, profile_hash: str | None = None) -> Path:
    normalized = normalize_pdb_id(pdb_id)
    return profile_root(store, profile_hash) / "entries" / pdb_archive_shard(normalized) / f"{normalized}.json.gz"


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _source_fingerprint(path: Path) -> dict[str, Any]:
    stat_result = path.stat()
    return {
        "name": path.name,
        "sha256": _sha256_file(path),
        "size_bytes": stat_result.st_size,
        "scope": SOURCE_SCOPE,
    }


def _assert_no_symlink_components(root: Path, target: Path) -> None:
    try:
        relative = target.relative_to(root)
    except ValueError as exc:
        raise InputValidationError("UNSAFE_CACHE_PATH", "Cache path escaped its profile root.") from exc
    current = root
    for part in relative.parts:
        current = current / part
        if current.is_symlink():
            raise InputValidationError("UNSAFE_CACHE_PATH", "Cache path contains a symlink.")


def _atomic_write(path: Path, data: bytes, *, safe_root: Path) -> None:
    _assert_no_symlink_components(safe_root, path.parent)
    path.parent.mkdir(parents=True, exist_ok=True)
    _assert_no_symlink_components(safe_root, path.parent)
    if path.is_symlink():
        raise InputValidationError("UNSAFE_CACHE_PATH", f"Cache target is a symlink: {path.name}")
    fd, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(fd, "wb") as handle:
            if hasattr(os, "fchmod"):
                os.fchmod(handle.fileno(), 0o640)
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        directory_flags = os.O_RDONLY
        if hasattr(os, "O_DIRECTORY"):
            directory_flags |= os.O_DIRECTORY
        directory_fd = os.open(path.parent, directory_flags)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        if temporary.exists():
            temporary.unlink()


def ensure_profile_manifest(store: Path | str) -> tuple[str, Path]:
    """Create or verify the immutable manifest for the active profile."""
    profile = scientific_profile()
    selected = profile_id(profile)
    root = profile_root(store, selected)
    manifest_path = root / "manifest.json"
    expected = {
        "artifact_schema_version": ARTIFACT_SCHEMA_VERSION,
        "profile_id": selected,
        "profile": profile,
    }
    if manifest_path.exists():
        if manifest_path.is_symlink() or manifest_path.stat().st_size > 1_000_000:
            raise InputValidationError("UNSAFE_CACHE_MANIFEST", "Cache profile manifest is unsafe.")
        try:
            current = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            raise InputValidationError("CORRUPT_CACHE_MANIFEST", "Cache profile manifest is invalid.") from exc
        if current != expected:
            raise InputValidationError(
                "CACHE_PROFILE_MISMATCH",
                "The cache profile manifest does not match the active scientific profile.",
            )
    else:
        _atomic_write(
            manifest_path,
            json.dumps(expected, indent=2, sort_keys=True).encode("utf-8"),
            safe_root=root,
        )
    return selected, root


def _bounded_gzip_json(
    path: Path,
    *,
    expanded_budget: int = MAX_DECOMPRESSED_ENTRY_BYTES,
) -> tuple[Any, int]:
    if path.is_symlink() or not path.is_file():
        raise InputValidationError("UNSAFE_CACHE_ENTRY", f"Cache entry is missing or a symlink: {path.name}")
    if path.stat().st_size > MAX_COMPRESSED_ENTRY_BYTES:
        raise InputValidationError("CACHE_ENTRY_TOO_LARGE", f"Compressed cache entry exceeds the size limit: {path.name}")
    read_limit = min(MAX_DECOMPRESSED_ENTRY_BYTES, max(0, int(expanded_budget)))
    if read_limit == 0:
        raise InputValidationError(
            "PRECOMPUTED_REQUEST_LIMIT_EXCEEDED",
            "Requested cache entries exceed the aggregate expanded-size limit.",
        )
    try:
        with gzip.open(path, "rb") as handle:
            raw = handle.read(read_limit + 1)
    except (OSError, EOFError) as exc:
        raise InputValidationError("CORRUPT_CACHE_ENTRY", f"Cannot decompress cache entry: {path.name}") from exc
    if len(raw) > read_limit:
        code = (
            "CACHE_ENTRY_TOO_LARGE"
            if read_limit == MAX_DECOMPRESSED_ENTRY_BYTES
            else "PRECOMPUTED_REQUEST_LIMIT_EXCEEDED"
        )
        raise InputValidationError(code, f"Expanded cache data exceeds the size limit at {path.name}.")
    try:
        return json.loads(raw.decode("utf-8")), len(raw)
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise InputValidationError("CORRUPT_CACHE_ENTRY", f"Cache entry is not valid UTF-8 JSON: {path.name}") from exc


def _valid_nonnegative_int(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value >= 0


def _bounded_optional_text(value: object, *, maximum: int) -> bool:
    return value is None or (isinstance(value, str) and len(value) <= maximum)


def validate_entry(
    payload: object,
    *,
    expected_pdb_id: str,
    expected_profile_id: str,
    require_current_annotations: bool = True,
) -> dict[str, Any]:
    """Validate an untrusted schema-2 entry and return a materialized view."""
    if not isinstance(payload, dict):
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache entry must be a JSON object.")
    expected_top_level = {
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
    }
    if set(payload) != expected_top_level:
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache entry top-level fields are invalid.")
    if not isinstance(payload.get("pdb2net_version"), str) or not payload["pdb2net_version"]:
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache package version is invalid.")
    if not isinstance(payload.get("created_at"), str) or len(payload["created_at"]) > 100:
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache creation timestamp is invalid.")
    if payload.get("artifact_schema_version") != ARTIFACT_SCHEMA_VERSION:
        raise InputValidationError("CACHE_SCHEMA_MISMATCH", "Cache entry schema is unsupported.")
    if payload.get("scientific_pipeline_version") != SCIENTIFIC_PIPELINE_VERSION:
        raise InputValidationError("CACHE_PIPELINE_MISMATCH", "Cache entry scientific pipeline is stale.")
    if payload.get("profile_id") != expected_profile_id:
        raise InputValidationError("CACHE_PROFILE_MISMATCH", "Cache entry geometry profile is stale.")
    if payload.get("pdb_id") != expected_pdb_id:
        raise InputValidationError("CACHE_PDB_ID_MISMATCH", "Cache entry PDB ID does not match its requested path.")

    identity_raw = payload.get("structure_identity")
    if not isinstance(identity_raw, dict):
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache structure identity is missing.")
    try:
        identity = StructureIdentity.from_mapping(identity_raw)
    except (TypeError, ValueError) as exc:
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache structure identity is invalid.") from exc
    if identity.source != "pdb" or identity.canonical_id != expected_pdb_id:
        raise InputValidationError("CACHE_PDB_ID_MISMATCH", "Cache identity does not match its requested path.")
    if identity_raw != identity.as_dict():
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache identity fields are invalid.")

    annotations = payload.get("annotations")
    if not isinstance(annotations, dict) or set(annotations) != {"profile_id", "references", "chains"}:
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache annotation section is invalid.")
    stored_annotation_profile = str(annotations.get("profile_id") or "")
    if not re.fullmatch(r"[a-f0-9]{64}", stored_annotation_profile):
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache annotation profile is malformed.")
    references = annotations.get("references")
    if (
        not isinstance(references, dict)
        or set(references) != {"manifest_id"}
        or not isinstance(references.get("manifest_id"), str)
    ):
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache annotation references are invalid.")
    if require_current_annotations and stored_annotation_profile != annotation_profile_id():
        raise InputValidationError("CACHE_ANNOTATION_PROFILE_MISMATCH", "Cache annotations are stale.")
    if require_current_annotations and references.get("manifest_id") != _reference_manifest_id():
        raise InputValidationError("CACHE_REFERENCE_MISMATCH", "Cache annotations use stale references.")

    source = payload.get("source")
    if (
        not isinstance(source, dict)
        or set(source) != {"name", "sha256", "size_bytes", "scope"}
        or not re.fullmatch(r"[a-f0-9]{64}", str(source.get("sha256") or ""))
    ):
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache entry source fingerprint is invalid.")
    if source.get("scope") != SOURCE_SCOPE or not _valid_nonnegative_int(source.get("size_bytes")):
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache entry source metadata is invalid.")
    source_name = source.get("name")
    if not isinstance(source_name, str) or not source_name or len(source_name) > 512 or Path(source_name).name != source_name:
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache source name is not portable.")

    geometry = payload.get("geometry")
    if not isinstance(geometry, dict) or set(geometry) != {"structure", "interactions"}:
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache geometry section is invalid.")
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
    if not isinstance(structure, dict) or not set(structure).issubset(allowed_structure_fields):
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache structure metadata is invalid.")
    if structure.get("pdb_id") != identity.display_id or structure.get("file_path") != f"pdb:{expected_pdb_id}":
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache structure label is invalid.")
    if structure.get("structure_identity") != identity.as_dict():
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache structure identity copies disagree.")

    geometry_chains = structure.get("atom_data")
    annotation_chains = annotations.get("chains")
    if not isinstance(geometry_chains, list) or not geometry_chains or len(geometry_chains) > MAX_CHAINS_PER_ENTRY:
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache chain collection is invalid.")
    if not isinstance(annotation_chains, list) or len(annotation_chains) != len(geometry_chains):
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache chain annotations are invalid.")
    annotations_by_chain: dict[str, dict[str, Any]] = {}
    for annotation in annotation_chains:
        if not isinstance(annotation, dict) or not set(annotation).issubset(ANNOTATION_CHAIN_FIELDS | {"unique_chain_id"}):
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache chain annotation contains unsupported fields.")
        unique_id = annotation.get("unique_chain_id")
        if not isinstance(unique_id, str) or unique_id in annotations_by_chain:
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache chain annotation IDs are invalid.")
        annotations_by_chain[unique_id] = annotation

    chains: list[dict[str, Any]] = []
    chain_ids: set[str] = set()
    geometry_allowed = CHAIN_FIELDS - ANNOTATION_CHAIN_FIELDS
    for raw_chain in geometry_chains:
        if not isinstance(raw_chain, dict) or not set(raw_chain).issubset(geometry_allowed):
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache chains must be compact geometry objects.")
        chain_id = raw_chain.get("chain_id")
        unique_id = raw_chain.get("unique_chain_id")
        if not isinstance(chain_id, str) or len(chain_id) > 256 or unique_id != f"{identity.display_id}:{chain_id}":
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache chain ID is invalid.")
        if unique_id in chain_ids or unique_id not in annotations_by_chain:
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache chain IDs must be unique and annotated.")
        if raw_chain.get("model_index") != 1:
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache schema 2 accepts first-model chains only.")
        try:
            structured_chain = ChainIdentity.from_mapping(raw_chain.get("chain_identity", {}))
        except (TypeError, ValueError) as exc:
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache structured chain identity is invalid.") from exc
        if (
            structured_chain.structure_key != expected_pdb_id
            or structured_chain.structure_display_id != identity.display_id
            or structured_chain.chain_id != chain_id
            or structured_chain.model_index != 1
            or raw_chain.get("chain_identity") != structured_chain.as_dict()
            or raw_chain.get("structure_key") != expected_pdb_id
            or raw_chain.get("structure_display_id") != identity.display_id
        ):
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache structured chain identity is inconsistent.")
        sequence = raw_chain.get("sequence")
        if not isinstance(sequence, str) or len(sequence) > 1_000_000:
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache chain sequence is invalid.")
        if not isinstance(raw_chain.get("is_hetatm"), bool):
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache chain hetero-atom flag is invalid.")
        if not _valid_nonnegative_int(raw_chain.get("aa_len")) or not _valid_nonnegative_int(raw_chain.get("nt_len")):
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache chain lengths are invalid.")
        embedded = raw_chain.get("embedded_annotations", {})
        if embedded:
            if not isinstance(embedded, dict) or any(key not in {"uniprot", "pfam", "cath", "scop2"} for key in embedded):
                raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache embedded annotations are invalid.")
            for database, segments in embedded.items():
                if not isinstance(segments, list) or len(segments) > MAX_EMBEDDED_SEGMENTS_PER_CHAIN:
                    raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache embedded annotation segments are invalid.")
                for segment in segments:
                    if not isinstance(segment, dict) or not set(segment).issubset(EMBEDDED_SEGMENT_FIELDS):
                        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache embedded annotation segment is invalid.")
                    for value in segment.values():
                        if not isinstance(value, (str, int, float, bool)) or (
                            isinstance(value, str) and len(value) > 2_000
                        ) or (
                            isinstance(value, float) and not math.isfinite(value)
                        ):
                            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache embedded annotation value is invalid.")
                    if segment.get("database") != database:
                        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache embedded annotation database is invalid.")
                    if segment.get("chain_id") != chain_id or not isinstance(segment.get("accession"), str):
                        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache embedded annotation chain is invalid.")
                    pdb_start, pdb_end = segment.get("pdb_start"), segment.get("pdb_end")
                    if (
                        not isinstance(pdb_start, int)
                        or isinstance(pdb_start, bool)
                        or not isinstance(pdb_end, int)
                        or isinstance(pdb_end, bool)
                        or pdb_start < 1
                        or pdb_end < pdb_start
                    ):
                        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache embedded annotation range is invalid.")
                    if database == "uniprot":
                        from .reference_data import is_uniprot_accession

                        database_start = segment.get("database_start")
                        database_end = segment.get("database_end")
                        identity_value = segment.get("identity")
                        if (
                            not is_uniprot_accession(segment.get("accession"))
                            or not isinstance(segment.get("best_mapping"), bool)
                            or not isinstance(identity_value, (int, float))
                            or isinstance(identity_value, bool)
                            or not math.isfinite(float(identity_value))
                            or not 0 <= float(identity_value) <= 1
                            or not isinstance(database_start, int)
                            or isinstance(database_start, bool)
                            or not isinstance(database_end, int)
                            or isinstance(database_end, bool)
                            or database_start < 1
                            or database_end < database_start
                        ):
                            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache embedded UniProt segment is invalid.")
        accessions = raw_chain.get("embedded_uniprot_accessions", [])
        if accessions:
            if (
                not isinstance(accessions, list)
                or len(accessions) > 1_000
                or any(not isinstance(value, str) or len(value) > 100 for value in accessions)
            ):
                raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache embedded UniProt accessions are invalid.")

        chain = dict(raw_chain)
        annotation = dict(annotations_by_chain[unique_id])
        annotation.pop("unique_chain_id", None)
        chain.update(annotation)
        if chain.get("molecule_type") not in ALLOWED_MOLECULE_TYPES:
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache molecule type is invalid.")
        if not _bounded_optional_text(chain.get("molecule_name"), maximum=10_000):
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache molecule name is invalid.")
        uniprot_id = chain.get("uniprot_id")
        if uniprot_id not in (None, ""):
            from .uniprot_matcher import _is_uniprotkb_accession
            if not isinstance(uniprot_id, str) or not _is_uniprotkb_accession(uniprot_id):
                raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache UniProt accession is invalid.")
        external_accessions = chain.get("external_sifts_accessions", [])
        if external_accessions and (
            not isinstance(external_accessions, list)
            or len(external_accessions) > 1_000
            or any(not isinstance(value, str) or len(value) > 100 for value in external_accessions)
        ):
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache external SIFTS accessions are invalid.")
        annotation_warnings = chain.get("annotation_warnings", [])
        if annotation_warnings and (
            not isinstance(annotation_warnings, list)
            or len(annotation_warnings) > 100
            or any(
                not isinstance(value, dict)
                or set(value) != {"code", "message"}
                or not isinstance(value.get("code"), str)
                or len(value["code"]) > 100
                or not isinstance(value.get("message"), str)
                or len(value["message"]) > 2_000
                for value in annotation_warnings
            )
        ):
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache annotation warnings are invalid.")
        for field in ANNOTATION_CHAIN_FIELDS - {
            "molecule_name",
            "molecule_type",
            "uniprot_id",
            "external_sifts_accessions",
            "annotation_warnings",
        }:
            if not _bounded_optional_text(chain.get(field), maximum=1_000):
                raise InputValidationError("CORRUPT_CACHE_ENTRY", f"Cache chain field {field} is invalid.")
        chains.append(chain)
        chain_ids.add(unique_id)

    interactions = geometry.get("interactions")
    if not isinstance(interactions, list) or len(interactions) > MAX_EDGES_PER_ENTRY:
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache interaction collection is invalid.")
    seen_edges: set[tuple[str, str]] = set()
    chains_by_id = {chain["unique_chain_id"]: chain for chain in chains}
    profile_filters = _profile_interaction_filters()
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
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache interactions must be objects.")
        chain_a, chain_b = edge.get("chain_a"), edge.get("chain_b")
        if chain_a not in chain_ids or chain_b not in chain_ids or chain_a == chain_b:
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache interaction endpoints are invalid.")
        unordered = tuple(sorted((chain_a, chain_b)))
        if unordered in seen_edges:
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache interactions contain duplicate chain pairs.")
        seen_edges.add(unordered)
        if edge.get("interaction_type") not in ALLOWED_INTERACTION_TYPES:
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache interaction type is invalid.")
        from .distances import determine_interaction_type
        expected_type = determine_interaction_type(
            chains_by_id[chain_a].get("molecule_type"),
            chains_by_id[chain_b].get("molecule_type"),
        )
        if edge.get("interaction_type") != expected_type:
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache interaction type disagrees with its chains.")
        if not _valid_nonnegative_int(edge.get("all_atoms_count")):
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache all-atom contact count is invalid.")
        if expected_type == "Protein-Protein":
            if (
                not _valid_nonnegative_int(edge.get("ca_neighbors"))
                or edge["ca_neighbors"] < profile_filters["protein_protein_min_ca_neighbors"]
                or edge["all_atoms_count"] < profile_filters["protein_protein_min_all_atom_contacts"]
            ):
                raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache protein contact counts are below profile filters.")
        else:
            minimum_key = (
                "protein_nucleic_acid_min_all_atom_contacts"
                if str(expected_type).startswith("Protein-")
                else "nucleic_acid_min_all_atom_contacts"
            )
            if "ca_neighbors" in edge or edge["all_atoms_count"] < profile_filters[minimum_key]:
                raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache contact counts disagree with profile filters.")
        if edge.get("file_path") != f"pdb:{expected_pdb_id}":
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache interaction source label is invalid.")
        if edge.get("structure_key") not in (None, expected_pdb_id):
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache interaction structure key is invalid.")
        if edge.get("model_index", 1) != 1:
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache interaction model index is invalid.")

    counts = payload.get("counts")
    if (
        not isinstance(counts, dict)
        or set(counts) != {"chains", "interactions"}
        or counts.get("chains") != len(chains)
        or counts.get("interactions") != len(interactions)
    ):
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache entry counts do not match its content.")
    materialized_structure = dict(structure)
    materialized_structure["atom_data"] = chains
    result = dict(payload)
    result["structure"] = materialized_structure
    result["interactions"] = interactions
    result["references"] = references
    return result


@contextmanager
def _entry_population_lock(
    store: Path | str,
    pdb_id: object,
    *,
    selected_profile: str | None = None,
    profile_directory: Path | None = None,
):
    """Serialize lazy population for one profile/PDB and release on process exit."""
    normalized = normalize_pdb_id(pdb_id)
    if selected_profile is None or profile_directory is None:
        selected_profile, profile_directory = ensure_profile_manifest(store)
    root = profile_directory
    lock_path = root / "locks" / pdb_archive_shard(normalized) / f"{normalized}.lock"
    _assert_no_symlink_components(root, lock_path.parent)
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    _assert_no_symlink_components(root, lock_path.parent)
    if lock_path.is_symlink():
        raise InputValidationError("UNSAFE_CACHE_PATH", "Cache population lock path is unsafe.")

    flags = os.O_CREAT | os.O_RDWR
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(lock_path, flags, 0o600)
    except OSError as exc:
        raise InputValidationError("CACHE_LOCK_UNAVAILABLE", "Cannot open the cache population lock.") from exc
    if not stat.S_ISREG(os.fstat(descriptor).st_mode):
        os.close(descriptor)
        raise InputValidationError("UNSAFE_CACHE_PATH", "Cache population lock is not a regular file.")
    if hasattr(os, "fchmod"):
        os.fchmod(descriptor, 0o660)

    acquired = False
    deadline = time.monotonic() + POPULATION_LOCK_TIMEOUT_SECONDS
    try:
        if os.name == "nt":
            import msvcrt

            if os.fstat(descriptor).st_size == 0:
                os.write(descriptor, b"\0")
            while True:
                try:
                    os.lseek(descriptor, 0, os.SEEK_SET)
                    msvcrt.locking(descriptor, msvcrt.LK_NBLCK, 1)
                    acquired = True
                    break
                except OSError:
                    if time.monotonic() >= deadline:
                        raise InputValidationError(
                            "PRECOMPUTED_ENTRY_BUSY",
                            f"Timed out waiting to populate {normalized.upper()}.",
                        )
                    time.sleep(0.1)
        else:
            import fcntl

            while True:
                try:
                    fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    acquired = True
                    break
                except BlockingIOError:
                    if time.monotonic() >= deadline:
                        raise InputValidationError(
                            "PRECOMPUTED_ENTRY_BUSY",
                            f"Timed out waiting to populate {normalized.upper()}.",
                        )
                    time.sleep(0.1)
        yield
    finally:
        if acquired:
            if os.name == "nt":
                import msvcrt

                os.lseek(descriptor, 0, os.SEEK_SET)
                msvcrt.locking(descriptor, msvcrt.LK_UNLCK, 1)
            else:
                import fcntl

                fcntl.flock(descriptor, fcntl.LOCK_UN)
        os.close(descriptor)


def load_entry(store: Path | str, pdb_id: object) -> dict[str, Any]:
    normalized = normalize_pdb_id(pdb_id)
    selected, _ = ensure_profile_manifest(store)
    return _load_entry_for_profile(store, normalized, selected)


def _load_entry_for_profile(
    store: Path | str,
    pdb_id: object,
    selected_profile: str,
) -> dict[str, Any]:
    payload, _expanded_bytes = _read_entry_for_profile(
        store,
        pdb_id,
        selected_profile,
    )
    return payload


def _read_entry_for_profile(
    store: Path | str,
    pdb_id: object,
    selected_profile: str,
    *,
    expanded_budget: int = MAX_DECOMPRESSED_ENTRY_BYTES,
    require_current_annotations: bool = True,
) -> tuple[dict[str, Any], int]:
    normalized = normalize_pdb_id(pdb_id)
    path = entry_path(store, normalized, selected_profile)
    _assert_no_symlink_components(profile_root(store, selected_profile), path)
    if not path.exists():
        raise InputValidationError("PRECOMPUTED_ENTRY_MISSING", f"No cache entry exists for PDB ID {normalized.upper()}.")
    raw_payload, expanded_bytes = _bounded_gzip_json(path, expanded_budget=expanded_budget)
    payload = validate_entry(
        raw_payload,
        expected_pdb_id=normalized,
        expected_profile_id=selected_profile,
        require_current_annotations=require_current_annotations,
    )
    return payload, expanded_bytes


def _entry_payload(
    source_path: Path,
    structure: Mapping[str, Any],
    interactions: Sequence[Mapping[str, Any]],
    *,
    selected_profile: str,
    source_fingerprint: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    from .pipeline import _compact_structure_summaries

    identity_raw = structure.get("structure_identity")
    if isinstance(identity_raw, Mapping) and identity_raw.get("canonical_id"):
        try:
            identity = StructureIdentity.from_mapping(identity_raw)
        except ValueError as exc:
            raise InputValidationError("INVALID_PDB_ID", "Structure identity is invalid.") from exc
    else:
        identity = identity_from_official_id(str(structure.get("pdb_id") or ""))
    if identity.source != "pdb":
        raise InputValidationError("INVALID_PDB_ID", "Precomputed archive entries require an official PDB identity.")
    normalized = identity.canonical_id
    label = f"pdb:{normalized}"
    compact = _compact_structure_summaries([dict(structure)])[0]
    compact.pop("identity_claims", None)
    compact.pop("input_warnings", None)
    compact["pdb_id"] = identity.display_id
    compact["file_path"] = label
    compact["structure_identity"] = identity.as_dict()
    geometry_chains: list[dict[str, Any]] = []
    annotation_chains: list[dict[str, Any]] = []
    for chain in compact.get("atom_data", []):
        chain_id = str(chain.get("chain_id") or "")
        structured_chain = ChainIdentity(
            structure_key=normalized,
            structure_display_id=identity.display_id,
            chain_id=chain_id,
            model_index=int(chain.get("model_index", 1)),
        )
        chain["unique_chain_id"] = structured_chain.node_id(include_model=False)
        chain["chain_identity"] = structured_chain.as_dict()
        chain["structure_key"] = normalized
        chain["structure_display_id"] = identity.display_id
        geometry_chains.append({
            key: value for key, value in chain.items()
            if key not in ANNOTATION_CHAIN_FIELDS
        })
        annotation = {"unique_chain_id": chain.get("unique_chain_id")}
        for field in ANNOTATION_CHAIN_FIELDS:
            if field in chain:
                annotation[field] = chain.get(field)
        annotation.setdefault("molecule_name", "Unknown")
        annotation.setdefault("molecule_type", "Unknown")
        annotation_chains.append(annotation)
    compact["atom_data"] = geometry_chains
    compact_edges: list[dict[str, Any]] = []
    for interaction in interactions:
        edge = {
            "file_path": label,
            "structure_key": normalized,
            "model_index": 1,
            "chain_a": interaction.get("chain_a"),
            "chain_b": interaction.get("chain_b"),
            "all_atoms_count": interaction.get("all_atoms_count"),
            "interaction_type": interaction.get("interaction_type"),
        }
        if "ca_neighbors" in interaction:
            edge["ca_neighbors"] = interaction.get("ca_neighbors")
        compact_edges.append(edge)
    return {
        "artifact_schema_version": ARTIFACT_SCHEMA_VERSION,
        "scientific_pipeline_version": SCIENTIFIC_PIPELINE_VERSION,
        "pdb2net_version": __version__,
        "profile_id": selected_profile,
        "pdb_id": normalized,
        "structure_identity": identity.as_dict(),
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source": dict(source_fingerprint or _source_fingerprint(source_path)),
        "geometry": {"structure": compact, "interactions": compact_edges},
        "annotations": {
            "profile_id": annotation_profile_id(),
            "references": {"manifest_id": _reference_manifest_id()},
            "chains": annotation_chains,
        },
        "counts": {"chains": len(compact.get("atom_data", [])), "interactions": len(compact_edges)},
    }


def write_entry(
    store: Path | str,
    source_path: Path,
    structure: Mapping[str, Any],
    interactions: Sequence[Mapping[str, Any]],
    *,
    source_fingerprint: Mapping[str, Any] | None = None,
    _selected_profile: str | None = None,
) -> tuple[Path, bool]:
    """Atomically publish one entry. Return ``(path, changed)``."""
    selected = _selected_profile
    if selected is None:
        selected, _ = ensure_profile_manifest(store)
    payload = _entry_payload(
        source_path,
        structure,
        interactions,
        selected_profile=selected,
        source_fingerprint=source_fingerprint,
    )
    normalized = payload["pdb_id"]
    validate_entry(payload, expected_pdb_id=normalized, expected_profile_id=selected)
    path = entry_path(store, normalized, selected)
    if path.exists():
        try:
            current = validate_entry(
                _bounded_gzip_json(path)[0], expected_pdb_id=normalized, expected_profile_id=selected
            )
            if (
                current.get("source", {}).get("sha256") == payload["source"]["sha256"]
                and current.get("annotations", {}).get("profile_id") == payload["annotations"]["profile_id"]
            ):
                return path, False
        except InputValidationError:
            pass
    serialized = _canonical_json(payload)
    if len(serialized) > MAX_DECOMPRESSED_ENTRY_BYTES:
        raise InputValidationError("CACHE_ENTRY_TOO_LARGE", f"Expanded cache entry for {normalized.upper()} is too large.")
    compressed = gzip.compress(serialized, compresslevel=6, mtime=0)
    if len(compressed) > MAX_COMPRESSED_ENTRY_BYTES:
        raise InputValidationError("CACHE_ENTRY_TOO_LARGE", f"Cache entry for {normalized.upper()} is too large.")
    _atomic_write(path, compressed, safe_root=profile_root(store, selected))
    return path, True


def discover_source_files(input_dir: Path | str, *, recursive: bool) -> list[Path]:
    raw_root = Path(input_dir).expanduser()
    if raw_root.is_symlink():
        raise InputValidationError(
            "SYMLINK_INPUT_NOT_ALLOWED", "Precompute input directory must not be a symlink."
        )
    root = raw_root.resolve()
    if not root.is_dir():
        raise InputValidationError("INPUT_PATH_NOT_DIRECTORY", f"Precompute input is not a directory: {root}")
    iterator = root.rglob("*") if recursive else root.iterdir()
    files: list[Path] = []
    for path in iterator:
        if path.is_symlink():
            if is_valid_file(str(path)):
                raise InputValidationError(
                    "SYMLINK_INPUT_NOT_ALLOWED",
                    f"Precompute structure inputs must not be symlinks: {path.name}",
                )
            continue
        if path.is_file() and is_valid_file(str(path)):
            files.append(path)
    files.sort()
    if not files:
        raise InputValidationError("NO_VALID_INPUT_FILES", f"No supported structure files were found in {root}.")
    return files


def resolve_archive_source(source_dir: Path | str, pdb_id: object) -> Path:
    """Resolve one PDB ID using fixed flat/wwPDB mirror conventions."""
    normalized = normalize_pdb_id(pdb_id)
    identity = identity_from_official_id(normalized)
    legacy = identity.legacy_id
    raw_root = Path(source_dir).expanduser()
    if raw_root.is_symlink():
        raise InputValidationError(
            "UNSAFE_PDB_ARCHIVE_PATH", "PDB archive root must not be a symlink."
        )
    root = raw_root.resolve()
    if not root.is_dir():
        raise InputValidationError("PDB_ARCHIVE_NOT_FOUND", f"PDB archive is not a directory: {root}")
    shard = pdb_archive_shard(normalized)
    basenames: list[str] = []
    for extension in ("cif", "mmcif"):
        basenames.extend((f"{normalized}.{extension}", f"{normalized}.{extension}.gz"))
    if legacy:
        for extension in ("cif", "mmcif", "pdb", "ent"):
            basenames.extend((f"{legacy}.{extension}", f"{legacy}.{extension}.gz"))
        basenames.extend((f"pdb{legacy}.ent", f"pdb{legacy}.ent.gz"))
    # Support a deliberately small flat archive as well as the fixed wwPDB
    # divided layouts.  This avoids an unbounded recursive scan on cache misses.
    archive_parents = (
        root,
        root / shard,
        root / "entries" / shard / normalized / "structures",
        root / "divided" / "mmCIF" / shard,
        root / "divided" / "pdb" / shard,
        root / "structures" / "divided" / "mmCIF" / shard,
        root / "structures" / "divided" / "pdb" / shard,
        root / "data" / "structures" / "divided" / "mmCIF" / shard,
        root / "data" / "structures" / "divided" / "pdb" / shard,
    )
    candidates: list[Path] = []
    for parent in archive_parents:
        for basename in basenames:
            candidate = parent / basename
            current = root
            has_symlink_component = False
            for component in candidate.relative_to(root).parts:
                current /= component
                if current.is_symlink():
                    has_symlink_component = True
                    break
            if has_symlink_component:
                if candidate.exists() or candidate.is_symlink():
                    raise InputValidationError(
                        "SYMLINK_INPUT_NOT_ALLOWED",
                        f"PDB archive source must not contain symlinks: {candidate.name}",
                    )
                continue
            if candidate.is_file():
                candidates.append(candidate.resolve())
    unique = sorted(set(candidates))
    if not unique:
        raise InputValidationError(
            "PDB_ARCHIVE_ENTRY_MISSING",
            f"PDB archive has no supported source file for {normalized.upper()}.",
        )
    if len(unique) > 1:
        raise InputValidationError(
            "PDB_ARCHIVE_ENTRY_AMBIGUOUS",
            f"PDB archive contains multiple source files for {normalized.upper()}.",
        )
    if root != unique[0] and root not in unique[0].parents:
        raise InputValidationError("UNSAFE_PDB_ARCHIVE_PATH", "Resolved archive file escaped its root.")
    return unique[0]


def _source_map(paths: Sequence[Path]) -> dict[str, Path]:
    mapped: dict[str, Path] = {}
    for path in paths:
        identity = get_structure_identity(str(path))
        if identity.source != "pdb":
            raise InputValidationError(
                "INVALID_PRECOMPUTE_SOURCE",
                f"Precomputed entries require an official PDB identity: {path.name}",
            )
        pdb_id = identity.canonical_id
        previous = mapped.get(pdb_id)
        if previous is not None and previous != path:
            raise InputValidationError(
                "DUPLICATE_STRUCTURE_IDENTITY",
                f"Multiple input files resolve to PDB ID {pdb_id.upper()}: {previous.name}, {path.name}",
            )
        mapped[pdb_id] = path
    return mapped


def _existing_entry_state(
    store: Path | str,
    pdb_id: str,
    selected_profile: str,
    source_fingerprint: Mapping[str, Any],
) -> tuple[str, dict[str, Any] | None]:
    try:
        payload, _ = _read_entry_for_profile(
            store,
            pdb_id,
            selected_profile,
            require_current_annotations=False,
        )
    except InputValidationError:
        return "miss", None
    if payload.get("source", {}).get("sha256") != source_fingerprint.get("sha256"):
        return "miss", None
    if payload.get("annotations", {}).get("profile_id") != annotation_profile_id():
        return "annotations_stale", payload
    return "hit", payload


def _existing_entry_matches(
    store: Path | str,
    pdb_id: str,
    selected_profile: str,
    source_fingerprint: Mapping[str, Any],
) -> bool:
    """Compatibility predicate used by older callers/tests."""
    return _existing_entry_state(store, pdb_id, selected_profile, source_fingerprint)[0] == "hit"


def precompute_sources(
    store: Path | str,
    paths: Sequence[Path],
    *,
    _locks_held: bool = False,
    _selected_profile: str | None = None,
) -> dict[str, Any]:
    """Precompute canonical graph entries for explicit source paths.

    Expensive annotation/search work is performed in the same bounded batches
    as the normal pipeline, preserving multi-FASTA BLAST/DIAMOND batching.
    """
    from . import pipeline
    from .unknown_molecule_uniprot import process_molecule_info

    _require_first_model()
    pipeline._validate_required_reference_files()
    if _selected_profile is None:
        selected, selected_root = ensure_profile_manifest(store)
    else:
        selected = _selected_profile
        selected_root = profile_root(store, selected)
    raw_paths = [Path(path).expanduser() for path in paths]
    for path in raw_paths:
        if path.is_symlink():
            raise InputValidationError("INVALID_PRECOMPUTE_SOURCE", f"Unsupported or unsafe source file: {path}")
    canonical_paths = [path.resolve() for path in raw_paths]
    for path in canonical_paths:
        if not path.is_file() or not is_valid_file(str(path)):
            raise InputValidationError("INVALID_PRECOMPUTE_SOURCE", f"Unsupported or unsafe source file: {path}")
    sources = _source_map(canonical_paths)

    if not _locks_held:
        aggregate: dict[str, Any] = {
            "artifact_schema_version": ARTIFACT_SCHEMA_VERSION,
            "profile_id": selected,
            "written": 0,
            "annotations_refreshed": 0,
            "cache_hits": 0,
            "failed": 0,
            "entries": [],
            "errors": [],
        }
        ordered = sorted(sources.items())
        for offset in range(0, len(ordered), MAX_SIMULTANEOUS_ENTRY_LOCKS):
            chunk = ordered[offset : offset + MAX_SIMULTANEOUS_ENTRY_LOCKS]
            with ExitStack() as locks:
                for pdb_id, _path in chunk:
                    locks.enter_context(
                        _entry_population_lock(
                            store,
                            pdb_id,
                            selected_profile=selected,
                            profile_directory=selected_root,
                        )
                    )
                # State and source fingerprints are recomputed only after all
                # locks for this deterministically ordered chunk are held.
                report = precompute_sources(
                    store,
                    [path for _pdb_id, path in chunk],
                    _locks_held=True,
                    _selected_profile=selected,
                )
            for field in ("written", "annotations_refreshed", "cache_hits", "failed"):
                aggregate[field] += int(report.get(field, 0))
            aggregate["entries"].extend(report.get("entries", []))
            aggregate["errors"].extend(report.get("errors", []))
        return aggregate

    initial_fingerprints = {
        pdb_id: _source_fingerprint(path)
        for pdb_id, path in sources.items()
    }

    report: dict[str, Any] = {
        "artifact_schema_version": ARTIFACT_SCHEMA_VERSION,
        "profile_id": selected,
        "written": 0,
        "annotations_refreshed": 0,
        "cache_hits": 0,
        "failed": 0,
        "entries": [],
        "errors": [],
    }
    pending: list[Path] = []
    stale_entries: dict[str, dict[str, Any]] = {}
    for pdb_id, path in sources.items():
        state, existing = _existing_entry_state(store, pdb_id, selected, initial_fingerprints[pdb_id])
        if state == "hit":
            report["cache_hits"] += 1
            report["entries"].append({"pdb_id": pdb_id, "status": "cache_hit"})
        else:
            pending.append(path)
            if state == "annotations_stale" and existing is not None:
                stale_entries[pdb_id] = existing
    if not pending:
        return report

    from .reference_data import load_pdb_fasta_headers

    pending_set = set(pending)
    pending_aliases: set[str] = set()
    for pdb_id, path in sources.items():
        if path not in pending_set:
            continue
        pending_aliases.add(pdb_id)
        legacy = identity_from_official_id(pdb_id).legacy_id
        if legacy:
            pending_aliases.add(legacy)
    pending_ids = tuple(sorted(pending_aliases))
    # The full archive build must scan pdb_seqres once, not once per processing
    # batch.  Only header metadata for requested IDs is retained; sequences are
    # still excluded from this index.
    pdb_fasta_headers = load_pdb_fasta_headers(
        str(config.get("pdb_fasta_path") or ""),
        pending_ids,
    )
    inventory = pipeline.inspect_input_files([str(path) for path in pending])
    for batch_paths_raw in pipeline.create_processing_batches([str(path) for path in pending], inventory):
        batch_paths = [Path(path).resolve() for path in batch_paths_raw]
        expected = {get_structure_identity(str(path)).canonical_id: path for path in batch_paths}
        try:
            batch_data = pipeline._parse_input_files([str(path) for path in batch_paths])
            process_molecule_info(batch_data, pdb_fasta_headers=pdb_fasta_headers)
            pipeline._run_blast_annotation(batch_data)
            geometry_batch: list[dict[str, Any]] = []
            for structure in batch_data:
                identity_raw = structure.get("structure_identity", {})
                canonical = (
                    StructureIdentity.from_mapping(identity_raw).canonical_id
                    if isinstance(identity_raw, Mapping) and identity_raw.get("canonical_id")
                    else normalize_pdb_id(structure.get("pdb_id"))
                )
                cached = stale_entries.get(canonical)
                if cached is None:
                    geometry_batch.append(structure)
                    continue
                cached_types = {
                    str(chain.get("unique_chain_id")): str(chain.get("molecule_type"))
                    for chain in cached.get("structure", {}).get("atom_data", [])
                }
                current_types = {
                    str(chain.get("unique_chain_id")): str(chain.get("molecule_type"))
                    for chain in structure.get("atom_data", [])
                }
                if cached_types != current_types:
                    geometry_batch.append(structure)
            batch_results = calculate_distances_with_ckdtree(geometry_batch) if geometry_batch else []
        except Exception as exc:
            for pdb_id in expected:
                report["failed"] += 1
                report["errors"].append({
                    "pdb_id": pdb_id,
                    "code": getattr(exc, "code", exc.__class__.__name__),
                    "message": str(exc)[:500],
                })
            tree_cache.clear()
            coords_cache.clear()
            continue

        produced: set[str] = set()
        for structure in batch_data:
            pdb_id: str | None = None
            try:
                identity_raw = structure.get("structure_identity", {})
                if isinstance(identity_raw, Mapping) and identity_raw.get("canonical_id"):
                    pdb_id = StructureIdentity.from_mapping(identity_raw).canonical_id
                else:
                    pdb_id = normalize_pdb_id(structure.get("pdb_id"))
                source = expected.get(pdb_id)
                if source is None:
                    raise InputValidationError(
                        "CACHE_PDB_ID_MISMATCH",
                        f"Parsed PDB ID {pdb_id.upper()} did not match the source inventory.",
                    )
                computed_edges = [
                    edge for edge in batch_results
                    if str(edge.get("structure_key") or "") == pdb_id
                ]
                cached = stale_entries.get(pdb_id)
                current_was_recomputed = any(
                    str(edge.get("structure_key") or "") == pdb_id
                    for edge in batch_results
                ) or any(
                    item is structure for item in geometry_batch
                )
                edges = computed_edges if current_was_recomputed or cached is None else list(cached["interactions"])
                original_fingerprint = initial_fingerprints[pdb_id]
                if _source_fingerprint(source) != original_fingerprint:
                    raise InputValidationError(
                        "INPUT_CHANGED_DURING_PROCESSING",
                        f"Source {source.name} changed while {pdb_id.upper()} was being precomputed.",
                    )
                path, changed = write_entry(
                    store,
                    source,
                    structure,
                    edges,
                    source_fingerprint=original_fingerprint,
                    _selected_profile=selected,
                )
                produced.add(pdb_id)
                report["written" if changed else "cache_hits"] += 1
                status = "written" if changed else "cache_hit"
                if cached is not None and changed and not current_was_recomputed:
                    report["annotations_refreshed"] += 1
                    status = "annotations_refreshed"
                report["entries"].append({
                    "pdb_id": pdb_id,
                    "status": status,
                    "artifact": str(path.relative_to(Path(store).expanduser().resolve())),
                })
            except Exception as exc:
                report["failed"] += 1
                report["errors"].append({
                    "pdb_id": pdb_id or str(structure.get("pdb_id") or "UNKNOWN"),
                    "code": getattr(exc, "code", exc.__class__.__name__),
                    "message": str(exc)[:500],
                })
        for missing in sorted(set(expected) - produced):
            if any(error.get("pdb_id") == missing for error in report["errors"]):
                continue
            report["failed"] += 1
            report["errors"].append({
                "pdb_id": missing,
                "code": "NO_PARSEABLE_STRUCTURES",
                "message": "The source file did not produce a parseable structure.",
            })
        tree_cache.clear()
        coords_cache.clear()
    return report


def precompute_directory(store: Path | str, input_dir: Path | str, *, recursive: bool) -> dict[str, Any]:
    return precompute_sources(store, discover_source_files(input_dir, recursive=recursive))


def _load_or_populate_entry(
    store: Path | str,
    pdb_id: str,
    *,
    source_dir: Path | str | None,
    populate_missing: bool,
) -> tuple[dict[str, Any], bool]:
    return _load_or_populate_entries(
        store,
        [pdb_id],
        source_dir=source_dir,
        populate_missing=populate_missing,
    )[0]


def _load_or_populate_entries(
    store: Path | str,
    pdb_ids: Sequence[str],
    *,
    source_dir: Path | str | None,
    populate_missing: bool,
    selected_profile: str | None = None,
    profile_directory: Path | None = None,
    resource_limits: Mapping[str, int | None] | None = None,
) -> list[tuple[dict[str, Any], bool]]:
    """Load entries in order and batch all cold misses under ordered locks."""
    if selected_profile is None or profile_directory is None:
        selected_profile, profile_directory = ensure_profile_manifest(store)
    hits: dict[str, bool] = {}
    missing: list[str] = []
    for pdb_id in pdb_ids:
        path = entry_path(store, pdb_id, selected_profile)
        if path.exists() or path.is_symlink():
            hits[pdb_id] = True
        else:
            if not populate_missing:
                raise InputValidationError(
                    "PRECOMPUTED_ENTRY_MISSING",
                    f"No cache entry exists for PDB ID {pdb_id.upper()}.",
                )
            missing.append(pdb_id)

    if missing:
        if source_dir is None:
            raise InputValidationError(
                "PDB_ARCHIVE_REQUIRED",
                "--populate-missing requires --source-dir.",
            )
        for offset in range(0, len(missing), MAX_SIMULTANEOUS_ENTRY_LOCKS):
            chunk = missing[offset : offset + MAX_SIMULTANEOUS_ENTRY_LOCKS]
            with ExitStack() as locks:
                # Stable ordering prevents deadlocks for overlapping jobs.
                for pdb_id in sorted(chunk):
                    locks.enter_context(
                        _entry_population_lock(
                            store,
                            pdb_id,
                            selected_profile=selected_profile,
                            profile_directory=profile_directory,
                        )
                    )

                sources: list[Path] = []
                populate_ids: list[str] = []
                for pdb_id in chunk:
                    path = entry_path(store, pdb_id, selected_profile)
                    if path.exists() or path.is_symlink():
                        hits[pdb_id] = True
                    else:
                        sources.append(resolve_archive_source(source_dir, pdb_id))
                        populate_ids.append(pdb_id)

                if sources:
                    report = precompute_sources(
                        store,
                        sources,
                        _locks_held=True,
                        _selected_profile=selected_profile,
                    )
                    if report["failed"]:
                        first = report["errors"][0]
                        raise InputValidationError(
                            str(first.get("code") or "PRECOMPUTE_FAILED"),
                            str(first.get("message") or "Could not populate requested cache entries."),
                        )
                    for pdb_id in populate_ids:
                        hits[pdb_id] = False

    limits = dict(resource_limits or {})
    compressed_total = 0
    for pdb_id in pdb_ids:
        path = entry_path(store, pdb_id, selected_profile)
        if path.is_symlink() or not path.is_file():
            raise InputValidationError("UNSAFE_CACHE_ENTRY", f"Cache entry is missing or a symlink: {path.name}")
        entry_bytes = path.stat().st_size
        maximum_single = limits.get("max_single_input_bytes")
        if maximum_single is not None and entry_bytes > maximum_single:
            raise InputValidationError(
                "INPUT_FILE_BYTES_LIMIT_EXCEEDED",
                f"Cache entry {pdb_id.upper()} is {entry_bytes} bytes; configured maximum is {maximum_single}.",
            )
        compressed_total += entry_bytes
        maximum_total = limits.get("max_total_input_bytes")
        if maximum_total is not None and compressed_total > maximum_total:
            raise InputValidationError(
                "INPUT_TOTAL_BYTES_LIMIT_EXCEEDED",
                f"Requested cache entries total {compressed_total} bytes; configured maximum is {maximum_total}.",
            )

    loaded: list[tuple[dict[str, Any], bool]] = []
    expanded_total = 0
    chain_total = 0
    edge_total = 0
    for pdb_id in pdb_ids:
        payload, expanded_bytes = _read_entry_for_profile(
            store,
            pdb_id,
            selected_profile,
            expanded_budget=MAX_REQUEST_EXPANDED_BYTES - expanded_total,
        )
        expanded_total += expanded_bytes
        chain_total += int(payload["counts"]["chains"])
        edge_total += int(payload["counts"]["interactions"])
        if chain_total > MAX_REQUEST_CHAINS or edge_total > MAX_REQUEST_EDGES:
            raise InputValidationError(
                "PRECOMPUTED_REQUEST_LIMIT_EXCEEDED",
                "Requested cache entries exceed aggregate chain or interaction limits.",
            )
        loaded.append((payload, hits[pdb_id]))
    return loaded


def _cache_resource_summary(
    requested_ids: Sequence[str],
    *,
    selected_profile: str,
    cache_hits: int,
    populated: int,
    compressed_bytes: int,
    expanded_bytes: int,
    largest_entry_bytes: int,
) -> dict[str, Any]:
    return {
        "input": {
            "files": len(requested_ids),
            "total_bytes": compressed_bytes,
            "largest_file_bytes": largest_entry_bytes,
        },
        "processing": {
            "mode": "precomputed",
            "cache_hits": cache_hits,
            "populated_missing": populated,
            "parsing_workers": 0 if populated == 0 else None,
        },
        "precomputed_store": {
            "artifact_schema_version": ARTIFACT_SCHEMA_VERSION,
            "profile_id": selected_profile,
            "entries": len(requested_ids),
            "compressed_bytes": compressed_bytes,
            "expanded_bytes": expanded_bytes,
        },
        "peak_rss_bytes": {"main_process": None, "child_processes": None},
    }


def run_assemble_pipeline(
    store: Path | str,
    pdb_ids: Sequence[object],
    *,
    web_output_dir: str | None = None,
    source_dir: Path | str | None = None,
    populate_missing: bool = False,
):
    """Assemble normal PDB2Net outputs from validated per-PDB cache entries."""
    from . import pipeline
    from .outputs import (
        collect_generated_outputs,
        collect_web_outputs,
        create_run_output_paths,
        write_failed_run_manifest,
        write_run_manifest,
        write_run_summary,
        write_runtime_analysis,
    )

    started_at = datetime.now().isoformat(timespec="seconds")
    start = time.monotonic()
    labels = []
    try:
        _require_first_model()
        network_annotation_config()
        if bool(config.get("export_detailed_interactions", False)):
            raise InputValidationError(
                "DETAILED_INTERACTIONS_UNAVAILABLE",
                "Compact precomputed entries do not contain atom pairs; use the normal file mode for detailed CSVs.",
            )
        pipeline._validate_output_root()
        pipeline._combined_graph_limits()
        requested = normalize_requested_ids(pdb_ids)
        resource_limits = pipeline._resource_limits()
        maximum = resource_limits.get("max_input_files")
        if maximum is not None and len(requested) > maximum:
            raise InputValidationError(
                "INPUT_FILE_COUNT_LIMIT_EXCEEDED",
                f"Input contains {len(requested)} PDB IDs; configured maximum is {maximum}.",
            )
        selected, selected_root = ensure_profile_manifest(store)
        output_paths = create_run_output_paths(config["output_path"])
        labels = [f"pdb:{pdb_id}" for pdb_id in requested]
        structures: list[dict[str, Any]] = []
        interactions: list[dict[str, Any]] = []
        cache_hits = 0
        populated = 0
        compressed_bytes = 0
        expanded_bytes = 0
        largest_entry_bytes = 0
        total_chains = 0
        total_edges = 0
        loaded_entries = _load_or_populate_entries(
            store,
            requested,
            source_dir=source_dir,
            populate_missing=populate_missing,
            selected_profile=selected,
            profile_directory=selected_root,
            resource_limits=resource_limits,
        )
        for pdb_id, (payload, hit) in zip(requested, loaded_entries):
            cache_hits += int(hit)
            populated += int(not hit)
            cache_path = entry_path(store, pdb_id, selected)
            entry_bytes = cache_path.stat().st_size
            maximum_single = resource_limits.get("max_single_input_bytes")
            if maximum_single is not None and entry_bytes > maximum_single:
                raise InputValidationError(
                    "INPUT_FILE_BYTES_LIMIT_EXCEEDED",
                    f"Cache entry {pdb_id.upper()} is {entry_bytes} bytes; configured maximum is {maximum_single}.",
                )
            compressed_bytes += entry_bytes
            largest_entry_bytes = max(largest_entry_bytes, entry_bytes)
            maximum_total = resource_limits.get("max_total_input_bytes")
            if maximum_total is not None and compressed_bytes > maximum_total:
                raise InputValidationError(
                    "INPUT_TOTAL_BYTES_LIMIT_EXCEEDED",
                    f"Requested cache entries total {compressed_bytes} bytes; configured maximum is {maximum_total}.",
                )
            expanded_bytes += len(_canonical_json(payload))
            total_chains += int(payload["counts"]["chains"])
            total_edges += int(payload["counts"]["interactions"])
            if (
                expanded_bytes > MAX_REQUEST_EXPANDED_BYTES
                or total_chains > MAX_REQUEST_CHAINS
                or total_edges > MAX_REQUEST_EDGES
            ):
                raise InputValidationError(
                    "PRECOMPUTED_REQUEST_LIMIT_EXCEEDED",
                    "Requested cache entries exceed aggregate expanded-size, chain, or interaction limits.",
                )
            structures.append(dict(payload["structure"]))
            interactions.extend(dict(edge) for edge in payload["interactions"])

        timings = pipeline.PipelineTimings()
        skipped = pipeline._export_network_outputs(
            structures,
            interactions,
            config["networks"],
            output_paths,
            timings,
        )
        warnings = [
            {
                "code": "COMBINED_GRAPH_LIMIT_EXCEEDED",
                "message": (
                    f"Skipped {entry['name']} with {entry['counts']['nodes']} nodes and "
                    f"{entry['counts']['edges']} edges because configured combined graph limits were exceeded."
                ),
                "output": entry["name"],
            }
            for entry in skipped
        ] + pipeline._embedded_annotation_warnings(structures)
        total = time.monotonic() - start
        write_runtime_analysis(output_paths.log_file, timings.as_dict(), total)
        generated = collect_generated_outputs(output_paths)
        write_run_manifest(
            output_paths.manifest_file,
            input_files=labels,
            output_paths=output_paths,
            config_snapshot=pipeline._config_snapshot(config["networks"]),
            status="success",
            started_at=started_at,
            finished_at=datetime.now().isoformat(timespec="seconds"),
            total_time=total,
            input_path=f"precomputed:{selected}",
            generated_outputs=generated,
            extra_counts={
                "structures": len(structures),
                "chains": sum(len(item.get("atom_data", [])) for item in structures),
                "interactions": len(interactions),
                "skipped_outputs": len(skipped),
            },
            annotations=pipeline._annotation_summary(structures),
            references={
                "manifest_id": _reference_manifest_id(),
                "precomputed_store": {
                    "artifact_schema_version": ARTIFACT_SCHEMA_VERSION,
                    "profile_id": selected,
                    "annotation_profile_id": annotation_profile_id(),
                    "source_scope": SOURCE_SCOPE,
                },
            },
            resources=_cache_resource_summary(
                requested,
                selected_profile=selected,
                cache_hits=cache_hits,
                populated=populated,
                compressed_bytes=compressed_bytes,
                expanded_bytes=expanded_bytes,
                largest_entry_bytes=largest_entry_bytes,
            ),
            skipped_outputs=skipped,
            warnings=warnings,
            identities=pipeline._identity_summary(structures),
            structure_inputs=pipeline._structure_input_summary(structures),
        )
        write_run_summary(output_paths)
        if web_output_dir:
            collect_web_outputs(output_paths, web_output_dir)
        return output_paths
    except Exception as exc:
        output_paths = write_failed_run_manifest(
            config.get("output_path", ""),
            input_path="precomputed:" + ",".join(labels or [str(value) for value in pdb_ids]),
            config_snapshot=pipeline._config_snapshot(config.get("networks", {})),
            error=exc,
            started_at=started_at,
        )
        if web_output_dir:
            collect_web_outputs(output_paths, web_output_dir)
        raise
    finally:
        tree_cache.clear()
        coords_cache.clear()
