"""Versioned, portable per-PDB graph cache and assembly pipeline.

The store is intentionally independent of the webserver.  It persists the
small scientific primitives needed to recreate per-PDB and combined networks:
annotated chain metadata and filtered chain-pair interaction edges.  Raw atom
coordinates and detailed atom-pair tables are deliberately not cached.
"""

from __future__ import annotations

import gzip
import hashlib
import json
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
from .file_parser import get_pdb_id, is_valid_file
from .input_contract import InputValidationError


ARTIFACT_SCHEMA_VERSION = "1"
SCIENTIFIC_PIPELINE_VERSION = "pdb2net-asu-first-standard-interactions-v1"
SOURCE_SCOPE = "asymmetric_unit"
MAX_COMPRESSED_ENTRY_BYTES = 64 * 1024 * 1024
MAX_DECOMPRESSED_ENTRY_BYTES = 128 * 1024 * 1024
MAX_CHAINS_PER_ENTRY = 50_000
MAX_EDGES_PER_ENTRY = 500_000
MAX_REQUEST_EXPANDED_BYTES = 256 * 1024 * 1024
MAX_REQUEST_CHAINS = 200_000
MAX_REQUEST_EDGES = 2_000_000
POPULATION_LOCK_TIMEOUT_SECONDS = 300.0

PDB_ID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
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
}


def normalize_pdb_id(value: object) -> str:
    """Return a canonical lowercase four-character PDB identifier."""
    raw = str(value or "").strip()
    if not PDB_ID_RE.fullmatch(raw):
        raise InputValidationError(
            "INVALID_PDB_ID",
            f"PDB IDs must contain one leading digit and three alphanumeric characters: {raw!r}",
        )
    return raw.lower()


def normalize_requested_ids(values: Iterable[object]) -> list[str]:
    """Normalize and stably de-duplicate requested PDB IDs."""
    normalized: list[str] = []
    seen: set[str] = set()
    for value in values:
        pdb_id = normalize_pdb_id(value)
        if pdb_id not in seen:
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
            "Precomputed store schema 1 supports structure_model_policy='first' only.",
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
    """Return the cache-relevant scientific profile."""
    import gemmi

    from .uniprot_matcher import _search_policy_payload

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
            "semantics": "gemmi-heavy-atoms-no-hydrogen-or-deuterium-v1",
            "gemmi_version": str(getattr(gemmi, "__version__", "unknown")),
        },
        "distance_thresholds": {
            "ca_radius": float(thresholds.get("ca_radius", 12.0)),
            "all_atoms_radius": float(thresholds.get("all_atoms_radius", 5.0)),
        },
        "interaction_filters": _profile_interaction_filters(),
        "annotation": {
            "reference_manifest_id": _reference_manifest_id(),
            # Reuse the sequence-search cache's result policy.  Chunk sizes,
            # thread counts, temporary paths, and executable locations are
            # intentionally excluded; DIAMOND block size is included because
            # it can affect sensitivity.
            "search_policy": _search_policy_payload(),
        },
    }
    # Search thresholds are naturally represented as tuples in Python.  Store
    # one canonical JSON-native shape so an in-memory profile compares exactly
    # with the same profile read back from manifest.json.
    return json.loads(_canonical_json(profile))


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
    return profile_root(store, profile_hash) / "entries" / normalized[1:3] / f"{normalized}.json.gz"


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


def validate_entry(payload: object, *, expected_pdb_id: str, expected_profile_id: str) -> dict[str, Any]:
    """Validate an untrusted cache entry before network generation."""
    if not isinstance(payload, dict):
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache entry must be a JSON object.")
    if payload.get("artifact_schema_version") != ARTIFACT_SCHEMA_VERSION:
        raise InputValidationError("CACHE_SCHEMA_MISMATCH", "Cache entry schema is unsupported.")
    if payload.get("scientific_pipeline_version") != SCIENTIFIC_PIPELINE_VERSION:
        raise InputValidationError("CACHE_PIPELINE_MISMATCH", "Cache entry scientific pipeline is stale.")
    if payload.get("profile_id") != expected_profile_id:
        raise InputValidationError("CACHE_PROFILE_MISMATCH", "Cache entry profile does not match the active profile.")
    if payload.get("pdb_id") != expected_pdb_id:
        raise InputValidationError("CACHE_PDB_ID_MISMATCH", "Cache entry PDB ID does not match its requested path.")

    references = payload.get("references")
    if not isinstance(references, dict) or references.get("manifest_id") != _reference_manifest_id():
        raise InputValidationError(
            "CACHE_REFERENCE_MISMATCH",
            "Cache entry annotations do not match the active reference manifest.",
        )

    source = payload.get("source")
    if not isinstance(source, dict) or not re.fullmatch(r"[a-f0-9]{64}", str(source.get("sha256") or "")):
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache entry source fingerprint is invalid.")
    if source.get("scope") != SOURCE_SCOPE or not _valid_nonnegative_int(source.get("size_bytes")):
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache entry source metadata is invalid.")
    source_name = source.get("name")
    if (
        not isinstance(source_name, str)
        or not source_name
        or len(source_name) > 512
        or Path(source_name).name != source_name
    ):
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache source name is not portable.")

    structure = payload.get("structure")
    if not isinstance(structure, dict) or structure.get("pdb_id") != expected_pdb_id.upper():
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache structure metadata is invalid.")
    if structure.get("file_path") != f"pdb:{expected_pdb_id.upper()}":
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache source label is not portable.")
    if set(structure) != {"file_path", "pdb_id", "atom_data"}:
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache structure contains unsupported fields.")
    chains = structure.get("atom_data")
    if not isinstance(chains, list) or not chains or len(chains) > MAX_CHAINS_PER_ENTRY:
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache chain collection is invalid.")
    chain_ids: set[str] = set()
    for chain in chains:
        if not isinstance(chain, dict) or not set(chain).issubset(CHAIN_FIELDS):
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache chains must be compact metadata objects.")
        required_chain_fields = {
            "chain_id",
            "unique_chain_id",
            "model_index",
            "molecule_name",
            "molecule_type",
            "sequence",
            "is_hetatm",
            "aa_len",
            "nt_len",
        }
        if not required_chain_fields.issubset(chain):
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache chain metadata is incomplete.")
        chain_id = chain.get("chain_id")
        unique_id = chain.get("unique_chain_id")
        if (
            not isinstance(chain_id, str)
            or len(chain_id) > 256
            or unique_id != f"{expected_pdb_id.upper()}:{chain_id}"
        ):
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache chain ID is invalid.")
        if unique_id in chain_ids:
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache chain IDs must be unique.")
        if chain.get("model_index") != 1:
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache schema 1 accepts first-model chains only.")
        if chain.get("molecule_type") not in ALLOWED_MOLECULE_TYPES:
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache molecule type is invalid.")
        if not _bounded_optional_text(chain.get("molecule_name"), maximum=10_000):
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache molecule name is invalid.")
        sequence = chain.get("sequence")
        if not isinstance(sequence, str) or len(sequence) > 1_000_000:
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache chain sequence is invalid.")
        if not isinstance(chain.get("is_hetatm"), bool):
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache chain hetero-atom flag is invalid.")
        if not _valid_nonnegative_int(chain.get("aa_len")) or not _valid_nonnegative_int(chain.get("nt_len")):
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache chain lengths are invalid.")
        uniprot_id = chain.get("uniprot_id")
        if uniprot_id not in (None, ""):
            from .uniprot_matcher import _is_uniprotkb_accession

            if not isinstance(uniprot_id, str) or not _is_uniprotkb_accession(uniprot_id):
                raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache UniProt accession is invalid.")
        for field in (
            "annotation_source",
            "matched_database",
            "matched_id",
            "representative_accession",
            "annotation_confidence",
        ):
            if not _bounded_optional_text(chain.get(field), maximum=1_000):
                raise InputValidationError("CORRUPT_CACHE_ENTRY", f"Cache chain field {field} is invalid.")
        chain_ids.add(unique_id)

    interactions = payload.get("interactions")
    if not isinstance(interactions, list) or len(interactions) > MAX_EDGES_PER_ENTRY:
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache interaction collection is invalid.")
    seen_edges: set[tuple[str, str]] = set()
    chains_by_id = {chain["unique_chain_id"]: chain for chain in chains}
    profile_filters = _profile_interaction_filters()
    for edge in interactions:
        if not isinstance(edge, dict) or not set(edge).issubset(
            {"file_path", "chain_a", "chain_b", "all_atoms_count", "ca_neighbors", "interaction_type"}
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
        if edge.get("file_path") != f"pdb:{expected_pdb_id.upper()}":
            raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache interaction source label is invalid.")

    counts = payload.get("counts")
    if not isinstance(counts, dict) or counts.get("chains") != len(chains) or counts.get("interactions") != len(interactions):
        raise InputValidationError("CORRUPT_CACHE_ENTRY", "Cache entry counts do not match its content.")
    return payload


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
    lock_path = root / "locks" / normalized[1:3] / f"{normalized}.lock"
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

    normalized = normalize_pdb_id(structure.get("pdb_id"))
    label = f"pdb:{normalized.upper()}"
    compact = _compact_structure_summaries([dict(structure)])[0]
    compact["pdb_id"] = normalized.upper()
    compact["file_path"] = label
    compact_edges: list[dict[str, Any]] = []
    for interaction in interactions:
        edge = {
            "file_path": label,
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
        "created_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source": dict(source_fingerprint or _source_fingerprint(source_path)),
        "references": {"manifest_id": _reference_manifest_id()},
        "structure": compact,
        "interactions": compact_edges,
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
            if current.get("source", {}).get("sha256") == payload["source"]["sha256"]:
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
    root = Path(input_dir).expanduser().resolve()
    if not root.is_dir():
        raise InputValidationError("INPUT_PATH_NOT_DIRECTORY", f"Precompute input is not a directory: {root}")
    iterator = root.rglob("*") if recursive else root.iterdir()
    files = sorted(path for path in iterator if path.is_file() and not path.is_symlink() and is_valid_file(str(path)))
    if not files:
        raise InputValidationError("NO_VALID_INPUT_FILES", f"No supported structure files were found in {root}.")
    return files


def resolve_archive_source(source_dir: Path | str, pdb_id: object) -> Path:
    """Resolve one PDB ID using fixed flat/wwPDB mirror conventions."""
    normalized = normalize_pdb_id(pdb_id)
    root = Path(source_dir).expanduser().resolve()
    if not root.is_dir():
        raise InputValidationError("PDB_ARCHIVE_NOT_FOUND", f"PDB archive is not a directory: {root}")
    shard = normalized[1:3]
    basenames: list[str] = []
    for extension in ("cif", "mmcif", "pdb", "ent"):
        basenames.extend((f"{normalized}.{extension}", f"{normalized}.{extension}.gz"))
    basenames.extend((f"pdb{normalized}.ent", f"pdb{normalized}.ent.gz"))
    basenames.extend((f"pdb_0000{normalized}.cif", f"pdb_0000{normalized}.cif.gz"))
    # Support a deliberately small flat archive as well as the fixed wwPDB
    # divided layouts.  This avoids an unbounded recursive scan on cache misses.
    archive_parents = (
        root,
        root / shard,
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
            if candidate.is_file() and not candidate.is_symlink():
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
        pdb_id = normalize_pdb_id(get_pdb_id(str(path)))
        previous = mapped.get(pdb_id)
        if previous is not None and previous != path:
            raise InputValidationError(
                "DUPLICATE_PDB_SOURCE",
                f"Multiple input files resolve to PDB ID {pdb_id.upper()}: {previous.name}, {path.name}",
            )
        mapped[pdb_id] = path
    return mapped


def _existing_entry_matches(
    store: Path | str,
    pdb_id: str,
    selected_profile: str,
    source_fingerprint: Mapping[str, Any],
) -> bool:
    try:
        payload = _load_entry_for_profile(store, pdb_id, selected_profile)
    except InputValidationError:
        return False
    return payload.get("source", {}).get("sha256") == source_fingerprint.get("sha256")


def precompute_sources(store: Path | str, paths: Sequence[Path]) -> dict[str, Any]:
    """Precompute canonical graph entries for explicit source paths.

    Expensive annotation/search work is performed in the same bounded batches
    as the normal pipeline, preserving multi-FASTA BLAST/DIAMOND batching.
    """
    from . import pipeline
    from .unknown_molecule_uniprot import process_molecule_info

    _require_first_model()
    pipeline._validate_required_reference_files()
    selected, _ = ensure_profile_manifest(store)
    raw_paths = [Path(path).expanduser() for path in paths]
    for path in raw_paths:
        if path.is_symlink():
            raise InputValidationError("INVALID_PRECOMPUTE_SOURCE", f"Unsupported or unsafe source file: {path}")
    canonical_paths = [path.resolve() for path in raw_paths]
    for path in canonical_paths:
        if not path.is_file() or not is_valid_file(str(path)):
            raise InputValidationError("INVALID_PRECOMPUTE_SOURCE", f"Unsupported or unsafe source file: {path}")
    sources = _source_map(canonical_paths)
    initial_fingerprints = {
        pdb_id: _source_fingerprint(path)
        for pdb_id, path in sources.items()
    }

    report: dict[str, Any] = {
        "artifact_schema_version": ARTIFACT_SCHEMA_VERSION,
        "profile_id": selected,
        "written": 0,
        "cache_hits": 0,
        "failed": 0,
        "entries": [],
        "errors": [],
    }
    pending: list[Path] = []
    for pdb_id, path in sources.items():
        if _existing_entry_matches(store, pdb_id, selected, initial_fingerprints[pdb_id]):
            report["cache_hits"] += 1
            report["entries"].append({"pdb_id": pdb_id.upper(), "status": "cache_hit"})
        else:
            pending.append(path)
    if not pending:
        return report

    from .reference_data import load_pdb_fasta_headers

    pending_set = set(pending)
    pending_ids = tuple(sorted(pdb_id for pdb_id, path in sources.items() if path in pending_set))
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
        expected = {normalize_pdb_id(get_pdb_id(str(path))): path for path in batch_paths}
        try:
            batch_data = pipeline._parse_input_files([str(path) for path in batch_paths])
            process_molecule_info(batch_data, pdb_fasta_headers=pdb_fasta_headers)
            pipeline._run_blast_annotation(batch_data)
            batch_results = calculate_distances_with_ckdtree(batch_data)
        except Exception as exc:
            for pdb_id in expected:
                report["failed"] += 1
                report["errors"].append({
                    "pdb_id": pdb_id.upper(),
                    "code": getattr(exc, "code", exc.__class__.__name__),
                    "message": str(exc)[:500],
                })
            tree_cache.clear()
            coords_cache.clear()
            continue

        produced: set[str] = set()
        for structure in batch_data:
            try:
                pdb_id = normalize_pdb_id(structure.get("pdb_id"))
                source = expected.get(pdb_id)
                if source is None:
                    raise InputValidationError(
                        "CACHE_PDB_ID_MISMATCH",
                        f"Parsed PDB ID {pdb_id.upper()} did not match the source inventory.",
                    )
                edges = [
                    edge
                    for edge in batch_results
                    if str(edge.get("chain_a") or "").startswith(pdb_id.upper() + ":")
                ]
                original_fingerprint = initial_fingerprints[pdb_id]
                if _source_fingerprint(source) != original_fingerprint:
                    raise InputValidationError(
                        "PRECOMPUTE_SOURCE_CHANGED",
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
                report["entries"].append({
                    "pdb_id": pdb_id.upper(),
                    "status": "written" if changed else "cache_hit",
                    "artifact": str(path.relative_to(Path(store).expanduser().resolve())),
                })
            except Exception as exc:
                report["failed"] += 1
                report["errors"].append({
                    "pdb_id": str(structure.get("pdb_id") or "UNKNOWN"),
                    "code": getattr(exc, "code", exc.__class__.__name__),
                    "message": str(exc)[:500],
                })
        for missing in sorted(set(expected) - produced):
            if any(error.get("pdb_id") == missing.upper() for error in report["errors"]):
                continue
            report["failed"] += 1
            report["errors"].append({
                "pdb_id": missing.upper(),
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
        with ExitStack() as locks:
            # Stable ordering prevents deadlocks for overlapping multi-ID jobs.
            for pdb_id in sorted(missing):
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
            for pdb_id in missing:
                path = entry_path(store, pdb_id, selected_profile)
                if path.exists() or path.is_symlink():
                    hits[pdb_id] = True
                else:
                    sources.append(resolve_archive_source(source_dir, pdb_id))
                    populate_ids.append(pdb_id)

            if sources:
                report = precompute_sources(store, sources)
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
        labels = [f"pdb:{pdb_id.upper()}" for pdb_id in requested]
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
        ]
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
