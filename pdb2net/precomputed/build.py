"""The sole offline writer for schema-3 precomputed stores."""

from __future__ import annotations

from pathlib import Path
import stat
from typing import Any, Mapping, Sequence

from ..config_loader import config
from ..distances import calculate_distances_with_ckdtree, coords_cache, tree_cache
from ..file_parser import get_structure_identity, is_valid_file
from ..input_contract import InputValidationError
from ..structure_identity import ChainIdentity, StructureIdentity, identity_from_official_id
from .io import (
    entry_path,
    manifest_path,
    read_entry_json,
    source_fingerprint,
    store_root,
    write_entry_json,
    write_manifest_json,
)
from .schema import (
    ANNOTATION_CHAIN_FIELDS,
    ARTIFACT_SCHEMA_VERSION,
    GEOMETRY_CHAIN_FIELDS,
    manifest_document,
    normalize_pdb_id,
    producer,
    profile_id,
    scientific_profile,
    utc_timestamp,
    validate_entry,
)


def discover_source_files(input_dir: Path | str, *, recursive: bool) -> list[Path]:
    """Discover supported regular structure files beneath one offline input root."""
    raw_root = Path(input_dir).expanduser()
    if raw_root.is_symlink():
        raise InputValidationError(
            "SYMLINK_INPUT_NOT_ALLOWED", "Precompute input directory must not be a symlink."
        )
    root = raw_root.resolve()
    if not root.is_dir():
        raise InputValidationError(
            "INPUT_PATH_NOT_DIRECTORY", f"Precompute input is not a directory: {root}"
        )
    iterator = root.rglob("*") if recursive else root.iterdir()
    files: list[Path] = []
    for path in iterator:
        if path.is_symlink():
            if is_valid_file(str(path)):
                raise InputValidationError(
                    "SYMLINK_INPUT_NOT_ALLOWED",
                    f"Precompute structure input must not be a symlink: {path.name}",
                )
            continue
        if path.is_file() and is_valid_file(str(path)):
            files.append(path.resolve())
    files.sort()
    if not files:
        raise InputValidationError(
            "NO_VALID_INPUT_FILES", f"No supported structure files were found in {root}."
        )
    return files


def _source_map(paths: Sequence[Path]) -> dict[str, Path]:
    mapped: dict[str, Path] = {}
    for path in paths:
        identity = get_structure_identity(str(path))
        if identity.source != "pdb":
            raise InputValidationError(
                "INVALID_PRECOMPUTE_SOURCE",
                f"Precomputed entries require an official PDB identity: {path.name}",
            )
        previous = mapped.get(identity.canonical_id)
        if previous is not None and previous != path:
            raise InputValidationError(
                "DUPLICATE_STRUCTURE_IDENTITY",
                "Multiple input files resolve to PDB ID "
                f"{identity.canonical_id.upper()}: {previous.name}, {path.name}",
            )
        mapped[identity.canonical_id] = path
    return mapped


def _entry_payload(
    source: Path,
    structure: Mapping[str, Any],
    interactions: Sequence[Mapping[str, Any]],
    *,
    profile: Mapping[str, Any],
    source_metadata: Mapping[str, Any],
) -> dict[str, Any]:
    from ..pipeline import _compact_structure_summaries

    identity_raw = structure.get("structure_identity")
    if isinstance(identity_raw, Mapping) and identity_raw.get("canonical_id"):
        try:
            identity = StructureIdentity.from_mapping(identity_raw)
        except ValueError as exc:
            raise InputValidationError(
                "INVALID_PDB_ID", "Precompute structure identity is invalid."
            ) from exc
    else:
        identity = identity_from_official_id(str(structure.get("pdb_id") or ""))
    if identity.source != "pdb":
        raise InputValidationError(
            "INVALID_PDB_ID", "Precomputed entries require an official PDB identity."
        )
    normalized = identity.canonical_id
    label = f"pdb:{normalized}"
    compact = _compact_structure_summaries([dict(structure)])[0]
    compact["pdb_id"] = identity.display_id
    compact["file_path"] = label
    compact["structure_identity"] = identity.as_dict()

    geometry_chains: list[dict[str, Any]] = []
    annotation_chains: list[dict[str, Any]] = []
    for chain_raw in compact.get("atom_data", []):
        chain = dict(chain_raw)
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
        geometry_chains.append(
            {key: chain[key] for key in GEOMETRY_CHAIN_FIELDS if key in chain}
        )
        annotation = {"unique_chain_id": chain["unique_chain_id"]}
        annotation.update(
            {key: chain[key] for key in ANNOTATION_CHAIN_FIELDS if key in chain}
        )
        annotation.setdefault("molecule_name", "Unknown")
        annotation.setdefault("molecule_type", "Unknown")
        annotation_chains.append(annotation)

    allowed_structure_fields = {
        "file_path",
        "pdb_id",
        "atom_data",
        "structure_identity",
        "input_format",
        "input_kind",
        "embedded_annotation_counts",
    }
    compact_structure = {
        key: compact[key]
        for key in allowed_structure_fields
        if key in compact
    }
    compact_structure["atom_data"] = geometry_chains
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
        "created_at": utc_timestamp(),
        "producer": producer(),
        "profile_id": profile_id(profile),
        "pdb_id": normalized,
        "structure_identity": identity.as_dict(),
        "source": dict(source_metadata),
        "geometry": {
            "structure": compact_structure,
            "interactions": compact_edges,
        },
        "annotations": {
            "references": {
                "manifest_id": profile["annotations"]["reference_manifest_id"]
            },
            "chains": annotation_chains,
        },
        "counts": {
            "chains": len(geometry_chains),
            "interactions": len(compact_edges),
        },
    }


def _existing_entry_matches(
    store: Path | str,
    pdb_id: str,
    *,
    profile: Mapping[str, Any],
    source_metadata: Mapping[str, Any],
) -> bool:
    path = entry_path(store, pdb_id)
    if not path.exists() and not path.is_symlink():
        return False
    try:
        raw, _expanded = read_entry_json(store, pdb_id)
        payload = validate_entry(
            raw,
            expected_pdb_id=pdb_id,
            expected_profile_id=profile_id(profile),
            expected_profile=profile,
        )
    except InputValidationError as exc:
        if exc.code.startswith("UNSAFE_"):
            raise
        return False
    return payload.get("source") == dict(source_metadata)


def _write_structure_entry(
    store: Path | str,
    source: Path,
    structure: Mapping[str, Any],
    interactions: Sequence[Mapping[str, Any]],
    *,
    profile: Mapping[str, Any],
    source_metadata: Mapping[str, Any],
) -> Path:
    payload = _entry_payload(
        source,
        structure,
        interactions,
        profile=profile,
        source_metadata=source_metadata,
    )
    normalized = payload["pdb_id"]
    validate_entry(
        payload,
        expected_pdb_id=normalized,
        expected_profile_id=profile_id(profile),
        expected_profile=profile,
    )
    return write_entry_json(store, normalized, payload)


def _entry_ids_on_disk(root: Path) -> set[str]:
    entries_root = root / "entries"
    if not entries_root.exists():
        return set()
    if entries_root.is_symlink() or not entries_root.is_dir():
        raise InputValidationError(
            "UNSAFE_PRECOMPUTED_PATH", "Precomputed entries directory is unsafe."
        )
    result: set[str] = set()
    for path in entries_root.rglob("*"):
        metadata = path.lstat()
        if stat.S_ISLNK(metadata.st_mode):
            raise InputValidationError(
                "UNSAFE_PRECOMPUTED_PATH", "Precomputed entries contain a symlink."
            )
        if not stat.S_ISREG(metadata.st_mode):
            continue
        if path.name.startswith(".") and path.name.endswith(".tmp"):
            raise InputValidationError(
                "UNSAFE_PRECOMPUTED_PATH", "Precomputed store contains a temporary file."
            )
        if not path.name.endswith(".json.gz"):
            raise InputValidationError(
                "UNSAFE_PRECOMPUTED_PATH", "Precomputed entries contain an unexpected file."
            )
        normalized = normalize_pdb_id(path.name.removesuffix(".json.gz"))
        expected = entry_path(root, normalized)
        if path != expected or normalized in result:
            raise InputValidationError(
                "PRECOMPUTED_ENTRY_LAYOUT_INVALID",
                "Precomputed entry is duplicated or outside its canonical shard.",
            )
        result.add(normalized)
    return result


def precompute_sources(store: Path | str, paths: Sequence[Path]) -> dict[str, Any]:
    """Build and publish one store from explicit source files."""
    from .. import pipeline
    from ..reference_data import load_pdb_fasta_headers
    from ..unknown_molecule_uniprot import process_molecule_info

    root = store_root(store, create=True)
    published_manifest = manifest_path(root)
    if published_manifest.exists() or published_manifest.is_symlink():
        raise InputValidationError(
            "PRECOMPUTED_STORE_PUBLISHED",
            "A published precomputed store is immutable; build into a new directory.",
        )
    profile = scientific_profile()
    selected_profile = profile_id(profile)
    pipeline._validate_required_reference_files()

    raw_paths = [Path(path).expanduser() for path in paths]
    canonical_paths: list[Path] = []
    for path in raw_paths:
        try:
            metadata = path.lstat()
        except OSError as exc:
            raise InputValidationError(
                "INVALID_PRECOMPUTE_SOURCE", f"Cannot inspect precompute source: {path}"
            ) from exc
        if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISREG(metadata.st_mode):
            raise InputValidationError(
                "INVALID_PRECOMPUTE_SOURCE", f"Unsupported precompute source: {path}"
            )
        resolved = path.resolve()
        if not is_valid_file(str(resolved)):
            raise InputValidationError(
                "INVALID_PRECOMPUTE_SOURCE", f"Unsupported precompute source: {path}"
            )
        canonical_paths.append(resolved)
    sources = _source_map(canonical_paths)
    if not sources:
        raise InputValidationError(
            "NO_VALID_INPUT_FILES", "At least one precompute source is required."
        )
    fingerprints = {
        pdb_id: source_fingerprint(path) for pdb_id, path in sources.items()
    }

    report: dict[str, Any] = {
        "artifact_schema_version": ARTIFACT_SCHEMA_VERSION,
        "profile_id": selected_profile,
        "written": 0,
        "cache_hits": 0,
        "failed": 0,
        "entries": [],
        "errors": [],
    }
    pending: list[Path] = []
    for pdb_id, path in sorted(sources.items()):
        if _existing_entry_matches(
            root,
            pdb_id,
            profile=profile,
            source_metadata=fingerprints[pdb_id],
        ):
            report["cache_hits"] += 1
            report["entries"].append({"pdb_id": pdb_id, "status": "reused"})
        else:
            pending.append(path)

    if pending:
        pending_ids: set[str] = set()
        for path in pending:
            pdb_id = get_structure_identity(str(path)).canonical_id
            pending_ids.add(pdb_id)
            legacy = identity_from_official_id(pdb_id).legacy_id
            if legacy:
                pending_ids.add(legacy)
        pdb_fasta_headers = load_pdb_fasta_headers(
            str(config.get("pdb_fasta_path") or ""), tuple(sorted(pending_ids))
        )
        inventory = pipeline.inspect_input_files([str(path) for path in pending])
        batches = pipeline.create_processing_batches(
            [str(path) for path in pending], inventory
        )
        for batch_paths_raw in batches:
            batch_paths = [Path(path).resolve() for path in batch_paths_raw]
            expected = {
                get_structure_identity(str(path)).canonical_id: path
                for path in batch_paths
            }
            try:
                batch_data = pipeline._parse_input_files([str(path) for path in batch_paths])
                process_molecule_info(batch_data, pdb_fasta_headers=pdb_fasta_headers)
                pipeline._run_blast_annotation(batch_data)
                batch_edges = calculate_distances_with_ckdtree(batch_data)
            except Exception as exc:
                for pdb_id in expected:
                    report["failed"] += 1
                    report["errors"].append(
                        {
                            "pdb_id": pdb_id,
                            "code": getattr(exc, "code", exc.__class__.__name__),
                            "message": str(exc)[:500],
                        }
                    )
                tree_cache.clear()
                coords_cache.clear()
                continue

            produced: set[str] = set()
            for structure in batch_data:
                pdb_id: str | None = None
                try:
                    identity_raw = structure.get("structure_identity", {})
                    if isinstance(identity_raw, Mapping) and identity_raw.get(
                        "canonical_id"
                    ):
                        pdb_id = StructureIdentity.from_mapping(
                            identity_raw
                        ).canonical_id
                    else:
                        pdb_id = normalize_pdb_id(structure.get("pdb_id"))
                    source = expected.get(pdb_id)
                    if source is None:
                        raise InputValidationError(
                            "PRECOMPUTED_PDB_ID_MISMATCH",
                            f"Parsed PDB ID {pdb_id.upper()} did not match the source inventory.",
                        )
                    edges = [
                        edge
                        for edge in batch_edges
                        if str(edge.get("structure_key") or "") == pdb_id
                    ]
                    if source_fingerprint(source) != fingerprints[pdb_id]:
                        raise InputValidationError(
                            "INPUT_CHANGED_DURING_PROCESSING",
                            f"Source {source.name} changed during precompute.",
                        )
                    artifact = _write_structure_entry(
                        root,
                        source,
                        structure,
                        edges,
                        profile=profile,
                        source_metadata=fingerprints[pdb_id],
                    )
                    produced.add(pdb_id)
                    report["written"] += 1
                    report["entries"].append(
                        {
                            "pdb_id": pdb_id,
                            "status": "written",
                            "artifact": artifact.relative_to(root).as_posix(),
                        }
                    )
                except Exception as exc:
                    report["failed"] += 1
                    report["errors"].append(
                        {
                            "pdb_id": pdb_id or str(structure.get("pdb_id") or "UNKNOWN"),
                            "code": getattr(exc, "code", exc.__class__.__name__),
                            "message": str(exc)[:500],
                        }
                    )
            for missing in sorted(set(expected) - produced):
                if any(error.get("pdb_id") == missing for error in report["errors"]):
                    continue
                report["failed"] += 1
                report["errors"].append(
                    {
                        "pdb_id": missing,
                        "code": "NO_PARSEABLE_STRUCTURES",
                        "message": "The source did not produce a parseable structure.",
                    }
                )
            tree_cache.clear()
            coords_cache.clear()

    if report["failed"]:
        return report

    try:
        disk_ids = _entry_ids_on_disk(root)
        if disk_ids != set(sources):
            raise InputValidationError(
                "PRECOMPUTED_ENTRY_SET_MISMATCH",
                "Unpublished store entries do not match the requested offline build.",
            )
        for pdb_id in sorted(sources):
            raw, _expanded = read_entry_json(root, pdb_id)
            validate_entry(
                raw,
                expected_pdb_id=pdb_id,
                expected_profile_id=selected_profile,
                expected_profile=profile,
            )
        write_manifest_json(root, manifest_document(profile, len(sources)))
    except Exception as exc:
        report["failed"] += 1
        report["errors"].append(
            {
                "pdb_id": "STORE",
                "code": getattr(exc, "code", exc.__class__.__name__),
                "message": str(exc)[:500],
            }
        )
        return report
    report["manifest"] = "manifest.json"
    return report


def precompute_directory(
    store: Path | str, input_dir: Path | str, *, recursive: bool
) -> dict[str, Any]:
    return precompute_sources(
        store, discover_source_files(input_dir, recursive=recursive)
    )


__all__ = ["discover_source_files", "precompute_directory", "precompute_sources"]
