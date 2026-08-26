"""Read-only assembly from a published schema-3 precomputed store."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
import time
from typing import Any, Mapping, Sequence

from ..config_loader import config
from ..distances import coords_cache, tree_cache
from ..input_contract import InputValidationError
from .io import entry_path, read_entry_json, read_manifest_json, store_root
from .schema import (
    ARTIFACT_SCHEMA_VERSION,
    MAX_REQUEST_CHAINS,
    MAX_REQUEST_EDGES,
    MAX_REQUEST_EXPANDED_BYTES,
    SOURCE_SCOPE,
    materialize_entry,
    normalize_pdb_id,
    normalize_requested_ids,
    scientific_profile,
    validate_entry,
    validate_manifest,
)


def load_manifest(store: Path | str) -> dict[str, Any]:
    """Load the published manifest and require the active fixed profile."""
    profile = scientific_profile()
    return validate_manifest(
        read_manifest_json(store),
        expected_profile=profile,
    )


def _read_validated_entry(
    store: Path | str,
    pdb_id: str,
    *,
    manifest: Mapping[str, Any],
    expanded_budget: int,
) -> tuple[dict[str, Any], int]:
    raw, expanded_bytes = read_entry_json(
        store, pdb_id, expanded_budget=expanded_budget
    )
    payload = validate_entry(
        raw,
        expected_pdb_id=pdb_id,
        expected_profile_id=str(manifest["profile_id"]),
        expected_profile=manifest["profile"],
    )
    return payload, expanded_bytes


def load_entry(store: Path | str, pdb_id: object) -> dict[str, Any]:
    """Load one exact-shape entry from a published store."""
    normalized = normalize_pdb_id(pdb_id)
    manifest = load_manifest(store)
    payload, _expanded = _read_validated_entry(
        store,
        normalized,
        manifest=manifest,
        expanded_budget=MAX_REQUEST_EXPANDED_BYTES,
    )
    return payload


def _load_requested_entries(
    store: Path | str,
    requested: Sequence[str],
    *,
    manifest: Mapping[str, Any],
    resource_limits: Mapping[str, int | None],
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    compressed_total = 0
    largest_entry = 0
    for pdb_id in requested:
        path = entry_path(store, pdb_id)
        if path.is_symlink() or not path.is_file():
            if not path.exists():
                raise InputValidationError(
                    "PRECOMPUTED_ENTRY_MISSING",
                    f"No precomputed entry exists for PDB ID {pdb_id.upper()}.",
                )
            raise InputValidationError(
                "UNSAFE_PRECOMPUTED_ENTRY",
                f"Precomputed entry is not a regular file: {path.name}",
            )
        entry_bytes = path.stat().st_size
        maximum_single = resource_limits.get("max_single_input_bytes")
        if maximum_single is not None and entry_bytes > maximum_single:
            raise InputValidationError(
                "INPUT_FILE_BYTES_LIMIT_EXCEEDED",
                f"Entry {pdb_id.upper()} exceeds the configured per-file limit.",
            )
        compressed_total += entry_bytes
        largest_entry = max(largest_entry, entry_bytes)
        maximum_total = resource_limits.get("max_total_input_bytes")
        if maximum_total is not None and compressed_total > maximum_total:
            raise InputValidationError(
                "INPUT_TOTAL_BYTES_LIMIT_EXCEEDED",
                "Requested entries exceed the configured total-input limit.",
            )

    loaded: list[dict[str, Any]] = []
    expanded_total = 0
    chain_total = 0
    edge_total = 0
    for pdb_id in requested:
        payload, expanded_bytes = _read_validated_entry(
            store,
            pdb_id,
            manifest=manifest,
            expanded_budget=MAX_REQUEST_EXPANDED_BYTES - expanded_total,
        )
        expanded_total += expanded_bytes
        chain_total += int(payload["counts"]["chains"])
        edge_total += int(payload["counts"]["interactions"])
        if chain_total > MAX_REQUEST_CHAINS or edge_total > MAX_REQUEST_EDGES:
            raise InputValidationError(
                "PRECOMPUTED_REQUEST_LIMIT_EXCEEDED",
                "Requested entries exceed aggregate chain or interaction limits.",
            )
        loaded.append(payload)
    return loaded, {
        "compressed_bytes": compressed_total,
        "expanded_bytes": expanded_total,
        "largest_entry_bytes": largest_entry,
    }


def _resource_summary(
    requested: Sequence[str],
    *,
    profile_id: str,
    compressed_bytes: int,
    expanded_bytes: int,
    largest_entry_bytes: int,
) -> dict[str, Any]:
    return {
        "input": {
            "files": len(requested),
            "total_bytes": compressed_bytes,
            "largest_file_bytes": largest_entry_bytes,
        },
        "processing": {
            "mode": "precomputed",
            "cache_hits": len(requested),
            "parsing_workers": 0,
        },
        "precomputed_store": {
            "artifact_schema_version": ARTIFACT_SCHEMA_VERSION,
            "profile_id": profile_id,
            "entries": len(requested),
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
):
    """Assemble normal outputs without ever modifying the published store."""
    from .. import pipeline
    from ..network_annotations import network_annotation_config
    from ..outputs import (
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
    labels: list[str] = []
    try:
        network_annotation_config()
        if bool(config.get("export_detailed_interactions", False)):
            raise InputValidationError(
                "DETAILED_INTERACTIONS_UNAVAILABLE",
                "Precomputed entries contain no atom pairs; use raw-file mode for detailed CSVs.",
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
        root = store_root(store)
        manifest = load_manifest(root)
        output_paths = create_run_output_paths(config["output_path"])
        labels = [f"pdb:{pdb_id}" for pdb_id in requested]
        payloads, sizes = _load_requested_entries(
            root,
            requested,
            manifest=manifest,
            resource_limits=resource_limits,
        )

        structures: list[dict[str, Any]] = []
        interactions: list[dict[str, Any]] = []
        for payload in payloads:
            structure, entry_interactions, _references = materialize_entry(payload)
            structures.append(structure)
            interactions.extend(entry_interactions)

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
                    f"{entry['counts']['edges']} edges because configured limits were exceeded."
                ),
                "output": entry["name"],
            }
            for entry in skipped
        ] + pipeline._embedded_annotation_warnings(structures)
        total = time.monotonic() - start
        write_runtime_analysis(output_paths.log_file, timings.as_dict(), total)
        generated = collect_generated_outputs(output_paths)
        profile = manifest["profile"]
        selected_profile = str(manifest["profile_id"])
        write_run_manifest(
            output_paths.manifest_file,
            input_files=labels,
            output_paths=output_paths,
            config_snapshot=pipeline._config_snapshot(config["networks"]),
            status="success",
            started_at=started_at,
            finished_at=datetime.now().isoformat(timespec="seconds"),
            total_time=total,
            input_path=f"precomputed:{selected_profile}",
            generated_outputs=generated,
            extra_counts={
                "structures": len(structures),
                "chains": sum(len(item.get("atom_data", [])) for item in structures),
                "interactions": len(interactions),
                "skipped_outputs": len(skipped),
            },
            annotations=pipeline._annotation_summary(structures),
            references={
                "manifest_id": profile["annotations"]["reference_manifest_id"],
                "precomputed_store": {
                    "artifact_schema_version": ARTIFACT_SCHEMA_VERSION,
                    "profile_id": selected_profile,
                    "source_scope": SOURCE_SCOPE,
                },
            },
            resources=_resource_summary(
                requested,
                profile_id=selected_profile,
                **sizes,
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
            input_path="precomputed:"
            + ",".join(labels or [str(value) for value in pdb_ids]),
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


__all__ = ["load_entry", "load_manifest", "run_assemble_pipeline"]
