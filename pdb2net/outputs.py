"""Helpers for per-run output directories and runtime logs."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Mapping

from . import __version__
from .contracts import OUTPUT_CONTRACT_VERSION
from .input_contract import InputValidationError
from .structure_identity import StructureIdentity

MAX_WEB_SUMMARY_BYTES = 1_000_000
MAX_WEB_ARTIFACTS = 1_000
MAX_WEB_ARTIFACT_BYTES = 512 * 1024 * 1024
MAX_WEB_ARTIFACT_TOTAL_BYTES = 2 * 1024 * 1024 * 1024
_FORBIDDEN_CX1_ASPECTS = {"nodeAttributes", "edgeAttributes", "cartesianLayout"}
_CX2_ASPECT_ORDER = (
    "metaData",
    "attributeDeclarations",
    "networkAttributes",
    "nodes",
    "edges",
    "visualProperties",
    "status",
)


@dataclass(frozen=True)
class RunOutputPaths:
    """Filesystem locations created for a single pipeline run."""

    run_output_path: str
    combined_dir: str
    protein_dir: str
    chain_dir: str
    distances_dir: str
    log_file: str
    manifest_file: str
    summary_file: str


def create_run_output_paths(base_output_path: str, timestamp: str | None = None) -> RunOutputPaths:
    """Create and return the standard output folders for one run."""
    run_timestamp = timestamp or datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    run_output_path = os.path.join(base_output_path, run_timestamp)
    os.makedirs(run_output_path, exist_ok=True)

    combined_dir = os.path.join(run_output_path, "combined")
    protein_dir = os.path.join(run_output_path, "protein")
    chain_dir = os.path.join(run_output_path, "chain")
    distances_dir = os.path.join(run_output_path, "distances")

    for directory in (combined_dir, protein_dir, chain_dir, distances_dir):
        os.makedirs(directory, exist_ok=True)

    return RunOutputPaths(
        run_output_path=run_output_path,
        combined_dir=combined_dir,
        protein_dir=protein_dir,
        chain_dir=chain_dir,
        distances_dir=distances_dir,
        log_file=os.path.join(run_output_path, "runtime_analysis.txt"),
        manifest_file=os.path.join(run_output_path, "manifest.json"),
        summary_file=os.path.join(run_output_path, "run_summary.json"),
    )


def _sorted_files(directory: str, suffix: str) -> list[str]:
    """Return generated files with stable ordering for manifests."""
    root = Path(directory)
    if not root.exists():
        return []
    return [str(path) for path in sorted(root.glob(f"*{suffix}")) if path.is_file()]


def collect_generated_outputs(output_paths: RunOutputPaths) -> dict[str, Any]:
    """Collect generated run artifacts without changing the internal layout."""
    network_files = (
        _sorted_files(output_paths.combined_dir, ".cx2")
        + _sorted_files(output_paths.protein_dir, ".cx2")
        + _sorted_files(output_paths.chain_dir, ".cx2")
    )
    interaction_csv_files = _sorted_files(output_paths.distances_dir, ".csv")
    runtime_analysis = output_paths.log_file if Path(output_paths.log_file).exists() else None

    return {
        "network_files": network_files,
        "interaction_csv_files": interaction_csv_files,
        "runtime_analysis": runtime_analysis,
        "counts": {
            "network_files": len(network_files),
            "interaction_csv_files": len(interaction_csv_files),
        },
    }


def write_runtime_analysis(
    log_file: str,
    sum_times: Mapping[str, float],
    total_time: float,
) -> None:
    """Write the existing runtime summary file for a run."""
    with open(log_file, "w", encoding="utf-8") as handle:
        handle.write("===== PDB2Net Batch Log =====\n\n")
        handle.write("Timing (total):\n")
        handle.write(f"- Parsing:                 {sum_times['parsing']:.1f} sec\n")
        handle.write(f"- Classification (SIFTS):  {sum_times['sifts']:.1f} sec\n")
        handle.write(f"- Classification (BLAST):  {sum_times['blast']:.1f} sec\n")
        handle.write(f"- Interaction:             {sum_times['interaction']:.1f} sec\n")
        handle.write(f"- Network export:          {sum_times['networks']:.1f} sec\n")
        handle.write(f"- Total:                   {total_time:.1f} sec\n")
        handle.write("\n===============================\n")


def write_run_manifest(
    manifest_file: str,
    *,
    input_files: list[str],
    output_paths: RunOutputPaths,
    config_snapshot: Mapping[str, Any],
    status: str,
    started_at: str,
    finished_at: str,
    total_time: float,
    errors: list[Any] | None = None,
    warnings: list[Any] | None = None,
    input_path: str | None = None,
    generated_outputs: Mapping[str, Any] | None = None,
    extra_counts: Mapping[str, Any] | None = None,
    annotations: Mapping[str, Any] | None = None,
    references: Mapping[str, Any] | None = None,
    resources: Mapping[str, Any] | None = None,
    skipped_outputs: list[Mapping[str, Any]] | None = None,
    identities: list[Mapping[str, Any]] | None = None,
    structure_inputs: list[Mapping[str, Any]] | None = None,
) -> None:
    """Write an additive per-run manifest for webserver/job tracking."""
    generated = dict(generated_outputs or collect_generated_outputs(output_paths))
    counts = {
        "input_files": len(input_files),
        **dict(generated.get("counts", {})),
        **dict(extra_counts or {}),
    }
    manifest = {
        "output_contract_version": OUTPUT_CONTRACT_VERSION,
        "pdb2net_version": __version__,
        "status": status,
        "started_at": started_at,
        "finished_at": finished_at,
        "total_seconds": round(float(total_time), 3),
        "input_files": list(input_files),
        "input_path": input_path,
        "identities": [dict(entry) for entry in (identities or [])],
        "structure_inputs": [dict(entry) for entry in (structure_inputs or [])],
        "outputs": {
            "run_output_path": output_paths.run_output_path,
            "combined_dir": output_paths.combined_dir,
            "protein_dir": output_paths.protein_dir,
            "chain_dir": output_paths.chain_dir,
            "distances_dir": output_paths.distances_dir,
            "runtime_analysis": output_paths.log_file,
            "network_files": list(generated.get("network_files", [])),
            "interaction_csv_files": list(generated.get("interaction_csv_files", [])),
        },
        "counts": counts,
        "annotations": dict(annotations or {}),
        "references": dict(references or {}),
        "resources": dict(resources or {}),
        "skipped_outputs": [dict(entry) for entry in (skipped_outputs or [])],
        "config": dict(config_snapshot),
        "errors": list(errors or []),
        "warnings": list(warnings or []),
    }

    with open(manifest_file, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)


def write_run_summary(output_paths: RunOutputPaths, manifest_file: str | None = None) -> None:
    """Mirror the run manifest to the stable run_summary.json filename."""
    source = Path(manifest_file or output_paths.manifest_file)
    if source.exists():
        shutil.copy2(source, output_paths.summary_file)


def _copy_unique(source: Path, destination_dir: Path) -> Path:
    """Copy a file into a flat directory, preserving names unless a collision occurs."""
    if source.is_symlink() or not source.is_file():
        raise InputValidationError(
            "UNSAFE_WEB_ARTIFACT",
            f"Generated artifact is missing or is a symlink: {source.name}",
        )
    candidate = destination_dir / source.name
    if not candidate.exists():
        shutil.copy2(source, candidate)
        return candidate

    stem = source.stem
    suffix = source.suffix
    index = 2
    while True:
        candidate = destination_dir / f"{stem}_{index}{suffix}"
        if not candidate.exists():
            shutil.copy2(source, candidate)
            return candidate
        index += 1


def _safe_input_label(value: object, *, fallback: str = "input") -> str:
    """Return a portable public label without exposing a server-side path."""
    raw = str(value or "").replace("\x00", "").strip()
    if raw.startswith("pdb:") and "/" not in raw and "\\" not in raw:
        label = raw
    else:
        label = raw.replace("\\", "/").rsplit("/", 1)[-1]
    label = "".join(character for character in label if character >= " " and character != "\x7f")
    encoded = label.encode("utf-8")
    if not label:
        return fallback
    if len(encoded) > 512:
        encoded = encoded[:512]
        label = encoded.decode("utf-8", errors="ignore") or fallback
    return label


def _redact_text(value: object, *, maximum: int = 1_000) -> str:
    """Remove path-shaped substrings from a public diagnostic string."""
    text = str(value or "").replace("\x00", "")
    text = re.sub(r"(?<![A-Za-z0-9_])/(?:[^\s,;:]+/)+([^\s,;:]+)", r"\1", text)
    text = re.sub(r"(?i)(?<![A-Za-z0-9_])[A-Z]:\\(?:[^\s,;:]+\\)+([^\s,;:]+)", r"\1", text)
    return text[:maximum]


def _public_diagnostics(entries: object) -> list[dict[str, Any]]:
    if not isinstance(entries, list):
        return []
    public: list[dict[str, Any]] = []
    for entry in entries[:100]:
        if not isinstance(entry, Mapping):
            continue
        code = str(entry.get("code") or "UNKNOWN")[:100]
        item: dict[str, Any] = {"code": code, "message": _redact_text(entry.get("message"))}
        if isinstance(entry.get("output"), str):
            item["output"] = _safe_input_label(entry["output"], fallback="output")
        public.append(item)
    return public


def _public_config(value: object) -> dict[str, Any]:
    """Expose only runtime choices a web submitter can select."""
    raw = value if isinstance(value, Mapping) else {}
    allowed = (
        "networks",
        "network_annotations",
        "distance_thresholds",
        "interaction_filters",
        "structure_model_policy",
        "export_detailed_interactions",
    )
    result: dict[str, Any] = {}
    for key in allowed:
        item = raw.get(key)
        if isinstance(item, Mapping):
            result[key] = dict(item)
        elif key in raw:
            result[key] = item
    return result


def _public_references(value: object) -> dict[str, Any]:
    raw = value if isinstance(value, Mapping) else {}
    result: dict[str, Any] = {"manifest_id": raw.get("manifest_id")}
    precomputed = raw.get("precomputed_store")
    if isinstance(precomputed, Mapping):
        result["precomputed_store"] = {
            key: precomputed[key]
            for key in (
                "artifact_schema_version",
                "profile_id",
                "source_scope",
            )
            if key in precomputed
        }
    return result


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _json_without_nonfinite(path: Path) -> Any:
    def reject_constant(value: str) -> None:
        raise ValueError(f"Non-finite JSON number {value} is not allowed.")

    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle, parse_constant=reject_constant)


def _cx2_counts(path: Path) -> tuple[int, int]:
    """Validate the native CX2 aspects written by Core and return node/edge counts."""
    try:
        document = _json_without_nonfinite(path)
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as exc:
        raise InputValidationError("INVALID_CX2_ARTIFACT", f"Invalid CX2 artifact {path.name}.") from exc
    if not isinstance(document, list) or len(document) < 5:
        raise InputValidationError("INVALID_CX2_ARTIFACT", f"Invalid CX2 document {path.name}.")
    header = document[0]
    if header != {"CXVersion": "2.0", "hasFragments": False}:
        raise InputValidationError("INVALID_CX2_ARTIFACT", f"CX2 version is invalid in {path.name}.")
    aspects: dict[str, Any] = {}
    aspect_order: list[str] = []
    for block in document[1:]:
        if not isinstance(block, dict) or len(block) != 1:
            raise InputValidationError(
                "INVALID_CX2_ARTIFACT", f"CX2 aspect blocks must have one field in {path.name}."
            )
        name, values = next(iter(block.items()))
        if name in aspects:
            raise InputValidationError("INVALID_CX2_ARTIFACT", f"Duplicate CX2 aspect {name} in {path.name}.")
        aspects[name] = values
        aspect_order.append(name)
    if tuple(aspect_order) != _CX2_ASPECT_ORDER:
        raise InputValidationError(
            "INVALID_CX2_ARTIFACT", f"CX2 native aspects are missing or out of order in {path.name}."
        )
    forbidden = _FORBIDDEN_CX1_ASPECTS.intersection(aspects)
    if forbidden:
        raise InputValidationError("INVALID_CX2_ARTIFACT", f"CX1-style aspects are forbidden in {path.name}.")
    nodes = aspects.get("nodes")
    edges = aspects.get("edges")
    if not isinstance(nodes, list) or not isinstance(edges, list):
        raise InputValidationError("INVALID_CX2_ARTIFACT", f"CX2 nodes/edges are missing in {path.name}.")
    node_ids: set[int] = set()
    declarations = aspects.get("attributeDeclarations")
    if (
        not isinstance(declarations, list)
        or len(declarations) != 1
        or not isinstance(declarations[0], dict)
        or set(declarations[0]) != {"networkAttributes", "nodes", "edges"}
        or any(not isinstance(declarations[0][scope], dict) for scope in declarations[0])
    ):
        raise InputValidationError(
            "INVALID_CX2_ARTIFACT", f"CX2 attribute declarations are invalid in {path.name}."
        )
    declared = declarations[0]
    supported_types = {"string", "long", "integer", "double", "boolean"}
    for scope in declared.values():
        if any(
            not isinstance(specification, dict)
            or set(specification) != {"d"}
            or specification.get("d") not in supported_types
            for specification in scope.values()
        ):
            raise InputValidationError(
                "INVALID_CX2_ARTIFACT", f"CX2 attribute declarations are invalid in {path.name}."
            )

    def declared_value(scope: str, key: str, value: object) -> bool:
        specification = declared[scope].get(key)
        if not isinstance(specification, dict):
            return False
        data_type = specification.get("d")
        if data_type == "string":
            return isinstance(value, str)
        if data_type == "boolean":
            return isinstance(value, bool)
        if data_type in {"integer", "long"}:
            return isinstance(value, int) and not isinstance(value, bool)
        if data_type == "double":
            return (
                isinstance(value, (int, float))
                and not isinstance(value, bool)
                and math.isfinite(float(value))
            )
        return False

    network_attributes = aspects.get("networkAttributes")
    visual_properties = aspects.get("visualProperties")
    if (
        not isinstance(network_attributes, list)
        or len(network_attributes) != 1
        or not isinstance(network_attributes[0], dict)
        or not set(network_attributes[0]).issubset(declared["networkAttributes"])
        or any(
            not declared_value("networkAttributes", key, value)
            for key, value in network_attributes[0].items()
        )
        or not isinstance(visual_properties, list)
        or len(visual_properties) != 1
        or not isinstance(visual_properties[0], dict)
    ):
        raise InputValidationError("INVALID_CX2_ARTIFACT", f"CX2 core aspects are invalid in {path.name}.")

    for expected_id, node in enumerate(nodes):
        if not isinstance(node, dict) or set(node) != {"id", "x", "y", "v"}:
            raise InputValidationError("INVALID_CX2_ARTIFACT", f"CX2 node fields are invalid in {path.name}.")
        node_id = node.get("id")
        if isinstance(node_id, bool) or not isinstance(node_id, int) or node_id != expected_id:
            raise InputValidationError("INVALID_CX2_ARTIFACT", f"CX2 node IDs are invalid in {path.name}.")
        if (
            not isinstance(node.get("v"), dict)
            or "id" in node["v"]
            or not set(node["v"]).issubset(declared["nodes"])
            or any(
                not declared_value("nodes", key, value)
                for key, value in node["v"].items()
            )
        ):
            raise InputValidationError("INVALID_CX2_ARTIFACT", f"CX2 node attributes are invalid in {path.name}.")
        if any(
            isinstance(node.get(axis), bool)
            or not isinstance(node.get(axis), (int, float))
            or not math.isfinite(float(node[axis]))
            for axis in ("x", "y")
        ):
            raise InputValidationError("INVALID_CX2_ARTIFACT", f"CX2 coordinates are invalid in {path.name}.")
        node_ids.add(node_id)
    edge_ids: set[int] = set()
    for expected_id, edge in enumerate(edges):
        if not isinstance(edge, dict) or set(edge) != {"id", "s", "t", "v"}:
            raise InputValidationError("INVALID_CX2_ARTIFACT", f"CX2 edge fields are invalid in {path.name}.")
        edge_id = edge.get("id")
        if isinstance(edge_id, bool) or not isinstance(edge_id, int) or edge_id != expected_id:
            raise InputValidationError("INVALID_CX2_ARTIFACT", f"CX2 edge IDs are invalid in {path.name}.")
        if (
            edge.get("s") not in node_ids
            or edge.get("t") not in node_ids
            or edge.get("s") == edge.get("t")
            or not isinstance(edge.get("v"), dict)
            or "id" in edge["v"]
            or not set(edge["v"]).issubset(declared["edges"])
            or any(
                not declared_value("edges", key, value)
                for key, value in edge["v"].items()
            )
        ):
            raise InputValidationError("INVALID_CX2_ARTIFACT", f"CX2 edge endpoints are invalid in {path.name}.")
        edge_ids.add(edge_id)
    metadata = aspects.get("metaData")
    if not isinstance(metadata, list):
        raise InputValidationError("INVALID_CX2_ARTIFACT", f"CX2 metadata is missing in {path.name}.")
    metadata_counts = {
        entry.get("name"): entry.get("elementCount")
        for entry in metadata
        if isinstance(entry, dict)
    }
    expected_metadata = {
        "attributeDeclarations": len(declarations),
        "networkAttributes": len(network_attributes),
        "nodes": len(nodes),
        "edges": len(edges),
        "visualProperties": len(visual_properties),
    }
    if metadata_counts != expected_metadata:
        raise InputValidationError("INVALID_CX2_ARTIFACT", f"CX2 metadata counts are invalid in {path.name}.")
    if document[-1] != {"status": [{"error": "", "success": True}]}:
        raise InputValidationError("INVALID_CX2_ARTIFACT", f"CX2 success status is invalid in {path.name}.")
    return len(nodes), len(edges)


def _csv_counts(path: Path) -> tuple[int, list[str]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle)
            header = next(reader)
            if not header or any(not column for column in header) or len(set(header)) != len(header):
                raise ValueError("invalid header")
            rows = 0
            for row in reader:
                if len(row) != len(header):
                    raise ValueError("invalid row width")
                rows += 1
    except (OSError, UnicodeError, csv.Error, StopIteration, ValueError) as exc:
        raise InputValidationError("INVALID_CSV_ARTIFACT", f"Invalid CSV artifact {path.name}.") from exc
    return rows, header


def _artifact_record(path: Path, web_root: Path, *, kind: str) -> dict[str, Any]:
    size = path.stat().st_size
    if size > MAX_WEB_ARTIFACT_BYTES:
        raise InputValidationError("WEB_ARTIFACT_TOO_LARGE", f"Artifact {path.name} exceeds the size limit.")
    record: dict[str, Any] = {
        "path": path.relative_to(web_root).as_posix(),
        "size_bytes": size,
        "sha256": _sha256_file(path),
    }
    if kind == "network":
        record["nodes"], record["edges"] = _cx2_counts(path)
    else:
        record["rows"], record["columns"] = _csv_counts(path)
    return record


def _public_input_metadata(internal_summary: Mapping[str, Any], *, success: bool) -> tuple[list[str], list[dict[str, Any]], list[dict[str, Any]]]:
    raw_files = internal_summary.get("input_files")
    raw_identities = internal_summary.get("identities")
    raw_structures = internal_summary.get("structure_inputs")
    if not isinstance(raw_files, list) or not isinstance(raw_identities, list) or not isinstance(raw_structures, list):
        raise InputValidationError("INVALID_WEB_OUTPUT_CONTRACT", "Input metadata arrays are invalid.")
    files = [_safe_input_label(value, fallback=f"input-{index + 1}") for index, value in enumerate(raw_files)]
    identities: list[dict[str, Any]] = []
    for value in raw_identities:
        if not isinstance(value, Mapping):
            continue
        try:
            identity = StructureIdentity.from_mapping(value)
        except (TypeError, ValueError) as exc:
            raise InputValidationError(
                "INVALID_WEB_OUTPUT_CONTRACT", "Input identity metadata is invalid."
            ) from exc
        normalized = identity.as_dict()
        if dict(value) != normalized:
            raise InputValidationError(
                "INVALID_WEB_OUTPUT_CONTRACT", "Input identity fields are not canonical."
            )
        identities.append(normalized)
    structures: list[dict[str, Any]] = []
    for index, value in enumerate(raw_structures):
        if not isinstance(value, Mapping):
            continue
        if set(value) != {"file", "identity", "format", "kind", "embedded_annotation_counts"}:
            raise InputValidationError(
                "INVALID_WEB_OUTPUT_CONTRACT", "Structure input descriptor fields are invalid."
            )
        counts = value.get("embedded_annotation_counts")
        if (
            not isinstance(counts, Mapping)
            or set(counts) != {"uniprot", "pfam", "cath", "scop2"}
            or any(
                isinstance(count, bool) or not isinstance(count, int) or count < 0
                for count in counts.values()
            )
        ):
            raise InputValidationError(
                "INVALID_WEB_OUTPUT_CONTRACT", "Embedded annotation counts are invalid."
            )
        input_format = value.get("format")
        input_kind = value.get("kind")
        if input_format not in {"pdb", "mmcif"} or input_kind not in {
            "pdb",
            "mmcif",
            "nextgen_enriched_mmcif",
        }:
            raise InputValidationError(
                "INVALID_WEB_OUTPUT_CONTRACT", "Structure input format metadata is invalid."
            )
        public = {
            "file": _safe_input_label(value.get("file"), fallback=f"input-{index + 1}"),
            "identity": dict(value.get("identity")) if isinstance(value.get("identity"), Mapping) else {},
            "format": input_format,
            "kind": input_kind,
            "embedded_annotation_counts": dict(counts),
        }
        structures.append(public)
    if success:
        if not files or len(files) != len(identities) or len(files) != len(structures):
            raise InputValidationError(
                "INVALID_WEB_OUTPUT_CONTRACT",
                "Successful web summaries require one ordered identity and structure descriptor per input.",
            )
        keys: list[str] = []
        for index, (identity, structure) in enumerate(zip(identities, structures)):
            key = identity.get("key")
            if (
                not isinstance(key, str)
                or not key
                or structure.get("identity") != identity
                or structure.get("file") != files[index]
            ):
                raise InputValidationError("INVALID_WEB_OUTPUT_CONTRACT", "Input identities are inconsistent.")
            keys.append(key)
        if len(set(keys)) != len(keys):
            raise InputValidationError("INVALID_WEB_OUTPUT_CONTRACT", "Input identities are not unique.")
    return files, identities, structures


def _write_public_summary(path: Path, summary: Mapping[str, Any]) -> None:
    serialized = json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False).encode("utf-8") + b"\n"
    if len(serialized) > MAX_WEB_SUMMARY_BYTES:
        raise InputValidationError("WEB_SUMMARY_TOO_LARGE", "Public summary exceeds the 1 MiB contract limit.")
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("wb") as handle:
        handle.write(serialized)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
    descriptor = os.open(path.parent, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def collect_web_outputs(output_paths: RunOutputPaths, web_output_dir: str) -> dict[str, Any]:
    """Copy run artifacts into the stable filesystem contract for web workers."""
    web_root = Path(web_output_dir)
    networks_dir = web_root / "networks"
    interactions_dir = web_root / "interactions"
    networks_dir.mkdir(parents=True, exist_ok=True)
    interactions_dir.mkdir(parents=True, exist_ok=True)

    source_summary = Path(output_paths.summary_file)
    if not source_summary.exists():
        source_summary = Path(output_paths.manifest_file)
    internal_summary = {}
    if source_summary.exists():
        internal_summary = json.loads(source_summary.read_text(encoding="utf-8"))

    status = internal_summary.get("status", "failed")
    success = status == "success"
    input_files, identities, structure_inputs = _public_input_metadata(internal_summary, success=success)
    if success and internal_summary.get("errors"):
        raise InputValidationError("INVALID_WEB_OUTPUT_CONTRACT", "Successful runs must not contain errors.")
    references = _public_references(internal_summary.get("references"))
    if success and (
        not isinstance(references.get("manifest_id"), str)
        or not str(references["manifest_id"]).strip()
    ):
        raise InputValidationError("INVALID_WEB_OUTPUT_CONTRACT", "Successful runs require a reference manifest ID.")

    # Failed/partial manifests are diagnostic-only.  Decide that before any
    # generated file is copied so a failed run cannot leave unlisted public
    # artifacts behind.
    generated = collect_generated_outputs(output_paths) if success else {
        "network_files": [],
        "interaction_csv_files": [],
        "runtime_analysis": None,
    }
    network_sources = [Path(path) for path in generated["network_files"]]
    interaction_sources = [Path(path) for path in generated["interaction_csv_files"]]
    artifact_sources = network_sources + interaction_sources
    if len(artifact_sources) > MAX_WEB_ARTIFACTS:
        raise InputValidationError("TOO_MANY_WEB_ARTIFACTS", "Web output exceeds the 1,000-artifact limit.")
    artifact_source_bytes = 0
    for source in artifact_sources:
        if source.is_symlink() or not source.is_file():
            raise InputValidationError(
                "UNSAFE_WEB_ARTIFACT",
                f"Generated artifact is missing or is a symlink: {source.name}",
            )
        size = source.stat().st_size
        if size > MAX_WEB_ARTIFACT_BYTES:
            raise InputValidationError(
                "WEB_ARTIFACT_TOO_LARGE", f"Artifact {source.name} exceeds the size limit."
            )
        artifact_source_bytes += size
    if artifact_source_bytes > MAX_WEB_ARTIFACT_TOTAL_BYTES:
        raise InputValidationError("WEB_ARTIFACT_TOTAL_TOO_LARGE", "Web artifacts exceed the aggregate size limit.")
    copied_network_paths = [_copy_unique(path, networks_dir) for path in network_sources]
    copied_interaction_paths = [_copy_unique(path, interactions_dir) for path in interaction_sources]
    runtime_analysis = None
    if generated.get("runtime_analysis"):
        runtime_analysis = _copy_unique(Path(generated["runtime_analysis"]), web_root).relative_to(web_root).as_posix()
    network_artifacts = [
        _artifact_record(path, web_root, kind="network") for path in copied_network_paths
    ]
    interaction_artifacts = [
        _artifact_record(path, web_root, kind="interaction") for path in copied_interaction_paths
    ]
    artifact_bytes = sum(
        entry["size_bytes"] for entry in network_artifacts + interaction_artifacts
    )
    if artifact_bytes > MAX_WEB_ARTIFACT_TOTAL_BYTES:
        raise InputValidationError("WEB_ARTIFACT_TOTAL_TOO_LARGE", "Web artifacts exceed the aggregate size limit.")
    copied_networks = [entry["path"] for entry in network_artifacts]
    copied_interactions = [entry["path"] for entry in interaction_artifacts]
    internal_counts = internal_summary.get("counts") if isinstance(internal_summary.get("counts"), Mapping) else {}
    summary = {
        "output_contract_version": internal_summary.get("output_contract_version", OUTPUT_CONTRACT_VERSION),
        "pdb2net_version": internal_summary.get("pdb2net_version", __version__),
        "status": status,
        "started_at": internal_summary.get("started_at"),
        "finished_at": internal_summary.get("finished_at"),
        "input_files": input_files,
        "identities": identities if success else [],
        "structure_inputs": structure_inputs if success else [],
        "networks": copied_networks,
        "interactions": copied_interactions,
        "artifacts": {
            "networks": network_artifacts,
            "interactions": interaction_artifacts,
        },
        "runtime_analysis": runtime_analysis,
        "counts": {
            "networks": len(copied_networks),
            "interactions": len(copied_interactions),
            "structures": int(internal_counts.get("structures", len(identities) if success else 0)),
            "chains": int(internal_counts.get("chains", 0)),
            "skipped_outputs": len(internal_summary.get("skipped_outputs", [])),
        },
        "annotations": internal_summary.get("annotations", {}),
        "references": references,
        "resources": internal_summary.get("resources", {}),
        "skipped_outputs": internal_summary.get("skipped_outputs", []),
        "config": _public_config(internal_summary.get("config")),
        "errors": [] if success else _public_diagnostics(internal_summary.get("errors")),
        "warnings": _public_diagnostics(internal_summary.get("warnings")),
    }

    summary_path = web_root / "summary.json"
    _write_public_summary(summary_path, summary)
    return summary


def write_failed_run_manifest(
    base_output_path: str,
    *,
    input_path: str,
    config_snapshot: Mapping[str, Any],
    error: Exception,
    started_at: str | None = None,
) -> RunOutputPaths:
    """Create a failed run folder and write a machine-readable manifest."""
    output_paths = create_run_output_paths(base_output_path)
    finished_at = datetime.now().isoformat(timespec="seconds")
    started = started_at or finished_at

    if isinstance(error, InputValidationError):
        error_entry: dict[str, Any] = {"code": error.code, "message": str(error)}
    else:
        error_entry = {"code": error.__class__.__name__, "message": str(error)}

    write_run_manifest(
        output_paths.manifest_file,
        input_files=[],
        output_paths=output_paths,
        config_snapshot=config_snapshot,
        status="failed",
        started_at=started,
        finished_at=finished_at,
        total_time=0.0,
        errors=[error_entry],
        warnings=[],
        input_path=input_path,
    )
    write_run_summary(output_paths)
    return output_paths
