"""Helpers for per-run output directories and runtime logs."""

from __future__ import annotations

import json
import os
import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Mapping

from . import __version__
from .input_contract import InputValidationError

OUTPUT_CONTRACT_VERSION = "1.1"


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


def collect_web_outputs(output_paths: RunOutputPaths, web_output_dir: str) -> dict[str, Any]:
    """Copy run artifacts into the stable filesystem contract for web workers."""
    web_root = Path(web_output_dir)
    networks_dir = web_root / "networks"
    interactions_dir = web_root / "interactions"
    networks_dir.mkdir(parents=True, exist_ok=True)
    interactions_dir.mkdir(parents=True, exist_ok=True)

    generated = collect_generated_outputs(output_paths)
    copied_networks = [
        str(_copy_unique(Path(path), networks_dir))
        for path in generated["network_files"]
        if Path(path).exists()
    ]
    copied_interactions = [
        str(_copy_unique(Path(path), interactions_dir))
        for path in generated["interaction_csv_files"]
        if Path(path).exists()
    ]

    runtime_analysis = None
    if generated.get("runtime_analysis"):
        runtime_analysis = str(_copy_unique(Path(generated["runtime_analysis"]), web_root))

    source_summary = Path(output_paths.summary_file)
    if not source_summary.exists():
        source_summary = Path(output_paths.manifest_file)
    internal_summary = {}
    if source_summary.exists():
        internal_summary = json.loads(source_summary.read_text(encoding="utf-8"))

    summary = {
        "output_contract_version": internal_summary.get("output_contract_version", OUTPUT_CONTRACT_VERSION),
        "pdb2net_version": internal_summary.get("pdb2net_version", __version__),
        "status": internal_summary.get("status", "completed"),
        "started_at": internal_summary.get("started_at"),
        "finished_at": internal_summary.get("finished_at"),
        "input_files": internal_summary.get("input_files", []),
        "input_path": internal_summary.get("input_path"),
        "run_summary": str(source_summary) if source_summary.exists() else None,
        "internal_run_output_path": output_paths.run_output_path,
        "networks": copied_networks,
        "interactions": copied_interactions,
        "runtime_analysis": runtime_analysis,
        "counts": {
            "networks": len(copied_networks),
            "interactions": len(copied_interactions),
            "skipped_outputs": len(internal_summary.get("skipped_outputs", [])),
        },
        "annotations": internal_summary.get("annotations", {}),
        "references": internal_summary.get("references", {}),
        "resources": internal_summary.get("resources", {}),
        "skipped_outputs": internal_summary.get("skipped_outputs", []),
        "config": internal_summary.get("config", {}),
        "errors": internal_summary.get("errors", []),
        "warnings": internal_summary.get("warnings", []),
    }

    summary_path = web_root / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
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
