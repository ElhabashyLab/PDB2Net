"""Helpers for per-run output directories and runtime logs."""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from datetime import datetime
from typing import Any, Mapping

from .input_contract import InputValidationError


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
    )


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
) -> None:
    """Write an additive per-run manifest for webserver/job tracking."""
    manifest = {
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
        },
        "config": dict(config_snapshot),
        "errors": list(errors or []),
        "warnings": list(warnings or []),
    }

    with open(manifest_file, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)


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
    return output_paths
