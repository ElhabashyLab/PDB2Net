"""Helpers for per-run output directories and runtime logs."""

from __future__ import annotations

import os
from dataclasses import dataclass
from datetime import datetime
from typing import Mapping


@dataclass(frozen=True)
class RunOutputPaths:
    """Filesystem locations created for a single pipeline run."""

    run_output_path: str
    combined_dir: str
    protein_dir: str
    chain_dir: str
    distances_dir: str
    log_file: str


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
