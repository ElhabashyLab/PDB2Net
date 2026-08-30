"""Public entry points for PDB2Net batch and single-run execution."""

from __future__ import annotations

import multiprocessing
import sys
from typing import List, Union


def run_main(batch_files: List[str]) -> None:
    """Helper wrapper to re-enter main() for a batch."""
    main(batch_files)


def main(input_path_or_filelist: Union[str, List[str]]) -> None:
    """Entry point for processing a folder or an explicit file list."""
    from .pipeline import run_pipeline

    run_pipeline(input_path_or_filelist)


def create_batches_streaming(file_paths, max_batch_kb: int):
    """Backward-compatible wrapper around the extracted batching helper."""
    from .batching import create_batches_streaming as create_batches

    return create_batches(file_paths, max_batch_kb)


def batch_run(input_folder: str, timeout_minutes: int = 10, size_limit_kb: int = 1000_100) -> None:
    """Backward-compatible wrapper around the extracted batch runner."""
    from .batching import batch_run as run_batches

    run_batches(
        input_folder,
        run_main,
        timeout_minutes=timeout_minutes,
        size_limit_kb=size_limit_kb,
        process_factory=multiprocessing.Process,
    )


if __name__ == "__main__":
    try:
        from .config_loader import ConfigError, get_config

        config = get_config()
        if config.get("open_in_cytoscape", True):
            from .cytoscape import ensure_cytoscape_running

            ensure_cytoscape_running()

        batch_run(config["input_folder_path"])
    except ConfigError as exc:
        print(f"PDB2Net configuration failed: {exc}", file=sys.stderr)
        raise SystemExit(1) from None
