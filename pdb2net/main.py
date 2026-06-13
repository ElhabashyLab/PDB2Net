"""Public entry points for PDB2Net batch and single-run execution."""

from __future__ import annotations

import multiprocessing
from typing import List, Union

from .batching import batch_run as _batch_run
from .batching import create_batches_streaming as _create_batches_streaming
from .config_loader import config
from .pipeline import run_pipeline


def run_main(batch_files: List[str]) -> None:
    """Helper wrapper to re-enter main() for a batch."""
    main(batch_files)


def main(input_path_or_filelist: Union[str, List[str]]) -> None:
    """Entry point for processing a folder or an explicit file list."""
    run_pipeline(input_path_or_filelist)


def create_batches_streaming(file_paths, max_batch_kb: int):
    """Backward-compatible wrapper around the extracted batching helper."""
    return _create_batches_streaming(file_paths, max_batch_kb)


def batch_run(input_folder: str, timeout_minutes: int = 10, size_limit_kb: int = 1000_100) -> None:
    """Backward-compatible wrapper around the extracted batch runner."""
    _batch_run(
        input_folder,
        run_main,
        timeout_minutes=timeout_minutes,
        size_limit_kb=size_limit_kb,
        process_factory=multiprocessing.Process,
    )


if __name__ == "__main__":
    from .cytoscape import ensure_cytoscape_running

    if config.get("open_in_cytoscape", True):
        ensure_cytoscape_running()

    INPUT_FOLDER_PATH = config["input_folder_path"]
    batch_run(INPUT_FOLDER_PATH)
