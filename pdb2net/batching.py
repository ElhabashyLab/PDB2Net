"""Helpers for streaming batch execution with timeout handling."""

from __future__ import annotations

import multiprocessing
import os
from pathlib import Path
import time
from typing import Callable, Iterable, Iterator, List, Protocol

from .config_loader import config
from .file_parser import get_pdb_id, is_valid_file
from .input_contract import InputValidationError
from .logging_utils import get_logger
from .outputs import write_failed_run_manifest


logger = get_logger(__name__)


class BatchProcess(Protocol):
    """Minimal process protocol needed by the batch runner."""

    exitcode: int | None

    def start(self) -> None: ...

    def join(self, timeout: float | None = None) -> None: ...

    def is_alive(self) -> bool: ...

    def terminate(self) -> None: ...


def stream_valid_files(folder: str) -> Iterator[str]:
    """Yield valid input files from a folder."""
    for entry in os.scandir(folder):
        if entry.is_file() and is_valid_file(entry.path):
            yield entry.path


def validate_input_folder(folder: str) -> None:
    """Validate the productive folder-based input contract before batching."""
    input_path = Path(folder)
    if not input_path.exists():
        raise InputValidationError(
            "INPUT_PATH_NOT_FOUND",
            f"input_folder_path does not exist: {input_path}",
        )
    if not input_path.is_dir():
        raise InputValidationError(
            "INPUT_PATH_NOT_DIRECTORY",
            f"input_folder_path must be a directory containing PDB/mmCIF files: {input_path}",
        )
    if not any(entry.is_file() and is_valid_file(entry.path) for entry in os.scandir(str(input_path))):
        raise InputValidationError(
            "NO_VALID_INPUT_FILES",
            f"input_folder_path contains no supported .pdb, .cif, or .mmcif files: {input_path}",
        )


def append_skipped_batch_log(log_path: str, batch_files: List[str], reason: str) -> None:
    """Append skipped batch details to the shared log file."""
    with open(log_path, "a", encoding="utf-8") as log_file:
        for file_path in batch_files:
            pdb_id = get_pdb_id(file_path)
            log_file.write(f"{pdb_id}\t{reason}\n")


def create_batches_streaming(file_paths: Iterable[str], max_batch_kb: int) -> Iterator[List[str]]:
    """Yield batches of files with a soft size cap, streaming over the input."""
    current_batch: List[str] = []
    current_size = 0

    for file_path in file_paths:
        try:
            size_kb = os.path.getsize(file_path) // 1024
        except Exception as exc:
            print(f"Error reading file size: {file_path} – {exc}")
            continue

        if size_kb > max_batch_kb:
            print(f"File too large for a single batch: {file_path} ({size_kb} KB)")
            continue

        if current_size + size_kb > max_batch_kb:
            yield current_batch
            current_batch = [file_path]
            current_size = size_kb
        else:
            current_batch.append(file_path)
            current_size += size_kb

    if current_batch:
        yield current_batch


def batch_run(
    input_folder: str,
    run_batch: Callable[[List[str]], None],
    *,
    timeout_minutes: int = 10,
    size_limit_kb: int = 1000_100,
    process_factory: Callable[..., BatchProcess] | None = None,
) -> None:
    """Run the pipeline repeatedly on streamed batches from a folder."""
    try:
        validate_input_folder(input_folder)
    except InputValidationError as exc:
        logger.error("Batch input validation failed: %s", exc)
        write_failed_run_manifest(
            config["output_path"],
            input_path=input_folder,
            config_snapshot={
                "workers": config.get("workers", {}),
                "open_in_cytoscape": config.get("open_in_cytoscape"),
                "export_detailed_interactions": config.get("export_detailed_interactions", False),
            },
            error=exc,
        )
        raise
    process_factory = process_factory or multiprocessing.Process

    timeout_seconds = timeout_minutes * 60
    logs_dir = os.path.join(config["output_path"], "error_in_batch_log")
    os.makedirs(logs_dir, exist_ok=True)
    log_path = os.path.join(logs_dir, "skipped_batches.txt")

    start_time_all = time.time()
    total_done = 0
    file_stream = stream_valid_files(input_folder)

    for index, batch_files in enumerate(create_batches_streaming(file_stream, size_limit_kb), start=1):
        logger.info("Processing batch %s (%s files)", index, len(batch_files))

        if config.get("open_in_cytoscape", False):
            start_time_batch = time.time()
            run_batch(batch_files)
            duration = time.time() - start_time_batch
            total_done += len(batch_files)
            logger.info("Batch %s finished in %.1f seconds", index, duration)
            logger.info("Processed so far: %s files", total_done)
            continue

        start_time_batch = time.time()
        process = process_factory(target=run_batch, args=(batch_files,))
        process.start()
        process.join(timeout=timeout_seconds)

        if process.is_alive():
            logger.error("Batch %s exceeded %s minutes. Terminating.", index, timeout_minutes)
            process.terminate()
            process.join()
            append_skipped_batch_log(log_path, batch_files, f"timeout>{timeout_minutes}min")
            continue

        duration = time.time() - start_time_batch
        if process.exitcode != 0:
            logger.error("Batch %s failed after %.1f seconds (exit code %s)", index, duration, process.exitcode)
            append_skipped_batch_log(log_path, batch_files, f"exitcode={process.exitcode}")
            continue

        total_done += len(batch_files)
        logger.info("Batch %s finished in %.1f seconds", index, duration)
        logger.info("Processed so far: %s files", total_done)

    total_time_all = time.time() - start_time_all
    logger.info("Total time for all batches: %.2f seconds", total_time_all)
