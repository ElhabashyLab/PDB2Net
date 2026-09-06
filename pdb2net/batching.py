"""Helpers for streaming batch execution with timeout handling."""

from __future__ import annotations

import multiprocessing
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
from typing import Callable, Iterable, Iterator, List, Protocol

from .config_loader import config
from .file_parser import FileSignature, get_pdb_id, input_file_signature, is_valid_file
from .input_contract import InputValidationError
from .logging_utils import get_logger
from .outputs import write_failed_run_manifest


logger = get_logger(__name__)
TERMINATION_GRACE_SECONDS = 5


class BatchProcess(Protocol):
    """Minimal process protocol needed by the batch runner."""

    exitcode: int | None
    pid: int | None

    def start(self) -> None: ...

    def join(self, timeout: float | None = None) -> None: ...

    def is_alive(self) -> bool: ...

    def terminate(self) -> None: ...

    def kill(self) -> None: ...


def _validate_batch_signatures(
    batch_files: List[str], expected_signatures: dict[str, FileSignature]
) -> None:
    """Reject inputs that changed after the parent inventory was captured."""
    for file_path in batch_files:
        if input_file_signature(file_path) != expected_signatures[file_path]:
            raise InputValidationError(
                "INPUT_CHANGED_DURING_PROCESSING",
                f"Input file changed before batch processing: {Path(file_path).name}",
            )


def _run_batch_process(
    run_batch: Callable[[List[str]], None],
    batch_files: List[str],
    expected_signatures: dict[str, FileSignature],
) -> None:
    """Run one batch in an isolated process group after checking its inventory."""
    try:
        if os.name == "posix":
            os.setsid()
        _validate_batch_signatures(batch_files, expected_signatures)
        run_batch(batch_files)
    except Exception as exc:
        print(f"PDB2Net batch failed: {exc}", file=sys.stderr)
        raise SystemExit(1) from None


def _stop_batch_process(process: BatchProcess, *, terminate_tree: bool) -> bool:
    """Stop a timed-out batch without an unbounded join."""
    if terminate_tree and process.pid is not None and os.name == "posix":
        group_terminated = True
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            group_terminated = False
        if not group_terminated and process.is_alive():
            process.terminate()
        process.join(timeout=TERMINATION_GRACE_SECONDS)
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        if process.is_alive():
            process.kill()
        process.join(timeout=TERMINATION_GRACE_SECONDS)
        return not process.is_alive()

    if terminate_tree and process.pid is not None and os.name == "nt":
        tree_stopped = False
        try:
            result = subprocess.run(
                ["taskkill", "/PID", str(process.pid), "/T", "/F"],
                check=False,
                capture_output=True,
                timeout=TERMINATION_GRACE_SECONDS,
            )
            tree_stopped = result.returncode == 0
        except (OSError, subprocess.TimeoutExpired):
            pass
        if process.is_alive():
            process.kill()
        process.join(timeout=TERMINATION_GRACE_SECONDS)
        if process.is_alive():
            process.kill()
            process.join(timeout=TERMINATION_GRACE_SECONDS)
        return tree_stopped and not process.is_alive()

    process.terminate()
    process.join(timeout=TERMINATION_GRACE_SECONDS)
    if process.is_alive():
        process.kill()
        process.join(timeout=TERMINATION_GRACE_SECONDS)
    return not process.is_alive()


def stream_valid_files(folder: str) -> Iterator[str]:
    """Yield valid input files from a folder."""
    for entry in os.scandir(folder):
        if entry.is_file() and is_valid_file(entry.path):
            yield entry.path


def validate_input_folder(folder: str) -> list[str]:
    """Validate and inventory the folder-based input contract before batching."""
    input_path = Path(folder)
    if input_path.is_symlink():
        raise InputValidationError(
            "SYMLINK_INPUT_NOT_ALLOWED",
            f"Input directories must not be symbolic links: {input_path}",
        )
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
    file_paths: list[str] = []
    with os.scandir(str(input_path)) as entries:
        for entry in entries:
            if entry.is_symlink() and is_valid_file(entry.path):
                raise InputValidationError(
                    "SYMLINK_INPUT_NOT_ALLOWED",
                    f"Structure inputs must not be symbolic links: {entry.name}",
                )
            if entry.is_file(follow_symlinks=False) and is_valid_file(entry.path):
                file_paths.append(entry.path)
    if not file_paths:
        raise InputValidationError(
            "NO_VALID_INPUT_FILES",
            f"input_folder_path contains no supported .pdb, .cif, or .mmcif files: {input_path}",
        )
    return sorted(file_paths)


def append_skipped_batch_log(log_path: str, batch_files: List[str], reason: str) -> None:
    """Append skipped batch details to the shared log file."""
    with open(log_path, "a", encoding="utf-8") as log_file:
        for file_path in batch_files:
            try:
                pdb_id = get_pdb_id(file_path)
            except Exception:
                pdb_id = Path(file_path).name
            log_file.write(f"{pdb_id}\t{reason}\n")


def create_batches_streaming(
    file_paths: Iterable[str],
    max_batch_kb: int,
    *,
    on_skipped: Callable[[str, str], None] | None = None,
) -> Iterator[List[str]]:
    """Yield batches of files with a soft size cap, streaming over the input."""
    current_batch: List[str] = []
    current_size = 0

    for file_path in file_paths:
        try:
            size_kb = os.path.getsize(file_path) // 1024
        except Exception as exc:
            print(f"Error reading file size: {file_path} – {exc}")
            if on_skipped is not None:
                on_skipped(file_path, f"size_error={exc.__class__.__name__}")
            continue

        if size_kb > max_batch_kb:
            print(f"File too large for a single batch: {file_path} ({size_kb} KB)")
            if on_skipped is not None:
                on_skipped(file_path, f"size_kb={size_kb}>{max_batch_kb}")
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
) -> bool:
    """Run the pipeline repeatedly on streamed batches from a folder."""
    try:
        input_files = validate_input_folder(input_folder)
        input_signatures = {
            file_path: input_file_signature(file_path) for file_path in input_files
        }
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
    terminate_tree = process_factory is None
    process_factory = process_factory or multiprocessing.Process

    timeout_seconds = timeout_minutes * 60
    logs_dir = os.path.join(config["output_path"], "error_in_batch_log")
    os.makedirs(logs_dir, exist_ok=True)
    log_path = os.path.join(logs_dir, "skipped_batches.txt")

    start_time_all = time.time()
    total_done = 0
    batches_seen = 0
    failed_batches = 0
    timed_out_batches = 0
    skipped_files = 0
    file_stream = iter(input_files)

    def record_skipped_file(file_path: str, reason: str) -> None:
        nonlocal skipped_files
        skipped_files += 1
        append_skipped_batch_log(log_path, [file_path], reason)

    batches = create_batches_streaming(
        file_stream,
        size_limit_kb,
        on_skipped=record_skipped_file,
    )
    for index, batch_files in enumerate(batches, start=1):
        batches_seen += 1
        logger.info("Processing batch %s (%s files)", index, len(batch_files))

        start_time_batch = time.time()
        process = process_factory(
            target=_run_batch_process,
            args=(run_batch, batch_files, input_signatures),
        )
        process.start()
        process.join(timeout=timeout_seconds)

        if process.is_alive():
            logger.error("Batch %s exceeded %s minutes. Terminating.", index, timeout_minutes)
            if not _stop_batch_process(process, terminate_tree=terminate_tree):
                logger.error("Batch %s remained alive after forced termination.", index)
            append_skipped_batch_log(log_path, batch_files, f"timeout>{timeout_minutes}min")
            timed_out_batches += 1
            continue

        duration = time.time() - start_time_batch
        if process.exitcode != 0:
            logger.error("Batch %s failed after %.1f seconds (exit code %s)", index, duration, process.exitcode)
            append_skipped_batch_log(log_path, batch_files, f"exitcode={process.exitcode}")
            failed_batches += 1
            continue

        total_done += len(batch_files)
        logger.info("Batch %s finished in %.1f seconds", index, duration)
        logger.info("Processed so far: %s files", total_done)

    total_time_all = time.time() - start_time_all
    logger.info("Total time for all batches: %.2f seconds", total_time_all)
    complete = (
        batches_seen > 0
        and failed_batches == 0
        and timed_out_batches == 0
        and skipped_files == 0
    )
    if not complete:
        if batches_seen == 0 and skipped_files == 0:
            logger.error("Batch run incomplete: no batches remained after input validation.")
        logger.error(
            "Batch run incomplete: %s failed, %s timed out, %s input file(s) skipped.",
            failed_batches,
            timed_out_batches,
            skipped_files,
        )
    return complete
