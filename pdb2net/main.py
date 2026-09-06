"""Public entry points for PDB2Net batch and single-run execution."""

from __future__ import annotations

import sys
from typing import List, Union


def run_main(batch_files: List[str]) -> None:
    """Helper wrapper to re-enter main() for a batch."""
    try:
        main(batch_files)
    except Exception as exc:
        print(f"PDB2Net batch failed: {exc}", file=sys.stderr)
        raise SystemExit(1) from None


def main(input_path_or_filelist: Union[str, List[str]]) -> None:
    """Entry point for processing a folder or an explicit file list."""
    from .pipeline import run_pipeline

    run_pipeline(input_path_or_filelist)


def create_batches_streaming(file_paths, max_batch_kb: int):
    """Backward-compatible wrapper around the extracted batching helper."""
    from .batching import create_batches_streaming as create_batches

    return create_batches(file_paths, max_batch_kb)


def batch_run(input_folder: str, timeout_minutes: int = 10, size_limit_kb: int = 1000_100) -> bool:
    """Backward-compatible wrapper around the extracted batch runner."""
    from .batching import batch_run as run_batches

    return run_batches(
        input_folder,
        run_main,
        timeout_minutes=timeout_minutes,
        size_limit_kb=size_limit_kb,
    )


def legacy_command() -> int:
    """Run the config-driven legacy batch entry point with a reliable exit code."""
    try:
        from .config_loader import ConfigError, get_config

        config = get_config()
        if config.get("open_in_cytoscape", True):
            from .cytoscape import ensure_cytoscape_running

            ensure_cytoscape_running()

        return 0 if batch_run(config["input_folder_path"]) else 1
    except ConfigError as exc:
        print(f"PDB2Net configuration failed: {exc}", file=sys.stderr)
        return 1
    except Exception as exc:
        print(f"PDB2Net legacy batch failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(legacy_command())
