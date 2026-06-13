"""Smoke-test that failed batch subprocesses are not counted as completed."""

from __future__ import annotations

import tempfile
from pathlib import Path
from unittest.mock import patch

from pdb2net import main as pipeline


class FailingProcess:
    exitcode = 1

    def __init__(self, target=None, args=()):
        self.target = target
        self.args = args

    def start(self) -> None:
        return None

    def join(self, timeout=None) -> None:
        return None

    def is_alive(self) -> bool:
        return False


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="pdb2net-batch-exitcode-") as temp_dir:
        root = Path(temp_dir)
        input_dir = root / "inputs"
        output_dir = root / "outputs"
        input_dir.mkdir()
        output_dir.mkdir()
        fixture = input_dir / "6m17.cif"
        fixture.write_text("data_6m17\n", encoding="utf-8")

        original_output = pipeline.config.get("output_path")
        original_open = pipeline.config.get("open_in_cytoscape")
        try:
            pipeline.config["output_path"] = str(output_dir)
            pipeline.config["open_in_cytoscape"] = False
            with patch.object(pipeline.multiprocessing, "Process", FailingProcess):
                pipeline.batch_run(str(input_dir), timeout_minutes=1, size_limit_kb=1000)
        finally:
            pipeline.config["output_path"] = original_output
            pipeline.config["open_in_cytoscape"] = original_open

        log_path = output_dir / "error_in_batch_log" / "skipped_batches.txt"
        log_text = log_path.read_text(encoding="utf-8")
        if "exitcode=1" not in log_text:
            raise AssertionError(f"expected exitcode log entry, got: {log_text!r}")

    print("batch_exitcode_failure: ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
