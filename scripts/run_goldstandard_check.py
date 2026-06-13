"""Run a small headless regression check against a CX2 goldstandard.

The script can either run the configured PDB2Net pipeline and compare the newest
output folder, or compare an already existing output folder. It delegates the
semantic CX2 comparison to scripts/compare_cx2_outputs.py.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]


def _load_output_root() -> Path:
    from pdb2net.config_loader import config

    return Path(str(config["output_path"])).expanduser()


def _timestamped_dirs(output_root: Path) -> list[Path]:
    if not output_root.exists():
        return []
    dirs = [path for path in output_root.iterdir() if path.is_dir()]
    return sorted(dirs, key=lambda path: path.stat().st_mtime, reverse=True)


def _run_pipeline() -> None:
    subprocess.run([sys.executable, "-m", "pdb2net.main"], cwd=ROOT, check=True)


def _find_new_output(output_root: Path, started_at: float) -> Path:
    candidates = [
        path
        for path in _timestamped_dirs(output_root)
        if path.stat().st_mtime >= started_at - 1.0
    ]
    if not candidates:
        raise SystemExit(f"No new output directory found in {output_root}")
    return candidates[0]


def _compare(actual: Path, expected: Path, markdown_out: Path, json_out: Path, overlap_threshold: float) -> dict[str, Any]:
    markdown_out.parent.mkdir(parents=True, exist_ok=True)
    json_out.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "compare_cx2_outputs.py"),
            "--overlap-threshold",
            str(overlap_threshold),
            "compare",
            "--actual",
            str(actual),
            "--expected",
            str(expected),
            "--out",
            str(markdown_out),
            "--json-out",
            str(json_out),
        ],
        cwd=ROOT,
        check=True,
    )
    with json_out.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expected", required=True, type=Path, help="Goldstandard directory containing expected .cx2 files.")
    parser.add_argument("--actual", type=Path, help="Existing actual output directory. If omitted, run the pipeline first.")
    parser.add_argument("--output-root", type=Path, help="PDB2Net output root. Defaults to configured output_path.")
    parser.add_argument("--report-dir", type=Path, default=Path("/tmp/pdb2net-goldstandard"), help="Directory for reports.")
    parser.add_argument("--overlap-threshold", type=float, default=20.0)
    parser.add_argument("--fail-on-warn", action="store_true", help="Exit non-zero for WARN results as well as FAIL.")
    return parser


def main() -> int:
    args = _build_parser().parse_args()
    expected = args.expected.expanduser().resolve()
    if not expected.exists():
        raise SystemExit(f"Expected goldstandard directory does not exist: {expected}")

    if args.actual:
        actual = args.actual.expanduser().resolve()
    else:
        output_root = (args.output_root.expanduser() if args.output_root else _load_output_root()).resolve()
        started_at = datetime.now().timestamp()
        _run_pipeline()
        actual = _find_new_output(output_root, started_at)

    if not actual.exists():
        raise SystemExit(f"Actual output directory does not exist: {actual}")

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    markdown_out = args.report_dir / f"goldstandard_compare_{stamp}.md"
    json_out = args.report_dir / f"goldstandard_compare_{stamp}.json"
    report = _compare(actual, expected, markdown_out, json_out, args.overlap_threshold)

    status = str(report.get("status", "FAIL"))
    print(f"Goldstandard status: {status}")
    print(f"Actual: {actual}")
    print(f"Expected: {expected}")
    print(f"Markdown report: {markdown_out}")
    print(f"JSON report: {json_out}")

    if status == "FAIL" or (args.fail_on_warn and status == "WARN"):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
