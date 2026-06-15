"""Report whether the local environment can run PDB2Net.

This script avoids importing PDB2Net modules so it can run even when project
dependencies are missing. It is intended as a quick pre-flight check before
running small fixtures or the full pipeline.
"""

from __future__ import annotations

import importlib.util
import json
import os
import platform
import shutil
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
CONFIG_DIR = ROOT / "pdb2net" / "configs"

PYTHON_IMPORTS = {
    "biopython": "Bio",
    "gemmi": "gemmi",
    "matplotlib": "matplotlib",
    "numpy": "numpy",
    "pandas": "pandas",
    "py4cytoscape": "py4cytoscape",
    "scipy": "scipy",
}

EXTERNAL_COMMANDS = ["blastp", "makeblastdb", "java", "cytoscape"]

PATH_KEYS = [
    "input_folder_path",
    "pdb_fasta_path",
    "uniprot_fasta_path",
    "sifts_tsv_path",
    "output_path",
    "blast_db_path",
]


def _read_json(path: Path) -> dict[str, Any]:
    try:
        with path.open("r", encoding="utf-8") as handle:
            return json.load(handle)
    except FileNotFoundError:
        return {}
    except json.JSONDecodeError as exc:
        print(f"[bad json] {path}: {exc}")
        return {}


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(base.get(key), dict):
            _deep_merge(base[key], value)
        else:
            base[key] = value
    return base


def _load_config_preview() -> dict[str, Any]:
    os_key = {
        "Windows": "windows",
        "Linux": "linux",
        "Darwin": "darwin",
    }.get(platform.system(), platform.system().lower())

    cfg: dict[str, Any] = {}
    for path in [
        CONFIG_DIR / "config.base.json",
        CONFIG_DIR / f"config.{os_key}.json",
        CONFIG_DIR / "config.local.json",
    ]:
        _deep_merge(cfg, _read_json(path))

    explicit_config = os.environ.get("PDB2NET_CONFIG_FILE")
    if explicit_config:
        _deep_merge(cfg, _read_json(Path(explicit_config)))

    env_overrides = {
        "PDB2NET_INPUT": "input_folder_path",
        "PDB2NET_OUTPUT": "output_path",
        "PDB2NET_PDB_FASTA": "pdb_fasta_path",
        "PDB2NET_UNIPROT_FASTA": "uniprot_fasta_path",
        "PDB2NET_SIFTS_TSV": "sifts_tsv_path",
        "PDB2NET_BLAST_DB": "blast_db_path",
        "PDB2NET_BLASTP": "blastp_executable",
    }
    for env_name, key in env_overrides.items():
        if os.environ.get(env_name):
            cfg[key] = os.environ[env_name]

    return cfg


def _expand_path(value: Any) -> Path | None:
    if not isinstance(value, str) or not value:
        return None
    return Path(os.path.expandvars(os.path.expanduser(value)))


def _status(ok: bool) -> str:
    return "ok" if ok else "missing"


def _path_exists(path: Path) -> tuple[bool, str | None]:
    try:
        return path.exists(), None
    except OSError as exc:
        return False, str(exc)


def main() -> int:
    print("PDB2Net environment check")
    print(f"Python: {sys.version.split()[0]} ({sys.executable})")
    print(f"Platform: {platform.system()} {platform.release()}")
    print()

    print("Python packages:")
    missing_python = []
    for package_name, import_name in PYTHON_IMPORTS.items():
        ok = importlib.util.find_spec(import_name) is not None
        if not ok:
            missing_python.append(package_name)
        print(f"  {package_name}: {_status(ok)}")
    print()

    print("External commands:")
    missing_external = []
    for command in EXTERNAL_COMMANDS:
        found = shutil.which(command)
        if not found and command not in {"cytoscape", "java"}:
            missing_external.append(command)
        detail = found if found else "missing"
        print(f"  {command}: {detail}")
    print()

    cfg = _load_config_preview()
    print("Configured paths:")
    missing_paths = []
    for key in PATH_KEYS:
        path = _expand_path(cfg.get(key))
        if path is None:
            missing_paths.append(key)
            print(f"  {key}: missing")
            continue
        exists, path_error = _path_exists(path)
        if not exists and key != "output_path":
            missing_paths.append(key)
        detail = _status(exists) if path_error is None else f"unavailable: {path_error}"
        print(f"  {key}: {path} ({detail})")

    blast_executable = cfg.get("blastp_executable", "blastp")
    print(f"  blastp_executable: {blast_executable}")
    print()

    if missing_python:
        print("Install missing Python packages with:")
        print("  python3 -m pip install -r requirements.txt")
    if missing_external:
        print("Install or configure missing BLAST commands before running BLAST-backed annotation.")
    if missing_paths:
        print("Create pdb2net/configs/config.local.json or set PDB2NET_* environment variables for local paths.")

    return 1 if missing_python else 0


if __name__ == "__main__":
    raise SystemExit(main())
