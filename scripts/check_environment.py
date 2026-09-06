"""Report whether the local environment can run PDB2Net.

This script imports only the standard-library configuration loader, so it can
still run when scientific project dependencies are missing. It is intended as
a quick pre-flight check before running small fixtures or the full pipeline.
"""

from __future__ import annotations

import importlib
import os
import platform
import shutil
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pdb2net.config_loader import ConfigError, load_config  # noqa: E402

PYTHON_IMPORTS = {
    "biopython": "Bio",
    "gemmi": "gemmi",
    "matplotlib": "matplotlib",
    "networkx": "networkx",
    "numpy": "numpy",
    "pandas": "pandas",
    "scipy": "scipy",
}
OPTIONAL_PYTHON_IMPORTS = {"py4cytoscape": "py4cytoscape"}

REQUIRED_REFERENCE_FILES = {
    "pdb_fasta_path": "PDB chain FASTA",
    "uniprot_fasta_path": "Swiss-Prot FASTA",
    "sifts_tsv_path": "SIFTS TSV",
}
STRUCTURE_SUFFIXES = (".pdb", ".ent", ".cif", ".mmcif")


def _expand_path(value: Any) -> Path | None:
    if not isinstance(value, str) or not value:
        return None
    return Path(os.path.expandvars(os.path.expanduser(value)))


def _status(ok: bool) -> str:
    return "ok" if ok else "missing"


def _find_executable(value: object) -> str | None:
    command = str(value or "").strip()
    if not command:
        return None
    return shutil.which(command)


def _import_available(import_name: str) -> tuple[bool, str | None]:
    try:
        importlib.import_module(import_name)
    except Exception as exc:
        return False, str(exc)
    return True, None


def _readable_file(path: Path | None) -> bool:
    if path is None or not path.is_file():
        return False
    try:
        with path.open("rb"):
            return True
    except OSError:
        return False


def _supported_structure_file(name: str) -> bool:
    filename = name.lower()
    if filename.endswith(".gz"):
        filename = filename[:-3]
    return filename.endswith(STRUCTURE_SUFFIXES)


def _input_directory_ready(path: Path | None) -> bool:
    if path is None or path.is_symlink() or not path.is_dir():
        return False
    try:
        with os.scandir(path) as entries:
            supported_found = False
            for entry in entries:
                if not _supported_structure_file(entry.name):
                    continue
                if entry.is_symlink():
                    return False
                if not entry.is_file(follow_symlinks=False):
                    continue
                if not _readable_file(Path(entry.path)):
                    return False
                supported_found = True
    except OSError:
        return False
    return supported_found


def _blast_database_files(cfg: dict[str, Any]) -> tuple[Path | None, list[Path]]:
    root = _expand_path(cfg.get("blast_db_path"))
    if root is None:
        return None, []
    prefix = root / "uniprot_db"
    files = [Path(str(prefix) + suffix) for suffix in (".pin", ".phr", ".psq")]
    return root, files


def _output_path_status(path: Path | None) -> tuple[bool, str]:
    if path is None:
        return False, "missing"
    access_mode = os.W_OK | (os.X_OK if os.name != "nt" else 0)
    try:
        if os.path.lexists(path) and not path.exists():
            return False, "dangling symbolic link"
        if path.exists():
            ok = path.is_dir() and os.access(path, access_mode)
            return ok, "ok" if ok else "not a writable directory"

        ancestor = path.parent
        while not ancestor.exists() and ancestor != ancestor.parent:
            if os.path.lexists(ancestor):
                return False, "parent is a dangling symbolic link"
            ancestor = ancestor.parent
        ok = ancestor.is_dir() and os.access(ancestor, access_mode)
        return ok, "not created yet" if ok else "parent is not writable"
    except OSError as exc:
        return False, f"unavailable: {exc}"


def _cytoscape_enabled(cfg: dict[str, Any]) -> bool:
    return bool(cfg["open_in_cytoscape"]) if "open_in_cytoscape" in cfg else True


def _cytoscape_running() -> bool:
    try:
        importlib.import_module("py4cytoscape").cytoscape_ping()
    except Exception:
        return False
    return True


def main() -> int:
    print("PDB2Net environment check")
    print(f"Python: {sys.version.split()[0]} ({sys.executable})")
    python_version_ok = sys.version_info >= (3, 11)
    print(f"Python >= 3.11: {_status(python_version_ok)}")
    print(f"Platform: {platform.system()} {platform.release()}")
    print()

    print("Python packages:")
    missing_python = []
    optional_python = {}
    for package_name, import_name in PYTHON_IMPORTS.items():
        ok, error = _import_available(import_name)
        if not ok:
            missing_python.append(package_name)
        detail = _status(ok) if error is None else f"unavailable: {error}"
        print(f"  {package_name}: {detail}")
    for package_name, import_name in OPTIONAL_PYTHON_IMPORTS.items():
        ok, error = _import_available(import_name)
        optional_python[package_name] = ok
        detail = _status(ok) if error is None else f"unavailable: {error}"
        print(f"  {package_name} (optional Cytoscape extra): {detail}")
    print()

    try:
        cfg = load_config()
    except ConfigError as exc:
        print(f"Configuration error: {exc}")
        return 1

    print("External commands:")
    blast_fallback_issues = []
    blast_setup_issues = []
    blastp = _find_executable(cfg.get("blastp_executable") or "blastp")
    if not blastp:
        blast_fallback_issues.append("blastp")
    print(f"  blastp (BLAST fallback): {blastp or 'missing'}")
    makeblastdb = _find_executable(
        cfg.get("makeblastdb_executable") or "makeblastdb"
    )
    if not makeblastdb:
        blast_setup_issues.append("makeblastdb")
    print(f"  makeblastdb (BLAST database setup): {makeblastdb or 'missing'}")
    cytoscape_path = cfg.get("cytoscape_path")
    cytoscape = _find_executable(cytoscape_path) if cytoscape_path else None
    print(f"  cytoscape (optional): {cytoscape or 'not configured'}")
    print()

    print("Configured paths:")
    missing_paths = []
    input_path = _expand_path(cfg.get("input_folder_path"))
    input_ok = _input_directory_ready(input_path)
    if not input_ok:
        missing_paths.append("input_folder_path")
    print(f"  input_folder_path: {input_path or 'missing'} ({_status(input_ok)})")

    for key, label in REQUIRED_REFERENCE_FILES.items():
        path = _expand_path(cfg.get(key))
        ok = _readable_file(path)
        if not ok:
            missing_paths.append(key)
        print(f"  {key}: {path or 'missing'} ({_status(ok)}; {label})")

    output_path = _expand_path(cfg.get("output_path"))
    output_ok, output_detail = _output_path_status(output_path)
    if not output_ok:
        missing_paths.append("output_path")
    print(f"  output_path: {output_path or 'missing'} ({output_detail})")

    blast_db_path, blast_db_files = _blast_database_files(cfg)
    blast_db_ok = bool(blast_db_files) and all(
        _readable_file(path) for path in blast_db_files
    )
    if not blast_db_ok:
        blast_fallback_issues.append("blast_db_path")
    print(f"  blast_db_path: {blast_db_path or 'missing'} ({_status(blast_db_ok)})")

    diamond_cfg = cfg.get("diamond", {}) if isinstance(cfg.get("diamond", {}), dict) else {}
    diamond_fallback_issues = []
    if diamond_cfg.get("enabled"):
        diamond_executable = str(diamond_cfg.get("executable") or "diamond")
        diamond_db = _expand_path(diamond_cfg.get("uniref90_db_path"))
        diamond_found = _find_executable(diamond_executable)
        diamond_db_ok = bool(diamond_db) and (
            _readable_file(diamond_db)
            or _readable_file(Path(str(diamond_db) + ".dmnd"))
        )
        if not diamond_found:
            diamond_fallback_issues.append("diamond")
        if not diamond_db_ok:
            diamond_fallback_issues.append("diamond.uniref90_db_path")
        print(f"  diamond_executable: {diamond_executable} ({diamond_found or 'missing'})")
        print(f"  diamond.uniref90_db_path: {diamond_db or 'missing'} ({_status(bool(diamond_db_ok))})")
    else:
        print("  diamond_uniref90: disabled")
    print()

    missing_active_optional = []
    if _cytoscape_enabled(cfg):
        if not optional_python.get("py4cytoscape", False):
            missing_active_optional.append("py4cytoscape")
        elif not cytoscape and not _cytoscape_running():
            missing_active_optional.append("running or executable Cytoscape")

    if missing_python:
        print("Install missing Python packages with:")
        print("  python3 -m pip install -r requirements.txt")
    if missing_paths:
        print("Create pdb2net/configs/config.local.json or set PDB2NET_* environment variables for local paths.")
    if blast_fallback_issues:
        print(
            "BLAST fallback is not ready; this is non-fatal for inputs that do not require fallback: "
            + ", ".join(blast_fallback_issues)
        )
    if blast_setup_issues:
        print(
            "BLAST database rebuild is unavailable; an existing complete database remains usable: "
            + ", ".join(blast_setup_issues)
        )
    if diamond_fallback_issues:
        print(
            "DIAMOND fallback is enabled but not ready; this is non-fatal for inputs that do not require fallback: "
            + ", ".join(diamond_fallback_issues)
        )
    if missing_active_optional:
        print(
            "Cytoscape live mode is enabled but unavailable: "
            + ", ".join(missing_active_optional)
        )
        if "py4cytoscape" in missing_active_optional:
            print("  python3 -m pip install '.[cytoscape]'")
        if "running or executable Cytoscape" in missing_active_optional:
            print("  Start Cytoscape or configure cytoscape_path.")

    return 1 if (
        not python_version_ok
        or missing_python
        or missing_paths
        or missing_active_optional
    ) else 0


if __name__ == "__main__":
    raise SystemExit(main())
