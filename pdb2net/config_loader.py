"""Configuration loader for PDB2Net.

This module assembles configuration values in ordered layers:
1) base file, 2) OS-specific file, 3) local overrides,
4) an explicit file via environment variable, and 5) environment overrides.

It also performs light post-processing (path expansion, defaults) and
detects headless/container environments to disable GUI-dependent features
(e.g., opening Cytoscape) unless explicitly enabled.
"""
from __future__ import annotations
import json, os, platform
from pathlib import Path
from typing import Any, Dict, Tuple

# === Minimal logging switch (default: quiet) ===
VERBOSE = os.environ.get("PDB2NET_VERBOSE", "").strip().lower() in {"1", "true", "yes", "on"}

# Global cache
_config_cache: Dict[str, Any] | None = None


def _log(msg: str) -> None:
    if VERBOSE:
        print(f"[config] {msg}")


def _read_json(p: Path, strict: bool = False) -> Dict[str, Any]:
    """Read a JSON file. If strict=True, a JSON error terminates the program."""
    try:
        with p.open("r", encoding="utf-8") as f:
            return json.load(f)
    except FileNotFoundError:
        if VERBOSE:
            _log(f"Info: {p} nicht gefunden.")
        return {}
    except json.JSONDecodeError as e:
        msg = f"Fehler in {p}: Ungültiges JSON ({e})."
        if strict:
            raise SystemExit(msg)
        _log("Warn: " + msg)
        return {}
    except Exception as e:
        if strict:
            raise SystemExit(f"Konnte {p} nicht lesen: {e}")
        _log(f"Warn: Konnte {p} nicht lesen: {e}")
        return {}


def _deep_merge(base: Dict[str, Any], override: Dict[str, Any]) -> Dict[str, Any]:
    """Recursive deep-merge strategy: values from `override` win; nested dicts are merged in place."""
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            _deep_merge(base[k], v)  # type: ignore[index]
        else:
            base[k] = v
    return base


def _bool_from_env(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "on", "y"}


def _normalize_path(value: Any) -> Any:
    if isinstance(value, str):
        # Expand ~ and $VARS; no OS-specific path conversion required
        return os.path.expanduser(os.path.expandvars(value))
    return value


def _is_headless_linux() -> bool:
    return platform.system() == "Linux" and not os.environ.get("DISPLAY")


def _is_container() -> bool:
    # simple heuristic for Docker/K8s
    if Path("/.dockerenv").exists():
        return True
    try:
        cgroup = Path("/proc/1/cgroup")
        if cgroup.exists():
            txt = cgroup.read_text(errors="ignore")
            return ("docker" in txt) or ("kubepods" in txt)
    except Exception:
        pass
    return False


def _apply_env_overrides(cfg: Dict[str, Any]) -> None:
    """Apply environment variable overrides to the in-memory config (flat and nested keys)."""
    flat: Dict[str, Tuple[str, Any]] = {
        "PDB2NET_INPUT": ("input_folder_path", str),
        "PDB2NET_OUTPUT": ("output_path", str),
        "PDB2NET_PDB_FASTA": ("pdb_fasta_path", str),
        "PDB2NET_UNIPROT_FASTA": ("uniprot_fasta_path", str),
        "PDB2NET_SIFTS_TSV": ("sifts_tsv_path", str),
        "PDB2NET_CYTO_PATH": ("cytoscape_path", str),
        "PDB2NET_BLAST_DB": ("blast_db_path", str),
        "PDB2NET_BLASTP": ("blastp_executable", str),
        "PDB2NET_OPEN_IN_CYTOSCAPE": ("open_in_cytoscape", _bool_from_env),
    }
    nested: Dict[str, Tuple[str, str]] = {
        "PDB2NET_WORKERS_PARSING": ("workers", "parsing"),
        "PDB2NET_WORKERS_BLAST": ("workers", "blast_threads"),
        "PDB2NET_CA_RADIUS": ("distance_thresholds", "ca_radius"),
        "PDB2NET_ALL_ATOMS_RADIUS": ("distance_thresholds", "all_atoms_radius"),
    }

    for env, (key, caster) in flat.items():
        raw = os.environ.get(env)
        if raw:
            cfg[key] = caster(raw) if callable(caster) else raw

    for env, (k1, k2) in nested.items():
        raw = os.environ.get(env)
        if not raw:
            continue
        cfg.setdefault(k1, {})
        # Allow numbers or strings (e.g., "auto")
        try:
            if raw.replace(".", "", 1).isdigit():
                cfg[k1][k2] = float(raw) if "." in raw else int(raw)
            else:
                cfg[k1][k2] = raw
        except Exception:
            cfg[k1][k2] = raw


def _postprocess(cfg: Dict[str, Any], os_key: str) -> None:
    """Normalize paths, set headless/container defaults, provide fallbacks, and ensure output directory."""
    for key in [
        "input_folder_path",
        "pdb_fasta_path",
        "uniprot_fasta_path",
        "sifts_tsv_path",
        "output_path",
        "cytoscape_path",
        "blast_db_path",
        "blastp_executable",
    ]:
        if key in cfg and cfg[key]:
            cfg[key] = _normalize_path(cfg[key])

    # Headless/Container => disable Cytoscape by default (unless explicitly set)
    if "open_in_cytoscape" not in cfg or cfg["open_in_cytoscape"] is None:
        if _is_headless_linux() or _is_container():
            cfg["open_in_cytoscape"] = False

    # BLAST fallback: Windows path vs. plain 'blastp' on Unix
    if not cfg.get("blastp_executable"):
        cfg["blastp_executable"] = (
            r"C:\Program Files\NCBI\blast-2.17.0+\bin\blastp.exe"
            if os_key == "windows"
            else "blastp"
        )

    # Keep worker defaults; "auto" is resolved later in the code
    cfg.setdefault("workers", {})
    cfg["workers"].setdefault("parsing", "auto")
    cfg["workers"].setdefault("blast_threads", "auto")

    # Create output folder proactively (if provided)
    out = cfg.get("output_path")
    if isinstance(out, str) and out:
        try:
            Path(out).mkdir(parents=True, exist_ok=True)
        except Exception as e:
            _log(f"Warn: Konnte output_path '{out}' nicht anlegen: {e}")


def load_config() -> Dict[str, Any]:
    """Load configuration in layers: base → OS-specific → local → explicit file → environment variables."""
    global _config_cache
    if _config_cache is not None:
        return _config_cache

    root = Path(__file__).resolve().parent
    os_key = {"Windows": "windows", "Linux": "linux", "Darwin": "darwin"}.get(
        platform.system(), platform.system().lower()
    )

    # Configuration directory (default: ./configs)
    cfg_dir = Path(os.environ.get("PDB2NET_CONFIG_DIR") or (root / "configs"))

    candidates: list[Tuple[Path, bool]] = [
        (cfg_dir / "config.base.json", False),
        (cfg_dir / f"config.{os_key}.json", False),
        (cfg_dir / "config.local.json", False),  # git-ignored, highest file-priority within configs/
    ]

    # Legacy fallback (backward compatibility) if nothing in configs/ applies
    legacy = root / "config.json"
    if legacy.exists():
        candidates.append((legacy, False))

    # Explicit file override via ENV
    env_cfg = os.environ.get("PDB2NET_CONFIG_FILE")
    if env_cfg:
        candidates.append((Path(env_cfg), True))

    # Merge in order
    cfg: Dict[str, Any] = {}
    for path, strict in candidates:
        part = _read_json(path, strict=strict)
        _deep_merge(cfg, part)

    # Post-processing and environment overrides
    _postprocess(cfg, os_key)
    _apply_env_overrides(cfg)

    if VERBOSE:
        _log(f"OS={os_key} | headless={_is_headless_linux()} | container={_is_container()}")
        safe = dict(cfg)
        if "cytoscape_path" in safe and isinstance(safe["cytoscape_path"], str) and safe["cytoscape_path"]:
            safe["cytoscape_path"] = "<set>"
        _log("Aktive Schlüssel: " + ", ".join(sorted(safe.keys())))

    _config_cache = cfg
    return cfg


# Load on import (lazy + cache)
config = load_config()

__all__ = ["config", "load_config"]
