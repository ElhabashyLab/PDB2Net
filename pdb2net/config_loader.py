"""Configuration loader for PDB2Net.

This module assembles configuration values in ordered layers:
1) base file, 2) OS-specific file, 3) local overrides,
4) an explicit file via environment variable, and 5) environment overrides.

It also performs light post-processing (path expansion, defaults) and
detects headless/container environments to disable GUI-dependent features
(e.g., opening Cytoscape) unless explicitly enabled.
"""
from __future__ import annotations
import json
import os
import platform
from pathlib import Path
from typing import Any, Dict, Tuple

SERVER_ENVIRONMENT = {
    "PDB2NET_PDB_FASTA": "pdb_fasta_path",
    "PDB2NET_SIFTS_TSV": "sifts_tsv_path",
    "PDB2NET_UNIPROT_FASTA": "uniprot_fasta_path",
    "PDB2NET_BLAST_DB": "blast_db_path",
    "PDB2NET_BLASTP": "blastp_executable",
    "PDB2NET_BLAST_CACHE_PATH": "blast_cache_path",
}

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
        **{
            name: (config_path, str)
            for name, config_path in SERVER_ENVIRONMENT.items()
        },
        "PDB2NET_CYTO_PATH": ("cytoscape_path", str),
        "PDB2NET_DIAMOND": ("diamond.executable", str),
        "PDB2NET_DIAMOND_UNIREF90_DB": ("diamond.uniref90_db_path", str),
        "PDB2NET_DIAMOND_TEMP_DIR": ("diamond.temp_dir", str),
        "PDB2NET_DIAMOND_THREADS": ("diamond.threads", int),
        "PDB2NET_DIAMOND_BLOCK_SIZE": ("diamond.block_size", float),
        "PDB2NET_DIAMOND_INDEX_CHUNKS": ("diamond.index_chunks", int),
        "PDB2NET_DIAMOND_MAX_TARGET_SEQS": ("diamond.max_target_seqs", int),
        "PDB2NET_DIAMOND_BATCH_MAX_SEQUENCES": ("diamond.batch_max_sequences", int),
        "PDB2NET_DIAMOND_BATCH_MAX_FASTA_BYTES": ("diamond.batch_max_fasta_bytes", int),
        "PDB2NET_DIAMOND_ITERATE": ("diamond.iterate", _bool_from_env),
        "PDB2NET_DIAMOND_SENSITIVITY": ("diamond.sensitivity", str),
        "PDB2NET_OPEN_IN_CYTOSCAPE": ("open_in_cytoscape", _bool_from_env),
        "PDB2NET_EXPORT_DETAILED_INTERACTIONS": ("export_detailed_interactions", _bool_from_env),
        "PDB2NET_DIAMOND_ENABLED": ("diamond.enabled", _bool_from_env),
        "PDB2NET_STRUCTURE_MODEL_POLICY": ("structure_model_policy", str),
        "PDB2NET_REFERENCE_MANIFEST_ID": ("reference_manifest_id", str),
        "PDB2NET_USE_EMBEDDED_SIFTS": ("network_annotations.use_embedded_sifts", _bool_from_env),
        "PDB2NET_ANNOTATION_TOOLTIP_FIELDS": (
            "network_annotations.tooltip_fields",
            lambda value: [field.strip() for field in value.split(",") if field.strip()],
        ),
        "PDB2NET_MAX_TOOLTIP_SEGMENTS_PER_DATABASE": (
            "network_annotations.max_tooltip_segments_per_database",
            int,
        ),
    }
    nested: Dict[str, Tuple[str, str]] = {
        "PDB2NET_WORKERS_PARSING": ("workers", "parsing"),
        "PDB2NET_WORKERS_BLAST": ("workers", "blast_threads"),
        "PDB2NET_CA_RADIUS": ("distance_thresholds", "ca_radius"),
        "PDB2NET_ALL_ATOMS_RADIUS": ("distance_thresholds", "all_atoms_radius"),
        "PDB2NET_PP_MIN_CA_NEIGHBORS": ("interaction_filters", "protein_protein_min_ca_neighbors"),
        "PDB2NET_PP_MIN_ALL_ATOM_CONTACTS": ("interaction_filters", "protein_protein_min_all_atom_contacts"),
        "PDB2NET_PNA_MIN_ALL_ATOM_CONTACTS": ("interaction_filters", "protein_nucleic_acid_min_all_atom_contacts"),
        "PDB2NET_NA_MIN_ALL_ATOM_CONTACTS": ("interaction_filters", "nucleic_acid_min_all_atom_contacts"),
        "PDB2NET_MAX_INPUT_FILES": ("resource_limits", "max_input_files"),
        "PDB2NET_MAX_TOTAL_INPUT_BYTES": ("resource_limits", "max_total_input_bytes"),
        "PDB2NET_MAX_SINGLE_INPUT_BYTES": ("resource_limits", "max_single_input_bytes"),
        "PDB2NET_MAX_PROCESSING_BATCH_BYTES": ("resource_limits", "max_processing_batch_bytes"),
        "PDB2NET_MAX_TOTAL_INPUT_EXPANDED_BYTES": ("resource_limits", "max_total_input_expanded_bytes"),
        "PDB2NET_MAX_SINGLE_INPUT_EXPANDED_BYTES": ("resource_limits", "max_single_input_expanded_bytes"),
        "PDB2NET_MAX_DETAILED_INTERACTION_ROWS": (
            "resource_limits",
            "max_detailed_interaction_rows",
        ),
        "PDB2NET_MAX_DETAILED_INTERACTION_BYTES": (
            "resource_limits",
            "max_detailed_interaction_bytes",
        ),
        "PDB2NET_MIN_OUTPUT_FREE_BYTES": (
            "resource_limits",
            "min_output_free_bytes",
        ),
        "PDB2NET_COMBINED_MAX_NODES": ("combined_graph_limits", "max_nodes"),
        "PDB2NET_COMBINED_MAX_EDGES": ("combined_graph_limits", "max_edges"),
    }

    for env, (key, caster) in flat.items():
        raw = os.environ.get(env)
        if raw:
            value = caster(raw) if callable(caster) else raw
            if "." in key:
                first, second = key.split(".", 1)
                cfg.setdefault(first, {})
                cfg[first][second] = value
            else:
                cfg[key] = value

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
    """Normalize paths, set headless/container defaults, and provide fallbacks."""
    for key in [
        "input_folder_path",
        "pdb_fasta_path",
        "uniprot_fasta_path",
        "sifts_tsv_path",
        "output_path",
        "cytoscape_path",
        "blast_db_path",
        "blast_cache_path",
        "blastp_executable",
        "diamond.uniref90_db_path",
        "diamond.executable",
        "diamond.temp_dir",
        "layout_engine_path",
    ]:
        if "." in key:
            first, second = key.split(".", 1)
            if isinstance(cfg.get(first), dict) and cfg[first].get(second):
                cfg[first][second] = _normalize_path(cfg[first][second])
        elif key in cfg and cfg[key]:
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
    cfg.setdefault("resource_limits", {})
    cfg["resource_limits"].setdefault("max_input_files", None)
    cfg["resource_limits"].setdefault("max_total_input_bytes", None)
    cfg["resource_limits"].setdefault("max_single_input_bytes", None)
    cfg["resource_limits"].setdefault("max_processing_batch_bytes", None)
    cfg["resource_limits"].setdefault("max_total_input_expanded_bytes", None)
    cfg["resource_limits"].setdefault("max_single_input_expanded_bytes", None)
    cfg["resource_limits"].setdefault("max_detailed_interaction_rows", None)
    cfg["resource_limits"].setdefault("max_detailed_interaction_bytes", None)
    cfg["resource_limits"].setdefault("min_output_free_bytes", None)
    cfg.setdefault("network_annotations", {})
    cfg["network_annotations"].setdefault("use_embedded_sifts", True)
    cfg["network_annotations"].setdefault("tooltip_fields", ["uniprot"])
    cfg["network_annotations"].setdefault("max_tooltip_segments_per_database", 20)
    cfg.setdefault("combined_graph_limits", {})
    cfg["combined_graph_limits"].setdefault("max_nodes", None)
    cfg["combined_graph_limits"].setdefault("max_edges", None)
    cfg.setdefault("reference_manifest_id", "")
    cfg.setdefault("layout_mode", "python_fast")
    cfg.setdefault("layout_engine_path", "")
    cfg.setdefault("layout_keep_temp_files", False)
    cfg.setdefault("structure_model_policy", "first")
    cfg.setdefault("diamond", {})
    cfg["diamond"].setdefault("enabled", False)
    cfg["diamond"].setdefault("executable", "diamond")
    cfg["diamond"].setdefault("uniref90_db_path", "")
    cfg["diamond"].setdefault("temp_dir", "")
    cfg["diamond"].setdefault("threads", 6)
    cfg["diamond"].setdefault("iterate", True)
    cfg["diamond"].setdefault("sensitivity", "sensitive")
    cfg["diamond"].setdefault("block_size", 1.0)
    cfg["diamond"].setdefault("index_chunks", 4)
    cfg["diamond"].setdefault("max_target_seqs", 50)
    cfg["diamond"].setdefault("batch_max_sequences", 5000)
    cfg["diamond"].setdefault("batch_max_fasta_bytes", 50 * 1024 * 1024)
    cfg["diamond"].setdefault("assign_uniprot_id", "never")
    if not isinstance(cfg.get("interaction_filters"), dict):
        cfg["interaction_filters"] = {}
    cfg["interaction_filters"].setdefault("protein_protein_min_ca_neighbors", 10)
    cfg["interaction_filters"].setdefault("protein_protein_min_all_atom_contacts", 1)
    cfg["interaction_filters"].setdefault("protein_nucleic_acid_min_all_atom_contacts", 1)
    cfg["interaction_filters"].setdefault("nucleic_acid_min_all_atom_contacts", 1)


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

__all__ = ["SERVER_ENVIRONMENT", "config", "load_config"]
