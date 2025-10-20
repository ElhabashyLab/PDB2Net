# config_loader.py
# OS-sensitiver, mehrstufiger Config-Loader mit ENV-Overrides und Headless-Erkennung.
from __future__ import annotations
import json, os, platform
from pathlib import Path
from typing import Any, Dict, Tuple

# === Minimaler Log-Schalter (Standard: leise) ===
VERBOSE = os.environ.get("PDB2NET_VERBOSE", "").strip().lower() in {"1", "true", "yes", "on"}

# Globaler Cache
_config_cache: Dict[str, Any] | None = None


def _log(msg: str) -> None:
    if VERBOSE:
        print(f"[config] {msg}")


def _read_json(p: Path, strict: bool = False) -> Dict[str, Any]:
    """Liest JSON-Datei. Bei strict=True führt JSON-Fehler zum Abbruch."""
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
    """rekursive Merge-Strategie: override gewinnt, Dictionaries werden tief gemerged"""
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
        # ~ und $VARS expandieren; keine OS-spezifische Pfad-Konvertierung nötig
        return os.path.expanduser(os.path.expandvars(value))
    return value


def _is_headless_linux() -> bool:
    return platform.system() == "Linux" and not os.environ.get("DISPLAY")


def _is_container() -> bool:
    # einfache Heuristik für Docker/K8s
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
    """ENV → Config (flache und verschachtelte Schlüssel)"""
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
        # Zahl vs. String erlauben (z. B. "auto")
        try:
            if raw.replace(".", "", 1).isdigit():
                cfg[k1][k2] = float(raw) if "." in raw else int(raw)
            else:
                cfg[k1][k2] = raw
        except Exception:
            cfg[k1][k2] = raw


def _postprocess(cfg: Dict[str, Any], os_key: str) -> None:
    """Pfad-Normalisierung, Headless-Default, Fallbacks, Ausgabeverzeichnis"""
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

    # Headless/Container => Cytoscape aus (falls nicht explizit gesetzt)
    if "open_in_cytoscape" not in cfg or cfg["open_in_cytoscape"] is None:
        if _is_headless_linux() or _is_container():
            cfg["open_in_cytoscape"] = False

    # BLAST-Fallback: Windows-Pfad vs. 'blastp' auf Unix
    if not cfg.get("blastp_executable"):
        cfg["blastp_executable"] = (
            r"C:\Program Files\NCBI\blast-2.17.0+\bin\blastp.exe"
            if os_key == "windows"
            else "blastp"
        )

    # Workers Defaults beibehalten; "auto" wird später im Code aufgelöst
    cfg.setdefault("workers", {})
    cfg["workers"].setdefault("parsing", "auto")
    cfg["workers"].setdefault("blast_threads", "auto")

    # Ausgabeordner vorsorglich anlegen (falls gesetzt)
    out = cfg.get("output_path")
    if isinstance(out, str) and out:
        try:
            Path(out).mkdir(parents=True, exist_ok=True)
        except Exception as e:
            _log(f"Warn: Konnte output_path '{out}' nicht anlegen: {e}")


def load_config() -> Dict[str, Any]:
    """Lädt Config in Layern: base → os-spezifisch → local → explizite Datei → ENV."""
    global _config_cache
    if _config_cache is not None:
        return _config_cache

    root = Path(__file__).resolve().parent
    os_key = {"Windows": "windows", "Linux": "linux", "Darwin": "darwin"}.get(
        platform.system(), platform.system().lower()
    )

    # Konfigurationsverzeichnis (standard: ./configs)
    cfg_dir = Path(os.environ.get("PDB2NET_CONFIG_DIR") or (root / "configs"))

    candidates: list[Tuple[Path, bool]] = [
        (cfg_dir / "config.base.json", False),
        (cfg_dir / f"config.{os_key}.json", False),
        (cfg_dir / "config.local.json", False),  # git-ignored, höchste Datei-Priorität unter configs/
    ]

    # Legacy-Fallback (Rückwärtskompatibilität), falls in configs/ nichts greift
    legacy = root / "config.json"

    # Explizite Datei per ENV (höchste Priorität & strict)
    explicit = os.environ.get("PDB2NET_CONFIG")
    if explicit:
        candidates.append((Path(explicit), True))  # strict=True

    cfg: Dict[str, Any] = {}
    loaded_any = False

    for path, strict in candidates:
        part = _read_json(path, strict=strict)
        if part:
            _deep_merge(cfg, part)
            loaded_any = True

    if not loaded_any:
        # Letzter Fallback auf alte Struktur
        if legacy.exists():
            _log("Nutze legacy config.json im Projekt-Root.")
            cfg = _read_json(legacy, strict=True)
        else:
            _log("Warn: Keine Konfigurationsdateien gefunden – starte mit leerer Config.")
            cfg = {}

    # ENV-Overrides anwenden
    _apply_env_overrides(cfg)

    # Post-Processing: Pfade normalisieren, Headless/Container, Defaults
    _postprocess(cfg, os_key=os_key)

    if VERBOSE:
        _log(f"OS={os_key} | headless={_is_headless_linux()} | container={_is_container()}")
        safe = dict(cfg)
        if "cytoscape_path" in safe and isinstance(safe["cytoscape_path"], str) and safe["cytoscape_path"]:
            safe["cytoscape_path"] = "<set>"
        _log("Aktive Schlüssel: " + ", ".join(sorted(safe.keys())))

    _config_cache = cfg
    return cfg


# Beim Import laden (lazy + Cache)
config = load_config()

__all__ = ["config", "load_config"]
