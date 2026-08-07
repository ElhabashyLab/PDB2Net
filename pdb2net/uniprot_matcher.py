"""UniProt matching via BLASTP.

This module:
- Optionally builds a local BLAST database from a UniProt FASTA.
- Extracts protein sequences from parsed chain data (3-letter → 1-letter).
- Runs BLASTP to find UniProt matches for eligible protein chains.
- Updates chain dictionaries in-place with:
    - 'uniprot_id'
    - 'molecule_name'
    - 'molecule_type'

Key behaviors
-------------
- Chains that already have a UniProt assignment are skipped:
    chain['uniprot_id'] not in {None, "Unknown"}
- Chains classified as nucleic acid are skipped for BLAST:
    'Nucleic Acid', 'DNA', 'RNA', 'DNA/RNA'
- For chains with unknown/empty molecule_type, a residue-based type inference is applied
  BEFORE BLAST selection (prevents nucleic-acid chains from being left as "Unknown").
- Chains with too many unknown residues are skipped for BLAST:
    X-fraction > 0.20
  If such a chain still has unknown/empty molecule_type, residue-based type inference is
  applied before skipping (prevents "Unknown" from leaking downstream).
- BLAST hit selection uses additional quality criteria:
    - query coverage (from BLAST-reported qcovs)
    - % identity
    - bitscore
    - e-value
- Ambiguity handling (softened, but guarded):
    - Dynamic bitscore margin based on best hit strength (floor/cap)
    - Small relaxations when best is reviewed (sp) or clearly better in qcov/pident
    - Escape hatch: accept "perfect ties" for qlen >= 80 when BOTH top hits are near-perfect
      (high identity + high coverage), choosing deterministically with preference for reviewed entries.

Debugging
---------
- PDB2NET_BLAST_DEBUG=1: verbose per-chain decisions.
- PDB2NET_BLAST_DIAG=1: aggregated counters printed at end of BLAST stage.
"""

from __future__ import annotations

import os
import re
import subprocess
import tempfile
import threading
import uuid
import hashlib
import json
import sqlite3
from datetime import datetime
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

from Bio.Data import IUPACData
from .config_loader import config
from .residue_types import (
    AMINO_ACIDS,
    DNA_RESIDUES,
    MODIFIED_RESIDUES_3TO3,
    NUCLEIC_ACID_TYPES,
    RNA_RESIDUES,
    normalize_residue_name,
)

# --- Paths from configuration ---
BLAST_DB_PATH: str = config["blast_db_path"]
BLAST_EXECUTABLE: str = config["blastp_executable"]
UNIPROT_FASTA_PATH: str = config["uniprot_fasta_path"]

_BLAST_CONFIG: Dict[str, Any] = config.get("blast", {}) if isinstance(config.get("blast", {}), dict) else {}
BLAST_MAX_TARGET_SEQS_DEFAULT: int = int(_BLAST_CONFIG.get("max_target_seqs_default", 50))
BLAST_MAX_TARGET_SEQS_SHORT: int = int(_BLAST_CONFIG.get("max_target_seqs_short", 100))
BLAST_SHORT_QUERY_LENGTH: int = int(_BLAST_CONFIG.get("short_query_length", 80))
BLAST_VERY_SHORT_QUERY_LENGTH: int = int(_BLAST_CONFIG.get("very_short_query_length", 30))
BLAST_USE_BLASTP_SHORT: bool = bool(_BLAST_CONFIG.get("use_blastp_short", True))
BLAST_MAX_HSPS: int = int(_BLAST_CONFIG.get("max_hsps", 1))
BLAST_MAX_X_FRACTION: float = float(_BLAST_CONFIG.get("max_x_fraction", 0.20))

# 3-letter → 1-letter amino-acid mapping (Biopython)
three_to_one: Dict[str, str] = IUPACData.protein_letters_3to1

# Common modified residues (PDB) → canonical amino acids.
# This reduces artificial 'X' inflation from frequent modifications
# (e.g., selenomethionine MSE), improving BLAST eligibility without
# changing the underlying polymer.
MODRES_3TO3: Dict[str, str] = MODIFIED_RESIDUES_3TO3

# Debug flag (opt-in)
_DEBUG_BLAST: bool = str(os.environ.get("PDB2NET_BLAST_DEBUG", "")).strip().lower() in {"1", "true", "yes", "on"}

# Diagnostics flag (opt-in). Prints an aggregated summary of where UniProt matching fails.
_DIAG_BLAST: bool = str(os.environ.get("PDB2NET_BLAST_DIAG", "")).strip().lower() in {"1", "true", "yes", "on"}
_DIAG_LOCK = threading.Lock()
_DIAG = Counter()

_UNIPROT_ACCESSION_PATTERN = (
    r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|"
    r"[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})"
)


def _is_uniprotkb_accession(value: Any) -> bool:
    """Return whether ``value`` is a complete UniProtKB accession or isoform ID."""
    accession = str(value or "").strip()
    return re.fullmatch(rf"{_UNIPROT_ACCESSION_PATTERN}(?:-[0-9]+)?", accession) is not None


def _diag_inc(key: str, n: int = 1) -> None:
    if not _DIAG_BLAST:
        return
    with _DIAG_LOCK:
        _DIAG[key] += n


def _diag_reset() -> None:
    if not _DIAG_BLAST:
        return
    with _DIAG_LOCK:
        _DIAG.clear()


def _diag_print_summary() -> None:
    if not _DIAG_BLAST:
        return
    with _DIAG_LOCK:
        if not _DIAG:
            print("\n[BLASTDIAG] Summary: (no diagnostic counters recorded)\n")
            return

        print("\n[BLASTDIAG] Summary (counts):")
        preferred_order = [
            "pre_chains_seen",
            "pre_skip_has_uniprot",
            "pre_skip_nucleic_acid",
            "pre_skip_empty_sequence",
            "pre_skip_xfrac",
            "pre_candidates_for_blast",
            "blast_calls",
            "accepted",
            "accepted_perfect_tie",
            "fail_ambiguous_margin",
            "fail_no_hits",
            "fail_evalue",
            "fail_qcov",
            "fail_pident",
            "fail_bitscore",
            "fail_thresholds_other",
            "blast_error",
            "blast_no_output",
            "blast_exception",
        ]

        shown = set()
        for k in preferred_order:
            if k in _DIAG:
                print(f"[BLASTDIAG]  {k:24s} {_DIAG[k]:8d}")
                shown.add(k)

        for k, v in sorted(((k, v) for k, v in _DIAG.items() if k not in shown), key=lambda kv: (-kv[1], kv[0])):
            print(f"[BLASTDIAG]  {k:24s} {v:8d}")
        print("")


# --- Persistent sequence-search cache (sequence → selected BlastHit / NO_HIT) ---
#
# This cache is used to avoid re-running BLAST for identical sequences across:
#   - multiple chains within one batch/run (dedupe is handled separately), and
#   - multiple batches / repeated runs (persistent SQLite cache).
#
# The cache is keyed by:
#   (db_signature, seq_key)
# where seq_key is sha1(sequence) + ":" + len(sequence).
#
# Cache can be disabled via: PDB2NET_BLAST_CACHE=0
_CACHE_ENABLED: bool = str(os.environ.get("PDB2NET_BLAST_CACHE", "1")).strip().lower() not in {"0", "false", "no", "off"}

BLAST_CACHE_PATH: str = str(
    config.get("blast_cache_path")
    or os.environ.get("PDB2NET_BLAST_CACHE_PATH", "").strip()
    or os.path.join(BLAST_DB_PATH, "blast_cache.sqlite3")
)

_CACHE_LOCK = threading.Lock()
_CACHE_CONN: sqlite3.Connection | None = None
_CACHE_DB_SIG: str | None = None


def _search_policy_payload() -> Dict[str, Any]:
    """Return result-relevant search settings used to namespace cache rows."""
    diamond_cfg = _diamond_config()
    return {
        "swissprot": {
            "max_target_seqs_default": BLAST_MAX_TARGET_SEQS_DEFAULT,
            "max_target_seqs_short": BLAST_MAX_TARGET_SEQS_SHORT,
            "short_query_length": BLAST_SHORT_QUERY_LENGTH,
            "very_short_query_length": BLAST_VERY_SHORT_QUERY_LENGTH,
            "use_blastp_short": BLAST_USE_BLASTP_SHORT,
            "max_hsps": BLAST_MAX_HSPS,
            "max_x_fraction": BLAST_MAX_X_FRACTION,
            "thresholds": {
                "short": _blast_search_parameters(20)[:4],
                "medium": _blast_search_parameters(100)[:4],
                "long": _blast_search_parameters(200)[:4],
            },
            "selection_policy": {
                "version": 1,
                "dynamic_margin_fraction": 0.05,
                "dynamic_margin_floor": 4.0,
                "dynamic_margin_cap": 10.0,
                "perfect_tie_min_length": 80,
                "perfect_tie_min_pident": 99.0,
                "perfect_tie_min_qcov": 0.95,
                "canonical_accession_validation": "uniprotkb-fullmatch-v1",
            },
        },
        "diamond": {
            "enabled": diamond_uniref90_enabled(),
            "sensitivity": str(diamond_cfg.get("sensitivity") or "sensitive").strip().lower(),
            "iterate": bool(diamond_cfg.get("iterate", True)),
            # DIAMOND documents a small sensitivity effect from block size,
            # so it belongs to the result policy rather than the chunk policy.
            "block_size": _positive_float_config(diamond_cfg.get("block_size"), 1.0),
            "max_target_seqs": _positive_int_config(diamond_cfg.get("max_target_seqs"), 50),
            "assign_uniprot_id": str(diamond_cfg.get("assign_uniprot_id") or "never").strip().lower(),
            "thresholds": {
                "short": _diamond_thresholds(20),
                "medium": _diamond_thresholds(100),
                "long": _diamond_thresholds(200),
            },
            "selection_policy": {
                "version": 1,
                "dynamic_margin_fraction": 0.05,
                "dynamic_margin_floor": 4.0,
                "dynamic_margin_cap": 10.0,
                "ambiguity_escape_min_pident": 99.0,
                "ambiguity_escape_min_qcov": 0.95,
                "high_confidence_min_pident": 95.0,
                "high_confidence_min_qcov": 0.90,
                "high_confidence_max_evalue": 1e-20,
                "canonical_accession_validation": "uniprotkb-fullmatch-v1",
            },
        },
    }


def _db_signature() -> str:
    """Compute a signature for active databases and result-affecting policies."""
    try:
        st = os.stat(UNIPROT_FASTA_PATH)
        swissprot_sig = f"uniprot_fasta:{st.st_size}:{st.st_mtime_ns}"
    except Exception:
        swissprot_sig = ""

    if not swissprot_sig:
        try:
            st = os.stat(os.path.join(BLAST_DB_PATH, "uniprot_db.pin"))
            swissprot_sig = f"blast_db:{st.st_size}:{st.st_mtime_ns}"
        except Exception:
            swissprot_sig = "unknown"

    diamond_sig = "diamond:off"
    if diamond_uniref90_enabled():
        db_path = get_diamond_uniref90_db_path()
        db_file = db_path if os.path.exists(db_path) else db_path + ".dmnd"
        try:
            st = os.stat(db_file)
            diamond_sig = f"diamond_uniref90:{db_file}:{st.st_size}:{st.st_mtime_ns}"
        except Exception:
            diamond_sig = f"diamond_uniref90:{db_path}:missing"

    policy_json = json.dumps(_search_policy_payload(), sort_keys=True, separators=(",", ":"))
    policy_sig = hashlib.sha256(policy_json.encode("utf-8")).hexdigest()
    return f"{swissprot_sig}|{diamond_sig}|policy:{policy_sig}"


def _seq_cache_key(sequence: str) -> str:
    """Stable cache key for a query sequence (robust across runs)."""
    h = hashlib.sha1(sequence.encode("utf-8")).hexdigest()
    return f"{h}:{len(sequence)}"


def _cache_init() -> None:
    global _CACHE_CONN, _CACHE_DB_SIG
    if not _CACHE_ENABLED:
        return
    if _CACHE_CONN is not None:
        return

    with _CACHE_LOCK:
        if _CACHE_CONN is not None:
            return
        try:
            cache_dir = os.path.dirname(BLAST_CACHE_PATH)
            if cache_dir:
                os.makedirs(cache_dir, exist_ok=True)

            conn = sqlite3.connect(BLAST_CACHE_PATH, timeout=30, check_same_thread=False)
            conn.execute("PRAGMA journal_mode=WAL;")
            conn.execute("PRAGMA synchronous=NORMAL;")
            conn.execute(
                "CREATE TABLE IF NOT EXISTS blast_cache ("
                "db_sig TEXT NOT NULL, "
                "seq_key TEXT NOT NULL, "
                "has_hit INTEGER NOT NULL, "
                "accession TEXT, "
                "reviewed INTEGER, "
                "bitscore REAL, "
                "evalue REAL, "
                "qcov REAL, "
                "pident REAL, "
                "title TEXT, "
                "source TEXT, "
                "database_name TEXT, "
                "matched_id TEXT, "
                "representative_accession TEXT, "
                "confidence TEXT, "
                "updated_at TEXT, "
                "PRIMARY KEY (db_sig, seq_key)"
                ")"
            )
            for ddl in [
                "ALTER TABLE blast_cache ADD COLUMN source TEXT",
                "ALTER TABLE blast_cache ADD COLUMN database_name TEXT",
                "ALTER TABLE blast_cache ADD COLUMN matched_id TEXT",
                "ALTER TABLE blast_cache ADD COLUMN representative_accession TEXT",
                "ALTER TABLE blast_cache ADD COLUMN confidence TEXT",
            ]:
                try:
                    conn.execute(ddl)
                except sqlite3.OperationalError:
                    pass
            conn.execute("CREATE INDEX IF NOT EXISTS idx_blast_cache_seq_key ON blast_cache(seq_key)")
            conn.commit()
            _CACHE_CONN = conn
            _CACHE_DB_SIG = _db_signature()
        except Exception:
            _CACHE_CONN = None
            _CACHE_DB_SIG = None
            _diag_inc("blast_cache_error")


def _cache_get(seq_key: str) -> Tuple[bool, Optional["BlastHit"]]:
    """Return (is_cached, BlastHit|None). is_cached=True includes negative cache (NO_HIT)."""
    if not _CACHE_ENABLED:
        return (False, None)
    try:
        _cache_init()
        if _CACHE_CONN is None or _CACHE_DB_SIG is None:
            return (False, None)

        with _CACHE_LOCK:
            if _CACHE_CONN is None or _CACHE_DB_SIG is None:
                return (False, None)
            row = _CACHE_CONN.execute(
                "SELECT has_hit, accession, reviewed, bitscore, evalue, qcov, pident, title, "
                "source, database_name, matched_id, representative_accession, confidence "
                "FROM blast_cache WHERE db_sig=? AND seq_key=?",
                (_CACHE_DB_SIG, seq_key),
            ).fetchone()

        if row is None:
            return (False, None)

        has_hit = int(row[0])
        if has_hit != 1:
            return (True, None)

        return (
            True,
            BlastHit(
                accession=str(row[1] or ""),
                reviewed=bool(int(row[2] or 0)),
                bitscore=float(row[3] or 0.0),
                evalue=float(row[4] or 0.0),
                qcov=float(row[5] or 0.0),
                pident=float(row[6] or 0.0),
                title=str(row[7] or ""),
                source=str(row[8] or "blastp_swissprot"),
                database=str(row[9] or "Swiss-Prot"),
                matched_id=str(row[10] or row[1] or ""),
                representative_accession=str(row[11] or row[1] or ""),
                confidence=str(row[12] or "high"),
            ),
        )
    except Exception:
        _diag_inc("blast_cache_error")
        return (False, None)


def _cache_put(seq_key: str, hit: Optional["BlastHit"]) -> None:
    """Store a BLAST result (including NO_HIT) in the persistent cache."""
    if not _CACHE_ENABLED:
        return
    try:
        _cache_init()
        if _CACHE_CONN is None or _CACHE_DB_SIG is None:
            return

        now = datetime.utcnow().isoformat(timespec="seconds") + "Z"

        with _CACHE_LOCK:
            if _CACHE_CONN is None or _CACHE_DB_SIG is None:
                return
            if hit is None:
                _CACHE_CONN.execute(
                    "INSERT OR REPLACE INTO blast_cache (db_sig, seq_key, has_hit, updated_at) VALUES (?,?,0,?)",
                    (_CACHE_DB_SIG, seq_key, now),
                )
            else:
                _CACHE_CONN.execute(
                    "INSERT OR REPLACE INTO blast_cache "
                    "(db_sig, seq_key, has_hit, accession, reviewed, bitscore, evalue, qcov, pident, title, "
                    "source, database_name, matched_id, representative_accession, confidence, updated_at) "
                    "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    (
                        _CACHE_DB_SIG,
                        seq_key,
                        1,
                        hit.accession,
                        int(bool(hit.reviewed)),
                        float(hit.bitscore),
                        float(hit.evalue),
                        float(hit.qcov),
                        float(hit.pident),
                        hit.title or "",
                        hit.source,
                        hit.database,
                        hit.matched_id or hit.accession,
                        hit.representative_accession or hit.accession,
                        hit.confidence,
                        now,
                    ),
                )
            _CACHE_CONN.commit()
        _diag_inc("blast_cache_store")
    except Exception:
        _diag_inc("blast_cache_error")


@dataclass(frozen=True)
class BlastHit:
    """Selected sequence-search hit used for chain annotation."""
    accession: str
    reviewed: bool
    bitscore: float
    evalue: float
    qcov: float
    pident: float
    title: str = ""
    source: str = "blastp_swissprot"
    database: str = "Swiss-Prot"
    matched_id: str = ""
    representative_accession: str = ""
    confidence: str = "high"


class SequenceSearchError(RuntimeError):
    """Raised when an external sequence-search process fails operationally."""


@dataclass(frozen=True)
class _SearchOutcome:
    """Internal tri-state result for a completed or failed search query."""

    status: str
    hit: Optional[BlastHit] = None
    error: str = ""

    @property
    def is_error(self) -> bool:
        return self.status == "error"


def _hit_outcome(hit: BlastHit) -> _SearchOutcome:
    return _SearchOutcome(status="hit", hit=hit)


def _no_hit_outcome() -> _SearchOutcome:
    return _SearchOutcome(status="no_hit")


def _error_outcome(message: str) -> _SearchOutcome:
    return _SearchOutcome(status="error", error=message)


def _coerce_search_outcome(value: Any) -> _SearchOutcome:
    """Normalize legacy test hooks returning ``BlastHit | None``."""
    if isinstance(value, _SearchOutcome):
        return value
    if isinstance(value, BlastHit):
        return _hit_outcome(value)
    return _no_hit_outcome()


def _outcome_hit_or_raise(outcome: _SearchOutcome, backend: str) -> Optional[BlastHit]:
    if outcome.is_error:
        raise SequenceSearchError(f"{backend} search failed: {outcome.error}")
    return outcome.hit


def _diamond_config() -> Dict[str, Any]:
    raw = config.get("diamond", {})
    return raw if isinstance(raw, dict) else {}


def diamond_uniref90_enabled() -> bool:
    """Return whether the optional DIAMOND/UniRef90 fallback is enabled."""
    return bool(_diamond_config().get("enabled", False))


def get_diamond_executable() -> str:
    """Return the configured DIAMOND executable name or path."""
    return str(_diamond_config().get("executable") or "diamond")


def get_diamond_uniref90_db_path() -> str:
    """Return the configured DIAMOND UniRef90 database path/prefix."""
    return str(_diamond_config().get("uniref90_db_path") or "")


def extract_direct_uniprot_accession(file_path_or_id: str) -> Optional[str]:
    """Return a UniProt accession encoded in an AlphaFold-style file name.

    This intentionally accepts only explicit AlphaFold/fragment-style names, e.g.
    ``AF-Q9BYF1-F1-model_v4.cif`` or ``Q8WZ42-F1.cif``. Generic custom names
    such as ``123456.cif`` are left for residue inference and BLAST fallback.
    """
    if not file_path_or_id:
        return None

    stem = os.path.splitext(os.path.basename(str(file_path_or_id)))[0].upper()
    acc = _UNIPROT_ACCESSION_PATTERN

    af_match = re.search(
        rf"(?:^|[^A-Z0-9])AF[-_](?P<accession>{acc}(?:-\d+)?)(?:[-_](?:F\d+|MODEL|V\d+)|$)",
        stem,
    )
    if af_match:
        return af_match.group("accession").replace("_", "-")

    fragment_match = re.search(rf"^(?P<accession>{acc}(?:-\d+)?)[-_]F\d+(?:[-_]|$)", stem)
    if fragment_match:
        return fragment_match.group("accession").replace("_", "-")

    return None


def _lookup_direct_uniprot_name(accession: str) -> Optional[str]:
    # process_molecule_info() runs before BLAST in the pipeline and loads this
    # cache already; reusing it avoids a second full UniProt FASTA parse.
    try:
        from . import unknown_molecule_uniprot

        names = unknown_molecule_uniprot.uniprot_dict
    except Exception:
        names = {}

    name = names.get(accession)
    if not name and "-" in accession:
        name = names.get(accession.split("-", 1)[0])
    if name and name != "Unknown Protein":
        return name
    return None


def _apply_direct_uniprot_to_structure(structure: Dict[str, Any], debug: bool = False) -> bool:
    """Annotate a single-chain AlphaFold-style structure from its file name."""
    file_path = str(structure.get("file_path") or "")
    pdb_id = str(structure.get("pdb_id") or "")
    accession = extract_direct_uniprot_accession(file_path) or extract_direct_uniprot_accession(pdb_id)
    if not accession:
        return False

    chains = list(structure.get("atom_data", []))
    if len(chains) != 1:
        _diag_inc("direct_uniprot_skip_multichain")
        if debug:
            print(
                f"[BLASTSEL][{os.path.basename(file_path) or pdb_id}] "
                f"SKIP_DIRECT_UNIPROT: accession={accession} chains={len(chains)}"
            )
        return False

    chain = chains[0]
    if chain.get("uniprot_id") not in [None, "Unknown"]:
        return False

    chain["uniprot_id"] = accession
    chain["molecule_type"] = "Protein"
    chain["molecule_name"] = _lookup_direct_uniprot_name(accession) or f"UniProt: {accession}"
    _diag_inc("direct_uniprot_assigned")

    if debug:
        chain_id = chain.get("chain_id", "?")
        uniq = chain.get("unique_chain_id", f"{pdb_id}:{chain_id}")
        print(f"[BLASTSEL][{uniq}] DIRECT_UNIPROT: {accession}")

    return True


def create_blast_database(force_rebuild: bool = False) -> None:
    """Create a BLAST database from the UniProt FASTA if needed."""
    os.makedirs(BLAST_DB_PATH, exist_ok=True)
    db_prefix = os.path.join(BLAST_DB_PATH, "uniprot_db")
    pin_file = db_prefix + ".pin"

    if os.path.exists(pin_file) and not force_rebuild:
        return

    makeblastdb = config.get("makeblastdb_executable", "makeblastdb")
    cmd = [
        makeblastdb,
        "-in", UNIPROT_FASTA_PATH,
        "-dbtype", "prot",
        "-out", db_prefix,
    ]
    subprocess.run(cmd, check=True)


def extract_sequence_from_parsed_data(chain_data: Dict[str, Any]) -> str:
    """Extract 1-letter protein-like sequence from parsed chain residues.

    Uses residue key 'residue_name' (project convention) with a fallback to 'resname'.
    Unknown / unmapped residues become 'X'.
    """
    residues = chain_data.get("residues", [])
    if not residues:
        return ""

    seq_chars: List[str] = []
    for res in residues:
        rn_raw = str(res.get("residue_name") or res.get("resname") or "").strip()
        if not rn_raw:
            continue

        rn = normalize_residue_name(rn_raw)

        # Biopython three_to_one uses capitalized 3-letter codes ("Ala", ...)
        seq_chars.append(three_to_one.get(rn.capitalize(), "X"))

    return "".join(seq_chars)


def _extract_dbtag_and_accession(sseqid: str) -> Tuple[str, str]:
    """Extract db tag + accession from BLAST sseqid.

    Supported common sseqid patterns:
      - sp|P12345|NAME
      - tr|Q9XXXX|NAME
      - P12345
      - gi|...|ref|P12345.1|
    """
    s = sseqid.strip()

    if "|" in s:
        parts = s.split("|")
        if len(parts) >= 2 and parts[0] in {"sp", "tr"}:
            dbtag = parts[0]
            accession = parts[1]
            return dbtag, accession

        for p in parts:
            p2 = p.strip()
            if p2 and p2[0].isalpha() and any(ch.isdigit() for ch in p2):
                accession = p2.split(".")[0]
                return "unk", accession

    accession = s.split()[0].split(".")[0]
    return "unk", accession


def _protein_name_from_uniprot_title(stitle: str) -> Optional[str]:
    """Try to extract a usable protein name from UniProt FASTA title."""
    if not stitle:
        return None

    left = stitle.strip()

    for token in ("OS=", "OX=", "GN=", "PE=", "SV="):
        if token in left:
            left = left.split(token, 1)[0].strip()

    if " " in left:
        _id, rest = left.split(" ", 1)
        name = rest.strip()
        return name or None

    return None


def _blast_search_parameters(qlen: int) -> Tuple[float, float, float, float, int, Optional[str]]:
    """Return thresholds and process options for one Swiss-Prot query."""
    if qlen < 80:
        thresholds = (1e-6, 0.85, 35.0, 60.0)
    elif qlen < 200:
        thresholds = (1e-10, 0.70, 30.0, 70.0)
    else:
        thresholds = (1e-20, 0.60, 25.0, 80.0)
    max_targets = BLAST_MAX_TARGET_SEQS_SHORT if qlen < BLAST_SHORT_QUERY_LENGTH else BLAST_MAX_TARGET_SEQS_DEFAULT
    blast_task = "blastp-short" if BLAST_USE_BLASTP_SHORT and qlen < BLAST_VERY_SHORT_QUERY_LENGTH else None
    return (*thresholds, max_targets, blast_task)


def _select_blastp_swissprot_hit(
    query_sequence: str,
    rows: List[str],
    *,
    use_qcovs: bool,
    label: str = "",
    debug: bool = False,
) -> Optional[BlastHit]:
    """Apply the established Swiss-Prot thresholds to one query's rows."""
    qlen0 = len(query_sequence)
    max_evalue, min_query_coverage, min_pident, min_bitscore, _max_targets, _task = (
        _blast_search_parameters(qlen0)
    )
    prefix = f"[BLASTDBG][{label or _seq_cache_key(query_sequence)[:8]}]"
    seen_any = False
    seen_pass_evalue = False
    seen_pass_qcov = False
    seen_pass_pident = False
    seen_pass_bitscore = False
    best_per_accession: Dict[str, Tuple[int, float, float, float, float, str, str]] = {}
    raw_top: List[Tuple[str, str, float, float, float, float, float]] = []

    for index, line in enumerate(rows):
        cols = line.rstrip("\n").split("\t")
        if (use_qcovs and len(cols) < 11) or (not use_qcovs and len(cols) < 10):
            continue
        dbtag, accession = _extract_dbtag_and_accession(cols[1].strip())
        try:
            pident = float(cols[2])
            qstart = int(float(cols[4]))
            qend = int(float(cols[5]))
            qlen = int(float(cols[6]))
            evalue = float(cols[8] if use_qcovs else cols[7])
            bitscore = float(cols[9] if use_qcovs else cols[8])
            title = cols[10] if use_qcovs else cols[9]
        except (TypeError, ValueError):
            continue
        if qlen <= 0:
            continue

        if use_qcovs:
            try:
                qcov = float(cols[7]) / 100.0
            except (TypeError, ValueError):
                qcov = max(0.0, min(1.0, abs(qend - qstart) + 1) / float(qlen))
        else:
            qcov = max(0.0, min(1.0, abs(qend - qstart) + 1) / float(qlen))

        seen_any = True
        reviewed = 1 if dbtag == "sp" else 0
        if index < 5:
            raw_top.append((dbtag, accession, pident, qcov, evalue, bitscore, min_query_coverage))
        seen_pass_evalue = seen_pass_evalue or evalue <= max_evalue
        seen_pass_qcov = seen_pass_qcov or qcov >= min_query_coverage
        seen_pass_pident = seen_pass_pident or pident >= min_pident
        seen_pass_bitscore = seen_pass_bitscore or bitscore >= min_bitscore

        if evalue > max_evalue:
            continue
        qcov_threshold = min_query_coverage
        if bitscore >= min_bitscore + 40.0 and pident >= min_pident + 5.0:
            qcov_threshold = max(0.0, qcov_threshold - 0.05)
        if qcov < qcov_threshold or pident < min_pident or bitscore < min_bitscore:
            continue

        current = (reviewed, bitscore, evalue, qcov, pident, accession, title)
        previous = best_per_accession.get(accession)
        if previous is None or current[:5] > previous[:5]:
            best_per_accession[accession] = current

    if not seen_any:
        _diag_inc("fail_no_hits")
        if debug:
            print(f"{prefix} NO_HITS: BLAST returned no rows.")
        return None
    if not best_per_accession:
        if not seen_pass_evalue:
            _diag_inc("fail_evalue")
        elif not seen_pass_qcov:
            _diag_inc("fail_qcov")
        elif not seen_pass_pident:
            _diag_inc("fail_pident")
        elif not seen_pass_bitscore:
            _diag_inc("fail_bitscore")
        else:
            _diag_inc("fail_thresholds_other")
        if debug:
            print(f"{prefix} NO_PASS: no accession passed thresholds.")
            for dbtag, accession, pident, qcov, evalue, bitscore, qcov_threshold in raw_top:
                print(
                    f"{prefix}   top: {dbtag}|{accession} pid={pident:.1f} qcov={qcov:.2f} "
                    f"ev={evalue:.2g} bs={bitscore:.1f} (qthr={qcov_threshold:.2f})"
                )
        return None

    ranked = sorted(
        best_per_accession.values(),
        key=lambda hit: (hit[0], hit[1], -hit[2], hit[3], hit[4]),
        reverse=True,
    )
    best = ranked[0]
    best_reviewed, best_bitscore, best_evalue, best_qcov, best_pident, best_acc, best_title = best
    dynamic_margin = max(4.0, min(10.0, best_bitscore * 0.05))
    relax = 2.0 if best_reviewed == 1 else 0.0
    if len(ranked) > 1:
        runner_up = ranked[1]
        if best_qcov - runner_up[3] >= 0.10 or best_pident - runner_up[4] >= 10.0:
            relax += 2.0
    effective_margin = max(2.0, dynamic_margin - relax)

    ties = []
    for candidate in ranked[1:]:
        _, candidate_bitscore, _, candidate_qcov, candidate_pident, candidate_acc, candidate_title = candidate
        if best_bitscore - candidate_bitscore <= effective_margin:
            ties.append((candidate_acc, candidate_bitscore, candidate_qcov, candidate_pident, candidate_title))
        else:
            break

    ambiguous = bool(ties)
    if ambiguous and qlen0 >= 80:
        contender = ranked[1]
        (
            contender_reviewed,
            contender_bitscore,
            _contender_evalue,
            contender_qcov,
            contender_pident,
            contender_acc,
            _contender_title,
        ) = contender
        if (
            best_pident >= 99.0
            and best_qcov >= 0.95
            and contender_pident >= 99.0
            and contender_qcov >= 0.95
        ):
            if best_reviewed != contender_reviewed:
                chosen = best if best_reviewed > contender_reviewed else contender
            elif best_bitscore != contender_bitscore:
                chosen = best if best_bitscore > contender_bitscore else contender
            else:
                chosen = best if best_acc <= contender_acc else contender
            best_reviewed, best_bitscore, best_evalue, best_qcov, best_pident, best_acc, best_title = chosen
            _diag_inc("accepted_perfect_tie")
            ambiguous = False

    if ambiguous:
        _diag_inc("fail_ambiguous_margin")
        if debug:
            print(
                f"{prefix} AMBIGUOUS: best={('sp' if best_reviewed else 'tr')}|{best_acc} "
                f"bitscore={best_bitscore:.1f} margin={effective_margin:.1f} ties={len(ties)}"
            )
        return None

    _diag_inc("accepted")
    return BlastHit(
        accession=best_acc,
        reviewed=bool(best_reviewed),
        bitscore=best_bitscore,
        evalue=best_evalue,
        qcov=best_qcov,
        pident=best_pident,
        title=best_title,
        source="blastp_swissprot",
        database="Swiss-Prot",
        matched_id=best_acc,
        representative_accession=best_acc,
        confidence="high",
    )


def _run_blastp_swissprot_batch(
    queries: List[Tuple[str, str, str]],
    *,
    max_workers: int = 1,
    debug: bool = False,
) -> Dict[str, _SearchOutcome]:
    """Run bounded Swiss-Prot Multi-FASTA chunks grouped by BLAST task."""
    results = {
        key: _error_outcome("Swiss-Prot search was not executed")
        for key, _sequence, _label in queries
    }
    if not queries:
        return results
    debug = bool(debug or _DEBUG_BLAST)
    blast_cfg = config.get("blast", {}) if isinstance(config.get("blast"), dict) else {}
    max_sequences = _positive_int_config(blast_cfg.get("batch_max_sequences"), 5000)
    max_fasta_bytes = _positive_int_config(blast_cfg.get("batch_max_fasta_bytes"), 50 * 1024 * 1024)
    groups: Dict[Tuple[float, int, Optional[str]], List[Tuple[str, str, str]]] = {}
    for query in queries:
        _qcov, _pident, _bitscore, max_targets, task = _blast_search_parameters(len(query[1]))[1:]
        max_evalue = _blast_search_parameters(len(query[1]))[0]
        groups.setdefault((max_evalue, max_targets, task), []).append(query)

    jobs: List[Tuple[float, int, Optional[str], List[Tuple[str, str, str]]]] = []
    for (max_evalue, max_targets, task), group_queries in groups.items():
        for chunk in _chunk_sequence_queries(
            group_queries,
            max_sequences=max_sequences,
            max_fasta_bytes=max_fasta_bytes,
        ):
            jobs.append((max_evalue, max_targets, task, chunk))

    def run_job(
        job: Tuple[float, int, Optional[str], List[Tuple[str, str, str]]]
    ) -> Dict[str, _SearchOutcome]:
        max_evalue, max_targets, task, chunk = job
        chunk_results = {
            key: _error_outcome("Swiss-Prot batch did not complete")
            for key, _sequence, _label in chunk
        }
        unique_id = str(uuid.uuid4())[:8]
        prefix = f"[BLASTDBG][batch {unique_id}]"
        _diag_inc("blast_calls")
        try:
            with tempfile.TemporaryDirectory(prefix="pdb2net-blast-") as temp_dir:
                query_file = os.path.join(temp_dir, "queries.fasta")
                output_file = os.path.join(temp_dir, "blast_results.tsv")
                query_ids: Dict[str, Tuple[str, str, str]] = {}
                with open(query_file, "w", encoding="utf-8") as handle:
                    for query_index, (key, sequence, label) in enumerate(chunk):
                        query_id = "query" if len(chunk) == 1 else f"q{query_index:06d}"
                        query_ids[query_id] = (key, sequence, label)
                        handle.write(f">{query_id}\n{sequence}\n")

                db_prefix = os.path.join(BLAST_DB_PATH, "uniprot_db")
                outfmt_qcovs = "6 qseqid sseqid pident length qstart qend qlen qcovs evalue bitscore stitle"
                outfmt_span = "6 qseqid sseqid pident length qstart qend qlen evalue bitscore stitle"

                def run_process(outfmt: str) -> subprocess.CompletedProcess[str]:
                    command = [
                        BLAST_EXECUTABLE,
                        "-query", query_file,
                        "-db", db_prefix,
                        "-out", output_file,
                        "-evalue", str(max_evalue),
                        "-max_target_seqs", str(max_targets),
                        "-outfmt", outfmt,
                    ]
                    if BLAST_MAX_HSPS > 0:
                        command.extend(["-max_hsps", str(BLAST_MAX_HSPS)])
                    if task:
                        command.extend(["-task", task])
                    return subprocess.run(command, capture_output=True, text=True)

                use_qcovs = True
                process_result = run_process(outfmt_qcovs)
                if process_result.returncode != 0:
                    stderr = (process_result.stderr or "").lower()
                    if any(token in stderr for token in ("qcovs", "outfmt", "unknown", "unrecognized")):
                        use_qcovs = False
                        process_result = run_process(outfmt_span)
                if process_result.returncode != 0:
                    _diag_inc("blast_error")
                    error = (process_result.stderr or "").strip()
                    for key in chunk_results:
                        chunk_results[key] = _error_outcome(
                            f"BLASTP exited with code {process_result.returncode}: {error[:400]}"
                        )
                    if debug:
                        print(f"{prefix} ERROR: {error[:400]}")
                    return chunk_results
                if not os.path.exists(output_file):
                    _diag_inc("blast_no_output")
                    for key in chunk_results:
                        chunk_results[key] = _error_outcome("BLASTP did not create its output file")
                    return chunk_results

                rows_by_query: Dict[str, List[str]] = {query_id: [] for query_id in query_ids}
                with open(output_file, "r", encoding="utf-8") as handle:
                    for line in handle:
                        query_id = line.split("\t", 1)[0]
                        if query_id in rows_by_query:
                            rows_by_query[query_id].append(line)
                for query_id, (key, sequence, label) in query_ids.items():
                    hit = _select_blastp_swissprot_hit(
                        sequence,
                        rows_by_query[query_id],
                        use_qcovs=use_qcovs,
                        label=label,
                        debug=debug,
                    )
                    chunk_results[key] = _hit_outcome(hit) if hit is not None else _no_hit_outcome()
        except Exception as exc:
            _diag_inc("blast_exception")
            for key in chunk_results:
                chunk_results[key] = _error_outcome(f"BLASTP batch exception: {exc}")
            if debug:
                print(f"{prefix} EXCEPTION: {exc}")
        return chunk_results

    worker_count = max(1, min(_positive_int_config(max_workers, 1), len(jobs)))
    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        for chunk_results in executor.map(run_job, jobs):
            results.update(chunk_results)
    return results


def _run_blastp_swissprot_search(query_sequence: str, label: str = "", debug: bool = False) -> Optional[BlastHit]:
    """Run Swiss-Prot BLASTP for one query via the bounded batch path."""
    outcome = _run_blastp_swissprot_batch(
        [("query", query_sequence, label)],
        max_workers=1,
        debug=debug,
    )["query"]
    return _outcome_hit_or_raise(outcome, "Swiss-Prot BLASTP")


def _extract_uniref90_representative(sseqid: str, title: str = "") -> Tuple[str, str]:
    """Return (matched UniRef90 ID, representative accession) from a DIAMOND subject ID."""
    matched_id = str(sseqid or "").strip().split()[0]
    representative = ""

    if matched_id.startswith("UniRef90_"):
        representative = matched_id[len("UniRef90_"):]
    else:
        match = re.search(r"(UniRef90_[A-Za-z0-9_.-]+)", f"{sseqid} {title}")
        if match:
            matched_id = match.group(1)
            representative = matched_id[len("UniRef90_"):]

    if not representative:
        accession_match = re.search(_UNIPROT_ACCESSION_PATTERN, f"{sseqid} {title}")
        if accession_match:
            representative = accession_match.group(0)

    return matched_id or representative, representative


def _diamond_thresholds(qlen: int) -> Tuple[float, float, float, float]:
    """Return e-value, qcov, pident, and bitscore thresholds for DIAMOND fallback."""
    if qlen < 80:
        return 1e-6, 0.85, 35.0, 60.0
    if qlen < 200:
        return 1e-10, 0.70, 30.0, 70.0
    return 1e-20, 0.60, 25.0, 80.0


def _positive_int_config(raw: Any, default: int) -> int:
    """Return a positive integer config value, falling back for invalid input."""
    try:
        value = int(raw)
    except (TypeError, ValueError):
        return default
    return value if value > 0 else default


def _positive_float_config(raw: Any, default: float) -> float:
    """Return a positive float config value, falling back for invalid input."""
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return default
    return value if value > 0 else default


def _chunk_sequence_queries(
    queries: List[Tuple[str, str, str]],
    *,
    max_sequences: int,
    max_fasta_bytes: int,
) -> List[List[Tuple[str, str, str]]]:
    """Split ``(key, sequence, label)`` queries at FASTA count/byte boundaries.

    A single sequence larger than the configured byte limit is kept as its own
    chunk; silently dropping it would change annotation semantics.
    """
    chunks: List[List[Tuple[str, str, str]]] = []
    current: List[Tuple[str, str, str]] = []
    current_bytes = 0

    for key, sequence, label in queries:
        record_bytes = len(f">{key}\n{sequence}\n".encode("utf-8"))
        exceeds_count = len(current) >= max_sequences
        exceeds_bytes = bool(current) and current_bytes + record_bytes > max_fasta_bytes
        if exceeds_count or exceeds_bytes:
            chunks.append(current)
            current = []
            current_bytes = 0
        current.append((key, sequence, label))
        current_bytes += record_bytes

    if current:
        chunks.append(current)
    return chunks


def _select_diamond_hit(
    query_sequence: str,
    rows: List[str],
    *,
    label: str = "",
    debug: bool = False,
) -> Optional[BlastHit]:
    """Select one conservative UniRef90 hit from rows for a single query."""
    qlen0 = len(query_sequence)
    max_evalue, min_query_coverage, min_pident, min_bitscore = _diamond_thresholds(qlen0)
    prefix = f"[DIAMONDDBG][{label or _seq_cache_key(query_sequence)[:8]}]"
    best_per_match: Dict[str, Tuple[float, float, float, float, str, str, str]] = {}
    seen_any = False

    for line in rows:
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 10:
            continue

        sseqid = cols[1].strip()
        try:
            pident = float(cols[2])
            qstart = int(float(cols[4]))
            qend = int(float(cols[5]))
            qlen = int(float(cols[6]))
            evalue = float(cols[7])
            bitscore = float(cols[8])
            title = cols[9]
        except (TypeError, ValueError):
            continue

        if qlen <= 0:
            continue

        seen_any = True
        qcov = max(0.0, min(1.0, (abs(qend - qstart) + 1) / float(qlen)))
        if (
            evalue > max_evalue
            or qcov < min_query_coverage
            or pident < min_pident
            or bitscore < min_bitscore
        ):
            continue

        matched_id, representative = _extract_uniref90_representative(sseqid, title)
        if not matched_id:
            continue

        cur = (bitscore, -evalue, qcov, pident, matched_id, representative, title)
        prev = best_per_match.get(matched_id)
        if prev is None or cur[:4] > prev[:4]:
            best_per_match[matched_id] = cur

    if not seen_any:
        _diag_inc("diamond_fail_no_hits")
        return None
    if not best_per_match:
        _diag_inc("diamond_fail_thresholds")
        return None

    ranked = sorted(best_per_match.values(), reverse=True)
    best_bitscore, neg_evalue, best_qcov, best_pident, best_id, representative, best_title = ranked[0]
    best_evalue = -neg_evalue

    dynamic_margin = max(4.0, min(10.0, best_bitscore * 0.05))
    if len(ranked) > 1 and (best_bitscore - ranked[1][0]) <= dynamic_margin:
        if not (best_pident >= 99.0 and best_qcov >= 0.95):
            _diag_inc("diamond_fail_ambiguous_margin")
            return None

    confidence = "high" if best_pident >= 95.0 and best_qcov >= 0.90 and best_evalue <= 1e-20 else "medium"
    _diag_inc("diamond_accepted")
    if debug:
        print(
            f"{prefix} HIT: {best_id} qcov={best_qcov:.2f} "
            f"pid={best_pident:.1f} conf={confidence}"
        )
    return BlastHit(
        accession=representative or best_id,
        reviewed=False,
        bitscore=best_bitscore,
        evalue=best_evalue,
        qcov=best_qcov,
        pident=best_pident,
        title=best_title,
        source="diamond_uniref90",
        database="UniRef90",
        matched_id=best_id,
        representative_accession=representative,
        confidence=confidence,
    )


def _run_diamond_uniref90_batch(
    queries: List[Tuple[str, str, str]],
    *,
    debug: bool = False,
) -> Dict[str, _SearchOutcome]:
    """Run one DIAMOND Multi-FASTA process per bounded query chunk."""
    results = {
        key: _error_outcome("DIAMOND search was not executed")
        for key, _sequence, _label in queries
    }
    if not queries:
        return results
    if not diamond_uniref90_enabled():
        return {key: _no_hit_outcome() for key in results}

    debug = bool(debug or _DEBUG_BLAST)
    diamond_cfg = _diamond_config()
    max_sequences = _positive_int_config(diamond_cfg.get("batch_max_sequences"), 5000)
    max_fasta_bytes = _positive_int_config(diamond_cfg.get("batch_max_fasta_bytes"), 50 * 1024 * 1024)
    threads = _positive_int_config(diamond_cfg.get("threads"), 6)
    max_targets = _positive_int_config(diamond_cfg.get("max_target_seqs"), 50)
    index_chunks = _positive_int_config(diamond_cfg.get("index_chunks"), 4)
    block_size = _positive_float_config(diamond_cfg.get("block_size"), 1.0)

    executable = get_diamond_executable()
    db_path = get_diamond_uniref90_db_path()
    sensitivity = str(diamond_cfg.get("sensitivity") or "sensitive").strip().lower()
    iterate = bool(diamond_cfg.get("iterate", True))
    temp_base = str(diamond_cfg.get("temp_dir") or "").strip() or None
    chunks = _chunk_sequence_queries(
        queries,
        max_sequences=max_sequences,
        max_fasta_bytes=max_fasta_bytes,
    )

    for chunk_index, chunk in enumerate(chunks, start=1):
        prefix = f"[DIAMONDDBG][batch {chunk_index}/{len(chunks)}]"
        try:
            with tempfile.TemporaryDirectory(prefix="pdb2net-diamond-", dir=temp_base) as temp_dir:
                query_file = os.path.join(temp_dir, "queries.fasta")
                output_file = os.path.join(temp_dir, "diamond_results.tsv")
                query_ids: Dict[str, Tuple[str, str, str]] = {}

                with open(query_file, "w", encoding="utf-8") as handle:
                    for query_index, (key, sequence, label) in enumerate(chunk):
                        query_id = "query" if len(chunk) == 1 else f"q{query_index:06d}"
                        query_ids[query_id] = (key, sequence, label)
                        handle.write(f">{query_id}\n{sequence}\n")

                # Use the loosest query-length threshold for the process and apply
                # exact per-query thresholds while parsing its rows below.
                process_evalue = max(_diamond_thresholds(len(sequence))[0] for _, sequence, _ in chunk)
                cmd = [
                    executable,
                    "blastp",
                    "--db", db_path,
                    "--query", query_file,
                    "--out", output_file,
                    "--evalue", str(process_evalue),
                    "--threads", str(threads),
                    "--block-size", str(block_size),
                    "--index-chunks", str(index_chunks),
                    "--max-target-seqs", str(max_targets),
                    "--tmpdir", temp_dir,
                    "--outfmt", "6", "qseqid", "sseqid", "pident", "length",
                    "qstart", "qend", "qlen", "evalue", "bitscore", "stitle",
                ]
                if iterate:
                    cmd.append("--iterate")
                if sensitivity and sensitivity != "default":
                    cmd.append(f"--{sensitivity}")

                result = subprocess.run(cmd, capture_output=True, text=True)
                if result.returncode != 0:
                    _diag_inc("diamond_error")
                    error = (result.stderr or "").strip()
                    for key, _sequence, _label in chunk:
                        results[key] = _error_outcome(
                            f"DIAMOND exited with code {result.returncode}: {error[:400]}"
                        )
                    if debug:
                        print(f"{prefix} ERROR returncode={result.returncode} stderr={error[:400]}")
                    continue
                if not os.path.exists(output_file):
                    _diag_inc("diamond_no_output")
                    for key, _sequence, _label in chunk:
                        results[key] = _error_outcome("DIAMOND did not create its output file")
                    continue

                rows_by_query: Dict[str, List[str]] = {query_id: [] for query_id in query_ids}
                with open(output_file, "r", encoding="utf-8") as handle:
                    for line in handle:
                        query_id = line.split("\t", 1)[0]
                        if query_id in rows_by_query:
                            rows_by_query[query_id].append(line)

                for query_id, (key, sequence, label) in query_ids.items():
                    hit = _select_diamond_hit(
                        sequence,
                        rows_by_query[query_id],
                        label=label,
                        debug=debug,
                    )
                    results[key] = _hit_outcome(hit) if hit is not None else _no_hit_outcome()
        except Exception as exc:
            _diag_inc("diamond_exception")
            for key, _sequence, _label in chunk:
                results[key] = _error_outcome(f"DIAMOND batch exception: {exc}")
            if debug:
                print(f"{prefix} EXCEPTION: {exc}")

    return results


def _run_diamond_uniref90_search(query_sequence: str, label: str = "", debug: bool = False) -> Optional[BlastHit]:
    """Run the optional DIAMOND fallback for one query via the batch path."""
    outcome = _run_diamond_uniref90_batch(
        [("query", query_sequence, label)],
        debug=debug,
    )["query"]
    return _outcome_hit_or_raise(outcome, "DIAMOND/UniRef90")


def run_blast_search(query_sequence: str, label: str = "", debug: bool = False) -> Optional[BlastHit]:
    """Run Swiss-Prot BLASTP, then optional DIAMOND/UniRef90 fallback."""
    swissprot = _run_blastp_swissprot_batch(
        [("query", query_sequence, label)],
        max_workers=1,
        debug=debug,
    )["query"]
    if swissprot.is_error:
        return _outcome_hit_or_raise(swissprot, "Swiss-Prot BLASTP")
    if swissprot.hit is not None:
        return swissprot.hit
    diamond = _run_diamond_uniref90_batch(
        [("query", query_sequence, label)],
        debug=debug,
    )["query"]
    return _outcome_hit_or_raise(diamond, "DIAMOND/UniRef90")


def classify_molecule_type(chain_data: Dict[str, Any], label: str = "", debug: bool = False) -> str:
    """Infer molecule type based on residue composition.

    Reads residue key 'residue_name' with fallback to 'resname'.
    """
    residues = chain_data.get("residues", [])
    if not residues:
        return "Unknown"

    aa_count = 0
    nt_count = 0
    unk_count = 0

    for r in residues:
        res = str(r.get("residue_name") or r.get("resname") or "").strip().upper()
        if not res:
            continue

        res = normalize_residue_name(res)

        if res in AMINO_ACIDS:
            aa_count += 1
        elif res in DNA_RESIDUES or res in RNA_RESIDUES:
            nt_count += 1
        else:
            unk_count += 1

    total = aa_count + nt_count + unk_count
    if total == 0:
        return "Unknown"

    aa_frac = aa_count / total
    nt_frac = nt_count / total

    if aa_frac >= 0.50:
        return "Protein"
    if nt_frac >= 0.50:
        if any(normalize_residue_name(rr.get("residue_name") or rr.get("resname")) in DNA_RESIDUES for rr in residues):
            if any(normalize_residue_name(rr.get("residue_name") or rr.get("resname")) in RNA_RESIDUES for rr in residues):
                return "DNA/RNA"
            return "DNA"
        if any(normalize_residue_name(rr.get("residue_name") or rr.get("resname")) in RNA_RESIDUES for rr in residues):
            return "RNA"
        return "Nucleic Acid"

    return "Unknown"


def parallel_blast_search(parsed_data: List[Dict[str, Any]], max_workers: int = 16) -> None:
    """Run BLASTP in parallel for all eligible chains and update them in-place."""
    debug = _DEBUG_BLAST

    _diag_reset()

    max_x_fraction: float = BLAST_MAX_X_FRACTION
    na_types = NUCLEIC_ACID_TYPES

    # Dedupe by sequence: seq_key -> { "qlen": int, "seq": str, "targets": [(chain_dict, label), ...] }
    seq_groups: Dict[str, Dict[str, Any]] = {}

    def _apply_result(chain: Dict[str, Any], label: str, hit: Optional[BlastHit]) -> None:
        if hit:
            chain["annotation_source"] = hit.source
            chain["matched_database"] = hit.database
            chain["matched_id"] = hit.matched_id or hit.accession
            chain["representative_accession"] = hit.representative_accession or hit.accession
            chain["annotation_confidence"] = hit.confidence

            assign_uniprot = True
            if hit.source == "diamond_uniref90":
                assign_policy = str(_diamond_config().get("assign_uniprot_id", "never")).strip().lower()
                assign_uniprot = assign_policy == "always" or (
                    assign_policy == "high_confidence" and hit.confidence == "high"
                )

            if assign_uniprot and _is_uniprotkb_accession(hit.accession):
                chain["uniprot_id"] = hit.accession
            elif assign_uniprot:
                _diag_inc("skip_non_uniprotkb_accession")
            chain["molecule_type"] = "Protein"

            name = _protein_name_from_uniprot_title(hit.title) if hit.title else None
            if name:
                chain["molecule_name"] = name
            else:
                chain["molecule_name"] = (
                    f"Matched {hit.database}: {hit.matched_id or hit.accession}"
                    if hit.source == "diamond_uniref90"
                    else f"Matched UniProt: {hit.accession}"
                )

            if debug:
                src = hit.source
                print(
                    f"[BLASTRES][{label}] HIT: {src}|{hit.accession} qcov={hit.qcov:.2f} "
                    f"pid={hit.pident:.1f} conf={hit.confidence} → molecule_type=Protein"
                )
        else:
            if debug:
                print(f"[BLASTRES][{label}] NO_HIT")
            mt = (chain.get("molecule_type") or "").strip()
            if mt in ("", "Unknown") or chain.get("molecule_type") is None:
                chain["molecule_type"] = classify_molecule_type(chain, label=label, debug=debug)

    for structure in parsed_data:
        pdb_id = structure.get("pdb_id", "?")
        file_path = structure.get("file_path", "?")
        _apply_direct_uniprot_to_structure(structure, debug=debug)

        for chain in structure.get("atom_data", []):
            _diag_inc("pre_chains_seen")
            chain_id = chain.get("chain_id", "?")
            uniq = chain.get("unique_chain_id", f"{pdb_id}:{chain_id}")
            label = f"{uniq} ({os.path.basename(str(file_path))})"

            if chain.get("uniprot_id") not in [None, "Unknown"]:
                _diag_inc("pre_skip_has_uniprot")
                if debug:
                    print(f"[BLASTSEL][{label}] SKIP: already has uniprot_id={chain.get('uniprot_id')}")
                continue

            mt0 = chain.get("molecule_type")
            mt = (mt0 or "").strip()
            if mt in ("", "Unknown") or mt0 is None:
                inferred = classify_molecule_type(chain, label=label, debug=debug)
                chain["molecule_type"] = inferred
                mt = inferred

            if mt in na_types:
                _diag_inc("pre_skip_nucleic_acid")
                if debug:
                    print(f"[BLASTSEL][{label}] SKIP: molecule_type={mt}")
                continue

            sequence = extract_sequence_from_parsed_data(chain)
            if not sequence:
                _diag_inc("pre_skip_empty_sequence")
                if debug:
                    print(f"[BLASTSEL][{label}] SKIP: empty sequence (residues={len(chain.get('residues', []))})")
                continue

            sequence = sequence.replace("O", "X")
            x_fraction = (sequence.count("X") / len(sequence)) if sequence else 1.0

            if debug:
                print(f"[BLASTSEL][{label}] candidate mol_type={mt or 'Unknown'} qlen={len(sequence)} Xfrac={x_fraction:.3f}")

            if x_fraction > max_x_fraction:
                _diag_inc("pre_skip_xfrac")
                mt1 = (chain.get("molecule_type") or "").strip()
                if mt1 in ("", "Unknown") or chain.get("molecule_type") is None:
                    chain["molecule_type"] = classify_molecule_type(chain, label=label, debug=debug)

                if debug:
                    print(
                        f"[BLASTSEL][{label}] SKIP_X: Xfrac={x_fraction:.3f} > {max_x_fraction:.2f} "
                        f"(molecule_type={chain.get('molecule_type')})"
                    )
                continue

            _diag_inc("pre_candidates_for_blast")

            seq_key = _seq_cache_key(sequence)
            entry = seq_groups.get(seq_key)
            if entry is None:
                seq_groups[seq_key] = {"qlen": len(sequence), "seq": sequence, "targets": [(chain, label)]}
            else:
                entry["targets"].append((chain, label))

    # Diagnostics: how much we saved by in-run dedupe
    if seq_groups:
        unique_sequences = len(seq_groups)
        candidate_chains = sum(len(e["targets"]) for e in seq_groups.values())
        saved_calls = candidate_chains - unique_sequences
        _diag_inc("pre_dedup_unique_sequences", unique_sequences)
        if saved_calls > 0:
            _diag_inc("pre_dedup_saved_blast_calls", saved_calls)

    # Split into cached vs. needs BLAST (persistent cache across runs/batches)
    cached_entries: List[Tuple[Optional[BlastHit], List[Tuple[Dict[str, Any], str]]]] = []
    to_blast: List[Tuple[int, str, str, List[Tuple[Dict[str, Any], str]]]] = []  # (qlen, seq_key, seq, targets)

    for seq_key, entry in seq_groups.items():
        is_cached, cached_hit = _cache_get(seq_key)
        if is_cached:
            _diag_inc("blast_cache_hit")
            cached_entries.append((cached_hit, entry["targets"]))
        else:
            _diag_inc("blast_cache_miss")
            to_blast.append((int(entry["qlen"]), seq_key, str(entry["seq"]), entry["targets"]))

    # Apply cached results immediately (hit or NO_HIT)
    for cached_hit, targets in cached_entries:
        for chain, label in targets:
            _apply_result(chain, label, cached_hit)

    # Keep the old behavior of sorting by query length for larger workloads
    if len(to_blast) > 40:
        to_blast.sort(key=lambda t: t[0])

    # Swiss-Prot remains the first fallback. DIAMOND receives only the unique,
    # uncached sequences for which Swiss-Prot did not produce an accepted hit.
    swissprot_queries: List[Tuple[str, str, str]] = []
    for _qlen, seq_key, seq, targets in to_blast:
        label0 = targets[0][1] if targets else ""
        label_call = f"{label0} (+{len(targets) - 1} dup)" if len(targets) > 1 else label0
        swissprot_queries.append((seq_key, seq, label_call))
    search_results = {
        key: _coerce_search_outcome(value)
        for key, value in _run_blastp_swissprot_batch(
            swissprot_queries,
            max_workers=max_workers,
            debug=debug,
        ).items()
    }

    diamond_queries: List[Tuple[str, str, str]] = []
    if diamond_uniref90_enabled():
        for _qlen, seq_key, seq, targets in to_blast:
            swissprot_outcome = search_results.get(
                seq_key,
                _error_outcome("Swiss-Prot batch omitted the query result"),
            )
            if swissprot_outcome.status != "no_hit":
                continue
            label0 = targets[0][1] if targets else ""
            label_call = f"{label0} (+{len(targets) - 1} dup)" if len(targets) > 1 else label0
            diamond_queries.append((seq_key, seq, label_call))

        if diamond_queries:
            search_results.update(
                {
                    key: _coerce_search_outcome(value)
                    for key, value in _run_diamond_uniref90_batch(
                        diamond_queries,
                        debug=debug,
                    ).items()
                }
            )

    search_errors: List[str] = []
    for _qlen, seq_key, _seq, targets in to_blast:
        outcome = search_results.get(
            seq_key,
            _error_outcome("Sequence-search batch omitted the query result"),
        )
        hit = outcome.hit
        if outcome.is_error:
            _diag_inc("search_error")
            search_errors.append(f"{seq_key[:12]}: {outcome.error}")
            if debug:
                print(f"[BLASTRES][{targets[0][1] if targets else seq_key}] SEARCH_ERROR: {outcome.error}")
        else:
            _cache_put(seq_key, hit)

        for chain, label in targets:
            _apply_result(chain, label, hit)

    if search_errors:
        raise SequenceSearchError(
            f"Sequence search failed for {len(search_errors)} unique query sequence(s); "
            f"first error: {search_errors[0]}"
        )

    _diag_print_summary()
