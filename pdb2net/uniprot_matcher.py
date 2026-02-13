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
import subprocess
import threading
import uuid
import hashlib
import sqlite3
from datetime import datetime
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

from Bio.Data import IUPACData
from config_loader import config

# --- Paths from configuration ---
BLAST_DB_PATH: str = config["blast_db_path"]
BLAST_EXECUTABLE: str = config["blastp_executable"]
UNIPROT_FASTA_PATH: str = config["uniprot_fasta_path"]

# 3-letter → 1-letter amino-acid mapping (Biopython)
three_to_one: Dict[str, str] = IUPACData.protein_letters_3to1

# Common modified residues (PDB) → canonical amino acids.
# This reduces artificial 'X' inflation from frequent modifications
# (e.g., selenomethionine MSE), improving BLAST eligibility without
# changing the underlying polymer.
MODRES_3TO3: Dict[str, str] = {
    # very common in crystal structures
    "MSE": "MET",  # selenomethionine
    # common phosphorylation variants
    "SEP": "SER",  # phosphoserine
    "TPO": "THR",  # phosphothreonine
    "PTR": "TYR",  # phosphotyrosine
    # a few frequent oxidations / special cases
    "CSO": "CYS",  # S-hydroxycysteine
    "HYP": "PRO",  # hydroxyproline
}

# Set of valid amino acid residue codes (3-letter)
AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "SEC", "PYL",
}

# Debug flag (opt-in)
_DEBUG_BLAST: bool = str(os.environ.get("PDB2NET_BLAST_DEBUG", "")).strip().lower() in {"1", "true", "yes", "on"}

# Diagnostics flag (opt-in). Prints an aggregated summary of where UniProt matching fails.
_DIAG_BLAST: bool = str(os.environ.get("PDB2NET_BLAST_DIAG", "")).strip().lower() in {"1", "true", "yes", "on"}
_DIAG_LOCK = threading.Lock()
_DIAG = Counter()


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


# --- Persistent BLAST cache (sequence → selected BlastHit / NO_HIT) ---
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


def _db_signature() -> str:
    """Compute a lightweight signature for the current UniProt/BLAST database."""
    try:
        st = os.stat(UNIPROT_FASTA_PATH)
        return f"uniprot_fasta:{st.st_size}:{st.st_mtime_ns}"
    except Exception:
        pass

    try:
        st = os.stat(os.path.join(BLAST_DB_PATH, "uniprot_db.pin"))
        return f"blast_db:{st.st_size}:{st.st_mtime_ns}"
    except Exception:
        pass

    return "unknown"


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
                "updated_at TEXT, "
                "PRIMARY KEY (db_sig, seq_key)"
                ")"
            )
            conn.execute("CREATE INDEX IF NOT EXISTS idx_blast_cache_seq_key ON blast_cache(seq_key)")
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

        row = _CACHE_CONN.execute(
            "SELECT has_hit, accession, reviewed, bitscore, evalue, qcov, pident, title "
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

        if hit is None:
            _CACHE_CONN.execute(
                "INSERT OR REPLACE INTO blast_cache (db_sig, seq_key, has_hit, updated_at) VALUES (?,?,0,?)",
                (_CACHE_DB_SIG, seq_key, now),
            )
        else:
            _CACHE_CONN.execute(
                "INSERT OR REPLACE INTO blast_cache "
                "(db_sig, seq_key, has_hit, accession, reviewed, bitscore, evalue, qcov, pident, title, updated_at) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?)",
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
                    now,
                ),
            )
        _diag_inc("blast_cache_store")
    except Exception:
        _diag_inc("blast_cache_error")


@dataclass(frozen=True)
class BlastHit:
    """Selected BLAST hit used for chain→UniProt assignment."""
    accession: str
    reviewed: bool
    bitscore: float
    evalue: float
    qcov: float
    pident: float
    title: str = ""


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

        rn = rn_raw.upper()
        rn = MODRES_3TO3.get(rn, rn)  # normalize common modified residues

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


def run_blast_search(query_sequence: str, label: str = "", debug: bool = False) -> Optional[BlastHit]:
    """Run BLASTP for a given protein sequence and return a robust UniProt match."""
    debug = bool(debug or _DEBUG_BLAST)

    _diag_inc("blast_calls")

    # --- Ambiguity handling (softened) ---
    # Old behavior: fixed margin of 10 bitscore points.
    # New behavior: dynamic margin based on best bitscore, with small relaxations
    # when best is reviewed (sp) or clearly better in qcov/pident.
    min_bitscore_margin_cap: float = 10.0
    margin_floor: float = 4.0
    margin_fraction_of_best: float = 0.05  # 5% of best bitscore (clamped to floor/cap)
    margin_relax: float = 2.0              # subtract when "safer" conditions apply
    min_margin_after_relax: float = 2.0

    # --- "Perfect-tie acceptance" guardrails ---
    allow_perfect_tie_if_qlen_at_least: int = 80
    perfect_min_pident: float = 99.0
    perfect_min_qcov: float = 0.95

    # --- Slight qcov slack only for clearly strong hits ---
    qcov_slack: float = 0.05
    strong_bitscore_bonus: float = 40.0
    strong_pident_bonus: float = 5.0

    qlen0 = len(query_sequence)

    # --- Length-dependent thresholds ---
    if qlen0 < 80:
        max_evalue = 1e-6
        min_query_coverage = 0.85
        min_pident = 35.0
        min_bitscore = 60.0
    elif qlen0 < 200:
        max_evalue = 1e-10
        min_query_coverage = 0.70
        min_pident = 30.0
        min_bitscore = 70.0
    else:
        max_evalue = 1e-20
        min_query_coverage = 0.60
        min_pident = 25.0
        min_bitscore = 80.0

    max_targets = 25  # give ambiguity rule enough candidates

    unique_id = str(uuid.uuid4())[:8]
    prefix = f"[BLASTDBG][{label or unique_id}]"

    query_file = f"query_{unique_id}.fasta"
    output_file = f"blast_results_{unique_id}.txt"

    try:
        if debug:
            print(
                f"{prefix} qlen={qlen0} thr: evalue<={max_evalue:g}, "
                f"qcov>={min_query_coverage:.2f}, pident>={min_pident:.1f}, bitscore>={min_bitscore:.1f}, "
                f"max_targets={max_targets}"
            )

        with open(query_file, "w", encoding="utf-8") as f:
            f.write(f">query\n{query_sequence}\n")

        db_prefix = os.path.join(BLAST_DB_PATH, "uniprot_db")
        if debug and not os.path.exists(db_prefix + ".pin"):
            print(f"{prefix} WARNING: BLAST DB seems missing: {db_prefix}.pin not found")

        # Preferred coverage signal: qcovs (query coverage per subject, reported by BLAST)
        # Fallback: compute qcov from qstart/qend span if qcovs is unsupported.
        outfmt_with_qcovs = "6 qseqid sseqid pident length qstart qend qlen qcovs evalue bitscore stitle"
        outfmt_span_only = "6 qseqid sseqid pident length qstart qend qlen evalue bitscore stitle"

        def _run_blast(outfmt: str) -> subprocess.CompletedProcess[str]:
            blast_cmd = [
                BLAST_EXECUTABLE,
                "-query", query_file,
                "-db", db_prefix,
                "-out", output_file,
                "-evalue", str(max_evalue),
                "-max_target_seqs", str(max_targets),
                "-outfmt", outfmt,
            ]
            return subprocess.run(blast_cmd, capture_output=True, text=True)

        use_qcovs = True
        result = _run_blast(outfmt_with_qcovs)
        if result.returncode != 0:
            stderr = (result.stderr or "").lower()
            if "qcovs" in stderr or "outfmt" in stderr or "unknown" in stderr or "unrecognized" in stderr:
                use_qcovs = False
                result = _run_blast(outfmt_span_only)

        if result.returncode != 0:
            _diag_inc("blast_error")
            if debug:
                stderr = (result.stderr or "").strip()
                print(f"{prefix} BLAST ERROR returncode={result.returncode} stderr={stderr[:400]}")
            return None

        if not os.path.exists(output_file):
            _diag_inc("blast_no_output")
            if debug:
                print(f"{prefix} No BLAST output file produced.")
            return None

        seen_any = False
        seen_pass_evalue = False
        seen_pass_qcov = False
        seen_pass_pident = False
        seen_pass_bitscore = False

        # Best hit per DISTINCT accession:
        # (reviewed, bitscore, evalue, qcov, pident, accession, title)
        best_per_accession: Dict[str, Tuple[int, float, float, float, float, str, str]] = {}

        # For debug on why nothing passed thresholds
        raw_top: List[Tuple[str, str, float, float, float, float, float]] = []  # (dbtag, acc, pident, qcov, evalue, bitscore, qcov_thr)

        with open(output_file, "r", encoding="utf-8") as f:
            for i, line in enumerate(f):
                cols = line.rstrip("\n").split("\t")
                if (use_qcovs and len(cols) < 11) or ((not use_qcovs) and len(cols) < 10):
                    continue

                sseqid = cols[1].strip()
                dbtag, accession = _extract_dbtag_and_accession(sseqid)

                try:
                    pident = float(cols[2])
                    qstart = int(float(cols[4]))
                    qend = int(float(cols[5]))
                    qlen = int(float(cols[6]))
                    evalue = float(cols[8] if use_qcovs else cols[7])
                    bitscore = float(cols[9] if use_qcovs else cols[8])
                    title = cols[10] if use_qcovs else cols[9]
                except Exception:
                    continue

                if qlen <= 0:
                    continue

                if use_qcovs:
                    try:
                        qcov = float(cols[7]) / 100.0
                    except Exception:
                        qcov = max(0.0, min(1.0, abs(qend - qstart) + 1) / float(qlen))
                else:
                    qcov = max(0.0, min(1.0, abs(qend - qstart) + 1) / float(qlen))

                seen_any = True
                reviewed = 1 if dbtag == "sp" else 0

                if i < 5:
                    raw_top.append((dbtag, accession, pident, qcov, evalue, bitscore, min_query_coverage))

                if evalue <= max_evalue:
                    seen_pass_evalue = True
                if qcov >= min_query_coverage:
                    seen_pass_qcov = True
                if pident >= min_pident:
                    seen_pass_pident = True
                if bitscore >= min_bitscore:
                    seen_pass_bitscore = True

                if evalue > max_evalue:
                    continue

                # qcov slack for clearly strong hits
                qcov_thr = min_query_coverage
                if bitscore >= (min_bitscore + strong_bitscore_bonus) and pident >= (min_pident + strong_pident_bonus):
                    qcov_thr = max(0.0, qcov_thr - qcov_slack)

                if qcov < qcov_thr:
                    continue
                if pident < min_pident:
                    continue
                if bitscore < min_bitscore:
                    continue

                prev = best_per_accession.get(accession)
                cur = (reviewed, bitscore, evalue, qcov, pident, accession, title)
                if prev is None:
                    best_per_accession[accession] = cur
                else:
                    # Keep the better one (prefer reviewed, higher bitscore, lower evalue, higher qcov/pident)
                    if cur[:5] > prev[:5]:
                        best_per_accession[accession] = cur

        if not seen_any:
            _diag_inc("fail_no_hits")
            if debug:
                print(f"{prefix} NO_HITS: BLAST returned no rows.")
            return None

        if not best_per_accession:
            # Diagnostics: understand what failed
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
                for dbtag, acc, pid, qcov, ev, bs, qthr in raw_top:
                    print(f"{prefix}   top: {dbtag}|{acc} pid={pid:.1f} qcov={qcov:.2f} ev={ev:.2g} bs={bs:.1f} (qthr={qthr:.2f})")
            return None

        # Rank hits by (reviewed, bitscore, -evalue, qcov, pident)
        ranked = sorted(best_per_accession.values(), key=lambda x: (x[0], x[1], -x[2], x[3], x[4]), reverse=True)
        best = ranked[0]

        best_reviewed, best_bitscore, best_evalue, best_qcov, best_pident, best_acc, best_title = best

        # Compute dynamic bitscore margin
        dynamic_margin = best_bitscore * margin_fraction_of_best
        dynamic_margin = max(margin_floor, min(min_bitscore_margin_cap, dynamic_margin))

        # Relax margin slightly if best is reviewed or clearly better in coverage/identity
        relax = 0.0
        if best_reviewed == 1:
            relax += margin_relax
        # If best is significantly higher in coverage/identity than runner-up, relax a bit
        if len(ranked) > 1:
            rb = ranked[1]
            if (best_qcov - rb[3]) >= 0.10 or (best_pident - rb[4]) >= 10.0:
                relax += margin_relax

        effective_margin = max(min_margin_after_relax, dynamic_margin - relax)

        # Ambiguity rule: reject if there is another accession close in bitscore
        ambiguous = False
        ties = []
        for cand in ranked[1:]:
            _, cand_bitscore, _, cand_qcov, cand_pident, cand_acc, cand_title = cand
            if (best_bitscore - cand_bitscore) <= effective_margin:
                ambiguous = True
                ties.append((cand_acc, cand_bitscore, cand_qcov, cand_pident, cand_title))
            else:
                break

        # Escape hatch: accept perfect ties for longer queries if both are near-perfect
        if ambiguous and qlen0 >= allow_perfect_tie_if_qlen_at_least and ties:
            # Identify top contender among ties (prefer reviewed)
            contender = ranked[1]
            cont_reviewed, cont_bitscore, cont_evalue, cont_qcov, cont_pident, cont_acc, cont_title = contender

            if (
                best_pident >= perfect_min_pident
                and best_qcov >= perfect_min_qcov
                and cont_pident >= perfect_min_pident
                and cont_qcov >= perfect_min_qcov
            ):
                # deterministically pick reviewed if different
                if best_reviewed != cont_reviewed:
                    chosen = best if best_reviewed > cont_reviewed else contender
                else:
                    # otherwise pick best bitscore, then lexicographic accession
                    if best_bitscore != cont_bitscore:
                        chosen = best if best_bitscore > cont_bitscore else contender
                    else:
                        chosen = best if best_acc <= cont_acc else contender

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
                for acc, bs, qc, pid, tt in ties[:5]:
                    print(f"{prefix}   tie: {acc} bs={bs:.1f} qcov={qc:.2f} pid={pid:.1f}")
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
        )

    except Exception as e:
        _diag_inc("blast_exception")
        if debug:
            print(f"{prefix} EXCEPTION: {e}")
        return None

    finally:
        for fp in (query_file, output_file):
            if os.path.exists(fp):
                try:
                    os.remove(fp)
                except Exception:
                    pass


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

    dna_res = {"DA", "DT", "DG", "DC", "DI"}
    rna_res = {"A", "U", "G", "C", "I"}

    for r in residues:
        res = str(r.get("residue_name") or r.get("resname") or "").strip().upper()
        if not res:
            continue

        res = MODRES_3TO3.get(res, res)

        if res in AMINO_ACIDS:
            aa_count += 1
        elif res in dna_res or res in rna_res:
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
        if any(str(rr.get("residue_name") or rr.get("resname") or "").strip().upper() in dna_res for rr in residues):
            if any(str(rr.get("residue_name") or rr.get("resname") or "").strip().upper() in rna_res for rr in residues):
                return "DNA/RNA"
            return "DNA"
        if any(str(rr.get("residue_name") or rr.get("resname") or "").strip().upper() in rna_res for rr in residues):
            return "RNA"
        return "Nucleic Acid"

    return "Unknown"


def parallel_blast_search(parsed_data: List[Dict[str, Any]], max_workers: int = 16) -> None:
    """Run BLASTP in parallel for all eligible chains and update them in-place."""
    debug = _DEBUG_BLAST

    _diag_reset()

    max_x_fraction: float = 0.20
    na_types = {"Nucleic Acid", "DNA", "RNA", "DNA/RNA"}

    # Dedupe by sequence: seq_key -> { "qlen": int, "seq": str, "targets": [(chain_dict, label), ...] }
    seq_groups: Dict[str, Dict[str, Any]] = {}

    def _apply_result(chain: Dict[str, Any], label: str, hit: Optional[BlastHit]) -> None:
        if hit:
            chain["uniprot_id"] = hit.accession
            chain["molecule_type"] = "Protein"

            name = _protein_name_from_uniprot_title(hit.title) if hit.title else None
            if name:
                chain["molecule_name"] = name
            else:
                chain["molecule_name"] = f"Matched UniProt: {hit.accession}"

            if debug:
                src = "sp" if hit.reviewed else "tr/unk"
                print(
                    f"[BLASTRES][{label}] HIT: {src}|{hit.accession} qcov={hit.qcov:.2f} pid={hit.pident:.1f} → molecule_type=Protein"
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

    # Run BLAST only for unique, uncached sequences
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures: List[Tuple[Any, str, List[Tuple[Dict[str, Any], str]]]] = []

        for _qlen, seq_key, seq, targets in to_blast:
            label0 = targets[0][1] if targets else ""
            label_call = f"{label0} (+{len(targets) - 1} dup)" if len(targets) > 1 else label0
            futures.append((executor.submit(run_blast_search, seq, label_call, debug), seq_key, targets))

        for future, seq_key, targets in futures:
            hit = future.result()
            _cache_put(seq_key, hit)

            for chain, label in targets:
                _apply_result(chain, label, hit)

    _diag_print_summary()
