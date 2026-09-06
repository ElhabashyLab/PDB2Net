"""Smoke-test the BLAST SQLite cache persistence."""

from __future__ import annotations

import tempfile
from pathlib import Path

from pdb2net import uniprot_matcher as matcher


def _reset_cache(path: Path) -> None:
    if matcher._CACHE_CONN is not None:
        matcher._CACHE_CONN.close()
    matcher._CACHE_CONN = None
    matcher._CACHE_DB_SIG = None
    matcher.config["blast_cache_path"] = str(path)


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="pdb2net-blast-cache-") as temp_dir:
        cache_path = Path(temp_dir) / "blast_cache.sqlite3"
        matcher._CACHE_ENABLED = True

        _reset_cache(cache_path)
        hit = matcher.BlastHit(
            accession="P12345",
            reviewed=True,
            bitscore=123.4,
            evalue=1e-30,
            qcov=0.95,
            pident=88.8,
            title="sp|P12345| Example protein",
            matched_id="P12345",
            representative_accession="P12345",
        )
        matcher._cache_put("positive-seq", hit)
        matcher._cache_put("negative-seq", None)

        # Reopen to prove entries were committed to disk, not only kept in memory.
        _reset_cache(cache_path)
        positive_cached, positive_hit = matcher._cache_get("positive-seq")
        negative_cached, negative_hit = matcher._cache_get("negative-seq")

        if not positive_cached or positive_hit != hit:
            raise AssertionError(f"positive cache entry mismatch: {positive_cached=} {positive_hit=}")
        if not negative_cached or negative_hit is not None:
            raise AssertionError(f"negative cache entry mismatch: {negative_cached=} {negative_hit=}")

        if not cache_path.is_file():
            raise AssertionError("cache was not created in the temporary test directory")
        if matcher._CACHE_CONN is not None:
            matcher._CACHE_CONN.close()
            matcher._CACHE_CONN = None

    print("blast_cache_persistence: ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
