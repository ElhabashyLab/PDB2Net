"""Portable, deterministic filenames for generated artifacts."""

from __future__ import annotations

import hashlib
import re

PORTABLE_ARTIFACT_STEM_SEMANTICS_ID = "pdb2net-portable-artifact-stem-v1"
MAX_ARTIFACT_STEM_BYTES = 180

_PORTABLE_STEM = re.compile(r"[A-Za-z0-9._-]+")
_WINDOWS_RESERVED = {
    "CON",
    "PRN",
    "AUX",
    "NUL",
    *(f"COM{index}" for index in range(1, 10)),
    *(f"LPT{index}" for index in range(1, 10)),
}


def _is_portable(stem: str) -> bool:
    basename = stem.split(".", 1)[0].upper()
    return (
        stem not in {"", ".", ".."}
        and len(stem.encode("utf-8")) <= MAX_ARTIFACT_STEM_BYTES
        and _PORTABLE_STEM.fullmatch(stem) is not None
        and basename not in _WINDOWS_RESERVED
    )


def portable_artifact_stem(value: object) -> str:
    """Preserve ordinary stems and bound unsafe ones with a stable digest."""
    raw = str(value or "").strip()
    if _is_portable(raw):
        return raw

    readable = re.sub(r"[^A-Za-z0-9._-]+", "_", raw)
    readable = re.sub(r"_+", "_", readable).strip("._-") or "artifact"
    digest_suffix = f"--{hashlib.sha256(raw.encode('utf-8')).hexdigest()[:16]}"
    maximum_prefix = MAX_ARTIFACT_STEM_BYTES - len(digest_suffix)
    prefix = readable[:maximum_prefix].rstrip("._-") or "artifact"
    return f"{prefix}{digest_suffix}"


__all__ = [
    "MAX_ARTIFACT_STEM_BYTES",
    "PORTABLE_ARTIFACT_STEM_SEMANTICS_ID",
    "portable_artifact_stem",
]
