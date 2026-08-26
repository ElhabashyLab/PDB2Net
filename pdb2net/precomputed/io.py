"""Bounded and atomic filesystem I/O for schema-3 stores."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
from pathlib import Path
import stat
import tempfile
from typing import Any, Mapping

from ..input_contract import InputValidationError
from ..structure_identity import pdb_archive_shard
from .schema import (
    MAX_COMPRESSED_ENTRY_BYTES,
    MAX_DECOMPRESSED_ENTRY_BYTES,
    SOURCE_SCOPE,
    canonical_json,
    normalize_pdb_id,
)


MAX_MANIFEST_BYTES = 1_000_000


def store_root(store: Path | str, *, create: bool = False) -> Path:
    raw = Path(store).expanduser()
    if raw.is_symlink():
        raise InputValidationError(
            "UNSAFE_PRECOMPUTED_PATH", "Precomputed store root must not be a symlink."
        )
    if raw.exists() and not raw.is_dir():
        raise InputValidationError(
            "UNSAFE_PRECOMPUTED_PATH", "Precomputed store root must be a directory."
        )
    if create:
        raw.mkdir(parents=True, exist_ok=True)
    root = raw.resolve()
    if create and (root.is_symlink() or not root.is_dir()):
        raise InputValidationError(
            "UNSAFE_PRECOMPUTED_PATH", "Precomputed store root is unsafe."
        )
    return root


def _assert_safe_target(root: Path, target: Path) -> None:
    try:
        relative = target.relative_to(root)
    except ValueError as exc:
        raise InputValidationError(
            "UNSAFE_PRECOMPUTED_PATH", "Precomputed path escaped its store root."
        ) from exc
    current = root
    if current.is_symlink():
        raise InputValidationError(
            "UNSAFE_PRECOMPUTED_PATH", "Precomputed store root must not be a symlink."
        )
    for part in relative.parts:
        current /= part
        if current.is_symlink():
            raise InputValidationError(
                "UNSAFE_PRECOMPUTED_PATH", "Precomputed path contains a symlink."
            )


def entry_path(store: Path | str, pdb_id: object) -> Path:
    normalized = normalize_pdb_id(pdb_id)
    root = store_root(store)
    path = root / "entries" / pdb_archive_shard(normalized) / f"{normalized}.json.gz"
    _assert_safe_target(root, path)
    return path


def manifest_path(store: Path | str) -> Path:
    root = store_root(store)
    path = root / "manifest.json"
    _assert_safe_target(root, path)
    return path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def source_fingerprint(path: Path) -> dict[str, Any]:
    try:
        metadata = path.lstat()
    except OSError as exc:
        raise InputValidationError(
            "INVALID_PRECOMPUTE_SOURCE", f"Cannot inspect precompute source: {path}"
        ) from exc
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISREG(metadata.st_mode):
        raise InputValidationError(
            "INVALID_PRECOMPUTE_SOURCE", f"Precompute source is not a regular file: {path}"
        )
    return {
        "sha256": sha256_file(path),
        "size_bytes": metadata.st_size,
        "scope": SOURCE_SCOPE,
    }


def atomic_write(path: Path, data: bytes, *, root: Path) -> None:
    _assert_safe_target(root, path.parent)
    path.parent.mkdir(parents=True, exist_ok=True)
    _assert_safe_target(root, path.parent)
    if path.is_symlink() or (path.exists() and not path.is_file()):
        raise InputValidationError(
            "UNSAFE_PRECOMPUTED_PATH", f"Precomputed target is unsafe: {path.name}"
        )
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            if hasattr(os, "fchmod"):
                os.fchmod(handle.fileno(), 0o640)
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        flags = os.O_RDONLY
        if hasattr(os, "O_DIRECTORY"):
            flags |= os.O_DIRECTORY
        directory_descriptor = os.open(path.parent, flags)
        try:
            os.fsync(directory_descriptor)
        finally:
            os.close(directory_descriptor)
    finally:
        if temporary.exists():
            temporary.unlink()


def write_entry_json(store: Path | str, pdb_id: object, payload: Mapping[str, Any]) -> Path:
    root = store_root(store, create=True)
    serialized = canonical_json(payload)
    normalized = normalize_pdb_id(pdb_id)
    if len(serialized) > MAX_DECOMPRESSED_ENTRY_BYTES:
        raise InputValidationError(
            "PRECOMPUTED_ENTRY_TOO_LARGE",
            f"Expanded precomputed entry for {normalized.upper()} exceeds the limit.",
        )
    compressed = gzip.compress(serialized, compresslevel=6, mtime=0)
    if len(compressed) > MAX_COMPRESSED_ENTRY_BYTES:
        raise InputValidationError(
            "PRECOMPUTED_ENTRY_TOO_LARGE",
            f"Compressed precomputed entry for {normalized.upper()} exceeds the limit.",
        )
    path = entry_path(root, normalized)
    atomic_write(path, compressed, root=root)
    return path


def read_entry_json(
    store: Path | str,
    pdb_id: object,
    *,
    expanded_budget: int = MAX_DECOMPRESSED_ENTRY_BYTES,
) -> tuple[Any, int]:
    root = store_root(store)
    path = entry_path(root, pdb_id)
    _assert_safe_target(root, path)
    try:
        metadata = path.lstat()
    except FileNotFoundError as exc:
        normalized = normalize_pdb_id(pdb_id)
        raise InputValidationError(
            "PRECOMPUTED_ENTRY_MISSING",
            f"No precomputed entry exists for PDB ID {normalized.upper()}.",
        ) from exc
    except OSError as exc:
        raise InputValidationError(
            "UNSAFE_PRECOMPUTED_ENTRY", f"Cannot inspect precomputed entry: {path.name}"
        ) from exc
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISREG(metadata.st_mode):
        raise InputValidationError(
            "UNSAFE_PRECOMPUTED_ENTRY",
            f"Precomputed entry is not a regular file: {path.name}",
        )
    if metadata.st_size > MAX_COMPRESSED_ENTRY_BYTES:
        raise InputValidationError(
            "PRECOMPUTED_ENTRY_TOO_LARGE",
            f"Compressed precomputed entry exceeds the limit: {path.name}",
        )
    read_limit = min(MAX_DECOMPRESSED_ENTRY_BYTES, max(0, int(expanded_budget)))
    if read_limit == 0:
        raise InputValidationError(
            "PRECOMPUTED_REQUEST_LIMIT_EXCEEDED",
            "Requested entries exceed the aggregate expanded-size limit.",
        )
    try:
        with gzip.open(path, "rb") as handle:
            raw = handle.read(read_limit + 1)
    except (OSError, EOFError) as exc:
        raise InputValidationError(
            "CORRUPT_PRECOMPUTED_ENTRY", f"Cannot decompress entry: {path.name}"
        ) from exc
    if len(raw) > read_limit:
        code = (
            "PRECOMPUTED_ENTRY_TOO_LARGE"
            if read_limit == MAX_DECOMPRESSED_ENTRY_BYTES
            else "PRECOMPUTED_REQUEST_LIMIT_EXCEEDED"
        )
        raise InputValidationError(code, f"Expanded entry exceeds the limit: {path.name}")
    try:
        return json.loads(raw.decode("utf-8")), len(raw)
    except (UnicodeError, json.JSONDecodeError) as exc:
        raise InputValidationError(
            "CORRUPT_PRECOMPUTED_ENTRY", f"Entry is not valid UTF-8 JSON: {path.name}"
        ) from exc


def write_manifest_json(store: Path | str, payload: Mapping[str, Any]) -> Path:
    root = store_root(store, create=True)
    serialized = json.dumps(
        payload, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False
    ).encode("utf-8") + b"\n"
    if len(serialized) > MAX_MANIFEST_BYTES:
        raise InputValidationError(
            "PRECOMPUTED_MANIFEST_TOO_LARGE", "Precomputed manifest exceeds 1 MiB."
        )
    path = manifest_path(root)
    atomic_write(path, serialized, root=root)
    return path


def read_manifest_json(store: Path | str) -> Any:
    root = store_root(store)
    path = manifest_path(root)
    try:
        metadata = path.lstat()
    except FileNotFoundError as exc:
        raise InputValidationError(
            "PRECOMPUTED_STORE_UNPUBLISHED",
            "Precomputed store has no published manifest.",
        ) from exc
    except OSError as exc:
        raise InputValidationError(
            "CORRUPT_PRECOMPUTED_MANIFEST", "Cannot inspect precomputed manifest."
        ) from exc
    if (
        stat.S_ISLNK(metadata.st_mode)
        or not stat.S_ISREG(metadata.st_mode)
        or metadata.st_size > MAX_MANIFEST_BYTES
    ):
        raise InputValidationError(
            "UNSAFE_PRECOMPUTED_MANIFEST", "Precomputed manifest is unsafe."
        )
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise InputValidationError(
            "CORRUPT_PRECOMPUTED_MANIFEST", "Precomputed manifest is invalid JSON."
        ) from exc


__all__ = [
    "atomic_write",
    "entry_path",
    "manifest_path",
    "read_entry_json",
    "read_manifest_json",
    "sha256_file",
    "source_fingerprint",
    "store_root",
    "write_entry_json",
    "write_manifest_json",
]
