"""Validated PDB/mmCIF transport parsing and stable structure identities."""

from __future__ import annotations

import csv
from contextlib import contextmanager
import hashlib
import os
import re
import stat
import zlib
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional

import gemmi

from .config_loader import config
from .input_contract import InputValidationError
from .network_annotations import extract_embedded_annotations
from .reference_data import load_valid_pdb_ids as _load_valid_pdb_ids
from .structure_identity import (
    EXTENDED_PDB_ID_RE,
    LEGACY_PDB_ID_RE,
    StructureIdentity,
    identity_from_official_id,
)

CORE_STRUCTURE_EXTENSIONS = (".pdb", ".ent", ".cif", ".mmcif")
ALLOWED_EXTENSIONS = set(CORE_STRUCTURE_EXTENSIONS)
GZIP_MAGIC = b"\x1f\x8b"
IO_CHUNK_BYTES = 1024 * 1024
FileSignature = tuple[int, int, int, int, int]


def load_valid_pdb_ids() -> set[str]:
    """Compatibility loader; identity detection never depends on this list."""
    try:
        return _load_valid_pdb_ids(str(config.get("pdb_fasta_path") or ""))
    except Exception as exc:
        print(f"Error loading pdb_seqres.txt: {exc}")
        return set()


VALID_PDB_IDS: set[str] | None = None


def _valid_pdb_ids() -> set[str]:
    global VALID_PDB_IDS
    if VALID_PDB_IDS is None:
        VALID_PDB_IDS = load_valid_pdb_ids()
    return VALID_PDB_IDS


def _structure_extension(file_path: str) -> str:
    filename = os.path.basename(file_path)
    if filename.lower().endswith(".gz"):
        filename = filename[:-3]
    return os.path.splitext(filename)[1].lower()


def _structure_stem(file_path: str) -> str:
    filename = os.path.basename(file_path)
    if filename.lower().endswith(".gz"):
        filename = filename[:-3]
    extension = os.path.splitext(filename)[1]
    return filename[: -len(extension)] if extension else filename


def is_valid_file(file_path: str) -> bool:
    """Return whether a filename has one supported structure suffix."""
    return _structure_extension(file_path) in ALLOWED_EXTENSIONS


def input_file_signature(file_path: str | os.PathLike[str]) -> FileSignature:
    """Return the immutable fields used to detect symlinks and input mutation."""
    path = os.fspath(file_path)
    try:
        result = os.lstat(path)
    except OSError as exc:
        raise InputValidationError("INPUT_FILE_STAT_FAILED", f"Cannot inspect input file: {path}") from exc
    if stat.S_ISLNK(result.st_mode):
        raise InputValidationError(
            "SYMLINK_INPUT_NOT_ALLOWED",
            f"Structure inputs must not be symbolic links: {Path(path).name}",
        )
    if not stat.S_ISREG(result.st_mode):
        raise InputValidationError(
            "INPUT_FILE_NOT_REGULAR",
            f"Structure input is not a regular file: {Path(path).name}",
        )
    return (
        result.st_dev,
        result.st_ino,
        result.st_size,
        result.st_mtime_ns,
        result.st_ctime_ns,
    )


def _raise_changed(file_path: str) -> None:
    raise InputValidationError(
        "INPUT_CHANGED_DURING_PROCESSING",
        f"Input file changed while it was being processed: {Path(file_path).name}",
    )


def _stat_signature(result: os.stat_result) -> FileSignature:
    return (
        result.st_dev,
        result.st_ino,
        result.st_size,
        result.st_mtime_ns,
        result.st_ctime_ns,
    )


@contextmanager
def _open_regular_input(file_path: str, expected_signature: FileSignature):
    """Open an inspected input without following a replacement symlink."""
    flags = os.O_RDONLY | getattr(os, "O_BINARY", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(file_path, flags)
    except OSError as exc:
        try:
            current = os.lstat(file_path)
        except OSError:
            current = None
        if current is not None and stat.S_ISLNK(current.st_mode):
            raise InputValidationError(
                "SYMLINK_INPUT_NOT_ALLOWED",
                f"Structure inputs must not be symbolic links: {Path(file_path).name}",
            ) from exc
        raise InputValidationError(
            "INPUT_FILE_READ_FAILED", f"Cannot read structure input: {Path(file_path).name}"
        ) from exc
    try:
        opened = os.fstat(descriptor)
        if not stat.S_ISREG(opened.st_mode) or _stat_signature(opened) != expected_signature:
            _raise_changed(file_path)
        with os.fdopen(descriptor, "rb", closefd=False) as handle:
            yield handle
    finally:
        os.close(descriptor)


def _append_expanded(
    decoded: bytes,
    *,
    output: list[bytes] | None,
    prefix: bytearray,
    total: int,
    maximum: int | None,
    file_path: str,
) -> int:
    if decoded and len(prefix) < 2:
        prefix.extend(decoded[: 2 - len(prefix)])
    updated = total + len(decoded)
    if maximum is not None and updated > maximum:
        raise InputValidationError(
            "INPUT_FILE_EXPANDED_BYTES_LIMIT_EXCEEDED",
            f"Expanded input file {Path(file_path).name} exceeds the configured maximum of {maximum} bytes.",
        )
    if output is not None and decoded:
        output.append(decoded)
    return updated


def _read_validated_gzip(
    file_path: str,
    *,
    expected_signature: FileSignature,
    maximum: int | None,
    collect: bool,
) -> tuple[bytes, int, str]:
    """Inflate every gzip member, checking every trailer and a shared limit."""
    output: list[bytes] | None = [] if collect else None
    prefix = bytearray()
    total = 0
    pending = b""
    raw_eof = False
    member_count = 0
    raw_digest = hashlib.sha256()

    try:
        with _open_regular_input(file_path, expected_signature) as raw:
            while True:
                while len(pending) < 2 and not raw_eof:
                    chunk = raw.read(IO_CHUNK_BYTES)
                    if chunk:
                        raw_digest.update(chunk)
                        pending += chunk
                    else:
                        raw_eof = True
                if not pending:
                    break
                if len(pending) < 2 or not pending.startswith(GZIP_MAGIC):
                    raise InputValidationError(
                        "INVALID_GZIP_INPUT",
                        f"Gzip input contains trailing data or an invalid member: {Path(file_path).name}",
                    )

                member_count += 1
                inflater = zlib.decompressobj(16 + zlib.MAX_WBITS)
                while not inflater.eof:
                    if not pending:
                        chunk = raw.read(IO_CHUNK_BYTES)
                        if not chunk:
                            raw_eof = True
                            raise InputValidationError(
                                "INVALID_GZIP_INPUT",
                                f"Gzip input is truncated or has no valid trailer: {Path(file_path).name}",
                            )
                        raw_digest.update(chunk)
                        pending = chunk
                    before = pending
                    try:
                        decoded = inflater.decompress(pending, IO_CHUNK_BYTES)
                    except zlib.error as exc:
                        raise InputValidationError(
                            "INVALID_GZIP_INPUT",
                            f"Gzip input has an invalid header, payload, CRC, or trailer: {Path(file_path).name}",
                        ) from exc
                    if inflater.eof:
                        pending = inflater.unused_data
                    else:
                        pending = inflater.unconsumed_tail
                    total = _append_expanded(
                        decoded,
                        output=output,
                        prefix=prefix,
                        total=total,
                        maximum=maximum,
                        file_path=file_path,
                    )
                    if not decoded and pending == before:
                        raise InputValidationError(
                            "INVALID_GZIP_INPUT",
                            f"Gzip input could not be decoded: {Path(file_path).name}",
                        )

                try:
                    flushed = inflater.flush()
                except zlib.error as exc:
                    raise InputValidationError(
                        "INVALID_GZIP_INPUT",
                        f"Gzip input has an invalid trailer: {Path(file_path).name}",
                    ) from exc
                total = _append_expanded(
                    flushed,
                    output=output,
                    prefix=prefix,
                    total=total,
                    maximum=maximum,
                    file_path=file_path,
                )
    except InputValidationError:
        raise
    except OSError as exc:
        raise InputValidationError(
            "INPUT_FILE_READ_FAILED",
            f"Cannot read structure input: {Path(file_path).name}",
        ) from exc

    if member_count == 0 or total == 0:
        raise InputValidationError(
            "INVALID_GZIP_INPUT",
            f"Gzip input expands to an empty structure: {Path(file_path).name}",
        )
    if bytes(prefix).startswith(GZIP_MAGIC):
        raise InputValidationError(
            "NESTED_GZIP_INPUT",
            f"Nested gzip input is not supported: {Path(file_path).name}",
        )
    return (b"".join(output) if output is not None else b"", total, raw_digest.hexdigest())


def input_file_sha256(
    file_path: str | os.PathLike[str],
    *,
    expected_signature: FileSignature | None = None,
) -> str:
    """Hash one regular input while binding the read to its inspected metadata."""
    path = os.fspath(file_path)
    current_signature = input_file_signature(path)
    if expected_signature is not None and current_signature != expected_signature:
        _raise_changed(path)
    signature = expected_signature or current_signature
    digest = hashlib.sha256()
    try:
        with _open_regular_input(path, signature) as handle:
            while chunk := handle.read(IO_CHUNK_BYTES):
                digest.update(chunk)
    except OSError as exc:
        raise InputValidationError(
            "INPUT_FILE_READ_FAILED", f"Cannot read structure input: {Path(path).name}"
        ) from exc
    if input_file_signature(path) != signature:
        _raise_changed(path)
    return digest.hexdigest()


def read_validated_structure_bytes(
    file_path: str,
    *,
    expected_signature: FileSignature | None = None,
    expected_sha256: str | None = None,
    maximum_expanded_bytes: int | None = None,
    collect: bool = True,
    include_sha256: bool = False,
) -> tuple[bytes, int] | tuple[bytes, int, str]:
    """Validate suffix/magic and return the single expanded structure stream."""
    current_signature = input_file_signature(file_path)
    if expected_signature is not None and current_signature != expected_signature:
        _raise_changed(file_path)
    signature = expected_signature or current_signature
    if not is_valid_file(file_path):
        raise InputValidationError(
            "UNSUPPORTED_INPUT_SUFFIX",
            f"Unsupported structure filename: {Path(file_path).name}",
        )
    compressed_suffix = str(file_path).lower().endswith(".gz")
    try:
        with _open_regular_input(file_path, signature) as handle:
            magic = handle.read(2)
    except OSError as exc:
        raise InputValidationError(
            "INPUT_FILE_READ_FAILED", f"Cannot read structure input: {Path(file_path).name}"
        ) from exc
    if not magic:
        raise InputValidationError("EMPTY_INPUT_FILE", f"Structure input is empty: {Path(file_path).name}")
    has_gzip_magic = magic == GZIP_MAGIC
    if compressed_suffix != has_gzip_magic:
        raise InputValidationError(
            "INPUT_COMPRESSION_MISMATCH",
            f"Filename suffix and gzip magic bytes disagree: {Path(file_path).name}",
        )

    if compressed_suffix:
        payload, expanded, raw_sha256 = _read_validated_gzip(
            file_path,
            expected_signature=signature,
            maximum=maximum_expanded_bytes,
            collect=collect,
        )
    else:
        expanded = signature[2]
        if maximum_expanded_bytes is not None and expanded > maximum_expanded_bytes:
            raise InputValidationError(
                "INPUT_FILE_EXPANDED_BYTES_LIMIT_EXCEEDED",
                f"Expanded input file {Path(file_path).name} exceeds the configured maximum of "
                f"{maximum_expanded_bytes} bytes.",
            )
        if collect:
            try:
                with _open_regular_input(file_path, signature) as handle:
                    payload = handle.read()
            except OSError as exc:
                raise InputValidationError(
                    "INPUT_FILE_READ_FAILED", f"Cannot read structure input: {Path(file_path).name}"
                ) from exc
        else:
            payload = b""
        if collect:
            raw_sha256 = hashlib.sha256(payload).hexdigest()
        elif expected_sha256 is not None or include_sha256:
            raw_sha256 = input_file_sha256(file_path, expected_signature=signature)
        else:
            raw_sha256 = ""
    if expected_sha256 is not None and raw_sha256 != expected_sha256:
        _raise_changed(file_path)
    if input_file_signature(file_path) != signature:
        _raise_changed(file_path)
    if include_sha256:
        return payload, expanded, raw_sha256
    return payload, expanded


def extract_structure_identity_from_filename(file_path: str) -> StructureIdentity | None:
    """Recognize legacy, extended, archive, and NextGen PDB filenames."""
    stem = _structure_stem(file_path)
    extended_match = re.fullmatch(
        r"(pdb_[A-Za-z0-9]{8})(?:_xyz-enrich)?",
        stem,
        flags=re.IGNORECASE,
    )
    if extended_match:
        return identity_from_official_id(extended_match.group(1))
    legacy_match = re.fullmatch(r"(?:pdb)?([0-9][A-Za-z0-9]{3})", stem, flags=re.IGNORECASE)
    if legacy_match:
        return identity_from_official_id(legacy_match.group(1))
    return None


def extract_pdb_id_from_filename(file_path: str) -> Optional[str]:
    identity = extract_structure_identity_from_filename(file_path)
    return identity.display_id if identity else None


def _official_identity_candidate(value: object) -> StructureIdentity | None:
    candidate = str(value or "").strip().strip("'\"")
    if LEGACY_PDB_ID_RE.fullmatch(candidate) or EXTENDED_PDB_ID_RE.fullmatch(candidate):
        return identity_from_official_id(candidate)
    return None


def _clean_cif_value(value: object) -> str | None:
    text = str(value or "").strip().strip("'\"")
    return None if not text or text in {".", "?"} else text


def _content_identity_claims(
    *,
    extension: str,
    text: str,
    block: Any | None,
) -> list[dict[str, Any]]:
    claims: list[dict[str, Any]] = []

    def add(source: str, raw: object, *, require_official: bool = False) -> None:
        value = _clean_cif_value(raw)
        if value is None:
            return
        identity = _official_identity_candidate(value)
        if require_official and identity is None:
            raise InputValidationError(
                "INVALID_STRUCTURE_IDENTITY",
                f"{source} contains an invalid PDB identifier: {value!r}",
            )
        claims.append(
            {
                "source": source,
                "value": value,
                "identity": identity.as_dict() if identity else None,
            }
        )

    if extension in {".cif", ".mmcif"} and block is not None:
        add("data_block", getattr(block, "name", ""))
        try:
            add("_entry.id", block.find_value("_entry.id"))
        except Exception:
            pass
        category = block.get_mmcif_category("_database_2.") or {}
        row_count = max((len(values) for values in category.values()), default=0)
        for index in range(row_count):
            database_id = _clean_cif_value(
                category.get("database_id", [None] * row_count)[index]
                if index < len(category.get("database_id", []))
                else None
            )
            if str(database_id or "").upper() != "PDB":
                continue
            code_values = category.get("database_code", [])
            code = code_values[index] if index < len(code_values) else None
            add("_database_2[PDB].database_code", code, require_official=True)
    else:
        for line in text.splitlines():
            if not line.startswith("HEADER"):
                continue
            fixed = line[62:66].strip() if len(line) >= 66 else ""
            if fixed:
                add("PDB_HEADER", fixed, require_official=True)
                continue
            tokens = line.split()
            candidate = tokens[-1] if tokens and _official_identity_candidate(tokens[-1]) else ""
            if candidate:
                add("PDB_HEADER", candidate, require_official=True)
    return claims


def _local_identity(file_path: str, content_value: str | None = None) -> StructureIdentity:
    stem = (content_value or _structure_stem(file_path)).strip() or "STRUCTURE"
    source = "alphafold" if re.match(r"^AF[-_]", stem, flags=re.IGNORECASE) else "local"
    return StructureIdentity(source, stem.upper())


def resolve_structure_identity(
    file_path: str,
    *,
    text: str | None = None,
    block: Any | None = None,
    expected_signature: FileSignature | None = None,
    expected_sha256: str | None = None,
    maximum_expanded_bytes: int | None = None,
) -> tuple[StructureIdentity, list[dict[str, Any]], list[dict[str, Any]]]:
    """Resolve every content claim, rejecting conflicts and warning on filenames."""
    extension = _structure_extension(file_path)
    if text is None:
        payload, _ = read_validated_structure_bytes(
            file_path,
            expected_signature=expected_signature,
            expected_sha256=expected_sha256,
            maximum_expanded_bytes=maximum_expanded_bytes,
        )
        try:
            text = payload.decode("utf-8-sig")
        except UnicodeDecodeError as exc:
            raise InputValidationError(
                "INVALID_STRUCTURE_TEXT_ENCODING",
                f"Structure input is not valid UTF-8 text: {Path(file_path).name}",
            ) from exc
        if extension in {".cif", ".mmcif"}:
            try:
                document = gemmi.cif.read_string(text)
            except Exception:
                document = None
            if document is not None and len(document) != 1:
                raise InputValidationError(
                    "INVALID_MMCIF_DATA_BLOCK_COUNT",
                    f"mmCIF input must contain exactly one data block: {Path(file_path).name}",
                )
            if document is not None:
                block = document[0]

    claims = _content_identity_claims(extension=extension, text=text or "", block=block)
    official: dict[str, StructureIdentity] = {}
    for claim in claims:
        raw_identity = claim.get("identity")
        if isinstance(raw_identity, Mapping):
            identity = StructureIdentity.from_mapping(raw_identity)
            official[identity.key] = identity
    if len(official) > 1:
        details = ", ".join(f"{claim['source']}={claim['value']}" for claim in claims)
        raise InputValidationError(
            "CONFLICTING_STRUCTURE_IDENTITY",
            f"Structure content contains conflicting identity claims ({details}).",
        )

    content_identity = next(iter(official.values()), None)
    local_claims = [
        str(claim["value"])
        for claim in claims
        if claim.get("source") in {"data_block", "_entry.id"} and not claim.get("identity")
    ]
    normalized_local = {value.casefold() for value in local_claims}
    if content_identity is not None and normalized_local:
        details = ", ".join(f"{claim['source']}={claim['value']}" for claim in claims)
        raise InputValidationError(
            "CONFLICTING_STRUCTURE_IDENTITY",
            f"Structure content contains conflicting identity claims ({details}).",
        )
    if len(normalized_local) > 1:
        details = ", ".join(f"{claim['source']}={claim['value']}" for claim in claims)
        raise InputValidationError(
            "CONFLICTING_STRUCTURE_IDENTITY",
            f"Structure content contains conflicting identity claims ({details}).",
        )
    if content_identity is None and local_claims:
        content_identity = _local_identity(file_path, local_claims[0])

    filename_identity = extract_structure_identity_from_filename(file_path) or _local_identity(file_path)
    claims.append(
        {
            "source": "filename",
            "value": _structure_stem(file_path),
            "identity": filename_identity.as_dict(),
        }
    )
    warnings: list[dict[str, Any]] = []
    if content_identity is not None:
        identity = content_identity
        if filename_identity.key != identity.key:
            warnings.append(
                {
                    "code": "STRUCTURE_FILENAME_ID_MISMATCH",
                    "message": (
                        f"Content identity {identity.display_id} takes precedence over conflicting "
                        f"filename identity {filename_identity.display_id} for {Path(file_path).name}."
                    ),
                }
            )
    else:
        identity = filename_identity
    return identity, warnings, claims


def extract_structure_identity_from_file(file_path: str) -> StructureIdentity | None:
    """Return a content identity, without falling back to the filename."""
    identity, _warnings, claims = resolve_structure_identity(file_path)
    content_claims = [claim for claim in claims if claim.get("source") != "filename"]
    has_content_identity = bool(content_claims)
    return identity if has_content_identity else None


def extract_pdb_id_from_file(file_path: str) -> Optional[str]:
    identity = extract_structure_identity_from_file(file_path)
    return identity.display_id if identity and identity.source == "pdb" else None


def get_structure_identity(
    file_path: str,
    *,
    expected_signature: FileSignature | None = None,
    expected_sha256: str | None = None,
    maximum_expanded_bytes: int | None = None,
) -> StructureIdentity:
    identity, _warnings, _claims = resolve_structure_identity(
        file_path,
        expected_signature=expected_signature,
        expected_sha256=expected_sha256,
        maximum_expanded_bytes=maximum_expanded_bytes,
    )
    if identity.source != "pdb" and not extract_structure_identity_from_filename(file_path):
        print(f"Info: No canonical PDB ID found. Using filename as structure ID: {identity.display_id}")
    return identity


def get_pdb_id(file_path: str) -> str:
    return get_structure_identity(file_path).display_id


def parse_structure_input(
    file_path: str,
    *,
    expected_signature: FileSignature | None = None,
    expected_sha256: str | None = None,
    maximum_expanded_bytes: int | None = None,
) -> Dict[str, Any]:
    """Parse exactly one validated structure and retain enriched-mmCIF data."""
    before = input_file_signature(file_path)
    if expected_signature is not None and before != expected_signature:
        _raise_changed(file_path)
    payload, _expanded = read_validated_structure_bytes(
        file_path,
        expected_signature=expected_signature or before,
        expected_sha256=expected_sha256,
        maximum_expanded_bytes=maximum_expanded_bytes,
    )
    try:
        text = payload.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise InputValidationError(
            "INVALID_STRUCTURE_TEXT_ENCODING",
            f"Structure input is not valid UTF-8 text: {Path(file_path).name}",
        ) from exc

    extension = _structure_extension(file_path)
    embedded: Dict[str, Any] = {"is_enriched": False, "by_chain": {}, "counts": {}, "warnings": []}
    block = None
    try:
        if extension in {".cif", ".mmcif"}:
            document = gemmi.cif.read_string(text)
            if len(document) != 1:
                raise InputValidationError(
                    "INVALID_MMCIF_DATA_BLOCK_COUNT",
                    f"mmCIF input must contain exactly one data block: {Path(file_path).name}",
                )
            block = document[0]
            embedded = extract_embedded_annotations(block)
            structure = gemmi.make_structure_from_block(block)
            input_format = "mmcif"
        else:
            structure = gemmi.read_pdb_string(text)
            input_format = "pdb"
    except InputValidationError:
        raise
    except Exception as exc:
        code = "INVALID_MMCIF_INPUT" if extension in {".cif", ".mmcif"} else "INVALID_PDB_INPUT"
        raise InputValidationError(code, f"Cannot parse structure input: {Path(file_path).name}") from exc

    identity, warnings, claims = resolve_structure_identity(file_path, text=text, block=block)
    structure.setup_entities()
    after = input_file_signature(file_path)
    if after != before or (expected_signature is not None and after != expected_signature):
        _raise_changed(file_path)
    return {
        "structure": structure,
        "structure_identity": identity.as_dict(),
        "identity_warnings": warnings,
        "identity_claims": claims,
        "embedded_annotations": embedded,
        "input_format": input_format,
        "input_kind": "nextgen_enriched_mmcif" if embedded.get("is_enriched") else input_format,
    }


def parse_structure(file_path: str, pdb_id: str) -> Optional[gemmi.Structure]:
    """Compatibility wrapper returning ``None`` for an invalid structure."""
    del pdb_id
    try:
        return parse_structure_input(file_path)["structure"]
    except Exception as exc:
        print(f"Error parsing {file_path}: {exc}")
        return None


def read_files_from_csv(csv_path: str) -> List[Dict[str, Any]]:
    with open(csv_path, newline="", encoding="utf-8") as csvfile:
        reader = csv.DictReader(csvfile)
        if "file_path" not in (reader.fieldnames or []):
            raise ValueError("CSV file must contain a 'file_path' column.")
        file_paths = [row["file_path"] for row in reader]
    return _read_compatibility_files(file_paths)


def _read_compatibility_files(file_paths: List[str]) -> List[Dict[str, Any]]:
    structures: List[Dict[str, Any]] = []
    for file_path in file_paths:
        if not is_valid_file(file_path):
            continue
        try:
            parsed = parse_structure_input(file_path)
        except Exception as exc:
            print(f"Error parsing {file_path}: {exc}")
            continue
        identity = StructureIdentity.from_mapping(parsed["structure_identity"])
        structures.append(
            {"file_path": file_path, "pdb_id": identity.display_id, "structure": parsed["structure"]}
        )
    return structures


def read_files_from_folder(folder_path: str) -> List[Dict[str, Any]]:
    file_paths: List[str] = []
    with os.scandir(folder_path) as entries:
        for entry in entries:
            if entry.is_symlink():
                raise InputValidationError(
                    "SYMLINK_INPUT_NOT_ALLOWED", f"Structure inputs must not be symbolic links: {entry.name}"
                )
            if entry.is_file(follow_symlinks=False) and is_valid_file(entry.path):
                file_paths.append(entry.path)
    return _read_compatibility_files(sorted(file_paths))


__all__ = [
    "ALLOWED_EXTENSIONS",
    "FileSignature",
    "extract_pdb_id_from_file",
    "extract_pdb_id_from_filename",
    "extract_structure_identity_from_file",
    "extract_structure_identity_from_filename",
    "get_pdb_id",
    "get_structure_identity",
    "input_file_signature",
    "is_valid_file",
    "parse_structure",
    "parse_structure_input",
    "read_files_from_csv",
    "read_files_from_folder",
    "read_validated_structure_bytes",
    "resolve_structure_identity",
]
