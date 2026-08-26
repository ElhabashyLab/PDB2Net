"""Parse and render compact annotations embedded in enriched mmCIF files."""

from __future__ import annotations

from collections import Counter
import json
import math
from typing import Any, Iterable, Mapping, Sequence

from .config_loader import config
from .input_contract import InputValidationError
from .reference_data import is_uniprot_accession

ANNOTATION_DATABASES = ("uniprot", "pfam", "cath", "scop2")
ANNOTATION_PIPELINE_VERSION = "pdb2net-sifts-fasta-search-annotations-v3"
NETWORK_ANNOTATION_RULES = {
    "use_embedded_sifts": {"type": "boolean", "default": True},
    "tooltip_fields": {
        "type": "array",
        "default": ["uniprot"],
        "items": {"type": "string", "choices": list(ANNOTATION_DATABASES)},
        "unique_items": True,
    },
    "max_tooltip_segments_per_database": {
        "type": "integer",
        "default": 20,
        "minimum": 1,
        "maximum": 1_000,
    },
}


MAX_EMBEDDED_SEGMENT_ROWS = 200_000
MAX_EMBEDDED_FIELD_LENGTH = 2_000
MAX_ANNOTATION_BYTES_PER_DATABASE = 65_536
MAX_ANNOTATION_BYTES_PER_NODE = 262_144
MAX_INVALID_SEGMENT_WARNINGS = 20


def network_annotation_config(values: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Return a validated annotation/output configuration."""
    raw = values if values is not None else config.get("network_annotations", {})
    if not isinstance(raw, Mapping):
        raise InputValidationError(
            "INVALID_NETWORK_ANNOTATIONS_CONFIG",
            "network_annotations must be a JSON object.",
        )
    use_embedded = raw.get(
        "use_embedded_sifts",
        NETWORK_ANNOTATION_RULES["use_embedded_sifts"]["default"],
    )
    if not isinstance(use_embedded, bool):
        raise InputValidationError(
            "INVALID_NETWORK_ANNOTATIONS_CONFIG",
            "network_annotations.use_embedded_sifts must be true or false.",
        )
    fields = raw.get(
        "tooltip_fields",
        NETWORK_ANNOTATION_RULES["tooltip_fields"]["default"],
    )
    if not isinstance(fields, list) or any(not isinstance(value, str) for value in fields):
        raise InputValidationError(
            "INVALID_NETWORK_ANNOTATIONS_CONFIG",
            "network_annotations.tooltip_fields must be a JSON array of strings.",
        )
    normalized_fields: list[str] = []
    for field in fields:
        normalized = field.strip().lower()
        if normalized not in ANNOTATION_DATABASES:
            raise InputValidationError(
                "INVALID_NETWORK_ANNOTATIONS_CONFIG",
                "Unknown network annotation tooltip field "
                f"{field!r}; allowed values are {', '.join(ANNOTATION_DATABASES)}.",
            )
        if normalized not in normalized_fields:
            normalized_fields.append(normalized)
    maximum_rule = NETWORK_ANNOTATION_RULES["max_tooltip_segments_per_database"]
    maximum = raw.get("max_tooltip_segments_per_database", maximum_rule["default"])
    if (
        isinstance(maximum, bool)
        or not isinstance(maximum, int)
        or maximum < maximum_rule["minimum"]
        or maximum > maximum_rule["maximum"]
    ):
        raise InputValidationError(
            "INVALID_NETWORK_ANNOTATIONS_CONFIG",
            "network_annotations.max_tooltip_segments_per_database must be an integer "
            f"from {maximum_rule['minimum']} to {maximum_rule['maximum']}.",
        )
    return {
        "use_embedded_sifts": use_embedded,
        "tooltip_fields": normalized_fields,
        "max_tooltip_segments_per_database": maximum,
    }


def _clean(value: Any) -> str | None:
    if value is None:
        return None
    cleaned = str(value).strip()
    if not cleaned or cleaned in {".", "?"}:
        return None
    if len(cleaned) > MAX_EMBEDDED_FIELD_LENGTH:
        raise InputValidationError(
            "EMBEDDED_ANNOTATION_FIELD_TOO_LONG",
            "An embedded mmCIF annotation field exceeds the supported length.",
        )
    return cleaned


def _integer(value: Any) -> int | None:
    cleaned = _clean(value)
    if cleaned is None:
        return None
    try:
        return int(cleaned)
    except ValueError:
        return None


def _number(value: Any) -> float | None:
    cleaned = _clean(value)
    if cleaned is None:
        return None
    try:
        parsed = float(cleaned)
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) else None


def _truthy(value: Any) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "y"}


def _category_rows(category: Mapping[str, Sequence[Any]]) -> Iterable[dict[str, Any]]:
    row_count = max((len(values) for values in category.values()), default=0)
    if row_count > MAX_EMBEDDED_SEGMENT_ROWS:
        raise InputValidationError(
            "EMBEDDED_ANNOTATION_LIMIT_EXCEEDED",
            f"An embedded mmCIF annotation category contains {row_count} rows; "
            f"the supported maximum is {MAX_EMBEDDED_SEGMENT_ROWS}.",
        )
    keys = tuple(category)
    for index in range(row_count):
        yield {
            key: category[key][index] if index < len(category[key]) else None
            for key in keys
        }


def _label_to_author_chains(block: Any) -> dict[str, tuple[str, ...]]:
    scheme_mapped: dict[str, set[str]] = {}
    scheme = block.get_mmcif_category("_pdbx_poly_seq_scheme.")
    if scheme:
        for row in _category_rows(scheme):
            label = _clean(row.get("asym_id"))
            author = _clean(row.get("pdb_strand_id"))
            if label and author:
                scheme_mapped.setdefault(label, set()).add(author)
    atom_mapped: dict[str, set[str]] = {}
    atom_site = block.get_mmcif_category("_atom_site.")
    if atom_site:
        labels = atom_site.get("label_asym_id", [])
        authors = atom_site.get("auth_asym_id", [])
        if max(len(labels), len(authors)) > MAX_EMBEDDED_SEGMENT_ROWS * 100:
            raise InputValidationError(
                "STRUCTURE_ATOM_LIMIT_EXCEEDED",
                "The mmCIF atom table exceeds the supported annotation-mapping limit.",
            )
        for label_raw, author_raw in zip(labels, authors):
            label = _clean(label_raw)
            author = _clean(author_raw)
            if label and author:
                atom_mapped.setdefault(label, set()).add(author)
    # The polymer scheme is authoritative when it mentions a label. Atom-site
    # mappings are a fallback only, not a way to vote an ambiguous scheme away.
    mapped = dict(scheme_mapped)
    for label, authors in atom_mapped.items():
        if label not in mapped:
            mapped[label] = authors
    return {key: tuple(sorted(values)) for key, values in mapped.items()}


def _xref_database(value: Any) -> str | None:
    cleaned = (_clean(value) or "").lower()
    if cleaned == "pfam":
        return "pfam"
    if cleaned == "cath":
        return "cath"
    if cleaned.startswith("scop2"):
        return "scop2"
    return None


def _deduplicate_segments(segments: Iterable[Mapping[str, Any]]) -> list[dict[str, Any]]:
    unique: dict[str, dict[str, Any]] = {}
    for segment in segments:
        cleaned = {key: value for key, value in segment.items() if value is not None}
        encoded = json.dumps(cleaned, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
        unique[encoded] = cleaned
    return [unique[key] for key in sorted(unique)]


def extract_embedded_annotations(block: Any) -> dict[str, Any]:
    """Extract compact SIFTS segment categories from one Gemmi CIF block."""
    by_chain: dict[str, dict[str, list[dict[str, Any]]]] = {}
    label_map = _label_to_author_chains(block)
    invalid: Counter[str] = Counter()
    examples: dict[str, str] = {}

    def reject(reason: str, row: Mapping[str, Any]) -> None:
        invalid[reason] += 1
        if reason not in examples:
            label = _clean(row.get("asym_id"))
            examples[reason] = f"label_asym_id={label!r}" if label else "missing label_asym_id"

    def positive_range(start: int | None, end: int | None) -> bool:
        return start is not None and end is not None and start > 0 and end >= start

    uniprot_category = block.get_mmcif_category("_pdbx_sifts_unp_segments.")
    for row in _category_rows(uniprot_category or {}):
        label_chain = _clean(row.get("asym_id"))
        accession = _clean(row.get("unp_acc"))
        required_text = (
            _clean(row.get("entity_id")),
            label_chain,
            accession,
            _clean(row.get("segment_id")),
            _clean(row.get("instance_id")),
        )
        if any(value is None for value in required_text):
            reject("missing_required_field", row)
            continue
        if not is_uniprot_accession(accession):
            reject("invalid_uniprot_accession", row)
            continue
        best_raw = (_clean(row.get("best_mapping")) or "").upper()
        if best_raw not in {"Y", "N"}:
            reject("invalid_best_mapping", row)
            continue
        identity = _number(row.get("identity"))
        if identity is None or not 0.0 <= identity <= 1.0:
            reject("invalid_identity", row)
            continue
        pdb_start = _integer(row.get("seq_id_start"))
        pdb_end = _integer(row.get("seq_id_end"))
        database_start = _integer(row.get("unp_start"))
        database_end = _integer(row.get("unp_end"))
        if not positive_range(pdb_start, pdb_end) or not positive_range(database_start, database_end):
            reject("invalid_range", row)
            continue
        author_chains = label_map.get(str(label_chain), ())
        if len(author_chains) != 1:
            reject("unresolvable_label_asym_id", row)
            continue
        author_chain = author_chains[0]
        segment = {
            "database": "uniprot",
            "accession": accession,
            "chain_id": author_chain,
            "label_asym_id": label_chain,
            "entity_id": required_text[0],
            "segment_id": required_text[3],
            "instance_id": required_text[4],
            "pdb_start": pdb_start,
            "pdb_end": pdb_end,
            "database_start": database_start,
            "database_end": database_end,
            "best_mapping": best_raw == "Y",
            "identity": identity,
        }
        by_chain.setdefault(author_chain, {}).setdefault("uniprot", []).append(segment)

    xref_category = block.get_mmcif_category("_pdbx_sifts_xref_db_segments.")
    for row in _category_rows(xref_category or {}):
        database = _xref_database(row.get("xref_db"))
        label_chain = _clean(row.get("asym_id"))
        accession = _clean(row.get("xref_db_acc"))
        required_text = (
            _clean(row.get("entity_id")),
            label_chain,
            _clean(row.get("xref_db")),
            accession,
            _clean(row.get("segment_id")),
            _clean(row.get("instance_id")),
        )
        if not database or any(value is None for value in required_text):
            reject("missing_or_unknown_xref_field", row)
            continue
        pdb_start = _integer(row.get("seq_id_start"))
        pdb_end = _integer(row.get("seq_id_end"))
        if not positive_range(pdb_start, pdb_end):
            reject("invalid_range", row)
            continue
        author_chains = label_map.get(str(label_chain), ())
        if len(author_chains) != 1:
            reject("unresolvable_label_asym_id", row)
            continue
        author_chain = author_chains[0]
        segment = {
            "database": database,
            "source_database": required_text[2],
            "accession": accession,
            "domain_name": _clean(row.get("domain_name")),
            "chain_id": author_chain,
            "label_asym_id": label_chain,
            "entity_id": required_text[0],
            "segment_id": required_text[4],
            "instance_id": required_text[5],
            "pdb_start": pdb_start,
            "pdb_end": pdb_end,
        }
        by_chain.setdefault(author_chain, {}).setdefault(database, []).append(segment)

    counts: Counter[str] = Counter()
    for databases in by_chain.values():
        for database, segments in list(databases.items()):
            databases[database] = _deduplicate_segments(segments)
            counts[database] += len(databases[database])
    warnings = [
        {
            "code": "INVALID_EMBEDDED_SIFTS_SEGMENT",
            "message": (
                f"Discarded {count} embedded SIFTS segment row(s): {reason} "
                f"({examples.get(reason, 'no example')})."
            ),
            "count": count,
            "reason": reason,
        }
        for reason, count in sorted(invalid.items())[:MAX_INVALID_SEGMENT_WARNINGS]
    ]
    return {
        "is_enriched": bool(sum(counts.values())),
        "by_chain": by_chain,
        "counts": {database: counts.get(database, 0) for database in ANNOTATION_DATABASES},
        "warnings": warnings,
        "invalid_segment_count": sum(invalid.values()),
    }


def apply_embedded_annotations(
    chains: list[dict[str, Any]],
    extracted: Mapping[str, Any],
    *,
    use_embedded_sifts: bool,
) -> None:
    """Attach segments and, when unambiguous, apply embedded UniProt identity."""
    by_chain = extracted.get("by_chain", {})
    if not isinstance(by_chain, Mapping):
        return
    for chain in chains:
        chain_id = str(chain.get("chain_id") or "")
        databases = by_chain.get(chain_id, {})
        if not isinstance(databases, Mapping) or not databases:
            continue
        annotations = {
            database: _deduplicate_segments(segments)
            for database, segments in databases.items()
            if database in ANNOTATION_DATABASES and isinstance(segments, list) and segments
        }
        if not annotations:
            continue
        chain["embedded_annotations"] = annotations
        chain["embedded_annotation_source"] = "wwpdb_nextgen"

        best_accessions = sorted({
            str(segment.get("accession"))
            for segment in annotations.get("uniprot", [])
            if segment.get("best_mapping") and segment.get("accession")
        })
        if not best_accessions:
            if annotations.get("uniprot"):
                chain["embedded_uniprot_status"] = "no_best_mapping"
            continue
        if len(best_accessions) == 1:
            chain["embedded_uniprot_status"] = "unique_best_mapping"
            if use_embedded_sifts:
                accession = best_accessions[0]
                chain["uniprot_id"] = accession
                chain["molecule_type"] = "Protein"
                chain["annotation_source"] = "embedded_sifts"
                chain["matched_database"] = "UniProtKB"
                chain["matched_id"] = accession
                chain["representative_accession"] = accession
                chain["annotation_confidence"] = "high"
        else:
            chain["embedded_uniprot_status"] = "ambiguous_multi_mapping"
            chain["embedded_uniprot_accessions"] = best_accessions
            if use_embedded_sifts:
                chain["uniprot_id"] = None
                chain["molecule_type"] = "Protein"
                chain["annotation_source"] = "embedded_sifts"
                chain["matched_database"] = "UniProtKB"
                chain["matched_id"] = ",".join(best_accessions)
                chain["annotation_confidence"] = "ambiguous"


def _range(start: Any, end: Any) -> str | None:
    if start is None and end is None:
        return None
    if start == end or end is None:
        return str(start)
    if start is None:
        return str(end)
    return f"{start}–{end}"


def _segment_tooltip(database: str, segment: Mapping[str, Any]) -> str:
    labels = {"uniprot": "UniProt mapping", "pfam": "Pfam", "cath": "CATH", "scop2": "SCOP2"}
    accession = str(segment.get("accession") or "Unknown")
    domain = str(segment.get("domain_name") or "").strip()
    if domain and domain != accession:
        accession = f"{accession} ({domain})"
    details: list[str] = []
    pdb_range = _range(segment.get("pdb_start"), segment.get("pdb_end"))
    database_range = _range(segment.get("database_start"), segment.get("database_end"))
    if pdb_range:
        details.append(f"PDB seq. {pdb_range}")
    if database == "uniprot" and database_range:
        details.append(f"UniProt {database_range}")
    identity = segment.get("identity")
    if database == "uniprot" and isinstance(identity, (int, float)):
        details.append(f"identity {identity * 100:.1f}%")
    if database == "uniprot" and segment.get("best_mapping"):
        details.append("best mapping")
    suffix = f"; {'; '.join(details)}" if details else ""
    return f"{labels[database]}: {accession}{suffix}"


def annotation_node_metadata(
    chains: Sequence[Mapping[str, Any]],
    *,
    uniprot_accession: str | None = None,
    annotation_config: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build tooltip lines and deterministic node attributes for selected DBs."""
    selected = network_annotation_config(annotation_config)
    maximum = selected["max_tooltip_segments_per_database"]
    metadata: dict[str, Any] = {"tooltip_lines": [], "annotation_truncated": False}
    node_annotation_bytes = 0
    for database in selected["tooltip_fields"]:
        segments: list[Mapping[str, Any]] = []
        for chain in chains:
            annotations = chain.get("embedded_annotations", {})
            if not isinstance(annotations, Mapping):
                continue
            values = annotations.get(database, [])
            if not isinstance(values, list):
                continue
            for segment in values:
                if not isinstance(segment, Mapping):
                    continue
                if database == "uniprot" and uniprot_accession:
                    if str(segment.get("accession") or "") != uniprot_accession:
                        continue
                segments.append(segment)
        unique = _deduplicate_segments(segments)
        if not unique:
            continue
        tooltip_visible = unique[:maximum]
        metadata["tooltip_lines"].extend(
            _segment_tooltip(database, segment) for segment in tooltip_visible
        )
        tooltip_hidden = len(unique) - len(tooltip_visible)
        if tooltip_hidden:
            metadata["tooltip_lines"].append(
                f"{database.upper()}: +{tooltip_hidden} weitere Segmente"
            )
            metadata["annotation_truncated"] = True

        included: list[dict[str, Any]] = []
        encoded = "[]"
        for segment in unique:
            candidate = included + [dict(segment)]
            candidate_encoded = json.dumps(
                candidate,
                ensure_ascii=False,
                sort_keys=True,
                separators=(",", ":"),
                allow_nan=False,
            )
            candidate_bytes = len(candidate_encoded.encode("utf-8"))
            if candidate_bytes > MAX_ANNOTATION_BYTES_PER_DATABASE:
                break
            if node_annotation_bytes + candidate_bytes > MAX_ANNOTATION_BYTES_PER_NODE:
                break
            included = candidate
            encoded = candidate_encoded
        encoded_bytes = len(encoded.encode("utf-8"))
        raw_truncated = len(included) < len(unique)
        if raw_truncated:
            metadata["annotation_truncated"] = True
        if node_annotation_bytes + encoded_bytes <= MAX_ANNOTATION_BYTES_PER_NODE:
            metadata[f"annotation_{database}"] = encoded
            node_annotation_bytes += encoded_bytes
        metadata[f"annotation_{database}_total"] = len(unique)
        metadata[f"annotation_{database}_included"] = len(included)
        metadata[f"annotation_{database}_truncated"] = raw_truncated
    return metadata


def embedded_annotation_counts(chains: Sequence[Mapping[str, Any]]) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for chain in chains:
        annotations = chain.get("embedded_annotations", {})
        if not isinstance(annotations, Mapping):
            continue
        for database in ANNOTATION_DATABASES:
            values = annotations.get(database, [])
            if isinstance(values, list):
                counts[database] += len(values)
    return {database: counts.get(database, 0) for database in ANNOTATION_DATABASES}


__all__ = [
    "ANNOTATION_DATABASES",
    "annotation_node_metadata",
    "apply_embedded_annotations",
    "embedded_annotation_counts",
    "extract_embedded_annotations",
    "network_annotation_config",
]
