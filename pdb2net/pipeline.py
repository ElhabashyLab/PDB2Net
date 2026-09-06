"""Pipeline orchestration helpers for single-run processing."""

from __future__ import annotations

import gc
import math
import os
import shutil
import sys
import time
from collections import Counter
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterator, List, Mapping, MutableMapping, Optional, Union

from .config_loader import config
from .components import build_identity_edges, find_linked_components, make_component_title
from .cytoscape import create_cytoscape_network, generate_nodes_from_atom_data
from .data_processor import process_structure
from .detailed_results_exporter import (
    DetailedInteractionBudget,
    export_detailed_interactions,
)
from .distances import (
    DISTANCE_THRESHOLD_RULES,
    INTERACTION_FILTER_RULES,
    calculate_distances_with_ckdtree,
    coords_cache,
    tree_cache,
)
from .file_parser import (
    FileSignature,
    get_structure_identity,
    input_file_sha256,
    input_file_signature,
    is_valid_file,
    parse_structure_input,
    read_validated_structure_bytes,
)
from .graph_limits import combined_graph_skip, normalize_combined_graph_limits
from .input_contract import InputValidationError
from .logging_utils import get_logger
from .network_annotations import (
    apply_embedded_annotations,
    embedded_annotation_counts,
    network_annotation_config,
)
from .outputs import (
    RunOutputPaths,
    collect_generated_outputs,
    collect_web_outputs,
    create_run_output_paths,
    reserve_web_output_directory,
    write_failed_run_manifest,
    write_run_manifest,
    write_run_summary,
    write_runtime_analysis,
)
from .protein_network import create_protein_network
from .residue_types import NUCLEIC_ACID_TYPES, count_polymer_lengths
from .structure_identity import StructureIdentity, identity_from_official_id
from .uniprot_matcher import (
    diamond_uniref90_enabled,
    extract_direct_uniprot_accession,
    get_diamond_executable,
    get_diamond_uniref90_db_path,
    parallel_blast_search,
)
from .unknown_molecule_uniprot import process_molecule_info


logger = get_logger(__name__)

NETWORK_OUTPUT_FIELDS = (
    "chain_per_pdb",
    "protein_per_pdb",
    "combined_chain_network",
    "combined_protein_network",
)
STRUCTURE_MODEL_POLICIES = ("first", "all")
ALL_MODEL_FORBIDDEN_NETWORKS = (
    "protein_per_pdb",
    "combined_protein_network",
)
RESOURCE_LIMIT_FIELDS = (
    "max_input_files",
    "max_total_input_bytes",
    "max_single_input_bytes",
    "max_processing_batch_bytes",
    "max_total_input_expanded_bytes",
    "max_single_input_expanded_bytes",
    "max_detailed_interaction_rows",
    "max_detailed_interaction_bytes",
    "min_output_free_bytes",
)
CHAIN_NETWORK_TITLE = "Chain_Interaction_Network"
COMBINED_CHAIN_NETWORK_TITLE = "Combined_Network"


@dataclass
class PipelineTimings:
    """Runtime counters for the major pipeline stages."""

    parsing: float = 0.0
    sifts: float = 0.0
    blast: float = 0.0
    interaction: float = 0.0
    networks: float = 0.0

    def as_dict(self) -> Dict[str, float]:
        """Return the mapping shape expected by the existing runtime log writer."""
        return {
            "parsing": self.parsing,
            "sifts": self.sifts,
            "blast": self.blast,
            "interaction": self.interaction,
            "networks": self.networks,
        }


@dataclass(frozen=True)
class InputInventory:
    """Validated input sizes used for limits, batching, and run metadata."""

    file_sizes: tuple[int, ...]
    expanded_file_sizes: tuple[int, ...] | None = None
    file_signatures: tuple[FileSignature, ...] | None = None
    file_sha256: tuple[str, ...] | None = None

    @property
    def total_bytes(self) -> int:
        return sum(self.file_sizes)

    @property
    def largest_file_bytes(self) -> int:
        return max(self.file_sizes, default=0)

    @property
    def effective_expanded_file_sizes(self) -> tuple[int, ...]:
        return self.expanded_file_sizes or self.file_sizes

    @property
    def total_expanded_bytes(self) -> int:
        return sum(self.effective_expanded_file_sizes)

    @property
    def largest_expanded_file_bytes(self) -> int:
        return max(self.effective_expanded_file_sizes, default=0)


def resolve_workers(value: Optional[Union[int, str]], *, kind: str = "parsing") -> int:
    """Convert a config value to a sensible worker/thread count."""
    cpu = os.cpu_count() or 4

    if isinstance(value, int):
        return max(1, value)

    if kind == "parsing":
        return max(1, cpu - 1)
    if kind == "blast":
        return max(1, cpu - 2)

    return max(1, cpu - 1)


def process_single_file(
    file_path: str,
    pdb_id_override: str | Mapping[str, Any] | None = None,
    expected_signature: FileSignature | None = None,
    expected_sha256: str | None = None,
    maximum_expanded_bytes: int | None = None,
) -> Optional[Dict[str, Any]]:
    """Process a single PDB/mmCIF file end-to-end for the parsing stage."""
    if pdb_id_override:
        if isinstance(pdb_id_override, Mapping):
            identity = StructureIdentity.from_mapping(pdb_id_override)
        else:
            try:
                identity = identity_from_official_id(pdb_id_override)
            except ValueError:
                identity = StructureIdentity("local", str(pdb_id_override))
    else:
        identity = None
    parsed = parse_structure_input(
        file_path,
        expected_signature=expected_signature,
        expected_sha256=expected_sha256,
        maximum_expanded_bytes=maximum_expanded_bytes,
    )
    parsed_identity = StructureIdentity.from_mapping(parsed["structure_identity"])
    if identity is not None and identity.key != parsed_identity.key:
        raise InputValidationError(
            "INPUT_CHANGED_DURING_PROCESSING",
            f"Resolved identity changed while processing {Path(file_path).name}.",
        )
    identity = identity or parsed_identity
    processed = process_structure(
        {
            "file_path": file_path,
            "pdb_id": identity.display_id,
            "structure_identity": identity.as_dict(),
            "structure": parsed["structure"],
        }
    )
    annotation_cfg = network_annotation_config()
    apply_embedded_annotations(
        processed["atom_data"],
        parsed["embedded_annotations"],
        use_embedded_sifts=annotation_cfg["use_embedded_sifts"],
    )
    processed["structure_identity"] = identity.as_dict()
    processed["identity_claims"] = parsed.get("identity_claims", [])
    processed["input_warnings"] = list(parsed.get("identity_warnings", [])) + list(
        parsed.get("embedded_annotations", {}).get("warnings", [])
    )
    processed["input_format"] = parsed["input_format"]
    processed["input_kind"] = parsed["input_kind"]
    processed["embedded_annotation_counts"] = embedded_annotation_counts(processed["atom_data"])
    return processed


def _activate_parsing_worker(config_snapshot: Mapping[str, Any]) -> None:
    """Install the parent run's config once in a spawned parsing worker."""
    from .config_loader import activate_config

    activate_config(dict(config_snapshot))


def discover_input_files(input_path_or_filelist: Union[str, List[str]]) -> List[str]:
    """Resolve a folder or internal file list into valid structure files."""
    if isinstance(input_path_or_filelist, list):
        return [file_path for file_path in input_path_or_filelist if is_valid_file(file_path)]

    input_path = Path(input_path_or_filelist)
    if input_path.is_symlink():
        raise InputValidationError(
            "SYMLINK_INPUT_NOT_ALLOWED",
            f"Input directories must not be symbolic links: {input_path}",
        )
    if not input_path.exists():
        raise InputValidationError(
            "INPUT_PATH_NOT_FOUND",
            f"input_folder_path does not exist: {input_path}",
        )
    if not input_path.is_dir():
        raise InputValidationError(
            "INPUT_PATH_NOT_DIRECTORY",
            f"input_folder_path must be a directory containing PDB/mmCIF files: {input_path}",
        )

    file_paths: list[str] = []
    with os.scandir(str(input_path)) as entries:
        for entry in entries:
            if entry.is_symlink() and is_valid_file(entry.path):
                raise InputValidationError(
                    "SYMLINK_INPUT_NOT_ALLOWED",
                    f"Structure inputs must not be symbolic links: {entry.name}",
                )
            if entry.is_file(follow_symlinks=False) and is_valid_file(entry.path):
                file_paths.append(entry.path)
    if not file_paths:
        raise InputValidationError(
            "NO_VALID_INPUT_FILES",
            f"input_folder_path contains no supported .pdb, .cif, or .mmcif files: {input_path}",
        )

    return sorted(file_paths)


RESOURCE_LIMIT_KEYS = RESOURCE_LIMIT_FIELDS


def _resource_limits() -> Dict[str, Optional[int]]:
    """Return validated optional resource limits; ``None`` means unlimited."""
    raw = config.get("resource_limits", {})
    values = raw if isinstance(raw, Mapping) else {}
    normalized: Dict[str, Optional[int]] = {}
    for key in RESOURCE_LIMIT_KEYS:
        value = values.get(key)
        if value in (None, ""):
            normalized[key] = None
            continue
        if isinstance(value, bool):
            raise InputValidationError(
                "INVALID_RESOURCE_LIMIT",
                f"resource_limits.{key} must be a positive integer or null.",
            )
        try:
            parsed = int(value)
        except (TypeError, ValueError) as exc:
            raise InputValidationError(
                "INVALID_RESOURCE_LIMIT",
                f"resource_limits.{key} must be a positive integer or null.",
            ) from exc
        if parsed <= 0:
            raise InputValidationError(
                "INVALID_RESOURCE_LIMIT",
                f"resource_limits.{key} must be a positive integer or null.",
            )
        normalized[key] = parsed
    return normalized


def _combined_graph_limits() -> Dict[str, Optional[int]]:
    """Validate combined-network caps before expensive pipeline stages."""
    try:
        return normalize_combined_graph_limits(config.get("combined_graph_limits"))
    except ValueError as exc:
        raise InputValidationError("INVALID_COMBINED_GRAPH_LIMIT", str(exc)) from exc


def _validate_analysis_config() -> None:
    """Validate scientific settings before references, parsing, or Cytoscape."""
    networks = config.get("networks")
    if not isinstance(networks, Mapping):
        raise InputValidationError("INVALID_NETWORK_CONFIG", "networks must be a JSON object.")
    for field in NETWORK_OUTPUT_FIELDS:
        if not isinstance(networks.get(field), bool):
            raise InputValidationError(
                "INVALID_NETWORK_CONFIG", f"networks.{field} must be true or false."
            )

    if config.get("layout_mode", "python_fast") != "python_fast":
        raise InputValidationError(
            "INVALID_LAYOUT_MODE",
            "layout_mode must be 'python_fast'.",
        )

    policy = config.get("structure_model_policy", "first")
    if not isinstance(policy, str) or policy.strip().lower() not in STRUCTURE_MODEL_POLICIES:
        raise InputValidationError(
            "INVALID_STRUCTURE_MODEL_POLICY",
            "structure_model_policy must be 'first' or 'all'.",
        )
    policy = policy.strip().lower()
    if policy == "all" and any(networks.get(field) for field in ALL_MODEL_FORBIDDEN_NETWORKS):
        raise InputValidationError(
            "PROTEIN_NETWORKS_UNSUPPORTED_FOR_ALL_MODELS",
            "Protein networks cannot be generated when structure_model_policy='all'; "
            "select model-separated chain networks or use the first model.",
        )

    thresholds = config.get("distance_thresholds")
    if not isinstance(thresholds, Mapping):
        raise InputValidationError(
            "INVALID_DISTANCE_THRESHOLDS", "distance_thresholds must be a JSON object."
        )
    for field, rule in DISTANCE_THRESHOLD_RULES.items():
        minimum = float(rule["minimum"])
        maximum = float(rule["maximum"])
        value = thresholds.get(field)
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise InputValidationError(
                "INVALID_DISTANCE_THRESHOLDS",
                f"distance_thresholds.{field} must be a number from {minimum:g} to {maximum:g}.",
            )
        number = float(value)
        if not math.isfinite(number) or number < minimum or number > maximum:
            raise InputValidationError(
                "INVALID_DISTANCE_THRESHOLDS",
                f"distance_thresholds.{field} must be a number from {minimum:g} to {maximum:g}.",
            )

    filters = config.get("interaction_filters")
    if not isinstance(filters, Mapping):
        raise InputValidationError(
            "INVALID_INTERACTION_FILTERS", "interaction_filters must be a JSON object."
        )
    for field, rule in INTERACTION_FILTER_RULES.items():
        value = filters.get(field)
        minimum = int(rule["minimum"])
        maximum = int(rule["maximum"])
        if (
            isinstance(value, bool)
            or not isinstance(value, int)
            or not minimum <= value <= maximum
        ):
            raise InputValidationError(
                "INVALID_INTERACTION_FILTERS",
                f"interaction_filters.{field} must be an integer from {minimum} to {maximum}.",
            )
    if not isinstance(config.get("export_detailed_interactions", False), bool):
        raise InputValidationError(
            "INVALID_EXPORT_CONFIG", "export_detailed_interactions must be true or false."
        )


def _preflight_structure_identities(
    file_paths: List[str], inventory: InputInventory
) -> dict[str, StructureIdentity]:
    """Resolve unique run identities before expensive parsing and annotation."""
    signatures = inventory.file_signatures or tuple(input_file_signature(path) for path in file_paths)
    digests = inventory.file_sha256 or tuple(
        input_file_sha256(path, expected_signature=signature)
        for path, signature in zip(file_paths, signatures)
    )
    maximum_expanded_bytes = _resource_limits()["max_single_input_expanded_bytes"]
    resolved: dict[str, StructureIdentity] = {}
    seen: dict[str, str] = {}
    for file_path, expected_signature, expected_sha256 in zip(file_paths, signatures, digests):
        if input_file_signature(file_path) != expected_signature:
            raise InputValidationError(
                "INPUT_CHANGED_DURING_PROCESSING",
                f"Input file changed before identity resolution: {Path(file_path).name}",
            )
        identity = get_structure_identity(
            file_path,
            expected_signature=expected_signature,
            expected_sha256=expected_sha256,
            maximum_expanded_bytes=maximum_expanded_bytes,
        )
        if input_file_signature(file_path) != expected_signature:
            raise InputValidationError(
                "INPUT_CHANGED_DURING_PROCESSING",
                f"Input file changed during identity resolution: {Path(file_path).name}",
            )
        previous = seen.get(identity.key)
        if previous is not None:
            raise InputValidationError(
                "DUPLICATE_STRUCTURE_IDENTITY",
                f"Inputs {Path(previous).name} and {Path(file_path).name} resolve to the same "
                f"structure identity {identity.key!r}.",
            )
        seen[identity.key] = file_path
        resolved[file_path] = identity
    return resolved


def _validate_input_file_count(
    file_paths: List[str], limits: Mapping[str, Optional[int]]
) -> None:
    maximum = limits["max_input_files"]
    if maximum is not None and len(file_paths) > maximum:
        raise InputValidationError(
            "INPUT_FILE_COUNT_LIMIT_EXCEEDED",
            f"Input contains {len(file_paths)} files; configured maximum is {maximum}.",
        )


def _validate_input_total_bytes(
    total_bytes: int, limits: Mapping[str, Optional[int]]
) -> None:
    maximum = limits["max_total_input_bytes"]
    if maximum is not None and total_bytes > maximum:
        raise InputValidationError(
            "INPUT_TOTAL_BYTES_LIMIT_EXCEEDED",
            f"Input totals {total_bytes} bytes; configured maximum is {maximum}.",
        )


def _validate_input_total_expanded_bytes(
    total_expanded_bytes: int, limits: Mapping[str, Optional[int]]
) -> None:
    maximum_expanded = limits["max_total_input_expanded_bytes"]
    if maximum_expanded is not None and total_expanded_bytes > maximum_expanded:
        raise InputValidationError(
            "INPUT_TOTAL_EXPANDED_BYTES_LIMIT_EXCEEDED",
            f"Expanded input totals {total_expanded_bytes} bytes; "
            f"configured maximum is {maximum_expanded}.",
        )


def inspect_input_files(
    file_paths: List[str], *, enforce_aggregate_limits: bool = True
) -> InputInventory:
    """Inspect inputs, enforcing single-file and optionally aggregate limits."""
    limits = _resource_limits()
    if enforce_aggregate_limits:
        _validate_input_file_count(file_paths, limits)

    file_sizes: List[int] = []
    expanded_file_sizes: List[int] = []
    file_signatures: List[FileSignature] = []
    file_sha256: List[str] = []
    for file_path in file_paths:
        path = Path(file_path)
        try:
            signature = input_file_signature(path)
        except InputValidationError as exc:
            if exc.code == "INPUT_FILE_STAT_FAILED" and not path.exists():
                raise InputValidationError("INPUT_FILE_NOT_FOUND", f"Input file does not exist: {path}") from exc
            raise
        size = signature[2]

        max_single = limits["max_single_input_bytes"]
        if max_single is not None and size > max_single:
            raise InputValidationError(
                "INPUT_FILE_BYTES_LIMIT_EXCEEDED",
                f"Input file {path.name} is {size} bytes; configured maximum is {max_single}.",
            )
        file_sizes.append(size)

        _payload, expanded_size, raw_sha256 = read_validated_structure_bytes(
            str(path),
            expected_signature=signature,
            maximum_expanded_bytes=limits["max_single_input_expanded_bytes"],
            collect=False,
            include_sha256=True,
        )
        maximum_expanded = limits["max_single_input_expanded_bytes"]
        if maximum_expanded is not None and expanded_size > maximum_expanded:
            raise InputValidationError(
                "INPUT_FILE_EXPANDED_BYTES_LIMIT_EXCEEDED",
                f"Expanded input file {path.name} is {expanded_size} bytes; "
                f"configured maximum is {maximum_expanded}.",
            )
        max_batch = limits["max_processing_batch_bytes"]
        if max_batch is not None and expanded_size > max_batch:
            raise InputValidationError(
                "INPUT_PROCESSING_BATCH_LIMIT_EXCEEDED",
                f"Expanded input file {path.name} is {expanded_size} bytes and cannot fit within the "
                f"configured processing batch limit of {max_batch} bytes.",
            )
        expanded_file_sizes.append(expanded_size)
        file_sha256.append(raw_sha256)
        if input_file_signature(path) != signature:
            raise InputValidationError(
                "INPUT_CHANGED_DURING_PROCESSING",
                f"Input file changed while it was being inspected: {path.name}",
            )
        file_signatures.append(signature)

    inventory = InputInventory(
        tuple(file_sizes),
        tuple(expanded_file_sizes),
        tuple(file_signatures),
        tuple(file_sha256),
    )
    if enforce_aggregate_limits:
        _validate_input_total_bytes(inventory.total_bytes, limits)
        _validate_input_total_expanded_bytes(inventory.total_expanded_bytes, limits)
    return inventory


def create_processing_batches(
    file_paths: List[str],
    inventory: InputInventory,
) -> Iterator[List[str]]:
    """Yield stable batches bounded by expanded input size."""
    max_batch = _resource_limits()["max_processing_batch_bytes"]
    if max_batch is None:
        if file_paths:
            yield list(file_paths)
        return

    current: List[str] = []
    current_bytes = 0
    for file_path, size in zip(file_paths, inventory.effective_expanded_file_sizes):
        if current and current_bytes + size > max_batch:
            yield current
            current = []
            current_bytes = 0
        current.append(file_path)
        current_bytes += size
    if current:
        yield current


def _validate_output_root() -> None:
    """Ensure the configured output folder can hold timestamped run outputs."""
    configured_output = str(config.get("output_path") or "").strip()
    if not configured_output:
        raise InputValidationError("OUTPUT_PATH_MISSING", "output_path is not configured.")
    output_root = Path(configured_output)
    try:
        output_root.mkdir(parents=True, exist_ok=True)
    except Exception as exc:
        raise InputValidationError("OUTPUT_PATH_NOT_WRITABLE", f"output_path cannot be created: {output_root}") from exc
    if not output_root.is_dir():
        raise InputValidationError("OUTPUT_PATH_NOT_DIRECTORY", f"output_path is not a directory: {output_root}")


def _validate_required_reference_files() -> None:
    """Fail early with actionable messages for reference files needed by normal runs."""
    required = {
        "pdb_fasta_path": "pdb_seqres.txt",
        "sifts_tsv_path": "pdb_chain_uniprot.tsv",
        "uniprot_fasta_path": "uniprot_sprot.fasta",
    }
    missing = []
    for config_key, label in required.items():
        path = Path(str(config.get(config_key) or ""))
        if not str(path) or not path.is_file():
            missing.append(f"{label} ({config_key}): {path}")
    if missing:
        raise InputValidationError(
            "REFERENCE_FILE_MISSING",
            "Required reference file(s) are missing. Configure them with config.local.json, "
            "PDB2NET_CONFIG_FILE, or CLI flags. Missing: " + "; ".join(missing),
        )


def _parse_input_files(
    file_paths: List[str],
    *,
    expected_identities: Mapping[str, StructureIdentity] | None = None,
    expected_signatures: Mapping[str, FileSignature] | None = None,
    expected_sha256: Mapping[str, str] | None = None,
    maximum_expanded_bytes: int | None = None,
    errors_by_path: MutableMapping[str, Exception] | None = None,
) -> List[Dict[str, Any]]:
    """Parse files with bounded scheduling and optional per-file error collection."""
    configured_workers = resolve_workers(config.get("workers", {}).get("parsing"), kind="parsing")
    parsing_workers = min(configured_workers, max(1, len(file_paths)))
    logger.info("[Workers] Parsing processes: %s", parsing_workers)
    indexed_results: Dict[int, Optional[Dict[str, Any]]] = {}
    with ProcessPoolExecutor(
        max_workers=parsing_workers,
        initializer=_activate_parsing_worker,
        initargs=(config,),
    ) as executor:
        pending: Dict[Any, int] = {}
        next_index = 0

        while next_index < len(file_paths) or pending:
            while next_index < len(file_paths) and len(pending) < parsing_workers:
                file_path = file_paths[next_index]
                expected_identity = (
                    expected_identities.get(file_path) if expected_identities is not None else None
                )
                expected_signature = (
                    expected_signatures.get(file_path) if expected_signatures is not None else None
                )
                expected_digest = (
                    expected_sha256.get(file_path) if expected_sha256 is not None else None
                )
                future = executor.submit(
                    process_single_file,
                    file_path,
                    expected_identity.as_dict() if expected_identity is not None else None,
                    expected_signature,
                    expected_digest,
                    maximum_expanded_bytes,
                )
                pending[future] = next_index
                next_index += 1

            completed, _ = wait(pending, return_when=FIRST_COMPLETED)
            for future in completed:
                index = pending.pop(future)
                try:
                    indexed_results[index] = future.result()
                except Exception as exc:
                    if errors_by_path is None:
                        raise
                    errors_by_path[file_paths[index]] = exc
                    indexed_results[index] = None

    return [
        result
        for index in range(len(file_paths))
        if (result := indexed_results.get(index)) is not None
    ]


def _register_unique_chain_ids(
    structures: List[Dict[str, Any]],
    seen: dict[str, str],
) -> None:
    """Reject public chain-node IDs that would alias across a run."""
    for structure in structures:
        file_label = Path(str(structure.get("file_path") or "input")).name
        for chain in structure.get("atom_data", []):
            unique_id = str(chain.get("unique_chain_id") or "")
            if not unique_id:
                raise InputValidationError(
                    "INVALID_CHAIN_NODE_ID",
                    f"Parsed chain in {file_label} has no public node ID.",
                )
            previous = seen.get(unique_id)
            if previous is not None:
                raise InputValidationError(
                    "DUPLICATE_CHAIN_NODE_ID",
                    f"Inputs {previous} and {file_label} produce the same chain node ID {unique_id!r}.",
                )
            seen[unique_id] = file_label


def _validate_web_reference_manifest(web_output_dir: str | None) -> None:
    if web_output_dir and not str(config.get("reference_manifest_id") or "").strip():
        raise InputValidationError(
            "REFERENCE_MANIFEST_ID_REQUIRED",
            "Web-output runs require a non-empty reference_manifest_id before processing starts.",
        )


def _run_blast_annotation(combined_data: List[Dict[str, Any]]) -> None:
    """Run BLAST fallback for unresolved chain annotations."""
    _validate_blast_ready_if_needed(combined_data)
    blast_workers = resolve_workers(config.get("workers", {}).get("blast_threads"), kind="blast")
    logger.info("[Workers] BLAST threads: %s", blast_workers)
    parallel_blast_search(combined_data, max_workers=blast_workers)


def _blast_fallback_needed(combined_data: List[Dict[str, Any]]) -> bool:
    """Return whether at least one chain may require BLAST fallback annotation."""
    for structure in combined_data:
        file_path = str(structure.get("file_path") or "")
        pdb_id = str(structure.get("pdb_id") or "")
        atom_data = list(structure.get("atom_data", []))
        has_direct_single_chain_uniprot = (
            len(atom_data) == 1
            and bool(extract_direct_uniprot_accession(file_path) or extract_direct_uniprot_accession(pdb_id))
        )
        for chain in structure.get("atom_data", []):
            if (
                chain.get("embedded_uniprot_status") == "ambiguous_multi_mapping"
                and chain.get("annotation_source") == "embedded_sifts"
            ):
                continue
            if chain.get("external_sifts_status") == "ambiguous_external_mapping":
                continue
            if chain.get("uniprot_id") not in (None, "Unknown"):
                continue
            if has_direct_single_chain_uniprot:
                continue
            molecule_type = (chain.get("molecule_type") or "").strip()
            if molecule_type not in NUCLEIC_ACID_TYPES:
                return True
    return False


def _validate_blast_ready_if_needed(combined_data: List[Dict[str, Any]]) -> None:
    """Fail early with clear messages when BLAST fallback is required but unavailable."""
    if not _blast_fallback_needed(combined_data):
        return

    blastp = str(config.get("blastp_executable") or "blastp")
    if os.path.sep in blastp or (os.path.altsep and os.path.altsep in blastp):
        blastp_ok = os.path.exists(blastp)
    else:
        blastp_ok = shutil.which(blastp) is not None
    if not blastp_ok:
        raise InputValidationError("BLASTP_NOT_FOUND", f"blastp executable is not available: {blastp}")

    blast_db_path = str(config.get("blast_db_path") or "")
    db_prefix = os.path.join(blast_db_path, "uniprot_db")
    missing = [path for path in (db_prefix + ".pin", db_prefix + ".phr", db_prefix + ".psq") if not os.path.exists(path)]
    if missing:
        raise InputValidationError(
            "BLAST_DATABASE_MISSING",
            "BLAST fallback is required but the UniProt BLAST database is incomplete. "
            f"Missing: {', '.join(missing)}",
        )

    if diamond_uniref90_enabled():
        diamond = get_diamond_executable()
        if os.path.sep in diamond or (os.path.altsep and os.path.altsep in diamond):
            diamond_ok = os.path.exists(diamond)
        else:
            diamond_ok = shutil.which(diamond) is not None
        if not diamond_ok:
            raise InputValidationError("DIAMOND_NOT_FOUND", f"DIAMOND executable is not available: {diamond}")

        diamond_db_path = get_diamond_uniref90_db_path()
        diamond_db_ok = bool(diamond_db_path) and (
            os.path.exists(diamond_db_path) or os.path.exists(diamond_db_path + ".dmnd")
        )
        if not diamond_db_ok:
            raise InputValidationError(
                "DIAMOND_DATABASE_MISSING",
                "DIAMOND/UniRef90 fallback is enabled but the DIAMOND database is missing. "
                f"Configure diamond.uniref90_db_path or PDB2NET_DIAMOND_UNIREF90_DB. Path: {diamond_db_path}",
            )


def _compact_structure_summaries(structures: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Drop atom coordinates after distance/detail work while preserving network metadata."""
    compacted: List[Dict[str, Any]] = []
    for structure in structures:
        structure_summary = {key: value for key, value in structure.items() if key != "atom_data"}
        chain_summaries: List[Dict[str, Any]] = []
        for chain in structure.get("atom_data", []):
            chain_summary = {key: value for key, value in chain.items() if key != "residues"}
            aa_len, nt_len = count_polymer_lengths(
                residue.get("residue_name") for residue in chain.get("residues", [])
            )
            chain_summary["aa_len"] = aa_len
            chain_summary["nt_len"] = nt_len
            chain_summaries.append(chain_summary)
        structure_summary["atom_data"] = chain_summaries
        compacted.append(structure_summary)
    return compacted


def _annotation_summary(combined_data: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Summarize final chain annotations without changing annotation decisions."""
    by_source: Counter[str] = Counter()
    by_database: Counter[str] = Counter()
    by_confidence: Counter[str] = Counter()
    chains_total = 0
    chains_annotated = 0
    embedded_counts: Counter[str] = Counter()
    enriched_structures = 0
    enriched_chains = 0

    for structure in combined_data:
        chains = list(structure.get("atom_data", []))
        if structure.get("input_kind") == "nextgen_enriched_mmcif":
            enriched_structures += 1
        direct_accession = None
        if len(chains) == 1:
            direct_accession = (
                extract_direct_uniprot_accession(str(structure.get("file_path") or ""))
                or extract_direct_uniprot_accession(str(structure.get("pdb_id") or ""))
            )

        for chain in chains:
            chains_total += 1
            chain_embedded = embedded_annotation_counts([chain])
            if sum(chain_embedded.values()):
                enriched_chains += 1
                embedded_counts.update(chain_embedded)
            source = str(chain.get("annotation_source") or "").strip()
            database = str(chain.get("matched_database") or "").strip()
            confidence = str(chain.get("annotation_confidence") or "").strip()
            uniprot_id = str(chain.get("uniprot_id") or "").strip()
            matched_id = str(chain.get("matched_id") or "").strip()

            if not source and uniprot_id:
                source = "direct_uniprot" if direct_accession == uniprot_id else "sifts"
            if source in {"direct_uniprot", "sifts"} and not database:
                database = "UniProtKB"

            if source or matched_id or uniprot_id:
                chains_annotated += 1
                by_source[source or "unspecified"] += 1
                by_database[database or "unspecified"] += 1
                by_confidence[confidence or "not_reported"] += 1

    annotation_cfg = network_annotation_config()
    return {
        "chains_total": chains_total,
        "chains_annotated": chains_annotated,
        "chains_unannotated": chains_total - chains_annotated,
        "by_source": dict(sorted(by_source.items())),
        "by_database": dict(sorted(by_database.items())),
        "by_confidence": dict(sorted(by_confidence.items())),
        "embedded_sifts": {
            "structures": enriched_structures,
            "chains": enriched_chains,
            "segments_by_database": {
                database: embedded_counts.get(database, 0)
                for database in ("uniprot", "pfam", "cath", "scop2")
            },
            "used_for_identity": annotation_cfg["use_embedded_sifts"],
        },
        "tooltip_fields": list(annotation_cfg["tooltip_fields"]),
    }


def _structure_input_summary(combined_data: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Return stable per-structure identity and input metadata for contract 2.0."""
    summaries: List[Dict[str, Any]] = []
    for structure in combined_data:
        identity = structure.get("structure_identity", {})
        summaries.append(
            {
                "file": str(structure.get("file_path") or ""),
                "identity": dict(identity) if isinstance(identity, Mapping) else {},
                "format": str(structure.get("input_format") or "unknown"),
                "kind": str(structure.get("input_kind") or "unknown"),
                "embedded_annotation_counts": dict(
                    structure.get("embedded_annotation_counts")
                    if isinstance(structure.get("embedded_annotation_counts"), Mapping)
                    else embedded_annotation_counts(structure.get("atom_data", []))
                ),
            }
        )
    return summaries


def _identity_summary(combined_data: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    identities: List[Dict[str, Any]] = []
    for structure in combined_data:
        identity = structure.get("structure_identity", {})
        if not isinstance(identity, Mapping):
            continue
        key = str(identity.get("key") or identity.get("canonical_id") or "")
        if key:
            identities.append(dict(identity))
    return identities


def _embedded_annotation_warnings(combined_data: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    cfg = network_annotation_config()
    maximum = cfg["max_tooltip_segments_per_database"]
    selected = set(cfg["tooltip_fields"])
    warnings: List[Dict[str, Any]] = []
    for structure in combined_data:
        for warning in structure.get("input_warnings", []):
            if isinstance(warning, Mapping):
                warnings.append(dict(warning))
        for chain in structure.get("atom_data", []):
            for warning in chain.get("annotation_warnings", []):
                if isinstance(warning, Mapping):
                    warnings.append(dict(warning))
            if chain.get("embedded_uniprot_status") == "ambiguous_multi_mapping" and cfg["use_embedded_sifts"]:
                warnings.append(
                    {
                        "code": "AMBIGUOUS_EMBEDDED_UNIPROT_MAPPING",
                        "message": (
                            f"{chain.get('unique_chain_id')} has multiple best embedded UniProt mappings; "
                            "no protein identity node was assigned."
                        ),
                    }
                )
            annotations = chain.get("embedded_annotations", {})
            if not isinstance(annotations, Mapping):
                continue
            for database in sorted(selected):
                segments = annotations.get(database, [])
                if isinstance(segments, list) and len(segments) > maximum:
                    warnings.append(
                        {
                            "code": "ANNOTATION_TOOLTIP_TRUNCATED",
                            "message": (
                                f"{chain.get('unique_chain_id')} has {len(segments)} {database} segments; "
                                f"the tooltip displays the first {maximum}."
                            ),
                        }
                    )
    return warnings


def _reference_file_metadata(path_value: object) -> Dict[str, Any]:
    raw_path = str(path_value or "")
    path = Path(raw_path) if raw_path else None
    try:
        size_bytes = path.stat().st_size if path and path.is_file() else None
    except OSError:
        size_bytes = None
    return {
        "path": str(path) if path else None,
        "size_bytes": size_bytes,
    }


def _reference_summary() -> Dict[str, Any]:
    """Record reference identity and lightweight file metadata without hashing large data."""
    blast_db_path = str(config.get("blast_db_path") or "")
    blast_prefix = Path(blast_db_path) / "uniprot_db" if blast_db_path else None
    diamond_path_value = get_diamond_uniref90_db_path() if diamond_uniref90_enabled() else ""
    diamond_path = Path(diamond_path_value)
    if diamond_path_value and not diamond_path.exists() and Path(diamond_path_value + ".dmnd").exists():
        diamond_path = Path(diamond_path_value + ".dmnd")
    try:
        diamond_size = diamond_path.stat().st_size if diamond_path_value and diamond_path.is_file() else None
    except OSError:
        diamond_size = None

    return {
        "manifest_id": str(config.get("reference_manifest_id") or "") or None,
        "pdb_fasta": _reference_file_metadata(config.get("pdb_fasta_path")),
        "sifts_tsv": _reference_file_metadata(config.get("sifts_tsv_path")),
        "uniprot_fasta": _reference_file_metadata(config.get("uniprot_fasta_path")),
        "blast_database": {
            "prefix": str(blast_prefix) if blast_prefix else None,
            "files_present": sum(
                Path(str(blast_prefix) + suffix).is_file()
                for suffix in (".pin", ".phr", ".psq")
            ) if blast_prefix else 0,
        },
        "diamond_database": {
            "enabled": diamond_uniref90_enabled(),
            "path": str(diamond_path) if diamond_path_value else None,
            "size_bytes": diamond_size,
        },
    }


def _peak_rss_bytes(who: int) -> Optional[int]:
    try:
        import resource

        value = int(resource.getrusage(who).ru_maxrss)
    except (ImportError, OSError, ValueError):
        return None
    return value if sys.platform == "darwin" else value * 1024


def _resource_summary(
    inventory: InputInventory,
    *,
    processing_batches: int,
    parsing_workers_used: int,
    detailed_interactions: DetailedInteractionBudget | None = None,
) -> Dict[str, Any]:
    limits = _resource_limits()
    try:
        import resource

        main_rss = _peak_rss_bytes(resource.RUSAGE_SELF)
        child_rss = _peak_rss_bytes(resource.RUSAGE_CHILDREN)
    except ImportError:
        main_rss = None
        child_rss = None
    detailed_summary = (
        detailed_interactions.as_dict()
        if detailed_interactions is not None
        else {
            "enabled": False,
            "rows": 0,
            "bytes": 0,
            "max_rows": limits["max_detailed_interaction_rows"],
            "max_bytes": limits["max_detailed_interaction_bytes"],
            "min_free_bytes": limits["min_output_free_bytes"],
        }
    )
    return {
        "input": {
            "files": len(inventory.file_sizes),
            "total_bytes": inventory.total_bytes,
            "largest_file_bytes": inventory.largest_file_bytes,
            "total_expanded_bytes": inventory.total_expanded_bytes,
            "largest_expanded_file_bytes": inventory.largest_expanded_file_bytes,
        },
        "processing": {
            "batches": processing_batches,
            "parsing_workers": parsing_workers_used,
            "max_batch_bytes": limits["max_processing_batch_bytes"],
            "batch_size_basis": "expanded_bytes",
        },
        "peak_rss_bytes": {
            "main_process": main_rss,
            "child_processes": child_rss,
        },
        "detailed_interactions": detailed_summary,
    }


def _export_detailed_interaction_tables(
    combined_data: List[Dict[str, Any]],
    results: List[Dict[str, Any]],
    distances_dir: str,
    budget: DetailedInteractionBudget,
) -> None:
    """Write optional atom-level interaction CSVs."""
    for structure_data in combined_data:
        structure_key = str(structure_data.get("structure_identity", {}).get("key") or "")
        chain_ids = {str(chain.get("unique_chain_id")) for chain in structure_data.get("atom_data", [])}
        pdb_interactions = [
            result
            for result in results
            if (
                str(result.get("structure_key") or "") == structure_key
                or (
                    not result.get("structure_key")
                    and result.get("chain_a") in chain_ids
                    and result.get("chain_b") in chain_ids
                )
            )
        ]
        export_detailed_interactions(
            structure_data,
            pdb_interactions,
            distances_dir,
            budget=budget,
        )


def _export_chain_networks(
    combined_data: List[Dict[str, Any]],
    results: List[Dict[str, Any]],
    chain_dir: str,
) -> None:
    """Export chain interaction networks per PDB, including nodes-only cases."""
    include_models = str(config.get("structure_model_policy", "first")).strip().lower() == "all"
    for structure in combined_data:
        pdb_id = structure["pdb_id"]
        structure_key = str(structure.get("structure_identity", {}).get("key") or "")
        _annotate_chain_parent_metadata(structure)
        model_indices = sorted(
            {int(chain.get("model_index", 1)) for chain in structure.get("atom_data", [])}
        )
        for model_index in model_indices:
            model_chains = [
                chain
                for chain in structure.get("atom_data", [])
                if int(chain.get("model_index", 1)) == model_index
            ]
            chain_ids = {str(chain.get("unique_chain_id")) for chain in model_chains}
            model_results = [
                entry
                for entry in results
                if int(entry.get("model_index", 1)) == model_index
                and (
                    str(entry.get("structure_key") or "") == structure_key
                    or (
                        not entry.get("structure_key")
                        and entry.get("chain_a") in chain_ids
                        and entry.get("chain_b") in chain_ids
                    )
                )
            ]
            nodes_data = generate_nodes_from_atom_data(model_chains, pdb_id)
            network_title = f"{CHAIN_NETWORK_TITLE}_{pdb_id}"
            if include_models:
                network_title += f"_model{model_index}"
            create_cytoscape_network(
                model_results,
                network_title,
                chain_dir,
                nodes_data=nodes_data,
            )


def _export_network_outputs(
    combined_data: List[Dict[str, Any]],
    results: List[Dict[str, Any]],
    network_config: Dict[str, Any],
    output_paths: RunOutputPaths,
    timings: PipelineTimings,
) -> List[Dict[str, Any]]:
    """Run the configured network export stages and accumulate timings."""
    graph_limits = _combined_graph_limits()
    skipped_outputs: List[Dict[str, Any]] = []

    if network_config["chain_per_pdb"]:
        start_time = time.time()
        _export_chain_networks(combined_data, results, output_paths.chain_dir)
        timings.networks += time.time() - start_time

    if network_config["combined_chain_network"]:
        start_time = time.time()
        skipped_outputs.extend(
            _create_linked_identity_network(
                results,
                combined_data,
                output_paths.combined_dir,
                combined_graph_limits=graph_limits,
            )
        )
        timings.networks += time.time() - start_time

    if network_config["protein_per_pdb"] or network_config["combined_protein_network"]:
        start_time = time.time()
        skipped_outputs.extend(
            create_protein_network(
                results,
                combined_data,
                output_paths.run_output_path,
                network_config,
                per_pdb_output_path=output_paths.protein_dir,
                combined_output_path=output_paths.combined_dir,
                combined_graph_limits=graph_limits,
            )
        )
        timings.networks += time.time() - start_time

    return skipped_outputs


def _config_snapshot(network_config: Dict[str, Any]) -> Dict[str, Any]:
    """Return the small config subset useful for manifests and webserver jobs."""
    return {
        "networks": network_config,
        "network_annotations": network_annotation_config(),
        "distance_thresholds": config.get("distance_thresholds", {}),
        "interaction_filters": config.get("interaction_filters", {}),
        "structure_model_policy": config.get("structure_model_policy", "first"),
        "workers": config.get("workers", {}),
        "resource_limits": config.get("resource_limits", {}),
        "combined_graph_limits": config.get("combined_graph_limits", {}),
        "layout_mode": config.get("layout_mode", "python_fast"),
        "open_in_cytoscape": config.get("open_in_cytoscape"),
        "export_detailed_interactions": config.get("export_detailed_interactions", False),
        "blast_cache_path": config.get("blast_cache_path", ""),
        "diamond": config.get("diamond", {}),
    }


def _get_border_color_for_uniprot(uniprot_id: Optional[str]) -> str:
    """Generate a consistent hex color from a UniProt ID string (or return gray)."""
    if not uniprot_id:
        return "#555555"
    import hashlib

    hash_val = int(hashlib.md5(uniprot_id.encode("utf-8")).hexdigest(), 16)
    return f"#{hash_val & 0xFFFFFF:06x}"


def _annotate_chain_parent_metadata(structure: Dict[str, Any]) -> None:
    """Attach source metadata used by chain-level Cytoscape node tables."""
    pdb_id = structure["pdb_id"]
    file_path = structure.get("file_path", "")
    file_label = os.path.basename(file_path) if file_path else pdb_id
    structure_key = str(structure.get("structure_identity", {}).get("key") or "")

    for chain in structure["atom_data"]:
        chain["_parent_pdb_id"] = pdb_id
        chain["_parent_structure_key"] = structure_key
        chain["_parent_file_path"] = file_path
        chain["_parent_file_label"] = file_label


def _create_linked_identity_network(
    results: List[Dict[str, Any]],
    combined_data: List[Dict[str, Any]],
    run_output_path: str,
    *,
    combined_graph_limits: Mapping[str, int | None] | None = None,
) -> List[Dict[str, Any]]:
    """Build separate linked-identity combined networks."""
    graph_limits = normalize_combined_graph_limits(combined_graph_limits)
    skipped_outputs: List[Dict[str, Any]] = []
    logger.info("Building Linked Identity Network (bridging different PDBs)")

    all_chains_flat = []
    uniprot_groups: Dict[str, List[Dict[str, Any]]] = {}

    for structure in combined_data:
        _annotate_chain_parent_metadata(structure)

        for chain in structure["atom_data"]:
            up_id = chain.get("uniprot_id")
            chain["uniprot_border_color"] = _get_border_color_for_uniprot(up_id)
            all_chains_flat.append(chain)

            if up_id:
                uniprot_groups.setdefault(up_id, []).append(chain)

    chain_lookup = {chain["unique_chain_id"]: chain for chain in all_chains_flat}
    chain_to_pdb = {
        chain["unique_chain_id"]: chain.get("_parent_pdb_id", "")
        for chain in all_chains_flat
    }
    uniprot_to_chain_ids = {
        up_id: [chain["unique_chain_id"] for chain in chains]
        for up_id, chains in uniprot_groups.items()
    }

    identity_edges = build_identity_edges(uniprot_to_chain_ids, chain_to_pdb)
    logger.info("Added %s identity bridges between files", len(identity_edges))

    if not identity_edges:
        logger.info("No interlinked components found. Skipping Combined_Network export.")
        return skipped_outputs

    component_node_sets = find_linked_components(results, identity_edges, valid_nodes=chain_lookup.keys())

    logger.info("Found %s interlinked component(s)", len(component_node_sets))

    for index, component_nodes in enumerate(component_node_sets, start=1):
        filtered_chains = [
            chain_lookup[node_id]
            for node_id in sorted(component_nodes)
            if node_id in chain_lookup
        ]
        if not filtered_chains:
            continue

        nodes_data = generate_nodes_from_atom_data(filtered_chains)
        for node in nodes_data:
            original = chain_lookup.get(node["id"])
            if not original:
                continue

            node["uniprot_border_color"] = original.get("uniprot_border_color", "#555555")
            node["color_group"] = original.get("_parent_structure_key") or original.get("_parent_pdb_id") or "Unknown"

            up_id = original.get("uniprot_id")
            if up_id:
                node["uniprot_id"] = up_id
                node["name"] = f"{node['id']}\n{up_id}"
            else:
                node["name"] = node["id"]

        final_edges = []
        for result in results:
            if result["chain_a"] in component_nodes and result["chain_b"] in component_nodes:
                final_edges.append(result)
        for edge in identity_edges:
            if edge["chain_a"] in component_nodes and edge["chain_b"] in component_nodes:
                final_edges.append(edge)

        component_uniprots = {
            chain_lookup[node_id].get("uniprot_id")
            for node_id in component_nodes
            if node_id in chain_lookup and chain_lookup[node_id].get("uniprot_id")
        }
        network_title = make_component_title(
            COMBINED_CHAIN_NETWORK_TITLE, component_uniprots
        )
        skipped = combined_graph_skip(
            network_kind="combined_chain_network",
            name=network_title,
            node_count=len(nodes_data),
            edge_count=len(final_edges),
            limits=graph_limits,
        )
        if skipped:
            skipped_outputs.append(skipped)
            logger.warning(
                "Skipping %s: %s nodes/%s edges exceed configured combined graph limits",
                network_title,
                len(nodes_data),
                len(final_edges),
            )
            continue
        logger.info(
            "Exporting component %s/%s: %s (%s chains, %s edges)",
            index,
            len(component_node_sets),
            network_title,
            len(filtered_chains),
            len(final_edges),
        )
        create_cytoscape_network(final_edges, network_title, run_output_path, nodes_data=nodes_data)

    return skipped_outputs


def run_pipeline(input_path_or_filelist: Union[str, List[str]], web_output_dir: str | None = None) -> RunOutputPaths:
    """Run the full single-process pipeline for a folder or explicit file list."""
    network_config = config["networks"]
    start_time_total = time.time()
    started_at = datetime.now().isoformat(timespec="seconds")
    timings = PipelineTimings()
    output_paths: RunOutputPaths | None = None
    web_output_reservation_attempted = False
    web_output_reserved = False

    try:
        _validate_output_root()
        if isinstance(input_path_or_filelist, list) and web_output_dir:
            raise InputValidationError(
                "WEB_OUTPUT_UNAVAILABLE_FOR_FILE_LIST",
                "Web output requires a folder-based run.",
            )
        _validate_analysis_config()
        network_annotation_config()
        file_paths = discover_input_files(input_path_or_filelist)
        _validate_web_reference_manifest(web_output_dir)
        if web_output_dir:
            web_output_reservation_attempted = True
            reserve_web_output_directory(web_output_dir)
            web_output_reserved = True
        inventory = inspect_input_files(file_paths)
        _combined_graph_limits()
        preflight_identities = _preflight_structure_identities(file_paths, inventory)
        _validate_required_reference_files()
        output_paths = create_run_output_paths(config["output_path"])
        limits = _resource_limits()
        detailed_interaction_budget = (
            DetailedInteractionBudget(
                max_rows=limits["max_detailed_interaction_rows"],
                max_bytes=limits["max_detailed_interaction_bytes"],
                min_free_bytes=limits["min_output_free_bytes"],
            )
            if config.get("export_detailed_interactions", False)
            else None
        )

        logger.info(
            "Run started with %s input file(s), %s total bytes",
            len(file_paths),
            inventory.total_bytes,
        )

        combined_data: List[Dict[str, Any]] = []
        results: List[Dict[str, Any]] = []
        processing_batches = 0
        parsing_workers_used = 0
        seen_chain_ids: dict[str, str] = {}
        signatures_by_path = {
            path: signature
            for path, signature in zip(file_paths, inventory.file_signatures or ())
        }
        digests_by_path = {
            path: digest
            for path, digest in zip(file_paths, inventory.file_sha256 or ())
        }
        maximum_expanded_bytes = _resource_limits()["max_single_input_expanded_bytes"]
        for processing_batches, batch_paths in enumerate(
            create_processing_batches(file_paths, inventory),
            start=1,
        ):
            batch_bytes = sum(Path(path).stat().st_size for path in batch_paths)
            logger.info(
                "Processing bounded batch %s (%s files, %s bytes)",
                processing_batches,
                len(batch_paths),
                batch_bytes,
            )
            parsing_workers_used = max(
                parsing_workers_used,
                min(
                    resolve_workers(
                        config.get("workers", {}).get("parsing"),
                        kind="parsing",
                    ),
                    len(batch_paths),
                ),
            )

            start_time = time.time()
            batch_data = _parse_input_files(
                batch_paths,
                expected_identities=preflight_identities,
                expected_signatures=signatures_by_path,
                expected_sha256=digests_by_path,
                maximum_expanded_bytes=maximum_expanded_bytes,
            )
            timings.parsing += time.time() - start_time
            if not batch_data:
                continue
            if len(batch_data) != len(batch_paths):
                raise InputValidationError(
                    "NO_PARSEABLE_STRUCTURES",
                    "One or more structure inputs did not produce a parsed structure.",
                )
            for source_path, parsed_structure in zip(batch_paths, batch_data):
                expected_identity = preflight_identities[source_path]
                if isinstance(parsed_structure.get("structure_identity"), Mapping):
                    current_identity = StructureIdentity.from_mapping(
                        parsed_structure["structure_identity"]
                    )
                    if current_identity.key != expected_identity.key:
                        raise InputValidationError(
                            "INPUT_CHANGED_DURING_PROCESSING",
                            f"Input identity changed while parsing {Path(source_path).name}.",
                        )
                expected_signature = signatures_by_path.get(source_path)
                if expected_signature is not None and input_file_signature(source_path) != expected_signature:
                    raise InputValidationError(
                        "INPUT_CHANGED_DURING_PROCESSING",
                        f"Input file changed while parsing {Path(source_path).name}.",
                    )

            _register_unique_chain_ids(batch_data, seen_chain_ids)

            start_time = time.time()
            process_molecule_info(batch_data)
            timings.sifts += time.time() - start_time

            start_time = time.time()
            _run_blast_annotation(batch_data)
            timings.blast += time.time() - start_time

            start_time = time.time()
            batch_results = calculate_distances_with_ckdtree(batch_data)
            timings.interaction += time.time() - start_time

            if config.get("export_detailed_interactions", False):
                if detailed_interaction_budget is None:
                    raise RuntimeError("Detailed interaction budget was not initialized")
                _export_detailed_interaction_tables(
                    batch_data,
                    batch_results,
                    output_paths.distances_dir,
                    detailed_interaction_budget,
                )

            results.extend(batch_results)
            combined_data.extend(_compact_structure_summaries(batch_data))

            tree_cache.clear()
            coords_cache.clear()
            del batch_results
            del batch_data
            gc.collect()

        if not combined_data:
            raise InputValidationError(
                "NO_PARSEABLE_STRUCTURES",
                "No input structure could be parsed successfully.",
            )
        if not any(structure.get("atom_data") for structure in combined_data):
            raise InputValidationError(
                "NO_VALID_CHAINS",
                "Parsed structures did not contain supported protein or nucleic-acid polymer chains.",
            )

        skipped_outputs = _export_network_outputs(
            combined_data,
            results,
            network_config,
            output_paths,
            timings,
        )
        warnings = [
            {
                "code": "COMBINED_GRAPH_LIMIT_EXCEEDED",
                "message": (
                    f"Skipped {entry['name']} with {entry['counts']['nodes']} nodes and "
                    f"{entry['counts']['edges']} edges because configured combined graph limits were exceeded."
                ),
                "output": entry["name"],
            }
            for entry in skipped_outputs
        ] + _embedded_annotation_warnings(combined_data)

        total_time = time.time() - start_time_total
        write_runtime_analysis(output_paths.log_file, timings.as_dict(), total_time)
        finished_at = datetime.now().isoformat(timespec="seconds")
        generated_outputs = collect_generated_outputs(output_paths)
        extra_counts = {
            "structures": len(combined_data),
            "chains": sum(len(structure.get("atom_data", [])) for structure in combined_data),
            "interactions": len(results),
            "skipped_outputs": len(skipped_outputs),
        }
        write_run_manifest(
            output_paths.manifest_file,
            input_files=file_paths,
            output_paths=output_paths,
            config_snapshot=_config_snapshot(network_config),
            status="success",
            started_at=started_at,
            finished_at=finished_at,
            total_time=total_time,
            input_path=str(input_path_or_filelist) if not isinstance(input_path_or_filelist, list) else None,
            generated_outputs=generated_outputs,
            extra_counts=extra_counts,
            annotations=_annotation_summary(combined_data),
            references=_reference_summary(),
            resources=_resource_summary(
                inventory,
                processing_batches=processing_batches,
                parsing_workers_used=parsing_workers_used,
                detailed_interactions=detailed_interaction_budget,
            ),
            skipped_outputs=skipped_outputs,
            warnings=warnings,
            identities=_identity_summary(combined_data),
            structure_inputs=_structure_input_summary(combined_data),
        )
        write_run_summary(output_paths)
        if web_output_dir:
            collect_web_outputs(
                output_paths,
                web_output_dir,
                web_output_prepared=True,
            )
        logger.info("Run finished successfully in %.1f seconds", total_time)
        return output_paths
    except Exception as exc:
        logger.error("Run failed: %s", exc)
        base_output_path = str(config.get("output_path") or "").strip()
        if not isinstance(input_path_or_filelist, list) and base_output_path:
            output_paths = write_failed_run_manifest(
                base_output_path,
                input_path=str(input_path_or_filelist),
                config_snapshot=_config_snapshot(network_config),
                error=exc,
                started_at=started_at,
                output_paths=output_paths,
            )
            if web_output_dir and (
                web_output_reserved or not web_output_reservation_attempted
            ):
                collect_web_outputs(
                    output_paths,
                    web_output_dir,
                    web_output_prepared=web_output_reserved,
                )
        raise
    finally:
        tree_cache.clear()
        coords_cache.clear()
        gc.collect()
