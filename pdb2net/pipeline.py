"""Pipeline orchestration helpers for single-run processing."""

from __future__ import annotations

import gc
import os
import shutil
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from .config_loader import config
from .components import build_identity_edges, find_linked_components, make_component_title
from .cytoscape import create_cytoscape_network, generate_nodes_from_atom_data
from .data_processor import process_structure
from .detailed_results_exporter import export_detailed_interactions
from .distances import calculate_distances_with_ckdtree, coords_cache, tree_cache
from .file_parser import get_pdb_id, is_valid_file, parse_structure
from .input_contract import InputValidationError
from .logging_utils import get_logger
from .outputs import (
    RunOutputPaths,
    collect_generated_outputs,
    collect_web_outputs,
    create_run_output_paths,
    write_failed_run_manifest,
    write_run_manifest,
    write_run_summary,
    write_runtime_analysis,
)
from .protein_network import create_protein_network
from .residue_types import NUCLEIC_ACID_TYPES
from .uniprot_matcher import extract_direct_uniprot_accession, parallel_blast_search
from .unknown_molecule_uniprot import process_molecule_info


logger = get_logger(__name__)


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


def process_single_file(file_path: str) -> Optional[Dict[str, Any]]:
    """Process a single PDB/mmCIF file end-to-end for the parsing stage."""
    try:
        pdb_id = get_pdb_id(file_path)
        structure = parse_structure(file_path, pdb_id)
        if not structure:
            return None
        return process_structure({"file_path": file_path, "pdb_id": pdb_id, "structure": structure})
    except Exception as exc:
        print(f"Error while processing {file_path}: {exc}")
        return None


def discover_input_files(input_path_or_filelist: Union[str, List[str]]) -> List[str]:
    """Resolve a folder or internal file list into valid structure files."""
    if isinstance(input_path_or_filelist, list):
        return [file_path for file_path in input_path_or_filelist if is_valid_file(file_path)]

    input_path = Path(input_path_or_filelist)
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

    file_paths = [
        entry.path for entry in os.scandir(str(input_path))
        if entry.is_file() and is_valid_file(entry.path)
    ]
    if not file_paths:
        raise InputValidationError(
            "NO_VALID_INPUT_FILES",
            f"input_folder_path contains no supported .pdb, .cif, or .mmcif files: {input_path}",
        )

    return file_paths


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


def _parse_input_files(file_paths: List[str]) -> List[Dict[str, Any]]:
    """Parse and preprocess structure files in parallel."""
    parsing_workers = resolve_workers(config.get("workers", {}).get("parsing"), kind="parsing")
    logger.info("[Workers] Parsing processes: %s", parsing_workers)
    with ProcessPoolExecutor(max_workers=parsing_workers) as executor:
        return list(filter(None, executor.map(process_single_file, file_paths)))


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


def _export_detailed_interaction_tables(
    combined_data: List[Dict[str, Any]],
    results: List[Dict[str, Any]],
    distances_dir: str,
) -> None:
    """Write optional atom-level interaction CSVs."""
    for structure_data in combined_data:
        pdb_id = structure_data["pdb_id"]
        pdb_interactions = [result for result in results if result["chain_a"].startswith(pdb_id)]
        export_detailed_interactions(structure_data, pdb_interactions, distances_dir)


def _export_chain_networks(
    combined_data: List[Dict[str, Any]],
    results: List[Dict[str, Any]],
    chain_dir: str,
) -> None:
    """Export chain interaction networks per PDB, including nodes-only cases."""
    results_by_pdb: Dict[str, List[Dict[str, Any]]] = {}
    for entry in results:
        pdb_id = entry["chain_a"].split(":")[0]
        results_by_pdb.setdefault(pdb_id, []).append(entry)

    for pdb_id, pdb_results in results_by_pdb.items():
        structure = next((structure for structure in combined_data if structure["pdb_id"] == pdb_id), None)
        if not structure:
            continue

        _annotate_chain_parent_metadata(structure)
        nodes_data = generate_nodes_from_atom_data(structure["atom_data"], pdb_id)
        network_title = f"Chain_Interaction_Network_{pdb_id}"
        create_cytoscape_network(
            pdb_results,
            network_title,
            chain_dir,
            nodes_data=nodes_data,
        )

    for structure in combined_data:
        pdb_id = structure["pdb_id"]
        if pdb_id not in results_by_pdb:
            _annotate_chain_parent_metadata(structure)
            nodes_data = generate_nodes_from_atom_data(structure["atom_data"], pdb_id)
            network_title = f"Chain_Interaction_Network_{pdb_id}"
            create_cytoscape_network(
                [],
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
) -> None:
    """Run the configured network export stages and accumulate timings."""
    if network_config["chain_per_pdb"]:
        start_time = time.time()
        _export_chain_networks(combined_data, results, output_paths.chain_dir)
        timings.networks += time.time() - start_time

    if network_config["combined_chain_network"]:
        start_time = time.time()
        _create_linked_identity_network(results, combined_data, output_paths.combined_dir)
        timings.networks += time.time() - start_time

    if network_config["protein_per_pdb"] or network_config["combined_protein_network"]:
        start_time = time.time()
        create_protein_network(
            results,
            combined_data,
            output_paths.run_output_path,
            network_config,
            per_pdb_output_path=output_paths.protein_dir,
            combined_output_path=output_paths.combined_dir,
        )
        timings.networks += time.time() - start_time


def _config_snapshot(network_config: Dict[str, Any]) -> Dict[str, Any]:
    """Return the small config subset useful for manifests and webserver jobs."""
    return {
        "networks": network_config,
        "distance_thresholds": config.get("distance_thresholds", {}),
        "interaction_filters": config.get("interaction_filters", {}),
        "structure_model_policy": config.get("structure_model_policy", "first"),
        "workers": config.get("workers", {}),
        "layout_mode": config.get("layout_mode"),
        "open_in_cytoscape": config.get("open_in_cytoscape"),
        "export_detailed_interactions": config.get("export_detailed_interactions", False),
        "blast_cache_path": config.get("blast_cache_path", ""),
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

    for chain in structure["atom_data"]:
        chain["_parent_pdb_id"] = pdb_id
        chain["_parent_file_path"] = file_path
        chain["_parent_file_label"] = file_label


def _create_linked_identity_network(
    results: List[Dict[str, Any]],
    combined_data: List[Dict[str, Any]],
    run_output_path: str,
) -> None:
    """Build separate linked-identity combined networks."""
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
        return

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
            node["color_group"] = original.get("_parent_file_path") or original.get("_parent_pdb_id") or "Unknown"

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
        network_title = make_component_title("Combined_Network", component_uniprots)
        logger.info(
            "Exporting component %s/%s: %s (%s chains, %s edges)",
            index,
            len(component_node_sets),
            network_title,
            len(filtered_chains),
            len(final_edges),
        )
        create_cytoscape_network(final_edges, network_title, run_output_path, nodes_data=nodes_data)


def run_pipeline(input_path_or_filelist: Union[str, List[str]], web_output_dir: str | None = None) -> RunOutputPaths:
    """Run the full single-process pipeline for a folder or explicit file list."""
    network_config = config["networks"]
    start_time_total = time.time()
    started_at = datetime.now().isoformat(timespec="seconds")
    timings = PipelineTimings()

    try:
        _validate_output_root()
        file_paths = discover_input_files(input_path_or_filelist)
        _validate_required_reference_files()
        output_paths = create_run_output_paths(config["output_path"])

        logger.info("Run started with %s input file(s)", len(file_paths))

        start_time = time.time()
        combined_data = _parse_input_files(file_paths)
        timings.parsing = time.time() - start_time
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

        start_time = time.time()
        process_molecule_info(combined_data)
        timings.sifts = time.time() - start_time

        start_time = time.time()
        _run_blast_annotation(combined_data)
        timings.blast = time.time() - start_time

        start_time = time.time()
        results = calculate_distances_with_ckdtree(combined_data)
        timings.interaction = time.time() - start_time

        if config.get("export_detailed_interactions", False):
            _export_detailed_interaction_tables(combined_data, results, output_paths.distances_dir)

        _export_network_outputs(combined_data, results, network_config, output_paths, timings)

        total_time = time.time() - start_time_total
        write_runtime_analysis(output_paths.log_file, timings.as_dict(), total_time)
        finished_at = datetime.now().isoformat(timespec="seconds")
        generated_outputs = collect_generated_outputs(output_paths)
        extra_counts = {
            "structures": len(combined_data),
            "chains": sum(len(structure.get("atom_data", [])) for structure in combined_data),
            "interactions": len(results),
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
        )
        write_run_summary(output_paths)
        if web_output_dir:
            collect_web_outputs(output_paths, web_output_dir)
        logger.info("Run finished successfully in %.1f seconds", total_time)
        return output_paths
    except Exception as exc:
        logger.error("Run failed: %s", exc)
        if not isinstance(input_path_or_filelist, list):
            output_paths = write_failed_run_manifest(
                config["output_path"],
                input_path=str(input_path_or_filelist),
                config_snapshot=_config_snapshot(network_config),
                error=exc,
                started_at=started_at,
            )
            if web_output_dir:
                collect_web_outputs(output_paths, web_output_dir)
        raise
    finally:
        tree_cache.clear()
        coords_cache.clear()
        gc.collect()
