"""PDB2Net main orchestration.

This module coordinates the end-to-end pipeline:
1) Input discovery (files/folders) and per-file parsing with Gemmi.
2) Lightweight structure processing (atoms/residues per chain).
3) Molecule classification (SIFTS / pdb_seqres-driven; then BLAST fallback).
4) Interaction detection (distance-based with cKDTree).
5) Optional detailed CSV export of atom-level interactions.
6) Network export(s): chain-level per PDB, combined linked-identity network, and protein-level networks.

Notes
-----
- Behavior is preserved; only documentation, messages, and style were unified.
- Worker counts are derived via `_resolve_workers()` using simple CPU heuristics.
- Headless mode vs. Cytoscape UI is controlled by `config["open_in_cytoscape"]`.
"""

from __future__ import annotations

import gc
import hashlib
import multiprocessing
import os
import time
from datetime import datetime
from typing import Any, Dict, Iterable, Iterator, List, Optional, Tuple, Union

from config_loader import config
from cytoscape import create_cytoscape_network, generate_nodes_from_atom_data
from data_processor import process_structure
from detailed_results_exporter import export_detailed_interactions
from distances import calculate_distances_with_ckdtree, tree_cache, coords_cache
from file_parser import get_pdb_id, is_valid_file, parse_structure
from protein_network import create_protein_network
from uniprot_matcher import parallel_blast_search
from unknown_molecule_uniprot import process_molecule_info


def _resolve_workers(value: Optional[Union[int, str]], *, kind: str = "parsing") -> int:
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
    except Exception as e:
        print(f"Error while processing {file_path}: {e}")
        return None


def _get_border_color_for_uniprot(uniprot_id: Optional[str]) -> str:
    """Generate a consistent hex color from a UniProt ID string (or return gray)."""
    if not uniprot_id:
        return "#555555"  # Dark gray for unknown/no ID
    # Generate deterministic color via MD5 hash of the ID
    hash_val = int(hashlib.md5(uniprot_id.encode('utf-8')).hexdigest(), 16)
    # Use the lower 24 bits for RGB
    return f"#{hash_val & 0xFFFFFF:06x}"


def _create_linked_identity_network(results: List[Dict[str, Any]], combined_data: List[Dict[str, Any]], run_output_path: str) -> None:
    """
    Internal helper to build separate linked-identity combined networks.

    Updated behavior:
    - Builds identity bridges only between chains from different files.
    - Splits the dataset into connected interlinked components.
    - Exports one Combined_Network per interlinked component.
    - Names each network using the UniProt IDs present in that component.
    """
    print("  > Building Linked Identity Network (bridging different PDBs)...")

    all_chains_flat = []
    uniprot_groups: Dict[str, List[Dict[str, Any]]] = {}

    for structure in combined_data:
        pdb_id = structure["pdb_id"]
        file_path = structure.get("file_path", "")
        file_label = os.path.basename(file_path) if file_path else pdb_id

        for chain in structure["atom_data"]:
            chain["_parent_pdb_id"] = pdb_id
            chain["_parent_file_path"] = file_path
            chain["_parent_file_label"] = file_label

            up_id = chain.get("uniprot_id")
            chain["uniprot_border_color"] = _get_border_color_for_uniprot(up_id)

            all_chains_flat.append(chain)

            if up_id:
                uniprot_groups.setdefault(up_id, []).append(chain)

    chain_lookup = {c["unique_chain_id"]: c for c in all_chains_flat}

    identity_edges = []
    count_identity_edges = 0

    for up_id, chains in uniprot_groups.items():
        if len(chains) < 2:
            continue

        chains = sorted(chains, key=lambda c: c["unique_chain_id"])

        for i in range(len(chains) - 1):
            chain_a = chains[i]
            chain_b = chains[i + 1]

            pdb_a = chain_a["_parent_pdb_id"]
            pdb_b = chain_b["_parent_pdb_id"]

            if pdb_a != pdb_b:
                identity_edges.append({
                    "chain_a": chain_a["unique_chain_id"],
                    "chain_b": chain_b["unique_chain_id"],
                    "all_atoms_count": 1000,
                    "interaction_type": "identity",
                })
                count_identity_edges += 1

    print(f"    - Added {count_identity_edges} identity bridges between files.")

    if not identity_edges:
        print("    - No interlinked components found. Skipping Combined_Network export.")
        return

    adjacency: Dict[str, set[str]] = {}

    def add_edge_to_adjacency(a: str, b: str) -> None:
        adjacency.setdefault(a, set()).add(b)
        adjacency.setdefault(b, set()).add(a)

    for res in results:
        add_edge_to_adjacency(res["chain_a"], res["chain_b"])

    for edge in identity_edges:
        add_edge_to_adjacency(edge["chain_a"], edge["chain_b"])

    identity_nodes = set()
    for edge in identity_edges:
        identity_nodes.add(edge["chain_a"])
        identity_nodes.add(edge["chain_b"])

    visited = set()
    component_node_sets: List[set[str]] = []

    for start_node in sorted(identity_nodes):
        if start_node in visited:
            continue

        stack = [start_node]
        component_nodes = set()

        while stack:
            node = stack.pop()
            if node in visited:
                continue

            visited.add(node)
            component_nodes.add(node)

            for neighbor in adjacency.get(node, set()):
                if neighbor not in visited:
                    stack.append(neighbor)

        if component_nodes:
            component_node_sets.append(component_nodes)

    def _sanitize_filename_part(text: str) -> str:
        allowed = []
        for ch in str(text):
            if ch.isalnum() or ch in {"-", "_"}:
                allowed.append(ch)
            else:
                allowed.append("_")
        sanitized = "".join(allowed).strip("_")
        return sanitized or "Unknown"

    def _make_component_network_title(component_nodes: set[str]) -> str:
        component_uniprots = sorted({
            chain_lookup[node_id].get("uniprot_id")
            for node_id in component_nodes
            if node_id in chain_lookup and chain_lookup[node_id].get("uniprot_id")
        })

        if not component_uniprots:
            return "Combined_Network_Unknown"

        sanitized_ids = [_sanitize_filename_part(up_id) for up_id in component_uniprots]
        preview = "_".join(sanitized_ids[:5])
        if len(sanitized_ids) > 5:
            digest_source = "|".join(sanitized_ids).encode("utf-8")
            digest = hashlib.md5(digest_source).hexdigest()[:8]
            return f"Combined_Network_{preview}__{digest}"
        return f"Combined_Network_{preview}"

    print(f"    - Found {len(component_node_sets)} interlinked component(s).")

    for idx, component_nodes in enumerate(component_node_sets, start=1):
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
            if original:
                node["uniprot_border_color"] = original.get("uniprot_border_color", "#555555")
                node["color_group"] = original.get("_parent_file_path") or original.get("_parent_pdb_id") or "Unknown"

                up_id = original.get("uniprot_id")
                if up_id:
                    node["uniprot_id"] = up_id
                    node["name"] = f"{node['id']}\n{up_id}"
                else:
                    node["name"] = node["id"]

        final_edges = []

        for res in results:
            if res["chain_a"] in component_nodes and res["chain_b"] in component_nodes:
                final_edges.append(res)

        for edge in identity_edges:
            if edge["chain_a"] in component_nodes and edge["chain_b"] in component_nodes:
                final_edges.append(edge)

        network_title = _make_component_network_title(component_nodes)

        print(
            f"    - Exporting component {idx}/{len(component_node_sets)}: "
            f"{network_title} ({len(filtered_chains)} chains, {len(final_edges)} edges)"
        )

        create_cytoscape_network(final_edges, network_title, run_output_path, nodes_data=nodes_data)

def run_main(batch_files: List[str]) -> None:
    """Helper wrapper to re-enter main() for a batch (used by batch_run)."""
    from main import main  # local import to preserve original behavior
    main(batch_files)


def main(input_path_or_filelist: Union[str, List[str]]) -> None:
    """Entry point for processing a folder or an explicit file list."""
    network_config = config["networks"]

    if isinstance(input_path_or_filelist, list):
        file_paths = [f for f in input_path_or_filelist if is_valid_file(f)]
    else:
        file_paths = [
            entry.path for entry in os.scandir(input_path_or_filelist)
            if entry.is_file() and is_valid_file(entry.path)
        ]

    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    run_output_path = os.path.join(config["output_path"], timestamp)
    os.makedirs(run_output_path, exist_ok=True)

    # --- NEW: fixed subfolders per run ---
    combined_dir = os.path.join(run_output_path, "combined")
    protein_dir = os.path.join(run_output_path, "protein")
    chain_dir = os.path.join(run_output_path, "chain")
    distances_dir = os.path.join(run_output_path, "distances")
    for d in (combined_dir, protein_dir, chain_dir, distances_dir):
        os.makedirs(d, exist_ok=True)

    start_time_total = time.time()
    sum_times = {"parsing": 0.0, "sifts": 0.0, "blast": 0.0, "interaction": 0.0, "networks": 0.0}

    # === Parsing (parallel via process pool) ===
    from concurrent.futures import ProcessPoolExecutor

    start_time = time.time()
    t_parse = time.time()
    parsing_workers = _resolve_workers(config.get("workers", {}).get("parsing"), kind="parsing")
    print(f"[Workers] Parsing processes: {parsing_workers}")
    with ProcessPoolExecutor(max_workers=parsing_workers) as executor:
        combined_data = list(filter(None, executor.map(process_single_file, file_paths)))
    parsing_duration = time.time() - t_parse
    sum_times["parsing"] = time.time() - start_time

    # === Molecule typing via SIFTS / pdb_seqres.txt ===
    start_time = time.time()
    process_molecule_info(combined_data)
    sum_times["sifts"] = time.time() - start_time

    # === Classification via BLAST (for unresolved IDs) ===
    start_time = time.time()
    blast_workers = _resolve_workers(config.get("workers", {}).get("blast_threads"), kind="blast")
    print(f"[Workers] BLAST threads: {blast_workers}")
    parallel_blast_search(combined_data, max_workers=blast_workers)
    sum_times["blast"] = time.time() - start_time

    # === Interaction detection ===
    start_time = time.time()
    results = calculate_distances_with_ckdtree(combined_data)
    sum_times["interaction"] = time.time() - start_time

    # === Optional: atom-level detailed CSV export ===
    if config.get("export_detailed_interactions", False):
        for structure_data in combined_data:
            pdb_id = structure_data["pdb_id"]
            pdb_interactions = [res for res in results if res["chain_a"].startswith(pdb_id)]
            # NEW: write all distance outputs into one shared folder
            export_detailed_interactions(structure_data, pdb_interactions, distances_dir)

    # === Network exports: chain-level per PDB ===
    if network_config["chain_per_pdb"]:
        start_time = time.time()
        results_by_pdb: Dict[str, List[Dict[str, Any]]] = {}
        for entry in results:
            pdb_id = entry["chain_a"].split(":")[0]
            results_by_pdb.setdefault(pdb_id, []).append(entry)

        # 1) PDBs with interactions
        for pdb_id, pdb_results in results_by_pdb.items():
            structure = next((s for s in combined_data if s["pdb_id"] == pdb_id), None)
            if not structure:
                continue
            nodes_data = generate_nodes_from_atom_data(structure["atom_data"], pdb_id)
            network_title = f"Chain_Interaction_Network_{pdb_id}"
            # NEW: chain networks folder
            create_cytoscape_network(pdb_results, network_title, chain_dir, nodes_data=nodes_data)

        # 2) PDBs without interactions → nodes-only networks
        for structure in combined_data:
            pdb_id = structure["pdb_id"]
            if pdb_id not in results_by_pdb:
                nodes_data = generate_nodes_from_atom_data(structure["atom_data"], pdb_id)
                network_title = f"Chain_Interaction_Network_{pdb_id}"
                # NEW: chain networks folder
                create_cytoscape_network([], network_title, chain_dir, nodes_data=nodes_data)

        sum_times["networks"] += time.time() - start_time

    # === Network export: combined chain network (LINKED MODE) ===
    if network_config["combined_chain_network"]:
        start_time = time.time()
        # NEW: combined networks folder
        _create_linked_identity_network(results, combined_data, combined_dir)
        sum_times["networks"] += time.time() - start_time

    # === Network export: protein-based (per PDB and/or combined) ===
    if network_config["protein_per_pdb"] or network_config["combined_protein_network"]:
        start_time = time.time()
        # NEW: split per-PDB vs combined outputs
        create_protein_network(
            results,
            combined_data,
            run_output_path,  # kept for backward-compat, but no longer used for exports
            network_config,
            per_pdb_output_path=protein_dir,
            combined_output_path=combined_dir,
        )
        sum_times["networks"] += time.time() - start_time

    # === Timing/logging ===
    total_time = time.time() - start_time_total
    log_file = os.path.join(run_output_path, "runtime_analysis.txt")

    with open(log_file, "w", encoding="utf-8") as f:
        f.write("===== PDB2Net Batch Log =====\n\n")
        f.write("Timing (total):\n")
        f.write(f"- Parsing:                 {sum_times['parsing']:.1f} sec\n")
        f.write(f"- Classification (SIFTS):  {sum_times['sifts']:.1f} sec\n")
        f.write(f"- Classification (BLAST):  {sum_times['blast']:.1f} sec\n")
        f.write(f"- Interaction:             {sum_times['interaction']:.1f} sec\n")
        f.write(f"- Network export:          {sum_times['networks']:.1f} sec\n")
        f.write(f"- Total:                   {total_time:.1f} sec\n")
        f.write("\n===============================\n")

    # Free caches/memory used by interaction stage
    tree_cache.clear()
    coords_cache.clear()
    gc.collect()



def create_batches_streaming(file_paths: Iterable[str], max_batch_kb: int) -> Iterator[List[str]]:
    """Yield batches of files with a soft size cap, streaming over the input."""
    current_batch: List[str] = []
    current_size = 0

    for file_path in file_paths:
        try:
            size_kb = os.path.getsize(file_path) // 1024
        except Exception as e:
            print(f"Error reading file size: {file_path} – {e}")
            continue

        if size_kb > max_batch_kb:
            print(f"File too large for a single batch: {file_path} ({size_kb} KB)")
            continue

        if current_size + size_kb > max_batch_kb:
            yield current_batch
            current_batch = [file_path]
            current_size = size_kb
        else:
            current_batch.append(file_path)
            current_size += size_kb

    if current_batch:
        yield current_batch


def batch_run(input_folder: str, timeout_minutes: int = 10, size_limit_kb: int = 1000_100) -> None:
    """Run the pipeline repeatedly on streamed batches from a folder."""
    def stream_valid_files(folder: str) -> Iterator[str]:
        for entry in os.scandir(folder):
            if entry.is_file() and is_valid_file(entry.path):
                yield entry.path

    timeout_seconds = timeout_minutes * 60
    logs_dir = os.path.join(config["output_path"], "error_in_batch_log")
    os.makedirs(logs_dir, exist_ok=True)
    log_path = os.path.join(logs_dir, "skipped_batches.txt")

    start_time_all = time.time()
    total_done = 0

    file_stream = stream_valid_files(input_folder)

    for i, batch_files in enumerate(create_batches_streaming(file_stream, size_limit_kb), start=1):
        print(f"\n--- Processing batch {i} ({len(batch_files)} files) ---")

        if config.get("open_in_cytoscape", False):
            start_time_batch = time.time()
            run_main(batch_files)
            duration = time.time() - start_time_batch
            total_done += len(batch_files)
            print(f"Batch {i} finished in {duration:.1f} seconds.")
            print(f"Processed so far: {total_done} files.")
            continue

        start_time_batch = time.time()
        process = multiprocessing.Process(target=run_main, args=(batch_files,))
        process.start()
        process.join(timeout=timeout_seconds)

        if process.is_alive():
            print(f"Batch {i} exceeded {timeout_minutes} minutes. Terminating.")
            process.terminate()
            process.join()
            with open(log_path, "a", encoding="utf-8") as log_file:
                for file_path in batch_files:
                    pdb_id = get_pdb_id(file_path)
                    log_file.write(f"{pdb_id}\n")
        else:
            duration = time.time() - start_time_batch
            total_done += len(batch_files)
            print(f"Batch {i} finished in {duration:.1f} seconds.")
            print(f"Processed so far: {total_done} files.")

    total_time_all = time.time() - start_time_all
    print(f"\nTotal time for all batches: {total_time_all:.2f} seconds")


if __name__ == "__main__":
    from cytoscape import ensure_cytoscape_running

    if config.get("open_in_cytoscape", True):
        ensure_cytoscape_running()

    INPUT_FOLDER_PATH = config["input_folder_path"]
    batch_run(INPUT_FOLDER_PATH)