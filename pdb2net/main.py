"""PDB2Net main orchestration.

This module coordinates the end-to-end pipeline:
1) Input discovery (files/folders) and per-file parsing with Gemmi.
2) Lightweight structure processing (atoms/residues per chain).
3) Molecule classification (SIFTS / pdb_seqres-driven; then BLAST fallback).
4) Interaction detection (distance-based with cKDTree).
5) Optional detailed CSV export of atom-level interactions.
6) Network export(s): chain-level per PDB, combined chain network, and protein-level networks.

Notes
-----
- Behavior is preserved; only documentation, messages, and style were unified.
- Worker counts are derived via `_resolve_workers()` using simple CPU heuristics.
- Headless mode vs. Cytoscape UI is controlled by `config["open_in_cytoscape"]`.
"""

from __future__ import annotations

import gc
import multiprocessing
import os
import time
from datetime import datetime
from typing import Any, Dict, Iterable, Iterator, List, Optional

from config_loader import config
from cytoscape import create_cytoscape_network, generate_nodes_from_atom_data
from data_processor import process_structure
from detailed_results_exporter import export_detailed_interactions
from distances import calculate_distances_with_ckdtree, tree_cache, coords_cache
from file_parser import get_pdb_id, is_valid_file, parse_structure
from protein_network import create_protein_network
from uniprot_matcher import parallel_blast_search
from unknown_molecule_uniprot import process_molecule_info


def _resolve_workers(value: Optional[int | str], *, kind: str = "parsing") -> int:
    """Convert a config value to a sensible worker/thread count.

    Rules
    -----
    - If value is an int → use `max(1, value)`.
    - If value is "auto" / None / missing → derive from CPU count:
        * parsing: max(1, cpu - 1)
        * blast:   max(1, cpu - 2)
        * default: max(1, cpu - 1)

    Parameters
    ----------
    value : int | str | None
        Configured value, often "auto" or an integer.
    kind : {"parsing", "blast"}
        Hint to select the heuristic for "auto".

    Returns
    -------
    int
        Chosen worker count (>= 1).
    """
    cpu = os.cpu_count() or 4

    if isinstance(value, int):
        return max(1, value)

    if kind == "parsing":
        return max(1, cpu - 1)
    if kind == "blast":
        return max(1, cpu - 2)

    return max(1, cpu - 1)


def process_single_file(file_path: str) -> Optional[Dict[str, Any]]:
    """Process a single PDB/mmCIF file end-to-end for the parsing stage.

    Steps
    -----
    - Extract PDB ID.
    - Parse structure via Gemmi.
    - Extract atoms/residues using `process_structure()`.

    Returns
    -------
    dict | None
        Structure payload compatible with downstream steps on success, else None.
    """
    try:
        pdb_id = get_pdb_id(file_path)
        structure = parse_structure(file_path, pdb_id)
        if not structure:
            return None
        return process_structure({"file_path": file_path, "pdb_id": pdb_id, "structure": structure})
    except Exception as e:
        print(f"Error while processing {file_path}: {e}")
        return None


def run_main(batch_files: List[str]) -> None:
    """Helper wrapper to re-enter main() for a batch (used by batch_run)."""
    from main import main  # local import to preserve original behavior
    main(batch_files)


def main(input_path_or_filelist: str | List[str]) -> None:
    """Entry point for processing a folder or an explicit file list.

    Parameters
    ----------
    input_path_or_filelist : str | list[str]
        Either a directory path containing input files, or an explicit list of file paths.
    """
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
            export_detailed_interactions(structure_data, pdb_interactions, run_output_path)

    # === Network exports: chain-level per PDB ===
    if network_config["chain_per_pdb"]:
        start_time = time.time()
        results_by_pdb: Dict[str, List[Dict[str, Any]]] = {}
        for entry in results:
            pdb_id = entry["chain_a"].split(":")[0]
            results_by_pdb.setdefault(pdb_id, []).append(entry)

        # 1) PDBs with interactions (as before)
        for pdb_id, pdb_results in results_by_pdb.items():
            structure = next((s for s in combined_data if s["pdb_id"] == pdb_id), None)
            if not structure:
                continue
            nodes_data = generate_nodes_from_atom_data(structure["atom_data"], pdb_id)
            network_title = f"Chain_Interaction_Network_{pdb_id}"
            create_cytoscape_network(pdb_results, network_title, run_output_path, nodes_data=nodes_data)

        # 2) PDBs without interactions → nodes-only networks
        for structure in combined_data:
            pdb_id = structure["pdb_id"]
            if pdb_id not in results_by_pdb:
                nodes_data = generate_nodes_from_atom_data(structure["atom_data"], pdb_id)
                network_title = f"Chain_Interaction_Network_{pdb_id}"
                create_cytoscape_network([], network_title, run_output_path, nodes_data=nodes_data)

        sum_times["networks"] += time.time() - start_time

    # === Network export: combined chain network ===
    if network_config["combined_chain_network"]:
        start_time = time.time()
        all_chains = [chain for structure in combined_data for chain in structure["atom_data"]]
        nodes_data = generate_nodes_from_atom_data(all_chains)
        create_cytoscape_network(results, "Combined_Network", run_output_path, nodes_data=nodes_data)
        sum_times["networks"] += time.time() - start_time

    # === Network export: protein-based (per PDB and/or combined) ===
    if network_config["protein_per_pdb"] or network_config["combined_protein_network"]:
        start_time = time.time()
        create_protein_network(results, combined_data, run_output_path, network_config)
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
    """Yield batches of files with a soft size cap, streaming over the input.

    The function does not precompute sizes for the entire set; it yields
    batches as it scans the input sequence to keep memory usage low.

    Parameters
    ----------
    file_paths : Iterable[str]
        Stream of file paths to be grouped.
    max_batch_kb : int
        Maximum accumulated size (in KB) per batch.

    Yields
    ------
    list[str]
        A batch of file paths whose total size is within `max_batch_kb`.
    """
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
    """Run the pipeline repeatedly on streamed batches from a folder.

    Each batch is processed either:
    - Inline in the main process when `open_in_cytoscape` is set (to avoid
      concurrent access to Cytoscape/py4cytoscape logging), or
    - In a subprocess with a timeout (default) to protect against stuck runs.

    Parameters
    ----------
    input_folder : str
        Folder to scan for valid input files.
    timeout_minutes : int, default 10
        Per-batch timeout for the subprocess mode.
    size_limit_kb : int, default 900_800
        Soft limit for total batch size.
    """
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

        # If Cytoscape/py4cytoscape is used, do not fork a subprocess to
        # avoid concurrent rotation/access of the Cytoscape log files.
        if config.get("open_in_cytoscape", False):
            start_time_batch = time.time()
            run_main(batch_files)  # run in the main process
            duration = time.time() - start_time_batch
            total_done += len(batch_files)
            print(f"Batch {i} finished in {duration:.1f} seconds.")
            print(f"Processed so far: {total_done} files.")
            continue

        # Default path: run the batch in a subprocess with a timeout.
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
