import logging
logging.getLogger("py4cytoscape").disabled = True
logging.getLogger("py4cytoscape.detail").disabled = True
import time
import os
import gc
import multiprocessing
from datetime import datetime

from file_parser import is_valid_file, get_pdb_id, parse_structure
from data_processor import process_structure
from unknown_molecule_uniprot import process_molecule_info
from distances import calculate_distances_with_ckdtree, tree_cache, coords_cache
from detailed_results_exporter import export_detailed_interactions
from uniprot_matcher import parallel_blast_search
from config_loader import config

from cytoscape import create_cytoscape_network, generate_nodes_from_atom_data
from protein_network import create_protein_network

def process_single_file(file_path):
    """
    Verarbeitet eine einzelne PDB/mmCIF-Datei vollständig:
    - PDB-ID extrahieren
    - Struktur mit Gemmi parsen
    - Atome/Residuen extrahieren mit process_structure

    Gibt None zurück, wenn Datei fehlschlägt.
    """
    try:
        pdb_id = get_pdb_id(file_path)
        structure = parse_structure(file_path, pdb_id)
        if not structure:
            return None
        return process_structure({
            "file_path": file_path,
            "pdb_id": pdb_id,
            "structure": structure
        })
    except Exception as e:
        print(f"⚠️ Fehler beim Verarbeiten von {file_path}: {e}")
        return None


def run_main(batch_files):
    import logging
    logging.getLogger("py4cytoscape").disabled = True
    logging.getLogger("py4cytoscape.detail").disabled = True
    from main import main
    main(batch_files)

def main(input_path_or_filelist):
    """
    Hauptfunktion für PDB-Verarbeitung. Unterstützt entweder einen Ordnerpfad oder eine Liste von Dateien.
    """
    network_config = config["networks"]

    if isinstance(input_path_or_filelist, list):
        file_paths = [f for f in input_path_or_filelist if is_valid_file(f)]
    else:
        file_paths = [entry.path for entry in os.scandir(input_path_or_filelist)
                      if entry.is_file() and is_valid_file(entry.path)]

    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    run_output_path = os.path.join(config["output_path"], timestamp)
    os.makedirs(run_output_path, exist_ok=True)

    start_time_total = time.time()
    sum_times = {"parsing": 0, "sifts": 0, "blast": 0, "interaction": 0, "networks": 0}

    # === PARSING (parallel über process_single_file) ===
    from concurrent.futures import ProcessPoolExecutor

    start_time = time.time()
    t_parse = time.time()
    with ProcessPoolExecutor(max_workers=12) as executor:
        combined_data = list(filter(None, executor.map(process_single_file, file_paths)))
    parsing_duration = time.time() - t_parse
    sum_times["parsing"] = time.time() - start_time

    # === Molekülklassifikation via SIFTS/pdb_seqres.txt ===
    start_time = time.time()
    process_molecule_info(combined_data)
    sum_times["sifts"] = time.time() - start_time

    # === Klassifikation via BLAST (nur fehlende IDs) ===
    start_time = time.time()
    parallel_blast_search(combined_data, max_workers=4)
    sum_times["blast"] = time.time() - start_time

    # === Interaktionsanalyse ===
    start_time = time.time()
    results = calculate_distances_with_ckdtree(combined_data)
    sum_times["interaction"] = time.time() - start_time

    # === Optional: Export detaillierter Interaktionen ===
    if config.get("export_detailed_interactions", False):
        for structure_data in combined_data:
            pdb_id = structure_data["pdb_id"]
            pdb_interactions = [res for res in results if res["chain_a"].startswith(pdb_id)]
            export_detailed_interactions(structure_data, pdb_interactions, run_output_path)

    # === Netzwerkexporte (Kette pro PDB) ===
    if network_config["chain_per_pdb"]:
        start_time = time.time()
        results_by_pdb = {}
        for entry in results:
            pdb_id = entry["chain_a"].split(":")[0]
            results_by_pdb.setdefault(pdb_id, []).append(entry)

        for pdb_id, pdb_results in results_by_pdb.items():
            structure = next((s for s in combined_data if s["pdb_id"] == pdb_id), None)
            if not structure:
                continue
            nodes_data = generate_nodes_from_atom_data(structure["atom_data"], pdb_id)
            network_title = f"Chain_Interaction_Network_{pdb_id}"
            create_cytoscape_network(pdb_results, network_title, run_output_path, nodes_data=nodes_data)
        sum_times["networks"] += time.time() - start_time

    # === Netzwerkexport: alle Ketten kombiniert ===
    if network_config["combined_chain_network"]:
        start_time = time.time()
        all_chains = [chain for structure in combined_data for chain in structure["atom_data"]]
        nodes_data = generate_nodes_from_atom_data(all_chains)
        create_cytoscape_network(results, "Combined_Network", run_output_path, nodes_data=nodes_data)
        sum_times["networks"] += time.time() - start_time

    # === Netzwerkexport: Proteinbasiert (per PDB oder kombiniert) ===
    if network_config["protein_per_pdb"] or network_config["combined_protein_network"]:
        start_time = time.time()
        create_protein_network(results, combined_data, run_output_path, network_config)
        sum_times["networks"] += time.time() - start_time

    # === Gesamtzeit + Logging ===
    total_time = time.time() - start_time_total
    log_file = os.path.join(run_output_path, "log.txt")

    with open(log_file, "w") as f:
        f.write("===== PDB2Net Batch Log =====\n\n")
        f.write("Timing (gesamt):\n")
        f.write(f"- Parsing:                 {sum_times['parsing']:.1f} sec\n")
        f.write(f"- Classification (SIFTS):  {sum_times['sifts']:.1f} sec\n")
        f.write(f"- Classification (BLAST):  {sum_times['blast']:.1f} sec\n")
        f.write(f"- Interaction:             {sum_times['interaction']:.1f} sec\n")
        f.write(f"- Network export:          {sum_times['networks']:.1f} sec\n")
        f.write(f"- Total:                   {total_time:.1f} sec\n")

        f.write("\nDetaillierte Parsing-Zeit (parallel):\n")
        f.write(f"- process_single_file():   {parsing_duration:.2f} sec\n")

        f.write("\n===============================\n")

    # Speicher freigeben
    tree_cache.clear()
    coords_cache.clear()
    gc.collect()


def create_batches_streaming(file_paths, max_batch_kb):
    """
    Generator, der Dateien nacheinander in Batches liefert,
    ohne vorher alle Batchgrößen zu berechnen.
    """
    current_batch = []
    current_size = 0

    for file_path in file_paths:
        try:
            size_kb = os.path.getsize(file_path) // 1024
        except Exception as e:
            print(f"⚠️ Fehler beim Lesen der Dateigröße: {file_path} – {e}")
            continue

        if size_kb > max_batch_kb:
            print(f"⚠️ Datei zu groß für ein einzelnes Batch: {file_path} ({size_kb} KB)")
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


def batch_run(input_folder, timeout_minutes=10, size_limit_kb=716800):  # 700 MB
    def stream_valid_files(folder):
        for entry in os.scandir(folder):
            if entry.is_file() and is_valid_file(entry.path):
                yield entry.path

    timeout_seconds = timeout_minutes * 60
    logs_dir = os.path.join(config["output_path"], "logs")
    os.makedirs(logs_dir, exist_ok=True)
    log_path = os.path.join(logs_dir, "skipped_batches.txt")

    start_time_all = time.time()
    total_done = 0

    file_stream = stream_valid_files(input_folder)

    for i, batch_files in enumerate(create_batches_streaming(file_stream, size_limit_kb), start=1):
        print(f"\n--- Bearbeite Batch {i} ({len(batch_files)} Dateien) ---")

        start_time_batch = time.time()
        process = multiprocessing.Process(target=run_main, args=(batch_files,))
        process.start()
        process.join(timeout=timeout_seconds)

        if process.is_alive():
            print(f"⚠️ Batch {i} zu lange (> {timeout_minutes} min). Wird abgebrochen.")
            process.terminate()
            process.join()
            with open(log_path, "a") as log_file:
                for file_path in batch_files:
                    pdb_id = get_pdb_id(file_path)
                    log_file.write(f"{pdb_id}\n")
        else:
            duration = time.time() - start_time_batch
            total_done += len(batch_files)
            print(f"✅ Batch {i} abgeschlossen in {duration:.1f} Sekunden.")
            print(f"📈 Verarbeitet bisher: {total_done} Dateien.")

    total_time_all = time.time() - start_time_all
    print(f"\n⏱ Gesamtzeit aller Batches: {total_time_all:.2f} Sekunden")



if __name__ == "__main__":
    from cytoscape import ensure_cytoscape_running
    ensure_cytoscape_running()
    INPUT_FOLDER_PATH = config["input_folder_path"]
    batch_run(INPUT_FOLDER_PATH)


