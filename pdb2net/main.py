import time
from file_parser import read_files_from_folder, is_valid_file, get_pdb_id, parse_structure
from data_processor import process_structure
from unknown_molecule_uniprot import process_molecule_info
from distances import calculate_distances_with_ckdtree
from cytoscape import create_cytoscape_network, generate_nodes_from_atom_data
from protein_network import create_protein_network
from config_loader import config
from datetime import datetime
import os
import py4cytoscape as p4c
import subprocess
from detailed_results_exporter import export_detailed_interactions
from uniprot_matcher import parallel_blast_search

# Cytoscape-Startversuch (nur einmal am Anfang)
CYTOSCAPE_PATH = config["cytoscape_path"]
try:
    p4c.cytoscape_ping()
    print("Cytoscape is already running.")
except:
    print("Starting Cytoscape...")
    subprocess.Popen(CYTOSCAPE_PATH)
    time.sleep(40)
    try:
        p4c.cytoscape_ping()
        print("Cytoscape started successfully.")
    except:
        print("Error: Cytoscape could not be started. Check the path in config.json.")
        exit(1)

def main(input_path_or_filelist):
    """
    Hauptfunktion für PDB-Verarbeitung. Unterstützt entweder einen Ordnerpfad oder eine Liste von Dateien.
    """
    network_config = config["networks"]

    # 🔹 INPUT-MODUS: Liste oder Ordner
    if isinstance(input_path_or_filelist, list):
        file_paths = [f for f in input_path_or_filelist if is_valid_file(f)]
    else:
        file_paths = [os.path.join(input_path_or_filelist, f) for f in os.listdir(input_path_or_filelist)
                      if os.path.isfile(os.path.join(input_path_or_filelist, f)) and is_valid_file(f)]

    # 🔹 Einmaliger Output-Ordner pro Batch
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    run_output_path = os.path.join(config["output_path"], timestamp)
    os.makedirs(run_output_path, exist_ok=True)
    #print(f"\nCreated run output directory: {run_output_path}")

    start_time_total = time.time()
    sum_times = {"parsing": 0, "sifts": 0, "blast": 0, "interaction": 0, "networks": 0}

    # 🔹 Strukturen einlesen
    #print("\nReading structure files...")
    start_time = time.time()
    structures = []
    for f in file_paths:
        pdb_id = get_pdb_id(f)
        structure = parse_structure(f, pdb_id)
        if structure:
            structures.append({
                "file_path": f,
                "pdb_id": pdb_id,
                "structure": structure
            })
    combined_data = [process_structure(s) for s in structures]
    sum_times["parsing"] = time.time() - start_time
    #print(f"File parsing completed in {sum_times['parsing']:.2f} seconds.")

    # 🔹 Molekülinformationen
    #print("\nDetermining molecule names and types...")

    # ⏱ SIFTS
    start_time = time.time()
    process_molecule_info(combined_data)
    sum_times["sifts"] = time.time() - start_time
    #print(f"SIFTS mapping completed in {sum_times['sifts']:.2f} seconds.")

    # ⏱ BLAST
    start_time = time.time()
    parallel_blast_search(combined_data, max_workers=4)
    sum_times["blast"] = time.time() - start_time
    #print(f"BLAST search completed in {sum_times['blast']:.2f} seconds.")

    # 🔹 Interaktionen berechnen
    #print("\nComputing atomic distances...")
    start_time = time.time()
    results = calculate_distances_with_ckdtree(combined_data)
    sum_times["interaction"] = time.time() - start_time
    #print(f"Distance calculation completed in {sum_times['interaction']:.2f} seconds.")

    # 🔹 Detail-Export
    if config.get("export_detailed_interactions", False):
        #print("\nExporting detailed interaction data...")
        for structure_data in combined_data:
            pdb_id = structure_data["pdb_id"]
            pdb_interactions = [res for res in results if res["chain_a"].startswith(pdb_id)]
            export_detailed_interactions(structure_data, pdb_interactions, run_output_path)
        #print("Detailed export completed.")

    # 🔹 Chain-Netzwerke pro PDB
    if network_config["chain_per_pdb"]:
        #print("\nCreating separate chain networks per PDB...")
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
        #print(f"Chain networks completed in {sum_times['networks']:.2f} seconds.")

    # 🔹 Kombiniertes Chain-Netzwerk
    if network_config["combined_chain_network"]:
        #print("\nCreating combined chain network...")
        start_time = time.time()
        all_chains = [chain for structure in combined_data for chain in structure["atom_data"]]
        nodes_data = generate_nodes_from_atom_data(all_chains)
        create_cytoscape_network(results, "Combined_Network", run_output_path, nodes_data=nodes_data)
        sum_times["networks"] += time.time() - start_time
        #print(f"Combined chain network completed.")

    # 🔹 Protein-Netzwerke
    if network_config["protein_per_pdb"] or network_config["combined_protein_network"]:
        #print("\nCreating protein-level networks...")
        start_time = time.time()
        create_protein_network(results, combined_data, run_output_path, network_config)
        protein_time = time.time() - start_time
        sum_times["networks"] += protein_time
        #print(f"Protein networks completed in {protein_time:.2f} seconds.")

    total_time = time.time() - start_time_total
    #print(f"\n✅ Total execution time: {total_time:.2f} seconds.")

    # 🔹 Zeit-Log
    log_file = os.path.join(run_output_path, "log.txt")
    with open(log_file, "w") as f:
        f.write("===== PDB2Net Batch Log =====\n\n")
        f.write("Timing:\n")
        f.write(f"- Parsing:                 {sum_times['parsing']:.1f} sec\n")
        f.write(f"- Classification (SIFTS):  {sum_times['sifts']:.1f} sec\n")
        f.write(f"- Classification (BLAST):  {sum_times['blast']:.1f} sec\n")
        f.write(f"- Interaction:             {sum_times['interaction']:.1f} sec\n")
        f.write(f"- Network export:          {sum_times['networks']:.1f} sec\n")
        f.write(f"- Total:                   {total_time:.1f} sec\n")
        f.write("\n===============================\n")

def batch_run(input_folder, batch_size=500):
    """
    Läuft main() in Batches über alle Dateien im Input-Ordner, ohne Kopieren.
    """
    all_files = [os.path.join(input_folder, f) for f in os.listdir(input_folder)
                 if os.path.isfile(os.path.join(input_folder, f)) and is_valid_file(f)]

    total_batches = (len(all_files) + batch_size - 1) // batch_size
    print(f"\nStarting batch run: {total_batches} batches of size {batch_size}")

    start_time_all = time.time()  # ⏱ Gesamtlaufzeit starten

    for i in range(total_batches):
        batch_files = all_files[i * batch_size:(i + 1) * batch_size]
        print(f"\n--- Processing batch {i + 1}/{total_batches} ({len(batch_files)} files) ---")
        main(batch_files)

    total_time_all = time.time() - start_time_all
    print(f"\n⏱ Gesamtzeit aller Batches: {total_time_all:.2f} Sekunden")

if __name__ == "__main__":
    INPUT_FOLDER_PATH = config["input_folder_path"]
    batch_run(INPUT_FOLDER_PATH, batch_size=20)