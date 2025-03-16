import time
from file_parser import read_files_from_csv
from data_processor import process_structure
from unknown_molecule_uniprot import process_molecule_info
from uniprot_matcher import match_sequence_to_uniprot
from distances import calculate_distances_with_ckdtree
from cytoscape import create_cytoscape_network
from protein_network import create_protein_network
from config_loader import config
from datetime import datetime
import os
import py4cytoscape as p4c
import subprocess
from detailed_results_exporter import export_detailed_interactions

# Load Cytoscape path from config
CYTOSCAPE_PATH = config["cytoscape_path"]

# Check if Cytoscape is running; if not, try to start it
try:
    p4c.cytoscape_ping()
    print("Cytoscape is already running.")
except:
    print("Starting Cytoscape...")
    subprocess.Popen(CYTOSCAPE_PATH)
    time.sleep(30)
    try:
        p4c.cytoscape_ping()
        print("Cytoscape started successfully.")
    except:
        print("Error: Cytoscape could not be started. Check the path in config.json.")
        exit(1)

def main(csv_path):
    """
    Main function to process PDB structures, compute distances, and visualize networks.

    Args:
        csv_path (str): Path to the CSV file containing file paths.
    """
    network_config = config["networks"]

    # Create a unique output directory for the run
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    run_output_path = os.path.join(config["output_path"], timestamp)
    os.makedirs(run_output_path, exist_ok=True)
    print(f"\nCreated run output directory: {run_output_path}")

    # Measure total execution time
    start_time_total = time.time()

    print("\nLoading PDB files from CSV...")
    start_time = time.time()
    structures = read_files_from_csv(csv_path)
    combined_data = [process_structure(structure_data) for structure_data in structures]
    print(f"File parsing completed in {time.time() - start_time:.2f} seconds.")

    print("\nDetermining molecule names and types...")
    start_time = time.time()
    process_molecule_info(combined_data)
    match_sequence_to_uniprot(combined_data)
    print(f"Molecule identification completed in {time.time() - start_time:.2f} seconds.")

    print("\nComputing atomic distances...")
    start_time = time.time()
    results = calculate_distances_with_ckdtree(combined_data)
    print(f"Distance calculation completed in {time.time() - start_time:.2f} seconds.")

    if config.get("export_detailed_interactions", False):
        print("\nExporting detailed interaction data for each PDB file...")
        start_time = time.time()
        for structure_data in combined_data:
            pdb_id = structure_data["pdb_id"]
            pdb_interactions = [res for res in results if res["chain_a"].startswith(pdb_id)]
            export_detailed_interactions(structure_data, pdb_interactions, run_output_path)
        print(f"Export completed in {time.time() - start_time:.2f} seconds.")

    if network_config["chain_per_pdb"]:
        print("\nCreating separate networks for each PDB file...")
        start_time = time.time()
        results_by_pdb = {}
        for entry in results:
            pdb_id = entry["chain_a"].split(":")[0]
            results_by_pdb.setdefault(pdb_id, []).append(entry)
        for pdb_id, pdb_results in results_by_pdb.items():
            create_cytoscape_network(pdb_results, network_title=f"Chain_Interaction_Network_{pdb_id}", run_output_path=run_output_path)
        print(f"Network creation completed in {time.time() - start_time:.2f} seconds.")

    if network_config["combined_chain_network"]:
        print("\nCreating a single combined chain network...")
        start_time = time.time()
        create_cytoscape_network(results, network_title="Combined_Network", run_output_path=run_output_path)
        print(f"Combined network created in {time.time() - start_time:.2f} seconds.")

    if network_config["protein_per_pdb"] or network_config["combined_protein_network"]:
        print("\nProcessing protein-level interactions...")
        start_time = time.time()
        create_protein_network(results, combined_data, run_output_path, network_config)
        print(f"Protein network creation completed in {time.time() - start_time:.2f} seconds.")

    # Display total execution time
    print(f"\nTotal execution time: {time.time() - start_time_total:.2f} seconds.")

if __name__ == "__main__":
    INPUT_CSV_PATH = config["input_csv_path"]
    main(INPUT_CSV_PATH)
