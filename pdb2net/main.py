from file_parser import read_files_from_csv
from data_processor import process_structure
from unknown_molecule_uniprot import process_molecule_info
from distances import calculate_distances_with_ckdtree
from cytoscape import create_cytoscape_network
from config_loader import config
from datetime import datetime
import os
import py4cytoscape as p4c
import subprocess
import time

CYTOSCAPE_PATH = config["cytoscape_path"]  # Aus Config laden

# Prüfen, ob Cytoscape läuft
try:
    p4c.cytoscape_ping()
    print("🌐 Cytoscape läuft bereits!")
except:
    print("⚙️ Cytoscape wird gestartet...")
    subprocess.Popen(CYTOSCAPE_PATH)
    time.sleep(30)  # Wartezeit, bis Cytoscape vollständig gestartet ist
    # Prüfen, ob Cytoscape nun läuft
    try:
        p4c.cytoscape_ping()
        print("✅ Cytoscape erfolgreich gestartet!")
    except:
        print("❌ Cytoscape konnte nicht gestartet werden. Prüfe den Pfad in config.json!")
        exit(1)

def main(csv_path):
    """
    Main function: Loads PDB structures, computes distances, and visualizes networks.
    """

    # 🔹 Create a unique run output folder with timestamp
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    run_output_path = os.path.join(config["output_path"], timestamp)
    os.makedirs(run_output_path, exist_ok=True)
    print(f"\n📁 Created run output directory: {run_output_path}")

    print("\n📂 Loading PDB files from CSV...")
    structures = read_files_from_csv(csv_path)
    combined_data = [process_structure(structure_data) for structure_data in structures]

    print("\n🔍 Determining molecule names and types...")
    process_molecule_info(combined_data)

    print("\n📏 Computing atomic distances...")
    results = calculate_distances_with_ckdtree(combined_data)

    # Decide if creating separate networks or one combined network
    if config["create_separate_networks"]:
        print("\n🌐 Creating separate networks for each PDB file...")
        results_by_pdb = {}
        for entry in results:
            pdb_id = entry["chain_a"].split(":")[0]  # Extract PDB-ID
            if pdb_id not in results_by_pdb:
                results_by_pdb[pdb_id] = []
            results_by_pdb[pdb_id].append(entry)

        for pdb_id, pdb_results in results_by_pdb.items():
            create_cytoscape_network(pdb_results, network_title=pdb_id, run_output_path=run_output_path)

    else:
        print("\n🌐 Creating a single combined network...")
        create_cytoscape_network(results, network_title="Combined_Network", run_output_path=run_output_path)

if __name__ == "__main__":
    csv_path = "C:\\Users\\Gregor\\Documents\\Uni Bioinformatik\\9. Semester\\B.A\\PDBFiles\\PathsCSV.csv"
    main(csv_path)
