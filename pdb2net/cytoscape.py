import py4cytoscape as p4c
import time
import os
import pandas as pd
from config_loader import config

def create_cytoscape_network(results, network_title="Protein_Interaction_Network", run_output_path="."):
    """
    Creates a network in Cytoscape based on interaction data.

    Args:
        results (list): List of dictionaries containing interaction data.
        network_title (str): Name of the network (usually PDB-ID or "Combined_Network").
        run_output_path (str): The output directory for the current run.
    """

    # 🔹 Check if Cytoscape is running
    try:
        p4c.cytoscape_ping()
        print(f"🌐 Connected to Cytoscape! Creating network: {network_title}")
    except Exception as e:
        print(f"❌ Cytoscape is not running! Start Cytoscape and try again. Error: {e}")
        return

    # 🔹 Delete excess networks (keep only the last N)
    existing_networks = p4c.get_network_list()
    while len(existing_networks) > config["keep_last_n_networks"]:
        oldest_network = existing_networks.pop(0)  # Remove oldest network
        p4c.delete_network(oldest_network)

    # 🔹 Prepare nodes & edges DataFrame
    unique_nodes = set()
    edges = []

    for entry in results:
        chain_a = entry["chain_a"]
        chain_b = entry["chain_b"]
        unique_nodes.add(chain_a)
        unique_nodes.add(chain_b)

        edges.append({
            "chain_a": chain_a,
            "chain_b": chain_b,
            "ca_nn_count": entry["ca_nn_count"],
            "all_atoms_close_count": entry["all_atoms_close_count"]
        })

    nodes_df = pd.DataFrame({"id": list(unique_nodes), "label": list(unique_nodes)})
    edges_df = pd.DataFrame(edges)

    # 🔹 Wichtig: Prüfe, welche Spalten wirklich existieren!
    print(edges_df.head())  # Debugging: zeigt die aktuellen Spaltennamen

    # 🔹 Korrekte Umbenennung in source und target für Cytoscape
    edges_df.rename(columns={"chain_a": "source", "chain_b": "target"}, inplace=True)
    edges_df["interaction"] = "interacts_with"

    # 🔹 Prüfe die endgültigen Spalten
    print(edges_df.head())  # Debugging: Prüfe ob Spalten korrekt sind

    # 🔹 Create the network in Cytoscape
    network_suid = p4c.create_network_from_data_frames(nodes_df, edges_df, title=network_title)

    # 🔹 Apply layout
    p4c.layout_network(layout_name="circular")

    # 🔹 Save the network to a specific PDB subfolder inside the run output path
    pdb_output_path = os.path.join(run_output_path, network_title)
    os.makedirs(pdb_output_path, exist_ok=True)
    network_file = os.path.join(pdb_output_path, f"{network_title}.cyjs")
    p4c.export_network(network_file, type="cyjs")
    print(f"✅ Network saved to: {network_file}")
