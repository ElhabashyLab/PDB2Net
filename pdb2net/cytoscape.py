import py4cytoscape as p4c
import os
import pandas as pd
from config_loader import config


def create_cytoscape_network(results, network_title="Protein_Interaction_Network", run_output_path="."):
    """
    Creates a network visualization in Cytoscape based on interaction data.

    Args:
        results (list): A list of dictionaries containing interaction data.
        network_title (str): The name of the network in Cytoscape.
        run_output_path (str): The directory where the network output should be saved.
    """

    # Check if Cytoscape is running
    try:
        p4c.cytoscape_ping()
    except Exception as e:
        print(f"Error: Cytoscape is not running. Please start Cytoscape and try again. Details: {e}")
        return

    # Remove older networks if configured
    existing_networks = p4c.get_network_list()
    while len(existing_networks) > config["keep_last_n_networks"]:
        oldest_network = existing_networks.pop(0)
        p4c.delete_network(oldest_network)

    # Prepare nodes and edges
    unique_nodes = set()
    edges = []

    for entry in results:
        if entry.get("all_atoms_count", 0) > 0:  # Only interactions where all_atoms_count > 0
            chain_a = entry["chain_a"]
            chain_b = entry["chain_b"]
            unique_nodes.add(chain_a)
            unique_nodes.add(chain_b)

            edges.append({
                "chain_a": chain_a,
                "chain_b": chain_b,
                "all_atoms_count": entry["all_atoms_count"]
            })

    # Create DataFrames for Cytoscape
    nodes_df = pd.DataFrame({"id": list(unique_nodes), "label": list(unique_nodes)})
    edges_df = pd.DataFrame(edges)

    if edges_df.empty:
        print("Warning: No valid edges found. Network will not be created.")
        return

    edges_df.rename(columns={"chain_a": "source", "chain_b": "target"}, inplace=True)
    edges_df["interaction"] = "interacts_with"

    # Create network in Cytoscape
    network_suid = p4c.create_network_from_data_frames(nodes_df, edges_df, title=network_title)

    # Apply layout
    p4c.layout_network(layout_name="force-directed")

    # Apply visualization style
    available_styles = p4c.get_visual_style_names()

    if "BioPAX_SIF" in available_styles:
        p4c.set_visual_style("BioPAX_SIF")
    else:
        print("Warning: BioPAX_SIF style not available. Please check Cytoscape settings.")

    # Save network to output directory
    pdb_output_path = os.path.join(run_output_path, network_title)
    os.makedirs(pdb_output_path, exist_ok=True)
    network_file = os.path.join(pdb_output_path, f"{network_title}.cyjs")
    p4c.export_network(network_file, type="cyjs")
