import py4cytoscape as p4c
import time
import pandas as pd

def create_cytoscape_network(results):
    """
    Creates a network in Cytoscape based on interaction data.

    Args:
        results (list): A list of dictionaries containing interaction data, each with:
            - "file_path": The path to the PDB file.
            - "chain_a": Unique identifier for the first chain.
            - "chain_b": Unique identifier for the second chain.
            - "ca_nn_count": Number of Cα/Cβ contacts within 15 Å.
            - "all_atoms_close_count": Number of all-atom contacts within 5 Å.
    """

    # 🔹 Step 1: Check if Cytoscape is running
    try:
        status = p4c.cytoscape_ping()
        print("🌐 Successfully connected to Cytoscape!")
    except Exception as e:
        print(f"❌ Cytoscape is not running! Please start Cytoscape and try again. Error: {e}")
        return

    # 🔹 Step 2: Delete existing networks (if any)
    try:
        existing_networks = p4c.get_network_list()
        if existing_networks:
            print(f"🗑️ Deleting {len(existing_networks)} existing networks...")
            for net in existing_networks:
                p4c.delete_network(net)
            time.sleep(1)  # Allow Cytoscape time to process
    except Exception as e:
        print(f"⚠️ Warning: Could not delete existing networks: {e}")

    # 🔹 Step 3: Extract unique nodes and edges
    unique_nodes = set()
    edges = []

    for entry in results:
        chain_a = entry["chain_a"]
        chain_b = entry["chain_b"]
        ca_nn = entry["ca_nn_count"]
        all_atoms = entry["all_atoms_close_count"]

        unique_nodes.add(chain_a)
        unique_nodes.add(chain_b)

        edges.append({
            "source": chain_a,
            "target": chain_b,
            "ca_nn_count": ca_nn,
            "all_atoms_close_count": all_atoms
        })

    # 🔹 Step 4: Convert nodes and edges into DataFrames
    nodes_df = pd.DataFrame({"id": list(unique_nodes), "label": list(unique_nodes)})
    edges_df = pd.DataFrame(edges)

    # 🔹 Step 5: Create the network in Cytoscape
    print("📡 Creating a new network in Cytoscape...")
    network_suid = p4c.create_network_from_data_frames(nodes_df, edges_df, title="Protein Interaction Network")

    # 🔹 Step 6: Apply a layout and display the network
    print("🎨 Applying layout...")
    time.sleep(1)  # Ensure the network is fully loaded before applying layout
    p4c.layout_network(layout_name="circular")

    print("✅ Network successfully created in Cytoscape!")
