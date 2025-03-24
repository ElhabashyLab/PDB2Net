import py4cytoscape as p4c
import os
import pandas as pd
from config_loader import config
from matplotlib import cm
from matplotlib.colors import to_hex


def create_cytoscape_network(results, network_title="Protein_Interaction_Network", run_output_path=".",
                             nodes_data=None):
    """
    Creates a network visualization in Cytoscape based on interaction data.

    Args:
        results (list): A list of dictionaries containing interaction data.
        network_title (str): The name of the network in Cytoscape.
        run_output_path (str): The directory where the network output should be saved.
        nodes_data (list of dict, optional): Node metadata including color_group, molecule_name, etc.
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
        if entry.get("all_atoms_count", 0) > 0:
            chain_a = entry["chain_a"]
            chain_b = entry["chain_b"]
            unique_nodes.add(chain_a)
            unique_nodes.add(chain_b)
            edges.append({
                "chain_a": chain_a,
                "chain_b": chain_b,
                "all_atoms_count": entry["all_atoms_count"]
            })

    # Create nodes dataframe
    if nodes_data:
        nodes_df = pd.DataFrame(nodes_data)
        nodes_df["name"] = nodes_df["id"]

        if "molecule_name" in nodes_df.columns:
            nodes_df["tooltip"] = nodes_df["molecule_name"]
        elif "label" in nodes_df.columns:
            nodes_df["tooltip"] = nodes_df["label"]
        else:
            nodes_df["tooltip"] = "Unknown"
    else:
        nodes_df = pd.DataFrame({"id": list(unique_nodes)})
        nodes_df["name"] = nodes_df["id"]
        nodes_df["tooltip"] = "Unknown"

    edges_df = pd.DataFrame(edges)
    if edges_df.empty:
        print("Warning: No valid edges found. Network will not be created.")
        return

    edges_df.rename(columns={"chain_a": "source", "chain_b": "target"}, inplace=True)
    edges_df["interaction"] = "interacts_with"

    # Create network in Cytoscape
    network_suid = p4c.create_network_from_data_frames(nodes_df, edges_df, title=network_title)
    p4c.layout_network(layout_name="force-directed")

    # Default style fallback
    if "BioPAX_SIF" in p4c.get_visual_style_names():
        p4c.set_visual_style("BioPAX_SIF")

    # Load node table data
    try:
        if "color_group" in nodes_df.columns:
            p4c.load_table_data(
                data=nodes_df[["id", "color_group"]],
                data_key_column="id",
                table="node",
                table_key_column="name"
            )

        p4c.load_table_data(
            data=nodes_df[["id", "name"]],
            data_key_column="id",
            table="node",
            table_key_column="id"
        )

        if "tooltip" in nodes_df.columns:
            p4c.load_table_data(
                data=nodes_df[["id", "tooltip"]],
                data_key_column="id",
                table="node",
                table_key_column="name"
            )

    except Exception as e:
        print(f"❌ Fehler beim Laden von Knotendaten: {e}")

    # Apply visual style
    try:
        color_groups = sorted(nodes_df["color_group"].dropna().unique())
        base_color_groups = [g for g in color_groups if g != "Multi"]
        print(f"[DEBUG] Found {len(color_groups)} unique color groups: {color_groups}")
        cmap = cm.get_cmap('tab20', len(base_color_groups))
        color_map = {group: to_hex(cmap(i)) for i, group in enumerate(base_color_groups)}

        # 🔴 Manuelle Farbe für "Multi" (rot)
        if "Multi" in color_groups:
            color_map["Multi"] = "#d62728"  # schönes kräftiges Rot aus tab10

        # ✅ Individueller Style pro Netzwerk
        style_name = f"PDB2Net_Style_{network_title}"

        defaults = {
            "NODE_SHAPE": "ELLIPSE",
            "NODE_SIZE": 40,
            "NODE_LABEL_POSITION": "C,C,c,0.00,0.00",
            "EDGE_TRANSPARENCY": 120
        }

        mappings = [
            {
                "mappingType": "passthrough",
                "mappingColumn": "name",
                "mappingColumnType": "String",
                "visualProperty": "NODE_LABEL"
            },
            {
                "mappingType": "passthrough",
                "mappingColumn": "tooltip",
                "mappingColumnType": "String",
                "visualProperty": "NODE_TOOLTIP"
            }
        ]

        if "color_group" in nodes_df.columns:
            mappings.append({
                "mappingType": "discrete",
                "mappingColumn": "color_group",
                "mappingColumnType": "String",
                "visualProperty": "NODE_FILL_COLOR",
                "map": [{"key": k, "value": v} for k, v in color_map.items()]
            })

        # Style immer neu erstellen – für volle Kontrolle über Farblogik
        p4c.create_visual_style(style_name, mappings=mappings, defaults=defaults)
        print(f"🎨 Visual Style '{style_name}' mit Tooltip und Farben erstellt.")

        p4c.set_current_network(network_title)
        p4c.set_visual_style(style_name)
        p4c.map_visual_property("NODE_LABEL", "name", "p")
        print("✅ Visual Style wurde angewendet.")

    except Exception as e:
        print(f"❌ Fehler beim Anwenden des Styles: {e}")

    # Export network
    pdb_output_path = os.path.join(run_output_path, network_title)
    os.makedirs(pdb_output_path, exist_ok=True)
    network_file = os.path.join(pdb_output_path, f"{network_title}.cyjs")
    p4c.export_network(network_file, type="cyjs")


def generate_nodes_from_atom_data(atom_data, pdb_id=None):
    """
    Erstellt eine Liste von nodes_data-Einträgen für create_cytoscape_network().

    Args:
        atom_data (list): Liste von Chains mit atomaren Infos.
        pdb_id (str, optional): Wird hier NICHT mehr für color_group verwendet.

    Returns:
        list of dict: Node-Einträge mit id, color_group (molecule_type), molecule_name.
    """
    return [
        {
            "id": chain["unique_chain_id"],
            "color_group": chain.get("molecule_type", "Unknown"),
            "molecule_name": chain.get("molecule_name", "Unknown")
        }
        for chain in atom_data
    ]
