import os
import pandas as pd
from config_loader import config
from matplotlib import cm
from matplotlib.colors import to_hex
import json
import time
import py4cytoscape as p4c
import subprocess
from config_loader import config

def ensure_cytoscape_running():
    try:
        p4c.cytoscape_ping()
        print("Cytoscape is already running.")
    except:
        print("Starting Cytoscape...")
        subprocess.Popen(config["cytoscape_path"])
        time.sleep(40)
        try:
            p4c.cytoscape_ping()
            print("Cytoscape started successfully.")
        except:
            print("Error: Cytoscape could not be started. Check the path in config.json.")
            exit(1)


def export_custom_cyjs(nodes_df, edges_df, network_title, run_output_path):
    """
    Saves a Cytoscape-compatible .cyjs file without using Cytoscape.
    Includes molecule_name, tooltip, color_group etc. in the node attributes.
    """
    import json
    os.makedirs(run_output_path, exist_ok=True)
    network_file = os.path.join(run_output_path, f"{network_title}.cyjs")

    nodes = []
    for _, row in nodes_df.iterrows():
        node = {
            "data": {
                "id": row["id"],
                "name": row["name"],
                "molecule_name": row.get("molecule_name", "Unknown")
            }
        }
        if "color_group" in row:
            node["data"]["color_group"] = row["color_group"]
        if "tooltip" in row:
            node["data"]["tooltip"] = row["tooltip"]
        nodes.append(node)

    edges = []
    for idx, row in edges_df.iterrows():
        edge = {
            "data": {
                "id": f"edge_{idx}",
                "source": row["source"],
                "target": row["target"],
                "interaction": row.get("interaction", "interacts_with"),
                "all_atoms_count": row.get("all_atoms_count", 0)
            }
        }
        edges.append(edge)

    cyjs_data = {
        "data": {"name": network_title},
        "elements": {"nodes": nodes, "edges": edges}
    }

    with open(network_file, "w") as f:
        json.dump(cyjs_data, f, indent=2)


def create_cytoscape_network(results, network_title="Protein_Interaction_Network", run_output_path=".", nodes_data=None):
    import py4cytoscape as p4c
    import pandas as pd
    from matplotlib import cm
    from matplotlib.colors import to_hex

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

    if not edges:
        print("Warning: No valid edges found. Network will not be created.")
        return

    if nodes_data:
        nodes_df = pd.DataFrame(nodes_data)
        nodes_df["name"] = nodes_df["id"]
        nodes_df["tooltip"] = nodes_df["molecule_name"]
    else:
        nodes_df = pd.DataFrame({"id": list(unique_nodes)})
        nodes_df["name"] = nodes_df["id"]
        nodes_df["tooltip"] = "Unknown"

    edges_df = pd.DataFrame(edges)
    edges_df.rename(columns={"chain_a": "source", "chain_b": "target"}, inplace=True)
    edges_df["interaction"] = "interacts_with"

    if not config.get("open_in_cytoscape", True):
        export_custom_cyjs(nodes_df, edges_df, network_title, run_output_path)
        return

    try:
        p4c.cytoscape_ping()
    except Exception as e:
        print(f"Error: Cytoscape is not running. Details: {e}")
        return

    existing_networks = p4c.get_network_list()
    while len(existing_networks) > config["keep_last_n_networks"]:
        oldest_network = existing_networks.pop(0)
        p4c.delete_network(oldest_network)

    p4c.create_network_from_data_frames(nodes_df, edges_df, title=network_title)

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
        print(f"Error while loading node data: {e}")

    try:
        color_groups = sorted(nodes_df["color_group"].dropna().unique()) if "color_group" in nodes_df.columns else []
        is_chain_network = "Chain" in network_title or network_title == "Combined_Network"
        is_protein_network = "Protein" in network_title
        is_combined_protein = is_protein_network and "combined" in network_title.lower()

        style_name = (
            "PDB2Net_Chain_Style" if is_chain_network else
            "PDB2Net_Protein_Combined_Style" if is_combined_protein else
            "PDB2Net_Protein_Style"
        )

        # 🎨 Feste Farben für bestimmte Molekültypen
        fixed_colors = {
            "Protein": "#1f77b4",      # Blau
            "DNA": "#ff7f0e",          # Orange
            "RNA": "#2ca02c",          # Grün
            "DNA/RNA": "#d62728",      # Rot
            "Nucleic Acid": "#9467bd"  # Lila
        }

        # Automatische Farben für alle anderen Gruppen
        base_color_groups = [g for g in color_groups if g not in fixed_colors and g != "Multi"]
        cmap = cm.get_cmap('tab20', len(base_color_groups))
        auto_colors = {group: to_hex(cmap(i)) for i, group in enumerate(base_color_groups)}

        # Multi-Farbe im kombinierten Protein-Netzwerk
        if is_combined_protein and "Multi" in color_groups:
            auto_colors["Multi"] = "#8c564b"

        # Farb-Map kombinieren: feste Farben + automatisch generierte
        color_map = {**auto_colors, **fixed_colors}

        if style_name not in p4c.get_visual_style_names():
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

            p4c.create_visual_style(style_name, mappings=mappings, defaults=defaults)

        p4c.set_current_network(network_title)
        p4c.set_visual_style(style_name)
        p4c.layout_network(layout_name="force-directed")

    except Exception as e:
        print(f"Error while applying style: {e}")

    os.makedirs(run_output_path, exist_ok=True)
    network_file = os.path.join(run_output_path, f"{network_title}.cyjs")
    p4c.export_network(network_file, type="cyjs")




def generate_nodes_from_atom_data(atom_data, pdb_id=None):
    """
    Generates a list of nodes from atom data for Cytoscape.

    Args:
        atom_data (list): List of chains with atom/residue info.
        pdb_id (str): Optional PDB ID for node labeling.

    Returns:
        list: List of node dictionaries for Cytoscape export.
    """
    return [
        {
            "id": chain["unique_chain_id"],
            "color_group": chain.get("molecule_type", "Unknown"),
            "molecule_name": chain.get("molecule_name", "Unknown")
        }
        for chain in atom_data
    ]
