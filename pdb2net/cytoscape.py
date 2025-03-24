import py4cytoscape as p4c
import os
import pandas as pd
from config_loader import config
from matplotlib import cm
from matplotlib.colors import to_hex


def create_cytoscape_network(results, network_title="Protein_Interaction_Network", run_output_path=".",
                             nodes_data=None):
    try:
        p4c.cytoscape_ping()
    except Exception as e:
        print(f"Error: Cytoscape is not running. Details: {e}")
        return

    # Remove older networks if the configured maximum is exceeded
    existing_networks = p4c.get_network_list()
    while len(existing_networks) > config["keep_last_n_networks"]:
        oldest_network = existing_networks.pop(0)
        p4c.delete_network(oldest_network)

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

        # Determine the type of network from the title
        is_chain_network = "Chain" in network_title
        is_protein_network = "Protein" in network_title
        is_combined_protein = is_protein_network and "combined" in network_title.lower()

        if is_chain_network:
            style_name = "PDB2Net_Chain_Style"
        elif is_combined_protein:
            style_name = "PDB2Net_Protein_Combined_Style"
        else:
            style_name = "PDB2Net_Protein_Style"

        base_color_groups = [g for g in color_groups if g != "Multi"]
        cmap = cm.get_cmap('tab20', len(base_color_groups))
        color_map = {group: to_hex(cmap(i)) for i, group in enumerate(base_color_groups)}

        # Set "Multi" to red only for combined protein network
        if is_combined_protein and "Multi" in color_groups:
            color_map["Multi"] = "#d62728"

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
        p4c.map_visual_property("NODE_LABEL", "name", "p")
        p4c.layout_network(layout_name="force-directed")

    except Exception as e:
        print(f"Error while applying style: {e}")

    pdb_output_path = os.path.join(run_output_path, network_title)
    os.makedirs(pdb_output_path, exist_ok=True)
    network_file = os.path.join(pdb_output_path, f"{network_title}.cyjs")
    p4c.export_network(network_file, type="cyjs")


def generate_nodes_from_atom_data(atom_data, pdb_id=None):
    """
    Generates node metadata for Cytoscape from atomic chain data.

    Args:
        atom_data (list): List of chain dictionaries with IDs and annotation.

    Returns:
        list: List of Cytoscape-compatible node dictionaries.
    """
    return [
        {
            "id": chain["unique_chain_id"],
            "color_group": chain.get("molecule_type", "Unknown"),
            "molecule_name": chain.get("molecule_name", "Unknown")
        }
        for chain in atom_data
    ]
