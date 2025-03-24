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
        nodes_data (list of dict, optional): Node metadata including color_group.
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
    else:
        nodes_df = pd.DataFrame({"id": list(unique_nodes)})
        nodes_df["name"] = nodes_df["id"]

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

    # Default style fallback
    if "BioPAX_SIF" in p4c.get_visual_style_names():
        p4c.set_visual_style("BioPAX_SIF")

    # Apply color_group mapping and label passthrough
    if "color_group" in nodes_df.columns:
        try:
            print("\n🔍 Vorschau auf Node-Daten mit color_group:")
            print(nodes_df[['id', 'color_group']].head())

            # color_group in Cytoscape laden
            p4c.load_table_data(
                data=nodes_df[['id', 'color_group']],
                data_key_column="id",
                table="node",
                table_key_column="name"
            )

            # UniProt-ID als "name" setzen
            try:
                nodes_df_corrected = nodes_df.copy()
                nodes_df_corrected["name"] = nodes_df_corrected["id"]
                p4c.load_table_data(
                    data=nodes_df_corrected[['id', 'name']],
                    data_key_column="id",
                    table="node",
                    table_key_column="id"
                )
                print("✅ Cytoscape-Spalte 'name' wurde erfolgreich mit UniProt-IDs überschrieben.")
            except Exception as e:
                print(f"❌ Fehler beim Überschreiben der 'name'-Spalte: {e}")

            # Farb-Mapping vorbereiten
            color_groups = sorted(nodes_df["color_group"].dropna().unique())
            cmap = cm.get_cmap('tab20', len(color_groups))
            color_map = {group: to_hex(cmap(i)) for i, group in enumerate(color_groups)}

            style_name = "PDB2Net_Style"

            if "combined" in network_title.lower():
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
                            "mappingType": "discrete",
                            "mappingColumn": "color_group",
                            "mappingColumnType": "String",
                            "visualProperty": "NODE_FILL_COLOR",
                            "map": [{"key": k, "value": v} for k, v in color_map.items()]
                        }
                    ]
                    p4c.create_visual_style(style_name, mappings=mappings, defaults=defaults)
                    print(f"🎨 Visual Style '{style_name}' mit NODE_LABEL passthrough erstellt.")
                else:
                    print(f"🎨 Visual Style '{style_name}' bereits vorhanden.")

                # Aktuelles Netzwerk setzen und Style anwenden
                try:
                    p4c.set_current_network(network_title)
                    p4c.set_visual_style(style_name)
                    print("✅ Visual Style wurde aktiv auf das kombinierte Netzwerk angewendet.")

                    # 🔥 **Hier wird das Mapping sicher nochmal erzwungen**
                    p4c.map_visual_property("NODE_LABEL", "name", "p")
                    print("✅ Passthrough-Mapping für NODE_LABEL wurde erzwungen.")

                except Exception as e:
                    print(f"❌ Fehler beim Anwenden des Styles: {e}")

                # Finales Layout erneut anwenden
                p4c.layout_network(layout_name="force-directed")

                # Debug
                try:
                    node_table = p4c.get_table_columns('node', columns=['name'])
                    print("🔍 Finale Cytoscape Node-Labels (name):", node_table.head())
                except Exception as e:
                    print(f"❌ Fehler beim Abrufen der Node-Tabelle: {e}")
            else:
                print("🔄 Standard-Style für separates Netzwerk beibehalten.")

        except Exception as e:
            print(f"❌ Fehler beim Anwenden des finalen Label-Fixes: {e}")

    # Export network
    pdb_output_path = os.path.join(run_output_path, network_title)
    os.makedirs(pdb_output_path, exist_ok=True)
    network_file = os.path.join(pdb_output_path, f"{network_title}.cyjs")
    p4c.export_network(network_file, type="cyjs")
