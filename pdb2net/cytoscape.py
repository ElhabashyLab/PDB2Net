import py4cytoscape as p4c
import time
import os
import pandas as pd
import random
from config_loader import config

# 🔹 PDB Standard-Kettenfarben (offiziell aus PyMOL / Chimera)
PDB_CHAIN_COLORS = [
    "#FF0000", "#0000FF", "#00FF00", "#FFFF00", "#FF00FF", "#00FFFF",
    "#F08080", "#800080", "#808000", "#008080", "#B22222", "#FF8C00",
    "#4682B4", "#32CD32", "#9400D3", "#FFD700", "#DC143C", "#008000",
    "#1E90FF", "#FF4500"
]

def get_random_color():
    """Generiert eine zufällige, aber kontrastreiche Farbe für zusätzliche Ketten."""
    return "#{:06x}".format(random.randint(0, 0xFFFFFF))

def create_cytoscape_network(results, network_title="Protein_Interaction_Network", run_output_path="."):
    """
    Erstellt ein Netzwerk in Cytoscape basierend auf Interaktionsdaten.

    Args:
        results (list): Liste von Dictionaries mit Interaktionsdaten.
        network_title (str): Name des Netzwerks.
        run_output_path (str): Das Ausgabeverzeichnis für den aktuellen Lauf.
    """

    # 🔹 Prüfe, ob Cytoscape läuft
    try:
        p4c.cytoscape_ping()
        print(f"🌐 Verbunden mit Cytoscape! Erstelle Netzwerk: {network_title}")
    except Exception as e:
        print(f"❌ Cytoscape läuft nicht! Starte Cytoscape und versuche es erneut. Fehler: {e}")
        return

    # 🔹 Entferne ältere Netzwerke (falls konfiguriert)
    existing_networks = p4c.get_network_list()
    while len(existing_networks) > config["keep_last_n_networks"]:
        oldest_network = existing_networks.pop(0)
        p4c.delete_network(oldest_network)

    # 🔹 Knoten und Kanten vorbereiten
    unique_nodes = set()
    edges = []

    for entry in results:
        if entry.get("all_atoms_count", 0) > 0:  # 🔥 Nur Interaktionen mit `all_atoms_count > 0`
            chain_a = entry["chain_a"]
            chain_b = entry["chain_b"]
            unique_nodes.add(chain_a)
            unique_nodes.add(chain_b)

            edges.append({
                "chain_a": chain_a,
                "chain_b": chain_b,
                "all_atoms_count": entry["all_atoms_count"]
            })

    # 🔹 Erstelle DataFrames für Cytoscape
    nodes_df = pd.DataFrame({"id": list(unique_nodes), "label": list(unique_nodes)})
    edges_df = pd.DataFrame(edges)

    if edges_df.empty:
        print("⚠️ Keine gültigen Kanten gefunden. Das Netzwerk wird nicht erstellt.")
        return

    edges_df.rename(columns={"chain_a": "source", "chain_b": "target"}, inplace=True)
    edges_df["interaction"] = "interacts_with"

    # 🔹 Netzwerk in Cytoscape erstellen
    network_suid = p4c.create_network_from_data_frames(nodes_df, edges_df, title=network_title)

    # 🔹 Layout anwenden
    p4c.layout_network(layout_name="force-directed")

    # 🔹 Setze runde Knoten
    p4c.set_node_shape_default("ELLIPSE")

    # 🔹 Setze Knoten-Transparenz (Alpha-Wert: 0 = unsichtbar, 255 = voll sichtbar)
    p4c.set_node_fill_opacity_default(150)  # 150 gibt eine halbtransparente Optik

    # 🔹 Falls relevant, wende Kettenfarben an
    if "Chain_Interaction_Network" in network_title or "Combined_Network" in network_title:
        apply_pdb_chain_colors_pdb_style()

    # 🔹 Netzwerk speichern
    pdb_output_path = os.path.join(run_output_path, network_title)
    os.makedirs(pdb_output_path, exist_ok=True)
    network_file = os.path.join(pdb_output_path, f"{network_title}.cyjs")
    p4c.export_network(network_file, type="cyjs")

    print(f"✅ Netzwerk gespeichert unter: {network_file}")

def apply_pdb_chain_colors_pdb_style():
    """
    Setzt PDB-Standardfarben für Ketten in Cytoscape.
    """
    nodes = p4c.get_table_columns(columns=["id"])
    if nodes is None or "id" not in nodes:
        print("⚠️ Keine Nodes gefunden!")
        return

    chain_list = sorted(set(node_id.split(":")[-1] for node_id in nodes["id"]))  # Alle Ketten finden
    node_colors = {}

    for index, chain_id in enumerate(chain_list):
        color = PDB_CHAIN_COLORS[index % len(PDB_CHAIN_COLORS)]  # Zyklisches Schema für PDB-Kettenfarben
        for node_id in nodes["id"]:
            if node_id.endswith(chain_id):
                node_colors[node_id] = color

    p4c.set_node_color_mapping(
        table_column="id",
        table_column_values=list(node_colors.keys()),
        colors=list(node_colors.values()),
        mapping_type="d"
    )

    print(f"✅ {len(chain_list)} Kettenfarben nach PDB-Standard angewendet!")
