import py4cytoscape as p4c
import time
import os
import pandas as pd
from config_loader import config


def create_cytoscape_network(results, network_title="Protein_Interaction_Network", run_output_path="."):
    """
    Erstellt ein Netzwerk in Cytoscape basierend auf den Interaktionsdaten.

    Args:
        results (list): Liste von Dictionaries mit Interaktionsdaten.
        network_title (str): Name des Netzwerks (meistens PDB-ID oder "Combined_Network").
        run_output_path (str): Das Ausgabeverzeichnis für den aktuellen Lauf.
    """

    # 🔹 Prüfe, ob Cytoscape läuft
    try:
        p4c.cytoscape_ping()
        print(f"🌐 Verbunden mit Cytoscape! Erstelle Netzwerk: {network_title}")
    except Exception as e:
        print(f"❌ Cytoscape läuft nicht! Starte Cytoscape und versuche es erneut. Fehler: {e}")
        return

    # 🔹 Lösche ältere Netzwerke (behalte nur die letzten N)
    existing_networks = p4c.get_network_list()
    while len(existing_networks) > config["keep_last_n_networks"]:
        oldest_network = existing_networks.pop(0)  # Entferne das älteste Netzwerk
        p4c.delete_network(oldest_network)

    # 🔹 Nodes & Edges vorbereiten
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
                "all_atoms_count": entry["all_atoms_count"]  # 🔥 Nur dieses Kriterium bleibt!
            })

    # 🔹 Erstelle DataFrames für Cytoscape
    nodes_df = pd.DataFrame({"id": list(unique_nodes), "label": list(unique_nodes)})
    edges_df = pd.DataFrame(edges)

    # 🔹 Falls keine Kanten vorhanden sind, Netzwerk nicht erstellen
    if edges_df.empty:
        print("⚠️ Keine gültigen Kanten gefunden. Das Netzwerk wird nicht erstellt.")
        return

    # 🔹 Korrekte Umbenennung in 'source' und 'target' für Cytoscape
    edges_df.rename(columns={"chain_a": "source", "chain_b": "target"}, inplace=True)
    edges_df["interaction"] = "interacts_with"

    # 🔹 Debugging: Zeige die ersten 5 Zeilen der Edges-Tabelle
    print(edges_df.head())

    # 🔹 Erstelle das Netzwerk in Cytoscape
    network_suid = p4c.create_network_from_data_frames(nodes_df, edges_df, title=network_title)

    # 🔹 Wende Layout an
    p4c.layout_network(layout_name="circular")

    # 🔹 Speichere das Netzwerk in einem PDB-Unterordner des Run-Output-Pfads
    pdb_output_path = os.path.join(run_output_path, network_title)
    os.makedirs(pdb_output_path, exist_ok=True)
    network_file = os.path.join(pdb_output_path, f"{network_title}.cyjs")
    p4c.export_network(network_file, type="cyjs")

    print(f"✅ Netzwerk gespeichert unter: {network_file}")
