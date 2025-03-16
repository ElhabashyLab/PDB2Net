import py4cytoscape as p4c
import os
import pandas as pd
from config_loader import config


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

    # 🔹 Setze den "BioPAX_SIF"-Style
    available_styles = p4c.get_visual_style_names()
    print("🎨 Verfügbare Styles in Cytoscape:", available_styles)

    if "BioPAX_SIF" in available_styles:
        p4c.set_visual_style("BioPAX_SIF")
        print("✅ BioPAX_SIF-Style angewendet!")
    else:
        print("⚠️ BioPAX_SIF-Style nicht verfügbar. Bitte prüfen!")

    # 🔹 Netzwerk speichern
    pdb_output_path = os.path.join(run_output_path, network_title)
    os.makedirs(pdb_output_path, exist_ok=True)
    network_file = os.path.join(pdb_output_path, f"{network_title}.cyjs")
    p4c.export_network(network_file, type="cyjs")

    print(f"✅ Netzwerk gespeichert unter: {network_file}")
