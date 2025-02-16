import py4cytoscape as p4c
import time
import pandas as pd

def create_cytoscape_network(results):
    """
    Erstellt ein Netzwerk in Cytoscape basierend auf den Daten in `results`.
    """

    # 🔹 1. Überprüfen, ob Cytoscape läuft
    try:
        status = p4c.cytoscape_ping()
        print("🌐 Verbindung zu Cytoscape erfolgreich! Cytoscape ist online.")
    except Exception as e:
        print(f"❌ Cytoscape scheint nicht zu laufen! Bitte starte Cytoscape und versuche es erneut. Fehler: {e}")
        return

    # 🔹 2. Falls Netzwerke existieren, löschen
    try:
        existing_networks = p4c.get_network_list()
        if existing_networks:
            print(f"🗑️ Lösche {len(existing_networks)} vorhandene Netzwerke...")
            for net in existing_networks:
                p4c.delete_network(net)
            time.sleep(1)  # Warte kurz, damit Cytoscape Zeit hat
    except Exception as e:
        print(f"⚠️ Warnung: Konnte alte Netzwerke nicht löschen: {e}")

    # 🔹 3. Knoten (Nodes) und Kanten (Edges) aus `results` extrahieren
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

    # 🔹 4. DataFrames für Cytoscape vorbereiten
    nodes_df = pd.DataFrame({"id": list(unique_nodes), "label": list(unique_nodes)})
    edges_df = pd.DataFrame(edges)

    # 🔹 5. Netzwerk in Cytoscape erstellen
    print("📡 Erstelle ein neues Netzwerk in Cytoscape...")
    network_suid = p4c.create_network_from_data_frames(nodes_df, edges_df, title="Protein Interaction Network")

    # 🔹 6. Layout anwenden & Netzwerk anzeigen
    print("🎨 Wende Layout an...")
    time.sleep(1)  # Warte kurz, um sicherzustellen, dass das Netzwerk geladen ist
    p4c.layout_network(layout_name="circular")

    print("✅ Netzwerk wurde erfolgreich in Cytoscape erstellt!")
