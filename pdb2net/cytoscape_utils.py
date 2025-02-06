from py2cytoscape.data.cyrest_client import CyRestClient


def create_cytoscape_network(results):
    """
    Erstellt ein Netzwerk in Cytoscape basierend auf den Analyseergebnissen.
    """
    # Verbinde mit Cytoscape
    cy = CyRestClient()

    # Knoten und Kanten aus den Ergebnissen erstellen
    nodes = set()
    edges = []
    for result in results:
        chain_a = result["chain_a"]
        chain_b = result["chain_b"]
        nodes.add(chain_a)
        nodes.add(chain_b)
        edges.append({
            'source': chain_a,
            'target': chain_b,
            'interaction': 'proximity',
            'ca_cb_count': result["ca_cb_count"],
            'all_atoms_close_count': result["all_atoms_close_count"],
            'file': result["file_path"]
        })

    # Erstelle die Knotenliste
    nodes_data = [{'data': {'id': node, 'label': node}} for node in nodes]

    # Erstelle die Kantenliste
    edges_data = [{'data': edge} for edge in edges]

    # Netzwerkdaten vorbereiten
    network_data = {'elements': {'nodes': nodes_data, 'edges': edges_data}}

    # Netzwerk in Cytoscape erstellen
    network = cy.network.create_from(network_data)

    # Layout anwenden
    apply_layout(cy, network)

    # Stil anwenden
    apply_default_style(cy, network)

    print(f"Netzwerk mit {len(nodes)} Knoten und {len(edges)} Kanten wurde in Cytoscape erstellt.")


def apply_layout(cy, network, layout_name='force-directed'):
    """
    Wendet ein Layout auf das Netzwerk an.
    """
    try:
        cy.layout.apply(name=layout_name, network=network)
        print(f"Layout '{layout_name}' erfolgreich angewendet.")
    except Exception as e:
        print(f"Fehler beim Anwenden des Layouts: {e}")


def apply_default_style(cy, network):
    """
    Wendet den Standardstil auf das Netzwerk an.
    """
    try:
        cy.style.apply(name='default', network=network)
        print("Standardstil erfolgreich angewendet.")
    except Exception as e:
        print(f"Fehler beim Anwenden des Standardstils: {e}")


def apply_custom_style(cy, network, style_name='custom', edge_attribute='ca_cb_count', color_map=None):
    """
    Erstellt und wendet einen benutzerdefinierten Stil auf das Netzwerk an.
    """
    try:
        style = cy.style.create(style_name)

        # Definiere eine kontinuierliche Kantenbreite basierend auf 'ca_cb_count'
        style.create_mapping(
            column=edge_attribute,
            col_type='edge',
            visual_prop='EDGE_WIDTH',
            mapping_type='continuous',
            points=[{'value': 0, 'lesser': 1, 'equal': 2, 'greater': 3}]
        )

        # Definiere eine kontinuierliche Farbe basierend auf einem Attribut
        if color_map:
            style.create_mapping(
                column=edge_attribute,
                col_type='edge',
                visual_prop='EDGE_STROKE_UNSELECTED_PAINT',
                mapping_type='continuous',
                points=color_map
            )

        # Anwenden des Stils
        cy.style.apply(style=style, network=network)
        print(f"Benutzerdefinierter Stil '{style_name}' erfolgreich angewendet.")
    except Exception as e:
        print(f"Fehler beim Erstellen des benutzerdefinierten Stils: {e}")


def export_network(cy, network, output_file, export_format='png'):
    """
    Exportiert das Netzwerk als Bild oder Projektdatei.
    """
    try:
        cy.network.export(network=network, file=output_file, type=export_format)
        print(f"Netzwerk erfolgreich exportiert nach: {output_file}")
    except Exception as e:
        print(f"Fehler beim Exportieren des Netzwerks: {e}")
