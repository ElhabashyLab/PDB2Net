from file_parser import read_files_from_csv
from data_processor import process_structure
from unknown_molecule_uniprot import process_unknown_molecules, list_unknown_molecules
from distances import calculate_distances_with_ckdtree
from cytoscape_utils import create_cytoscape_network  # Neues File importiert


def main(csv_path):
    """
    Hauptfunktion, die die Datenstrukturen erstellt, unbekannte Moleküle verarbeitet, Distanzen berechnet
    und das Netzwerk in Cytoscape erstellt.
    """
    try:
        # Schritt 1: Dateien lesen und Strukturen verarbeiten
        structures = read_files_from_csv(csv_path)
        combined_data = []
        for structure_data in structures:
            processed_data = process_structure(structure_data)
            combined_data.append(processed_data)

        # Schritt 2: Liste aller UNKNOWN-Ketten vor UniProt-Verarbeitung
        unknown_molecules = list_unknown_molecules(combined_data)
        print("\nMolekülnamen vor UniProt-Verarbeitung:")
        for entry in unknown_molecules:
            print(f"Datei: {entry['file_path']}, Kette: {entry['chain_id']}")

        # Schritt 3: Unbekannte Moleküle verarbeiten
        print("\nÜberprüfung auf unbekannte Moleküle...")
        process_unknown_molecules(combined_data)

        # Schritt 4: Ergebnisse ausgeben
        print("\nKontrolle der extrahierten Daten:")
        for data in combined_data:
            file_name = data["file_path"].split("\\")[-1]  # Extrahiere nur den Dateinamen
            print(f"\nDatei: {file_name}")
            for chain in data["atom_data"]:
                print(
                    f"  Kette: {chain['chain_id']}, Molekültyp: {chain['molecule_type']}, Molekülname: {chain['molecule_name']}"
                )

        # Schritt 5: Distanzberechnungen durchführen
        print("\nBerechnung der Atomdistanzen...")
        results = calculate_distances_with_ckdtree(combined_data)

        # Ergebnisse der Distanzberechnung ausgeben
        print("\nErgebnisse der Distanzberechnung:")
        for result in results:
            print(f"Datei: {result['file_path']}")
            print(f"  Kettenpaar: {result['chain_a']} - {result['chain_b']}")
            print(f"  Cα/Cβ mit Abstand <15 Å: {result['ca_cb_count']}")
            print(f"  Atome mit Abstand <5 Å: {result['all_atoms_close_count']}")

        # Schritt 6: Netzwerk in Cytoscape erstellen
        print("\nErstelle Netzwerk in Cytoscape...")
        create_cytoscape_network(results)

    except Exception as e:
        print(f"Fehler: {e}")


if __name__ == "__main__":
    csv_path = "C:\\Users\\Gregor\\Documents\\Uni Bioinformatik\\9. Semester\\B.A\\PDBFiles\\PathsCSV.csv"
    main(csv_path)
