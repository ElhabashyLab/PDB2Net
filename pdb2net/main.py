from file_parser import read_files_from_csv
from data_processor import process_structure
from unknown_molecule_uniprot import process_unknown_molecules, list_unknown_molecules


def main(csv_path):
    """
    Hauptfunktion, die die Datenstrukturen erstellt, unbekannte Moleküle verarbeitet und ausgibt.
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

    except Exception as e:
        print(f"Fehler: {e}")


if __name__ == "__main__":
    csv_path = "C:\\Users\\Gregor\\Documents\\Uni Bioinformatik\\9. Semester\\B.A\\PDBFiles\\PathsCSV.csv"
    main(csv_path)
