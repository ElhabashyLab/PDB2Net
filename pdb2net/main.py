import os
from file_parser import parse_structure
from csv_reader import read_csv

def main():
    # Pfad zur Eingabe-CSV
    csv_path = "C:\\Users\\Gregor\\Documents\\Uni Bioinformatik\\9. Semester\\B.A\\PDBFiles\\PathsCSV.csv"

    # CSV-Datei lesen
    try:
        file_paths = read_csv(csv_path)
        print(f"Gefundene Dateien: {len(file_paths)}")
    except Exception as e:
        print(f"Fehler beim Lesen der CSV-Datei: {e}")
        return

    # Jede Datei parsen
    for file_path in file_paths:
        if os.path.exists(file_path):
            print(f"\nAnalysiere Datei: {file_path}")
            try:
                atom_chains, hetatm_molecules = parse_structure(file_path)

                # ATOM-Daten anzeigen
                print("\n--- ATOM-Daten ---")
                for chain in atom_chains:
                    print(f"Kette: {chain['chain_id']}, "
                          f"Residuen: {len(chain['residues'])}, "
                          f"Molekültyp: {chain['molecule_type']}, "
                          f"Molekülname: {chain['molecule_name']}")

                # HETATM-Daten anzeigen
                print("\n--- HETATM-Daten ---")
                for hetatm in hetatm_molecules:
                    print(f"HETATM-Kette: {hetatm['chain_id']}, "
                          f"Residuen: {len(hetatm['residues'])}, "
                          f"Molekültyp: {hetatm['molecule_type']}, "
                          f"Molekülname: {hetatm['molecule_name']}")

            except Exception as e:
                print(f"Fehler bei der Verarbeitung von {file_path}: {e}")
        else:
            print(f"Datei nicht gefunden: {file_path}")

if __name__ == "__main__":
    main()
