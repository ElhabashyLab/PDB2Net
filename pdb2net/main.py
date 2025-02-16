from file_parser import read_files_from_csv
from data_processor import process_structure
from unknown_molecule_uniprot import process_molecule_info
from distances import calculate_distances_with_ckdtree
from cytoscape import create_cytoscape_network  # 🔹 Cytoscape-Integration hinzufügen

# 🔹 Externe Datei-Pfade
PDB_FASTA_PATH = "C:\\Users\\Gregor\\Documents\\Uni Bioinformatik\\9. Semester\\B.A\\Neuer Ordner\\pdb_seqres.txt"
UNIPROT_FASTA_PATH = "C:\\Users\\Gregor\\Documents\\Uni Bioinformatik\\9. Semester\\B.A\\Neuer Ordner\\uniprot_sprot.fasta"


def main(csv_path):
    """
    Hauptfunktion, die PDB-Strukturen lädt, verarbeitet, Molekülinformationen ergänzt,
    Distanzberechnungen durchführt und das Netzwerk in Cytoscape erstellt.
    """
    try:
        # 🔹 Schritt 1: Dateien aus CSV einlesen & Strukturen parsen
        print("\n📂 Lade PDB-Dateien aus CSV...")
        structures = read_files_from_csv(csv_path)
        combined_data = []

        for structure_data in structures:
            processed_data = process_structure(structure_data)
            combined_data.append(processed_data)

        # 🔹 Schritt 2: Namen & Molekültyp mit SIFTS & UniProt bestimmen
        print("\n🔍 Bestimme Namen, Typen und Sequenzen aus PDB FASTA...")
        process_molecule_info(combined_data)

        # 🔹 Kontrollausgabe: Zeige die verarbeiteten Ketten
        print("\n📌 Kontrollausgabe: Alle Ketten mit Name, Typ und Sequenz:")
        for structure in combined_data:
            print(f"\n📄 Datei: {structure['file_path']} (PDB-ID: {structure['pdb_id']})")
            for chain in structure["atom_data"]:
                print(f"  🔹 Kette: {chain['chain_id']}")
                print(f"     🏷 Name: {chain['molecule_name']}")
                print(f"     🔬 Typ: {chain['molecule_type']}")
                print(f"     🧬 Sequenz: {chain['sequence'][:50]}... (gekürzt)")

        # 🔹 Schritt 3: Distanzberechnung
        print("\n📏 Berechnung der Atomdistanzen...")
        results = calculate_distances_with_ckdtree(combined_data)

        # 🔹 Ergebnisse der Distanzberechnung ausgeben
        print("\n📊 Ergebnisse der Distanzberechnung:")
        for result in results:
            print(f"\n📄 Datei: {result['file_path']}")
            print(f"  🔗 Kettenpaar: {result['chain_a']} - {result['chain_b']}")
            print(f"  ⚛ Cα/Cβ mit Abstand <15 Å: {result['ca_nn_count']}")
            print(f"  🔍 Atome mit Abstand <5 Å: {result['all_atoms_close_count']}")

        # 🔹 Schritt 4: Netzwerk mit Cytoscape erstellen
        print("\n🌐 Erstelle Netzwerk in Cytoscape...")
        create_cytoscape_network(results)

    except Exception as e:
        print(f"❌ Fehler: {e}")


if __name__ == "__main__":
    csv_path = "C:\\Users\\Gregor\\Documents\\Uni Bioinformatik\\9. Semester\\B.A\\PDBFiles\\PathsCSV.csv"
    main(csv_path)
