from file_parser import read_files_from_csv
from data_processor import process_structure
from unknown_molecule_uniprot import process_molecule_info
from distances import calculate_distances_with_ckdtree

PDB_FASTA_PATH = "C:\\Users\\Gregor\\Documents\\Uni Bioinformatik\\9. Semester\\B.A\\Neuer Ordner\\pdb_seqres.txt"
UNIPROT_FASTA_PATH = "C:\\Users\\Gregor\\Documents\\Uni Bioinformatik\\9. Semester\\B.A\\Neuer Ordner\\uniprot_sprot.fasta"


def main(csv_path):
    """Hauptfunktion, die die Daten verarbeitet."""
    try:
        print("\n📂 Lade PDB-Dateien aus CSV...")
        structures = read_files_from_csv(csv_path)
        combined_data = [process_structure(structure) for structure in structures]

        print("\n🔍 Bestimme Namen, Typen und Sequenzen aus PDB FASTA...")
        process_molecule_info(combined_data)

        print("\n📏 Berechnung der Atomdistanzen...")
        results = calculate_distances_with_ckdtree(combined_data)

        print("\n📊 Ergebnisse der Distanzberechnung:")
        for result in results:
            print(f"\n📄 Datei: {result['file_path']}")
            print(f"  🔗 Kettenpaar: {result['chain_a']} - {result['chain_b']}")
            print(f"  ⚛ Cα/CNN mit Abstand <15 Å: {result['ca_nn_count']}")
            print(f"  🔍 Atome mit Abstand <5 Å: {result['all_atoms_close_count']}")

    except Exception as e:
        print(f"❌ Fehler: {e}")


if __name__ == "__main__":
    csv_path = "C:\\Users\\Gregor\\Documents\\Uni Bioinformatik\\9. Semester\\B.A\\PDBFiles\\PathsCSV.csv"
    main(csv_path)
