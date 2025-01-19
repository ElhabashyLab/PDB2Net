import os
import requests
from rcsbsearchapi import TextQuery, AttributeQuery
from rcsbsearchapi import rcsb_attributes as attrs

# **1. Einzeldatei herunterladen**
def download_single_file():
    pdb_id = input("Geben Sie die PDB-ID der Datei ein (z. B. 1A4W): ").strip().upper()
    file_format = input("Wählen Sie das Dateiformat ('pdb' oder 'cif', default: pdb): ").strip() or "pdb"
    save_dir = input("Geben Sie das Zielverzeichnis ein (default: ./structures): ").strip() or "structures"
    os.makedirs(save_dir, exist_ok=True)
    url = f"https://files.rcsb.org/download/{pdb_id}.{file_format}"

    print(f"Lade Datei {pdb_id} herunter...")
    response = requests.get(url)
    if response.status_code == 200:
        file_path = os.path.join(save_dir, f"{pdb_id}.{file_format}")
        with open(file_path, "wb") as file:
            file.write(response.content)
        print(f"Datei gespeichert: {file_path}")
    else:
        print(f"Fehler beim Herunterladen von {pdb_id}: {response.status_code}")

# **2. Suche basierend auf Nutzerwahl**
def search_and_download():
    print("\nWählen Sie eine Suchoption aus:")
    print("1. Volltextsuche (TextQuery)")
    print("2. Attributbasierte Suche (AttributeQuery)")
    print("3. Sequenzähnlichkeit (SequenceQuery)")
    print("4. Strukturähnlichkeit (StructSimilarityQuery)")
    print("5. Strukturmotiv-Suche (StructMotifQuery)")
    print("6. Chemische Ähnlichkeit (ChemSimilarityQuery)")
    choice = input("Geben Sie Ihre Wahl ein (1-6): ").strip()

    # Initialisierung der Query
    query = None

    if choice == "1":
        text_query = input("Geben Sie einen Textsuchbegriff ein (z. B. 'ATP'): ").strip()
        query = TextQuery(text_query)

    elif choice == "2":
        print("\nVerfügbare Attribute:\n")
        print("1. Ligandensuche (rcsb_nonpolymer_entity.chem_comp.chem_comp_id)")
        print("2. Organismensuche (rcsb_entity_source_organism.taxonomy_lineage.name)")
        print("3. Auflösung (rcsb_entry_info.resolution_combined)")
        print("4. Experimenttyp (rcsb_entry_info.experimental_method)")

        attr_choice = input("Wählen Sie ein Attribut (1-4): ").strip()
        if attr_choice == "1":
            attribute = attrs.rcsb_nonpolymer_entity.chem_comp.chem_comp_id
        elif attr_choice == "2":
            attribute = attrs.rcsb_entity_source_organism.taxonomy_lineage.name
        elif attr_choice == "3":
            attribute = attrs.rcsb_entry_info.resolution_combined
        elif attr_choice == "4":
            attribute = attrs.rcsb_entry_info.experimental_method
        else:
            print("Ungültige Auswahl. Kehre zum Hauptmenü zurück.")
            return

        value = input(f"Geben Sie den Wert für das Attribut ein (z. B. 'ATP', 'Homo sapiens', '< 2.0'): ").strip()
        query = AttributeQuery(attribute == value)

    elif choice == "3":
        sequence = input("Geben Sie eine Protein- oder Nukleotidsequenz ein: ").strip()
        query_type = input("Ist dies eine Protein- oder DNA-Sequenz? ('protein' oder 'dna'): ").strip()
        query = SequenceQuery(sequence, query_type)

    elif choice == "4":
        pdb_id = input("Geben Sie die PDB-ID einer Struktur ein, nach der gesucht werden soll: ").strip()
        query = StructSimilarityQuery(pdb_id, "structure")

    elif choice == "5":
        motif = input("Geben Sie ein Strukturmotiv ein (z. B. 'C-alpha-helix'): ").strip()
        query = StructMotifQuery(motif)

    elif choice == "6":
        chemical = input("Geben Sie eine chemische Verbindung oder Formel ein (z. B. 'ATP'): ").strip()
        query = ChemSimilarityQuery(chemical, "exact_match")

    else:
        print("Ungültige Eingabe. Kehre zum Hauptmenü zurück.")
        return

    try:
        # Ergebnisse korrekt abrufen
        results = list(query)
        print(f"Gefundene PDB-IDs: {results}")

        if results:
            file_format = input("Wählen Sie das Dateiformat ('pdb' oder 'cif', default: pdb): ").strip() or "pdb"
            save_dir = input("Geben Sie das Zielverzeichnis ein (default: ./structures): ").strip() or "structures"
            download_pdb_files(results, save_dir, file_format)
        else:
            print("Keine Ergebnisse gefunden.")
    except Exception as e:
        print(f"Fehler: {e}")


# **3. Herunterladen von mehreren Dateien**
def download_pdb_files(pdb_ids, save_dir="structures", file_format="pdb"):
    os.makedirs(save_dir, exist_ok=True)
    for pdb_id in pdb_ids:
        url = f"https://files.rcsb.org/download/{pdb_id}.{file_format}"
        response = requests.get(url)
        if response.status_code == 200:
            file_path = os.path.join(save_dir, f"{pdb_id}.{file_format}")
            with open(file_path, "wb") as file:
                file.write(response.content)
            print(f"Datei gespeichert: {file_path}")
        else:
            print(f"Fehler beim Herunterladen von {pdb_id}: {response.status_code}")

# **4. Hauptfunktion**
if __name__ == "__main__":
    print("Willkommen zur flexiblen PDB-Datenverwaltung!")
    print("1. Einzelne Datei herunterladen")
    print("2. Mehrere Dateien suchen und herunterladen")
    choice = input("Geben Sie Ihre Wahl ein (1 oder 2): ").strip()

    if choice == "1":
        download_single_file()
    elif choice == "2":
        search_and_download()
    else:
        print("Ungültige Eingabe. Programm beendet.")
