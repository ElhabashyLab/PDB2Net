import os
import csv
from Bio import SeqIO
import re

# 🔹 Globale Pfade
PDB_FASTA_PATH = "C:\\Users\\Gregor\\Documents\\Uni Bioinformatik\\9. Semester\\B.A\\Neuer Ordner\\pdb_seqres.txt"
UNIPROT_FASTA_PATH = "C:\\Users\\Gregor\\Documents\\Uni Bioinformatik\\9. Semester\\B.A\\Neuer Ordner\\uniprot_sprot.fasta"
SIFTS_TSV_PATH = "C:\\Users\\Gregor\\Documents\\Uni Bioinformatik\\9. Semester\\B.A\\SIFTS\\pdb_chain_uniprot.tsv"

# 🔹 Dictionaries für schnelle Lookups
pdb_to_uniprot = {}  # PDB-Kette → UniProt-ID
uniprot_dict = {}  # UniProt-ID → Proteinname


def load_sifts_mapping(tsv_path):
    """ Lädt SIFTS pdb_chain_uniprot.tsv und speichert die Zuordnung PDB-Kette → UniProt-ID. """
    global pdb_to_uniprot
    print("\n📥 Lade SIFTS PDB-UniProt-Mapping...")

    with open(tsv_path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or len(row) < 3:  # Falls Zeile leer ist oder weniger als 3 Spalten hat
                print(f"⚠ Warnung: Ungültige Zeile in SIFTS: {row}")
                continue  # Überspringen

            try:
                pdb_id, chain, uniprot_id = row[0].strip().lower(), row[1].strip().upper(), row[2].strip()
                key = f"{pdb_id}_{chain}"
                pdb_to_uniprot[key] = uniprot_id
            except Exception as e:
                print(f"⚠ Fehler beim Verarbeiten der Zeile {row}: {e}")

    print(f"✅ {len(pdb_to_uniprot)} PDB-Ketten erfolgreich mit UniProt verknüpft.")


def load_uniprot_fasta(fasta_path):
    """ Lädt UniProt SwissProt FASTA und speichert UniProt-ID → Proteinname. """
    global uniprot_dict
    print("\n📥 Lade UniProt FASTA-Daten...")

    for record in SeqIO.parse(fasta_path, "fasta"):
        uniprot_id = record.id.split("|")[1]  # sp|P12345|HBB_HUMAN → "P12345"
        protein_name = record.description.split(" ", 1)[1]  # Proteinname nach ID
        uniprot_dict[uniprot_id] = protein_name

    print(f"✅ {len(uniprot_dict)} Einträge aus UniProt FASTA geladen.")


def load_pdb_fasta(pdb_fasta_path):
    """ Lädt PDB FASTA-Datei und speichert sie als Dictionary {pdb_id_chain: {info, sequence}}. """
    pdb_sequences = {}
    with open(pdb_fasta_path, "r") as f:
        current_key, current_seq = None, []
        for line in f:
            if line.startswith(">"):
                if current_key and current_seq:
                    pdb_sequences[current_key]["sequence"] = "".join(current_seq)
                parts = line.split()
                fasta_header = parts[0][1:]
                if "_" in fasta_header:
                    pdb_id, chain_id = fasta_header.split("_")
                    formatted_key = f"{pdb_id.lower()}_{chain_id.upper()}"
                    pdb_sequences[formatted_key] = {"info": " ".join(parts[1:]), "sequence": ""}
                    current_key = formatted_key
                    current_seq = []
            else:
                current_seq.append(line.strip())
        if current_key and current_seq:
            pdb_sequences[current_key]["sequence"] = "".join(current_seq)
    return pdb_sequences


def determine_molecule_info(pdb_id, chain_id, pdb_fasta):
    """
    Bestimmt den Namen & Molekültyp für eine PDB-Kette mit SIFTS, UniProt oder PDB FASTA.
    Entfernt alle Metadaten wie 'mol:protein', 'mol:na' und 'length:XYZ' aus dem Namen.
    """
    search_key = f"{pdb_id.lower()}_{chain_id.upper()}"
    print(f"🔍 Bestimme Name/Typ für {search_key}")  # Debugging

    # 1️⃣ Falls SIFTS einen UniProt-Match hat, nimm den Namen von UniProt
    if search_key in pdb_to_uniprot:
        uniprot_id = pdb_to_uniprot[search_key]
        if uniprot_id in uniprot_dict:
            print(f"✅ {search_key}: UniProt-Match gefunden: {uniprot_id} → {uniprot_dict[uniprot_id]}")  # Debugging
            return uniprot_dict[uniprot_id], "Protein"

    # 2️⃣ Falls kein UniProt-Match, nutze PDB FASTA
    if search_key in pdb_fasta:
        fasta_info = pdb_fasta[search_key]["info"]
        sequence = pdb_fasta[search_key]["sequence"]

        # **Entferne 'mol:protein', 'mol:na', 'length:XYZ' und extra Leerzeichen**
        cleaned_info = re.sub(r"mol:\w+\s*", "", fasta_info)  # Entfernt 'mol:protein' oder 'mol:na'
        cleaned_info = re.sub(r"length:\d+\s*", "", cleaned_info)  # Entfernt 'length:XYZ'
        cleaned_info = cleaned_info.strip()  # Entfernt überflüssige Leerzeichen

        # Bestimme den Molekültyp (Protein oder Nukleinsäure)
        molecule_type = "Protein" if "mol:protein" in fasta_info else "Nucleic Acid"

        print(f"✅ {search_key}: PDB FASTA Fallback → {cleaned_info}")  # Debugging
        return cleaned_info, molecule_type

    # 3️⃣ Falls nichts gefunden → "Unknown"
    print(f"⚠ {search_key}: KEINE Zuordnung gefunden")  # Debugging
    return "Unknown", "Unknown"


def process_molecule_info(combined_data):
    """
    Geht alle PDB-Strukturen durch und ergänzt Name & Typ.
    """
    print("\n🔍 Bestimme Namen und Typen für PDB-Ketten...")

    # Lade PDB FASTA für Fallback
    pdb_fasta = load_pdb_fasta(PDB_FASTA_PATH)

    for structure_data in combined_data:
        pdb_id = structure_data["pdb_id"].lower()  # ✅ NUR DIE EXTRAHIERTE ID VERWENDEN!
        for chain in structure_data["atom_data"]:
            chain_id = chain["chain_id"].upper()
            name, mol_type = determine_molecule_info(pdb_id, chain_id, pdb_fasta)
            chain["molecule_name"] = name
            chain["molecule_type"] = mol_type

            print(f"✅ {pdb_id}_{chain_id}: {name} ({mol_type})")



# 🔹 Initialisierung
def initialize():
    """
    Lädt SIFTS, UniProt FASTA und PDB FASTA einmalig beim Start.
    """
    load_sifts_mapping(SIFTS_TSV_PATH)
    load_uniprot_fasta(UNIPROT_FASTA_PATH)


# 🔹 Starte den Ladevorgang bei Import
initialize()