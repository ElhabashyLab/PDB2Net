import os
import csv
from Bio import SeqIO

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
    """
    search_key = f"{pdb_id.lower()}_{chain_id.upper()}"

    # 🔹 Falls SIFTS einen UniProt-Match hat → Namen aus UniProt nehmen
    if search_key in pdb_to_uniprot:
        uniprot_id = pdb_to_uniprot[search_key]
        if uniprot_id in uniprot_dict:
            return uniprot_dict[uniprot_id], "Protein"

    # 🔹 Falls kein UniProt-Match → Fallback auf PDB FASTA
    if search_key in pdb_fasta:
        fasta_info = pdb_fasta[search_key]["info"]
        sequence = pdb_fasta[search_key]["sequence"]
        molecule_type = "Protein" if "mol:protein" in fasta_info else "Nucleic Acid"
        return fasta_info.split("length:")[-1].strip(), molecule_type

    # 🔹 Falls keine Daten vorhanden sind
    return "Unknown", "Unknown"


def process_molecule_info(combined_data):
    """
    Geht alle PDB-Strukturen durch und ergänzt Name & Typ.
    """
    print("\n🔍 Bestimme Namen und Typen für PDB-Ketten...")

    # Lade PDB FASTA für Fallback
    pdb_fasta = load_pdb_fasta(PDB_FASTA_PATH)

    for structure_data in combined_data:
        pdb_id = os.path.basename(structure_data["file_path"]).split(".")[0].lower()
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