import os
import csv
from Bio import SeqIO
import re
from config_loader import config

# 🔹 Pfade aus Config laden
PDB_FASTA_PATH = config["pdb_fasta_path"]
UNIPROT_FASTA_PATH = config["uniprot_fasta_path"]
SIFTS_TSV_PATH = config["sifts_tsv_path"]

# 🔹 Dictionaries für schnelles Lookup
pdb_to_uniprot = {}  # Maps PDB-chain → UniProt ID
uniprot_dict = {}  # Maps UniProt ID → Protein name

def load_sifts_mapping(tsv_path):
    """
    Loads SIFTS mapping from 'pdb_chain_uniprot.tsv' to map PDB chains to UniProt IDs.
    """
    global pdb_to_uniprot
    print("\n📥 Loading SIFTS PDB-UniProt mapping...")

    with open(tsv_path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or len(row) < 3:
                continue  # Überspringe fehlerhafte Zeilen

            pdb_id, chain, uniprot_id = row[0].strip().lower(), row[1].strip().upper(), row[2].strip()
            key = f"{pdb_id}_{chain}"
            pdb_to_uniprot[key] = uniprot_id

    print(f"✅ {len(pdb_to_uniprot)} PDB chains successfully mapped to UniProt IDs.")

    # Debugging: Zeige nur die ersten 10 Einträge
    for i, (key, value) in enumerate(pdb_to_uniprot.items()):
        if i >= 10:
            break
        print(f"🔍 Mapping Check {i+1}: {key} → {value}")

def load_uniprot_fasta(fasta_path):
    """
    Loads UniProt FASTA and extracts UniProt ID → Protein name mappings.
    """
    global uniprot_dict
    print("\n📥 Loading UniProt FASTA data...")

    for record in SeqIO.parse(fasta_path, "fasta"):
        uniprot_id = record.id.split("|")[1]  # Extracts UniProt ID from FASTA header
        protein_name = record.description.split(" ", 1)[1]  # Extracts protein name
        uniprot_dict[uniprot_id] = protein_name

    print(f"✅ {len(uniprot_dict)} entries loaded from UniProt FASTA.")

def load_pdb_fasta(pdb_fasta_path):
    """
    Loads PDB FASTA file and stores it as a dictionary.

    Returns:
        dict: Dictionary mapping pdb_id_chain → {info, sequence}.
    """
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
    Bestimmt den Molekülnamen, Typ und die UniProt-ID für eine gegebene PDB-Kette.
    Falls die PDB-ID nicht echt ist, wird der SIFTS-Abgleich übersprungen.
    """
    search_key = f"{pdb_id.lower()}_{chain_id.upper()}"

    # Falls die PDB-ID nicht aus 4 Zeichen besteht → Kein SIFTS-Abgleich!
    if len(pdb_id) != 4:
        print(f"⚠ Keine echte PDB-ID ({pdb_id}). Überspringe SIFTS-Abgleich und nutze UniProt FASTA...")
        return determine_from_fasta(search_key, pdb_fasta)

    # Normaler SIFTS-Abgleich für echte PDB-IDs
    uniprot_id = pdb_to_uniprot.get(search_key)
    if uniprot_id:
        protein_name = uniprot_dict.get(uniprot_id, "Unknown Protein")
        return protein_name, "Protein", uniprot_id

    return determine_from_fasta(search_key, pdb_fasta)

def determine_from_fasta(search_key, pdb_fasta):
    """Falls SIFTS nicht verwendet wird, bestimme Namen und Typ nur aus der PDB FASTA."""
    if search_key in pdb_fasta:
        fasta_info = pdb_fasta[search_key]["info"]
        sequence = pdb_fasta[search_key]["sequence"]
        cleaned_info = re.sub(r"mol:\w+\s*", "", fasta_info)
        cleaned_info = re.sub(r"length:\d+\s*", "", cleaned_info).strip()
        molecule_type = "Protein" if "mol:protein" in fasta_info else "Nucleic Acid"
        return cleaned_info, molecule_type, None

    return "Unknown", "Unknown", None


def process_molecule_info(combined_data):
    """
    Assigns molecule names and types to chains.
    """
    print("\n🔍 Assigning molecule names and types...")
    pdb_fasta = load_pdb_fasta(PDB_FASTA_PATH)

    for structure_data in combined_data:
        pdb_id = structure_data["pdb_id"].lower()
        for chain in structure_data["atom_data"]:
            chain_id = chain["chain_id"].upper()
            name, mol_type, uniprot_id = determine_molecule_info(pdb_id, chain_id, pdb_fasta)

            chain["molecule_name"] = name
            chain["molecule_type"] = mol_type
            chain["uniprot_id"] = uniprot_id

    # Debugging: Zeige maximal 10 Ketten
    print("\n📌 UniProt Assignments (max. 10 examples):")
    for i, structure_data in enumerate(combined_data):
        if i >= 10:
            break
        pdb_id = structure_data["pdb_id"]
        for chain in structure_data["atom_data"]:
            print(f"  🔹 {pdb_id}_{chain['chain_id']}: {chain['molecule_name']} ({chain['molecule_type']}) UniProt-ID: {chain['uniprot_id']}")

# Lade die Daten bei Modulimport
load_sifts_mapping(SIFTS_TSV_PATH)
load_uniprot_fasta(UNIPROT_FASTA_PATH)
