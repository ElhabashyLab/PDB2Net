import os
import csv
from Bio import SeqIO
import re
from config_loader import config

# Load paths from configuration
PDB_FASTA_PATH = config["pdb_fasta_path"]
UNIPROT_FASTA_PATH = config["uniprot_fasta_path"]
SIFTS_TSV_PATH = config["sifts_tsv_path"]

# Dictionaries for fast lookups
pdb_to_uniprot = {}  # Maps PDB-chain → UniProt ID
uniprot_dict = {}    # Maps UniProt ID → Protein name

# Lazy-loading flags
_sifts_loaded = False
_uniprot_loaded = False


def load_sifts_mapping(tsv_path):
    global pdb_to_uniprot, _sifts_loaded
    if _sifts_loaded:
        return

    #print("\nLoading SIFTS PDB-UniProt mapping...")
    with open(tsv_path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or len(row) < 3:
                continue
            pdb_id, chain, uniprot_id = row[0].strip().lower(), row[1].strip().upper(), row[2].strip()
            key = f"{pdb_id}_{chain}"
            pdb_to_uniprot[key] = uniprot_id

    #print(f"{len(pdb_to_uniprot)} PDB chains successfully mapped to UniProt IDs.")
    _sifts_loaded = True


def load_uniprot_fasta(fasta_path):
    global uniprot_dict, _uniprot_loaded
    if _uniprot_loaded:
        return

    #print("\nLoading UniProt FASTA data...")
    for record in SeqIO.parse(fasta_path, "fasta"):
        uniprot_id = record.id.split("|")[1]  # Extract UniProt ID
        protein_name = record.description.split(" ", 1)[1]
        uniprot_dict[uniprot_id] = protein_name

    #print(f"{len(uniprot_dict)} entries loaded from UniProt FASTA.")
    _uniprot_loaded = True


def load_pdb_fasta(pdb_fasta_path):
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


def determine_from_fasta(search_key, pdb_fasta):
    if search_key in pdb_fasta:
        fasta_info = pdb_fasta[search_key]["info"]
        sequence = pdb_fasta[search_key]["sequence"]
        cleaned_info = re.sub(r"mol:\w+\s*", "", fasta_info)
        cleaned_info = re.sub(r"length:\d+\s*", "", cleaned_info).strip()

        if "mol:protein" in fasta_info.lower():
            molecule_type = "Protein"
        elif "mol:na" in fasta_info.lower():
            # Detailauswertung nach length:...
            if "dna/rna" in fasta_info.lower():
                molecule_type = "DNA/RNA"
            elif "dna" in fasta_info.lower():
                molecule_type = "DNA"
            elif "rna" in fasta_info.lower():
                molecule_type = "RNA"
            elif "polyribonucleotide" in fasta_info.lower():
                molecule_type = "RNA"
            elif "polydeoxyribonucleotide" in fasta_info.lower():
                molecule_type = "DNA"
            else:
                molecule_type = "Nucleic Acid"
        else:
            molecule_type = "Unknown"

        return cleaned_info, molecule_type, None
    return "Unknown", "Unknown", None


def determine_molecule_info(pdb_id, chain_id, pdb_fasta):
    search_key = f"{pdb_id.lower()}_{chain_id.upper()}"
    fasta_name, mol_type, _ = determine_from_fasta(search_key, pdb_fasta)

    name = fasta_name
    uniprot_id = None

    if len(pdb_id) == 4:
        uniprot_id = pdb_to_uniprot.get(search_key)
        if uniprot_id:
            better_name = uniprot_dict.get(uniprot_id)
            if better_name and better_name != "Unknown Protein":
                name = better_name
            mol_type = "Protein"

    return name, mol_type, uniprot_id


def process_molecule_info(combined_data):
    """
    Assigns molecule names and types to chains in the given dataset.
    Only loads mapping files once (lazy-loading).
    """
    load_sifts_mapping(SIFTS_TSV_PATH)
    load_uniprot_fasta(UNIPROT_FASTA_PATH)

    #print("\nAssigning molecule names and types...")
    pdb_fasta = load_pdb_fasta(PDB_FASTA_PATH)

    for structure_data in combined_data:
        pdb_id = structure_data["pdb_id"].lower()
        for chain in structure_data["atom_data"]:
            chain_id = chain["chain_id"].upper()
            name, mol_type, uniprot_id = determine_molecule_info(pdb_id, chain_id, pdb_fasta)

            chain["molecule_name"] = name
            chain["molecule_type"] = mol_type
            chain["uniprot_id"] = uniprot_id

    #Optional: Einmalige Ausgabe für Debugging
    print("\nUniProt Assignments (example for one file):")
    for i, structure_data in enumerate(combined_data):
        if i >= 20:
            break
        pdb_id = structure_data["pdb_id"]
        for chain in structure_data["atom_data"]:
            print(f"  {pdb_id}_{chain['chain_id']}: {chain['molecule_name']} ({chain['molecule_type']}) UniProt-ID: {chain['uniprot_id']}")
