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
uniprot_dict = {}  # Maps UniProt ID → Protein name

def load_sifts_mapping(tsv_path):
    """
    Loads the SIFTS mapping from 'pdb_chain_uniprot.tsv' to map PDB chains to UniProt IDs.

    Args:
        tsv_path (str): Path to the SIFTS TSV file.

    Updates:
        pdb_to_uniprot (dict): Dictionary mapping PDB-chain pairs to UniProt IDs.
    """
    global pdb_to_uniprot
    print("\nLoading SIFTS PDB-UniProt mapping...")

    with open(tsv_path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or len(row) < 3:
                continue  # Skip invalid rows

            pdb_id, chain, uniprot_id = row[0].strip().lower(), row[1].strip().upper(), row[2].strip()
            key = f"{pdb_id}_{chain}"
            pdb_to_uniprot[key] = uniprot_id

    print(f"{len(pdb_to_uniprot)} PDB chains successfully mapped to UniProt IDs.")

def load_uniprot_fasta(fasta_path):
    """
    Loads the UniProt FASTA file and extracts mappings of UniProt ID → Protein name.

    Args:
        fasta_path (str): Path to the UniProt FASTA file.

    Updates:
        uniprot_dict (dict): Dictionary mapping UniProt IDs to protein names.
    """
    global uniprot_dict
    print("\nLoading UniProt FASTA data...")

    for record in SeqIO.parse(fasta_path, "fasta"):
        uniprot_id = record.id.split("|")[1]  # Extract UniProt ID from FASTA header
        protein_name = record.description.split(" ", 1)[1]  # Extract protein name
        uniprot_dict[uniprot_id] = protein_name

    print(f"{len(uniprot_dict)} entries loaded from UniProt FASTA.")

def load_pdb_fasta(pdb_fasta_path):
    """
    Loads the PDB FASTA file and stores it as a dictionary.

    Args:
        pdb_fasta_path (str): Path to the PDB FASTA file.

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
    Determines the molecule name, type, and UniProt ID for a given PDB chain.

    Args:
        pdb_id (str): The PDB ID of the structure.
        chain_id (str): The chain identifier.
        pdb_fasta (dict): Dictionary containing PDB FASTA sequences.

    Returns:
        tuple: (molecule_name, molecule_type, uniprot_id)
    """
    search_key = f"{pdb_id.lower()}_{chain_id.upper()}"

    # If PDB ID is not 4 characters long, do not use SIFTS mapping
    if len(pdb_id) != 4:
        return determine_from_fasta(search_key, pdb_fasta)

    # Standard SIFTS lookup for valid PDB IDs
    uniprot_id = pdb_to_uniprot.get(search_key)
    if uniprot_id:
        protein_name = uniprot_dict.get(uniprot_id, "Unknown Protein")
        return protein_name, "Protein", uniprot_id

    return determine_from_fasta(search_key, pdb_fasta)

def determine_from_fasta(search_key, pdb_fasta):
    """
    If SIFTS is not used, determine molecule name and type from PDB FASTA.

    Args:
        search_key (str): Formatted key for PDB FASTA lookup.
        pdb_fasta (dict): Dictionary containing PDB FASTA sequences.

    Returns:
        tuple: (molecule_name, molecule_type, uniprot_id)
    """
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
    Assigns molecule names and types to chains in the given dataset.

    Args:
        combined_data (list): List of parsed structure data.
    """
    print("\nAssigning molecule names and types...")
    pdb_fasta = load_pdb_fasta(PDB_FASTA_PATH)

    for structure_data in combined_data:
        pdb_id = structure_data["pdb_id"].lower()
        for chain in structure_data["atom_data"]:
            chain_id = chain["chain_id"].upper()
            name, mol_type, uniprot_id = determine_molecule_info(pdb_id, chain_id, pdb_fasta)

            chain["molecule_name"] = name
            chain["molecule_type"] = mol_type
            chain["uniprot_id"] = uniprot_id

    # Debugging: Show up to 10 assignments
    print("\nUniProt Assignments (example for one file):")
    for i, structure_data in enumerate(combined_data):
        if i >= 1:
            break
        pdb_id = structure_data["pdb_id"]
        for chain in structure_data["atom_data"]:
            print(f"  {pdb_id}_{chain['chain_id']}: {chain['molecule_name']} ({chain['molecule_type']}) UniProt-ID: {chain['uniprot_id']}")

# Load data when module is imported
load_sifts_mapping(SIFTS_TSV_PATH)
load_uniprot_fasta(UNIPROT_FASTA_PATH)
