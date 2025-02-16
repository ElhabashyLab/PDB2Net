import os
import csv
from Bio import SeqIO
import re

# 🔹 Global paths for required data files
PDB_FASTA_PATH = "C:\\Users\\Gregor\\Documents\\Uni Bioinformatik\\9. Semester\\B.A\\Neuer Ordner\\pdb_seqres.txt"
UNIPROT_FASTA_PATH = "C:\\Users\\Gregor\\Documents\\Uni Bioinformatik\\9. Semester\\B.A\\Neuer Ordner\\uniprot_sprot.fasta"
SIFTS_TSV_PATH = "C:\\Users\\Gregor\\Documents\\Uni Bioinformatik\\9. Semester\\B.A\\SIFTS\\pdb_chain_uniprot.tsv"

# 🔹 Dictionaries for fast lookup
pdb_to_uniprot = {}  # Maps PDB-chain → UniProt ID
uniprot_dict = {}  # Maps UniProt ID → Protein name

def load_sifts_mapping(tsv_path):
    """
    Loads SIFTS mapping from 'pdb_chain_uniprot.tsv' to map PDB chains to UniProt IDs.

    Args:
        tsv_path (str): Path to the SIFTS TSV file.

    Populates:
        pdb_to_uniprot (dict): Dictionary mapping PDB chain IDs (pdb_id_chain) to UniProt IDs.
    """
    global pdb_to_uniprot
    print("\n📥 Loading SIFTS PDB-UniProt mapping...")

    with open(tsv_path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row or len(row) < 3:  # Skip invalid or empty rows
                print(f"⚠ Warning: Invalid row in SIFTS: {row}")
                continue

            try:
                pdb_id, chain, uniprot_id = row[0].strip().lower(), row[1].strip().upper(), row[2].strip()
                key = f"{pdb_id}_{chain}"
                pdb_to_uniprot[key] = uniprot_id
            except Exception as e:
                print(f"⚠ Error processing row {row}: {e}")

    print(f"✅ {len(pdb_to_uniprot)} PDB chains successfully mapped to UniProt IDs.")

def load_uniprot_fasta(fasta_path):
    """
    Loads UniProt SwissProt FASTA file and extracts UniProt ID → Protein name mappings.

    Args:
        fasta_path (str): Path to the UniProt FASTA file.

    Populates:
        uniprot_dict (dict): Dictionary mapping UniProt IDs to protein names.
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
                # Store previous entry before processing the new one
                if current_key and current_seq:
                    pdb_sequences[current_key]["sequence"] = "".join(current_seq)

                # Extract PDB ID and chain
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

        # Store the last sequence
        if current_key and current_seq:
            pdb_sequences[current_key]["sequence"] = "".join(current_seq)

    return pdb_sequences

def determine_molecule_info(pdb_id, chain_id, pdb_fasta):
    """
    Determines the molecule name and type for a given PDB chain.

    Args:
        pdb_id (str): PDB ID in lowercase.
        chain_id (str): Chain ID in uppercase.
        pdb_fasta (dict): Dictionary with PDB FASTA data.

    Returns:
        tuple: (Molecule name, Molecule type)
    """
    search_key = f"{pdb_id.lower()}_{chain_id.upper()}"
    print(f"🔍 Determining name/type for {search_key}")  # Debugging

    # Check SIFTS mapping for a UniProt ID
    if search_key in pdb_to_uniprot:
        uniprot_id = pdb_to_uniprot[search_key]
        if uniprot_id in uniprot_dict:
            print(f"✅ {search_key}: UniProt match found: {uniprot_id} → {uniprot_dict[uniprot_id]}")
            return uniprot_dict[uniprot_id], "Protein"

    # If no UniProt match, fallback to PDB FASTA
    if search_key in pdb_fasta:
        fasta_info = pdb_fasta[search_key]["info"]
        sequence = pdb_fasta[search_key]["sequence"]

        # **Remove metadata like 'mol:protein', 'mol:na', 'length:XYZ'**
        cleaned_info = re.sub(r"mol:\w+\s*", "", fasta_info)  # Remove 'mol:protein' or 'mol:na'
        cleaned_info = re.sub(r"length:\d+\s*", "", cleaned_info)  # Remove 'length:XYZ'
        cleaned_info = cleaned_info.strip()  # Trim whitespace

        # Determine molecule type (Protein or Nucleic Acid)
        molecule_type = "Protein" if "mol:protein" in fasta_info else "Nucleic Acid"

        print(f"✅ {search_key}: PDB FASTA fallback → {cleaned_info}")
        return cleaned_info, molecule_type

    # If no match found, return "Unknown"
    print(f"⚠ {search_key}: NO match found")  # Debugging
    return "Unknown", "Unknown"

def process_molecule_info(combined_data):
    """
    Processes molecular information for each chain in the dataset.

    Args:
        combined_data (list): List of dictionaries containing parsed PDB structures.

    Updates:
        Each chain in combined_data with:
            - molecule_name
            - molecule_type
    """
    print("\n🔍 Determining molecule names and types for PDB chains...")

    # Load PDB FASTA as fallback data
    pdb_fasta = load_pdb_fasta(PDB_FASTA_PATH)

    for structure_data in combined_data:
        pdb_id = structure_data["pdb_id"].lower()  # ✅ Always use extracted PDB ID
        for chain in structure_data["atom_data"]:
            chain_id = chain["chain_id"].upper()
            name, mol_type = determine_molecule_info(pdb_id, chain_id, pdb_fasta)
            chain["molecule_name"] = name
            chain["molecule_type"] = mol_type

            print(f"✅ {pdb_id}_{chain_id}: {name} ({mol_type})")

# 🔹 Initialization function to preload mappings
def initialize():
    """
    Loads SIFTS, UniProt FASTA, and PDB FASTA for fast lookups.
    """
    load_sifts_mapping(SIFTS_TSV_PATH)
    load_uniprot_fasta(UNIPROT_FASTA_PATH)

# 🔹 Run initialization on module import
initialize()
