import os
import subprocess
from Bio import SeqIO
from Bio.Data import IUPACData
from config_loader import config

# Load paths from configuration
UNIPROT_FASTA_PATH = config["uniprot_fasta_path"]
BLAST_DB_PATH = config["blast_db_path"]
BLAST_EXECUTABLE = config["blastp_executable"]

# Biopython mapping for three-letter to one-letter amino acid codes
three_to_one = IUPACData.protein_letters_3to1

# List of 20 canonical amino acids plus selenocysteine and pyrrolysine
AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "SEC", "PYL"
}

def create_blast_database():
    """
    Creates a BLAST database from the UniProt FASTA file if it does not already exist.
    """
    db_path = os.path.join(BLAST_DB_PATH, "uniprot_db")
    if not os.path.exists(db_path + ".pin"):  # Check if BLAST database exists
        print("Creating BLAST database...")
        subprocess.run([
            "makeblastdb", "-in", UNIPROT_FASTA_PATH, "-dbtype", "prot", "-out", db_path
        ], check=True)
        print("BLAST database created successfully.")
    else:
        print("BLAST database already exists.")

def extract_sequence_from_parsed_data(chain_data):
    """
    Converts residue names (three-letter code) into a protein sequence (one-letter code).

    Args:
        chain_data (dict): Dictionary containing residue information.

    Returns:
        str or None: The converted protein sequence, or None if no valid sequence is found.
    """
    sequence = "".join([three_to_one.get(res["residue_name"].capitalize(), "X") for res in chain_data["residues"]])
    return sequence if sequence else None

def run_blast_search(query_sequence):
    """
    Performs a BLAST search using the given protein sequence.

    Args:
        query_sequence (str): The amino acid sequence to be matched.

    Returns:
        str or None: The best UniProt ID match, or None if no match is found.
    """
    query_file = "query.fasta"
    output_file = "blast_results.txt"

    # Write sequence to a FASTA file
    with open(query_file, "w") as f:
        f.write(f">query\n{query_sequence}\n")

    # Run BLAST search
    subprocess.run([
        BLAST_EXECUTABLE, "-query", query_file, "-db", os.path.join(BLAST_DB_PATH, "uniprot_db"),
        "-out", output_file, "-evalue", "1e-5", "-max_target_seqs", "8", "-outfmt", "6"
    ], check=True)

    # Read first hit result
    best_match = None
    with open(output_file, "r") as f:
        for line in f:
            best_match = line.split("\t")[1].split("|")[1]  # Extract only the UniProt ID
            break

    # Delete temporary files
    os.remove(query_file)
    os.remove(output_file)

    return best_match if best_match else None

def classify_molecule_type(chain):
    """
    Determines if a chain is a protein based on its residues.

    Args:
        chain (dict): Chain dictionary containing residue information.

    Returns:
        str: "Protein" if the chain contains only amino acids, otherwise "Nucleic Acid".
    """
    residues = chain["residues"]
    if len(residues) == 0:
        return "Unknown"

    sample_names = {res["residue_name"].strip() for res in residues}

    # Check if all residues belong to known amino acids
    return "Protein" if all(res in AMINO_ACIDS for res in sample_names) else "Nucleic Acid"

def match_sequence_to_uniprot(parsed_data):
    """
    Matches extracted sequences with UniProt via BLAST and classifies molecules.

    Args:
        parsed_data (list): List of parsed structure data with atom and sequence information.
    """
    print("Starting UniProt BLAST search...")
    create_blast_database()

    for structure in parsed_data:
        for chain in structure["atom_data"]:
            # Skip chains that already have a UniProt ID or are nucleic acids
            if chain.get("uniprot_id") not in [None, "Unknown"] or chain.get("molecule_type") in ["Nucleic Acid"]:
                continue

            sequence = extract_sequence_from_parsed_data(chain)  # Convert to one-letter code
            if not sequence:
                print(f"Warning: No sequence extracted for {chain['unique_chain_id']}")
                continue

            sequence = sequence.replace("O", "X")  # Fix incorrect "O" characters

            uniprot_id = run_blast_search(sequence)
            if uniprot_id:
                chain["uniprot_id"] = uniprot_id
                chain["molecule_name"] = f"Matched UniProt: {uniprot_id}"
                chain["molecule_type"] = "Protein"
                print(f"Matched {chain['unique_chain_id']} → {uniprot_id}")
            else:
                chain["molecule_type"] = classify_molecule_type(chain)  # Use classification if no match
                print(f"No UniProt match for {chain['unique_chain_id']}, classified as {chain['molecule_type']}")

    print("UniProt BLAST search completed.")
