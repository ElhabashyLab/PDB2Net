import os
import subprocess
from Bio.Data import IUPACData
from concurrent.futures import ThreadPoolExecutor
from config_loader import config
import uuid

# Load paths from configuration
BLAST_DB_PATH = config["blast_db_path"]
BLAST_EXECUTABLE = config["blastp_executable"]
UNIPROT_FASTA_PATH = config["uniprot_fasta_path"]

# 3-letter to 1-letter amino acid mapping
three_to_one = IUPACData.protein_letters_3to1

# Set of valid amino acid residue codes
AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "SEC", "PYL"
}


def create_blast_database():
    """
    Creates the UniProt BLAST database if it does not already exist.
    """
    db_path = os.path.join(BLAST_DB_PATH, "uniprot_db")
    if not os.path.exists(db_path + ".pin"):
        result = subprocess.run([
            "makeblastdb", "-in", UNIPROT_FASTA_PATH, "-dbtype", "prot", "-out", db_path
        ], capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error while creating BLAST database:\n{result.stderr}")
        else:
            print("BLAST database successfully created.")
    else:
        print("BLAST database already exists.")


def extract_sequence_from_parsed_data(chain_data):
    """
    Converts residue names from 3-letter code to a protein sequence.

    Args:
        chain_data (dict): Chain with residue information.

    Returns:
        str: One-letter sequence string or None.
    """
    sequence = "".join([three_to_one.get(res["residue_name"].capitalize(), "X") for res in chain_data["residues"]])
    return sequence if sequence else None


def run_blast_search(query_sequence):
    """
    Performs a BLAST search for a given protein sequence and returns the best UniProt match.

    Args:
        query_sequence (str): Protein sequence to search.

    Returns:
        str or None: Matched UniProt ID or None if no match found.
    """
    unique_id = str(uuid.uuid4())[:8]
    query_file = f"query_{unique_id}.fasta"
    output_file = f"blast_results_{unique_id}.txt"

    try:
        with open(query_file, "w") as f:
            f.write(f">query\n{query_sequence}\n")

        blast_cmd = [
            BLAST_EXECUTABLE,
            "-query", query_file,
            "-db", os.path.join(BLAST_DB_PATH, "uniprot_db"),
            "-out", output_file,
            "-evalue", "1e-5",
            "-max_target_seqs", "8",
            "-outfmt", "6"
        ]

        result = subprocess.run(blast_cmd, capture_output=True, text=True)

        if result.returncode != 0:
            return None

        if not os.path.exists(output_file):
            return None

        best_match = None
        with open(output_file, "r") as f:
            for line in f:
                best_match = line.split("\t")[1].split("|")[1]
                break

        return best_match

    except Exception:
        return None

    finally:
        if os.path.exists(query_file):
            os.remove(query_file)
        if os.path.exists(output_file):
            os.remove(output_file)


def classify_molecule_type(chain):
    """
    Classifies the molecule type of a chain based on its residues.

    Args:
        chain (dict): Chain with residues.

    Returns:
        str: "Protein", "Nucleic Acid", or "Unknown".
    """
    residues = chain["residues"]
    if len(residues) == 0:
        return "Unknown"

    sample_names = {res["residue_name"].strip() for res in residues}
    return "Protein" if all(res in AMINO_ACIDS for res in sample_names) else "Nucleic Acid"


def parallel_blast_search(parsed_data, max_workers=8):
    """
    Performs parallelized BLAST searches for all eligible chains in the dataset.

    Args:
        parsed_data (list): List of structure dictionaries.
        max_workers (int): Number of parallel threads.
    """
    tasks = []

    # Prepare list of chains and sequences to BLAST
    for structure in parsed_data:
        for chain in structure["atom_data"]:
            if chain.get("uniprot_id") not in [None, "Unknown"] or chain.get("molecule_type") == "Nucleic Acid":
                continue

            sequence = extract_sequence_from_parsed_data(chain)
            if not sequence:
                continue

            sequence = sequence.replace("O", "X")
            tasks.append((len(sequence), sequence, chain))

    # Sort tasks by sequence length if many tasks
    if len(tasks) > 40:
        tasks.sort(key=lambda t: t[0])

    # Perform parallel BLAST searches
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [(executor.submit(run_blast_search, seq), chain) for _, seq, chain in tasks]

        for future, chain in futures:
            uniprot_id = future.result()
            if uniprot_id:
                chain["uniprot_id"] = uniprot_id
                chain["molecule_name"] = f"Matched UniProt: {uniprot_id}"
                chain["molecule_type"] = "Protein"
            else:
                chain["molecule_type"] = classify_molecule_type(chain)
