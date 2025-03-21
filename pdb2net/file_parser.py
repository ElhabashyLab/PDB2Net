import os
import re
from Bio import PDB
from config_loader import config
import csv
import gemmi

# Define allowed file extensions
ALLOWED_EXTENSIONS = {'.pdb', '.cif', '.mmcif'}

# Load the path to pdb_seqres.txt from the configuration
PDB_FASTA_PATH = config["pdb_fasta_path"]

def load_valid_pdb_ids():
    """
    Loads all valid PDB IDs from pdb_seqres.txt.

    Returns:
        set: A set of all valid PDB IDs.
    """
    valid_pdb_ids = set()
    try:
        with open(PDB_FASTA_PATH, "r") as f:
            for line in f:
                if line.startswith(">"):
                    parts = line.split()[0][1:].split("_")  # Extract PDB ID from header
                    if len(parts[0]) == 4:
                        valid_pdb_ids.add(parts[0].upper())  # Store IDs in uppercase
    except Exception as e:
        print(f"Error loading pdb_seqres.txt: {e}")

    return valid_pdb_ids

# Load PDB IDs once at startup
VALID_PDB_IDS = load_valid_pdb_ids()

def is_valid_file(file_path):
    """
    Checks if the given file has a valid extension.

    Args:
        file_path (str): Path to the file.

    Returns:
        bool: True if the file has a valid extension, otherwise False.
    """
    _, ext = os.path.splitext(file_path)
    return ext.lower() in ALLOWED_EXTENSIONS

def extract_pdb_id_from_filename(file_path):
    """
    Extracts the PDB ID from the filename if it follows the standard format.

    Args:
        file_path (str): Path to the file.

    Returns:
        str or None: The extracted PDB ID (uppercase) if valid, otherwise None.
    """
    filename = os.path.basename(file_path)
    match = re.match(r"^([0-9][A-Za-z0-9]{3})$", filename.split(".")[0])  # PDB ID is exactly 4 characters
    if match:
        pdb_id = match.group(1).upper()
        if pdb_id in VALID_PDB_IDS:
            return pdb_id
        print(f"Warning: File {filename} contains an invalid PDB ID ({pdb_id} not in pdb_seqres.txt).")
    return None

def extract_pdb_id_from_file(file_path):
    """
    Extracts the PDB ID from the file content if not present in the filename.

    Args:
        file_path (str): Path to the file.

    Returns:
        str or None: The extracted PDB ID (uppercase) if valid, otherwise None.
    """
    _, ext = os.path.splitext(file_path)
    ext = ext.lower()

    try:
        with open(file_path, "r") as f:
            for line in f:
                if ext == ".pdb" and line.startswith("HEADER"):
                    pdb_id = line.split()[-1].strip().upper()
                elif ext in {".cif", ".mmcif"}:
                    if line.lower().startswith("data_"):
                        pdb_id = line.split("_", 1)[1].strip().upper()
                    elif "_entry.id" in line:
                        pdb_id = line.split()[-1].strip().upper()
                else:
                    continue

                if len(pdb_id) == 4 and pdb_id in VALID_PDB_IDS:
                    return pdb_id
                print(f"Warning: File {file_path} contains an invalid PDB ID ({pdb_id} not in pdb_seqres.txt).")
                return None

    except Exception as e:
        print(f"Error reading {file_path}: {e}")

    return None

def get_pdb_id(file_path):
    """
    Determines the PDB ID, preferring extraction from the filename.
    If no valid PDB ID can be extracted, the filename is used as a fallback.

    Args:
        file_path (str): Path to the PDB/mmCIF file.

    Returns:
        str: The determined PDB ID.
    """
    pdb_id = extract_pdb_id_from_filename(file_path)
    if pdb_id:
        return pdb_id

    pdb_id = extract_pdb_id_from_file(file_path)
    if pdb_id:
        return pdb_id

    # If no valid PDB ID is found, use the filename as a fallback
    filename = os.path.basename(file_path).split('.')[0].upper()
    print(f"Warning: No valid PDB ID found. Using filename as PDB ID: {filename}")
    return filename  # Use filename as fallback PDB ID

def parse_structure(file_path, pdb_id):
    """
    Lädt eine PDB/mmCIF-Datei mit gemmi.
    """
    try:
        structure = gemmi.read_structure(file_path)
        structure.setup_entities()
        return structure
    except Exception as e:
        print(f"Error parsing {file_path}: {e}")
        return None

def read_files_from_csv(csv_path):
    """
    Reads file paths from a CSV file, extracts PDB IDs, and parses valid structures.

    Args:
        csv_path (str): Path to the CSV file.

    Returns:
        list: A list of dictionaries containing file path, PDB ID, and parsed structure.
    """
    file_paths = []

    # Open and read the CSV file
    with open(csv_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        if "file_path" not in reader.fieldnames:
            raise ValueError("CSV file must contain a 'file_path' column.")

        # Collect file paths
        for row in reader:
            file_paths.append(row["file_path"])

    structures = []
    for file_path in file_paths:
        if is_valid_file(file_path):
            pdb_id = get_pdb_id(file_path)
            if pdb_id:

                structure = parse_structure(file_path, pdb_id)
                if structure:
                    structures.append({"file_path": file_path, "pdb_id": pdb_id, "structure": structure})
            else:
                print(f"Skipping file (no valid PDB ID found): {file_path}")
        else:
            print(f"Skipping file (invalid extension): {file_path}")

    return structures
