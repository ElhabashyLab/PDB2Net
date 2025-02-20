import os
import csv
import re
from Bio import PDB
from config_loader import config

# Define allowed file extensions
ALLOWED_EXTENSIONS = {'.pdb', '.cif', '.mmcif'}


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
        str or None: The extracted PDB ID (uppercase) if found, otherwise None.
    """
    filename = os.path.basename(file_path)
    match = re.match(r"([0-9][A-Za-z0-9]{3})", filename)
    return match.group(1).upper() if match else None


def extract_pdb_id_from_file(file_path):
    """
    Extracts the PDB ID from the file content if it is not present in the filename.

    Args:
        file_path (str): Path to the file.

    Returns:
        str or None: The extracted PDB ID (uppercase) if found, otherwise None.
    """
    _, ext = os.path.splitext(file_path)
    ext = ext.lower()

    try:
        with open(file_path, "r") as f:
            for line in f:
                # PDB format: PDB ID appears in the HEADER line
                if ext == ".pdb" and line.startswith("HEADER"):
                    return line.split()[-1].strip().upper()

                # mmCIF format: PDB ID appears in 'data_' or '_entry.id' lines
                if ext in {".cif", ".mmcif"}:
                    if line.lower().startswith("data_"):
                        return line.split("_", 1)[1].strip().upper()
                    if "_entry.id" in line:
                        return line.split()[-1].strip().upper()

    except Exception as e:
        print(f"⚠ Error reading PDB ID from {file_path}: {e}")

    return None


def get_pdb_id(file_path):
    """
    Determines the PDB ID, preferring extraction from the filename.
    If not found in the filename, attempts extraction from the file content.

    Args:
        file_path (str): Path to the file.

    Returns:
        str: The determined PDB ID, or None if it could not be extracted.
    """
    pdb_id = extract_pdb_id_from_filename(file_path)
    if pdb_id:
        return pdb_id

    pdb_id = extract_pdb_id_from_file(file_path)
    if pdb_id:
        return pdb_id

    print(f"⚠ WARNING: No PDB ID found for file: {file_path}")
    return None


def parse_structure(file_path, pdb_id):
    """
    Parses a PDB or mmCIF file and returns the structure.

    Args:
        file_path (str): Path to the PDB or mmCIF file.
        pdb_id (str): The PDB ID associated with the structure.

    Returns:
        Bio.PDB.Structure.Structure or None: Parsed structure object, or None on failure.
    """
    _, ext = os.path.splitext(file_path)

    # Select the appropriate parser based on the file type
    parser = PDB.MMCIFParser(QUIET=True) if ext.lower() in {'.cif', '.mmcif'} else PDB.PDBParser(QUIET=True)

    try:
        structure = parser.get_structure(pdb_id, file_path)
        return structure
    except Exception as e:
        print(f"⚠ Error parsing {file_path}: {e}")
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
            pdb_id = get_pdb_id(file_path) or "UNKNOWN"
            print(f"📂 Processing file: {file_path} → PDB ID: {pdb_id}")

            structure = parse_structure(file_path, pdb_id)
            if structure:
                structures.append({"file_path": file_path, "pdb_id": pdb_id, "structure": structure})
        else:
            print(f"⚠ Skipping file (invalid extension): {file_path}")

    return structures
