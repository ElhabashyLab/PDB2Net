import os
import csv

# Supported file extensions
ALLOWED_EXTENSIONS = {'.pdb', '.cif', '.mmcif'}

def is_valid_file(file_path):
    """
    Checks if a file has a valid extension.
    """
    _, ext = os.path.splitext(file_path)
    return ext.lower() in ALLOWED_EXTENSIONS

def read_file(file_path):
    """
    Reads an individual file if it has a valid extension and returns its content.
    """
    if is_valid_file(file_path):
        print(f"Reading file: {file_path}")
        # Replace with actual file reading logic
        with open(file_path, 'r') as f:
            content = f.read()
        return content
    else:
        raise ValueError(f"Invalid file type: {file_path}")

def read_files_from_csv(csv_path):
    """
    Reads a CSV file containing paths to files and returns a list of valid file contents.
    """
    contents = []
    with open(csv_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            file_path = row[0]  # Assumes each row has a single file path in the first column
            if is_valid_file(file_path):
                contents.append(read_file(file_path))
            else:
                print(f"Skipped invalid file: {file_path}")
    return contents
