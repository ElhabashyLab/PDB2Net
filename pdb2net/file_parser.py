import os
import csv
from Bio import PDB

ALLOWED_EXTENSIONS = {'.pdb', '.cif', '.mmcif'}


def is_valid_file(file_path):
    """
    Prüft, ob die Datei eine gültige Erweiterung hat.
    """
    _, ext = os.path.splitext(file_path)
    return ext.lower() in ALLOWED_EXTENSIONS


def parse_structure(file_path):
    """
    Parst PDB- oder mmCIF-Dateien und gibt die Struktur zurück.
    """
    _, ext = os.path.splitext(file_path)
    parser = PDB.MMCIFParser(QUIET=True) if ext.lower() in {'.cif', '.mmcif'} else PDB.PDBParser(QUIET=True)
    structure_id = os.path.basename(file_path).split('.')[0]

    try:
        structure = parser.get_structure(structure_id, file_path)
        return structure
    except Exception as e:
        print(f"Fehler beim Parsen von {file_path}: {e}")
        return None


def read_csv(file_path):
    """
    Liest eine CSV-Datei und extrahiert Dateipfade.
    """
    paths = []
    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        if "file_path" not in reader.fieldnames:
            raise ValueError("Die CSV-Datei muss eine Spalte 'file_path' enthalten.")
        for row in reader:
            paths.append(row["file_path"])
    return paths


def read_files_from_csv(csv_path):
    """
    Liest Pfade aus der CSV-Datei und parst gültige Dateien.
    """
    file_paths = read_csv(csv_path)
    structures = []
    for file_path in file_paths:
        if is_valid_file(file_path):
            structure = parse_structure(file_path)
            if structure:
                structures.append({"file_path": file_path, "structure": structure})
        else:
            print(f"Skipped invalid file: {file_path}")
    return structures
