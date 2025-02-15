import os
import csv
import re
from Bio import PDB

ALLOWED_EXTENSIONS = {'.pdb', '.cif', '.mmcif'}


def is_valid_file(file_path):
    """Prüft, ob die Datei eine gültige Erweiterung hat."""
    _, ext = os.path.splitext(file_path)
    return ext.lower() in ALLOWED_EXTENSIONS


def extract_pdb_id_from_filename(file_path):
    """
    Extrahiert die PDB-ID aus dem Dateinamen, falls vorhanden.
    Erwartet eine 4-stellige alphanumerische ID.
    """
    filename = os.path.basename(file_path)
    match = re.match(r"([0-9][A-Za-z0-9]{3})", filename)  # PDB-IDs sind genau 4 Zeichen lang
    return match.group(1).upper() if match else None  # Großbuchstaben zurückgeben


def extract_pdb_id_from_file(file_path):
    """
    Versucht die PDB-ID direkt aus der Datei zu extrahieren, falls sie nicht im Dateinamen steht.
    - Für PDB-Dateien wird die ID aus der HEADER-Zeile gelesen.
    - Für mmCIF-Dateien wird die ID aus der 'data_' Zeile oder '_entry.id' gelesen.
    """
    _, ext = os.path.splitext(file_path)
    ext = ext.lower()

    try:
        with open(file_path, "r") as f:
            for line in f:
                if ext == ".pdb" and line.startswith("HEADER"):
                    return line.split()[-1].strip().upper()  # Letztes Element der Zeile ist die PDB-ID

                if ext in {".cif", ".mmcif"}:
                    if line.lower().startswith("data_"):
                        return line.split("_", 1)[1].strip().upper()  # 'data_1A4Y' → '1A4Y'
                    if "_entry.id" in line:
                        return line.split()[-1].strip().upper()  # Letztes Element der Zeile

    except Exception as e:
        print(f"⚠ Fehler beim Lesen der PDB-ID aus {file_path}: {e}")

    return None  # Falls keine ID gefunden wird


def get_pdb_id(file_path):
    """
    Bestimmt die PDB-ID, indem zuerst der Dateiname und dann der Dateiinhalt überprüft wird.
    """
    pdb_id = extract_pdb_id_from_filename(file_path)
    if pdb_id:
        return pdb_id

    pdb_id = extract_pdb_id_from_file(file_path)
    if pdb_id:
        return pdb_id

    print(f"⚠ Keine PDB-ID gefunden für Datei: {file_path}")
    return None


def parse_structure(file_path):
    """Parst PDB- oder mmCIF-Dateien und gibt die Struktur zurück."""
    _, ext = os.path.splitext(file_path)
    parser = PDB.MMCIFParser(QUIET=True) if ext.lower() in {'.cif', '.mmcif'} else PDB.PDBParser(QUIET=True)
    pdb_id = get_pdb_id(file_path) or "UNKNOWN"

    try:
        structure = parser.get_structure(pdb_id, file_path)
        return structure
    except Exception as e:
        print(f"Fehler beim Parsen von {file_path}: {e}")
        return None


def read_csv(file_path):
    """Liest eine CSV-Datei und extrahiert Dateipfade."""
    paths = []
    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        if "file_path" not in reader.fieldnames:
            raise ValueError("Die CSV-Datei muss eine Spalte 'file_path' enthalten.")
        for row in reader:
            paths.append(row["file_path"])
    return paths


def read_files_from_csv(csv_path):
    """Liest Pfade aus der CSV-Datei und parst gültige Dateien."""
    file_paths = read_csv(csv_path)
    structures = []
    for file_path in file_paths:
        if is_valid_file(file_path):
            structure = parse_structure(file_path)
            if structure:
                pdb_id = get_pdb_id(file_path) or "UNKNOWN"
                structures.append({"file_path": file_path, "pdb_id": pdb_id, "structure": structure})
        else:
            print(f"Skipped invalid file: {file_path}")
    return structures