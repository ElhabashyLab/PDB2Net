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
    """Extrahiert die PDB-ID aus dem Dateinamen, falls vorhanden."""
    filename = os.path.basename(file_path)
    match = re.match(r"([0-9][A-Za-z0-9]{3})", filename)
    return match.group(1).upper() if match else None


def extract_pdb_id_from_file(file_path):
    """Extrahiert die PDB-ID aus dem Dateiinhalt, falls nicht im Dateinamen enthalten."""
    _, ext = os.path.splitext(file_path)
    ext = ext.lower()

    try:
        with open(file_path, "r") as f:
            for line in f:
                if ext == ".pdb" and line.startswith("HEADER"):
                    print(f"🔎 HEADER gefunden in {file_path}: {line.strip()}")  # Debugging
                    return line.split()[-1].strip().upper()

                if ext in {".cif", ".mmcif"}:
                    if line.lower().startswith("data_"):
                        print(f"🔎 data_ Zeile gefunden in {file_path}: {line.strip()}")  # Debugging
                        return line.split("_", 1)[1].strip().upper()
                    if "_entry.id" in line:
                        print(f"🔎 _entry.id gefunden in {file_path}: {line.strip()}")  # Debugging
                        return line.split()[-1].strip().upper()

    except Exception as e:
        print(f"⚠ Fehler beim Lesen der PDB-ID aus {file_path}: {e}")

    return None


def get_pdb_id(file_path):
    """Bestimmt die PDB-ID nur einmal beim Einlesen der Datei."""
    pdb_id = extract_pdb_id_from_filename(file_path)
    if pdb_id:
        print(f"✅ PDB-ID aus Dateiname extrahiert: {pdb_id}")  # Debugging
        return pdb_id

    pdb_id = extract_pdb_id_from_file(file_path)
    if pdb_id:
        print(f"✅ PDB-ID aus Dateiinhalt extrahiert: {pdb_id}")  # Debugging
        return pdb_id

    print(f"⚠ WARNUNG: Keine PDB-ID gefunden für Datei: {file_path}")  # Debugging
    return None


def parse_structure(file_path, pdb_id):
    """Parst PDB- oder mmCIF-Dateien und gibt die Struktur zurück."""
    _, ext = os.path.splitext(file_path)
    parser = PDB.MMCIFParser(QUIET=True) if ext.lower() in {'.cif', '.mmcif'} else PDB.PDBParser(QUIET=True)

    try:
        structure = parser.get_structure(pdb_id, file_path)
        return structure
    except Exception as e:
        print(f"⚠ Fehler beim Parsen von {file_path}: {e}")
        return None


def read_files_from_csv(csv_path):
    """Liest Pfade aus der CSV-Datei, bestimmt die PDB-ID einmal und parst gültige Dateien."""
    file_paths = []
    with open(csv_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        if "file_path" not in reader.fieldnames:
            raise ValueError("Die CSV-Datei muss eine Spalte 'file_path' enthalten.")
        for row in reader:
            file_paths.append(row["file_path"])

    structures = []
    for file_path in file_paths:
        if is_valid_file(file_path):
            pdb_id = get_pdb_id(file_path) or "UNKNOWN"
            print(f"📂 Datei: {file_path} → PDB-ID: {pdb_id}")  # Debugging

            structure = parse_structure(file_path, pdb_id)
            if structure:
                structures.append({"file_path": file_path, "pdb_id": pdb_id, "structure": structure})
        else:
            print(f"⚠ Datei übersprungen (ungültige Erweiterung): {file_path}")  # Debugging

    return structures
