import os
import subprocess
from Bio import SeqIO
from Bio.Data import IUPACData
from concurrent.futures import ThreadPoolExecutor
from config_loader import config
import uuid

# Lade Pfade aus der Konfiguration
BLAST_DB_PATH = config["blast_db_path"]
BLAST_EXECUTABLE = config["blastp_executable"]
UNIPROT_FASTA_PATH = config["uniprot_fasta_path"]

# Aminosäure-Mapping
three_to_one = IUPACData.protein_letters_3to1

AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "SEC", "PYL"
}

def create_blast_database():
    """
    Erstellt eine BLAST-Datenbank, falls sie nicht existiert.
    """
    db_path = os.path.join(BLAST_DB_PATH, "uniprot_db")
    if not os.path.exists(db_path + ".pin"):
        print("Erstelle BLAST-Datenbank...")
        result = subprocess.run([
            "makeblastdb", "-in", UNIPROT_FASTA_PATH, "-dbtype", "prot", "-out", db_path
        ], capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Fehler beim Erstellen der BLAST-Datenbank:\n{result.stderr}")
        else:
            print("BLAST-Datenbank erfolgreich erstellt.")
    else:
        print("BLAST-Datenbank existiert bereits.")

def extract_sequence_from_parsed_data(chain_data):
    """
    Konvertiert Residuen-Namen (Dreiletter-Code) in eine Sequenz.
    """
    sequence = "".join([three_to_one.get(res["residue_name"].capitalize(), "X") for res in chain_data["residues"]])
    return sequence if sequence else None

def run_blast_search(query_sequence):
    """
    Führt eine BLAST-Suche für eine gegebene Proteinsequenz aus und gibt die beste UniProt-ID zurück.
    Nutzt individuelle temporäre Dateien und prüft nur im Fehlerfall, ob die Datenbank existiert.
    """
    unique_id = str(uuid.uuid4())[:8]  # Eindeutiger Bezeichner für temporäre Dateien
    query_file = f"query_{unique_id}.fasta"
    output_file = f"blast_results_{unique_id}.txt"

    try:
        # Schreibe die Sequenz in eine temporäre FASTA-Datei
        with open(query_file, "w") as f:
            f.write(f">query\n{query_sequence}\n")

        # BLAST-Befehl zusammenbauen
        blast_cmd = [
            BLAST_EXECUTABLE,
            "-query", query_file,
            "-db", os.path.join(BLAST_DB_PATH, "uniprot_db"),
            "-out", output_file,
            "-evalue", "1e-5",
            "-max_target_seqs", "8",
            "-outfmt", "6"
        ]

        print(f"Starte BLAST: {' '.join(blast_cmd)}")
        result = subprocess.run(blast_cmd, capture_output=True, text=True)

        # Fehlerbehandlung
        if result.returncode != 0:
            stderr = result.stderr.strip()
            if "No alias or index file found" in stderr:
                print("❌ Fehler: Die BLAST-Datenbank wurde nicht gefunden oder ist beschädigt!")
                print("→ Erwarteter Pfad:", os.path.join(BLAST_DB_PATH, "uniprot_db"))
            else:
                print(f"❌ Allgemeiner BLAST-Fehler:\n{stderr}")
            return None

        if not os.path.exists(output_file):
            print(f"❌ Fehler: Die BLAST-Ausgabedatei {output_file} wurde nicht erstellt!")
            return None

        # UniProt-ID aus erster Trefferzeile extrahieren
        best_match = None
        with open(output_file, "r") as f:
            for line in f:
                best_match = line.split("\t")[1].split("|")[1]
                break

        return best_match

    except Exception as e:
        print(f"❌ Unerwarteter Fehler bei BLAST-Suche: {e}")
        return None

    finally:
        # Temporäre Dateien aufräumen
        if os.path.exists(query_file):
            os.remove(query_file)
        if os.path.exists(output_file):
            os.remove(output_file)

def classify_molecule_type(chain):
    """
    Bestimmt, ob eine Kette ein Protein oder Nukleinsäure ist.
    """
    residues = chain["residues"]
    if len(residues) == 0:
        return "Unknown"

    sample_names = {res["residue_name"].strip() for res in residues}
    return "Protein" if all(res in AMINO_ACIDS for res in sample_names) else "Nucleic Acid"

def parallel_blast_search(parsed_data, max_workers=8):
    """
    Parallelisiert die BLAST-Suche für alle Ketten.
    Sortiert optional nach Sequenzlänge, wenn mehr als 30 Anfragen anstehen.
    """
    tasks = []

    # 1. Vorbereitung: extrahiere relevante Ketten und Sequenzen
    for structure in parsed_data:
        for chain in structure["atom_data"]:
            if chain.get("uniprot_id") not in [None, "Unknown"] or chain.get("molecule_type") == "Nucleic Acid":
                continue

            sequence = extract_sequence_from_parsed_data(chain)
            if not sequence:
                continue

            sequence = sequence.replace("O", "X")
            tasks.append((len(sequence), sequence, chain))

    print(f"🔍 Es stehen {len(tasks)} BLAST-Anfragen an.")

    # 2. Sortieren bei Bedarf (wenn viele und unterschiedlich lange)
    if len(tasks) > 40:
        tasks.sort(key=lambda t: t[0])  # sortiere nach Sequenzlänge (aufsteigend)
        print("⚙️ Sortiere Sequenzen nach Länge für effizientere Parallelisierung...")

    # 3. Parallel ausführen
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [(executor.submit(run_blast_search, seq), chain) for _, seq, chain in tasks]

        for future, chain in futures:
            uniprot_id = future.result()
            if uniprot_id:
                chain["uniprot_id"] = uniprot_id
                chain["molecule_name"] = f"Matched UniProt: {uniprot_id}"
                chain["molecule_type"] = "Protein"
                print(f"✅ Gefunden: {chain['unique_chain_id']} → {uniprot_id}")
            else:
                chain["molecule_type"] = classify_molecule_type(chain)
                print(f"❌ Kein UniProt-Match für {chain['unique_chain_id']}, klassifiziert als {chain['molecule_type']}")

    print("✅ Parallelisierte UniProt BLAST-Suche abgeschlossen.")

