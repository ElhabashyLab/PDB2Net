import os
import subprocess
from Bio import SeqIO
from Bio.Data import IUPACData
from config_loader import config

# Lade Pfade aus Config
UNIPROT_FASTA_PATH = config["uniprot_fasta_path"]
BLAST_DB_PATH = config["blast_db_path"]
BLAST_EXECUTABLE = config["blastp_executable"]

# Biopython Mapping für 3-Buchstaben-Code → 1-Buchstaben-Code
three_to_one = IUPACData.protein_letters_3to1

# Liste der 20 kanonischen Aminosäuren + Selenocystein & Pyrrolysin
AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "SEC", "PYL"
}


def create_blast_database():
    """Erstellt eine BLAST-Datenbank aus der UniProt FASTA-Datei, falls nicht vorhanden."""
    db_path = os.path.join(BLAST_DB_PATH, "uniprot_db")
    if not os.path.exists(db_path + ".pin"):  # Prüfe, ob BLAST-DB existiert
        print("🛠 Erstelle BLAST-Datenbank...")
        subprocess.run([
            "makeblastdb", "-in", UNIPROT_FASTA_PATH, "-dbtype", "prot", "-out", db_path
        ], check=True)
        print("✅ BLAST-Datenbank erstellt!")
    else:
        print("✅ BLAST-Datenbank bereits vorhanden.")


def extract_sequence_from_parsed_data(chain_data):
    """Konvertiert Residue-Namen (3-Buchstaben-Code) in eine Aminosäuresequenz (1-Buchstaben-Code)."""
    sequence = "".join([three_to_one.get(res["residue_name"].capitalize(), "X") for res in chain_data["residues"]])
    return sequence if sequence else None


def run_blast_search(query_sequence):
    """Führt eine BLAST-Suche mit der gegebenen Proteinsequenz durch."""
    query_file = "query.fasta"
    output_file = "blast_results.txt"

    # Speichere die Sequenz in einer FASTA-Datei
    with open(query_file, "w") as f:
        f.write(f">query\n{query_sequence}\n")

    # Führe BLAST-Suche durch
    subprocess.run([
        BLAST_EXECUTABLE, "-query", query_file, "-db", os.path.join(BLAST_DB_PATH, "uniprot_db"),
        "-out", output_file, "-evalue", "1e-5", "-max_target_seqs", "8", "-outfmt", "6"
    ], check=True)

    # Lese das erste Treffer-Ergebnis
    best_match = None
    with open(output_file, "r") as f:
        for line in f:
            best_match = line.split("\t")[1]  # Zweite Spalte enthält die UniProt-ID
            break

    # Lösche temporäre Dateien
    os.remove(query_file)
    os.remove(output_file)

    return best_match if best_match else None


def classify_molecule_type(chain):
    """Bestimmt, ob eine Kette ein Protein ist, basierend auf ihren Residuen."""
    residues = chain["residues"]
    if len(residues) == 0:
        return "Unknown"

    sample_names = {res["residue_name"].strip() for res in residues}

    # Prüft, ob alle Residuen in der bekannten Aminosäure-Liste sind
    return "Protein" if all(res in AMINO_ACIDS for res in sample_names) else "Unknown"


def match_sequence_to_uniprot(parsed_data):
    """Vergleicht extrahierte Sequenzen mit UniProt via BLAST und klassifiziert Moleküle."""
    print("🔍 Starte UniProt BLAST-Suche...")
    create_blast_database()

    for structure in parsed_data:
        if len(structure["pdb_id"]) == 4:  # ❌ PDB-Dateien werden übersprungen
            continue

        for chain in structure["atom_data"]:
            if chain["is_hetatm"]:
                continue  # HETATM-Ketten ignorieren

            sequence = extract_sequence_from_parsed_data(chain)  # ✅ Umwandlung in 1-Buchstaben-Code
            if not sequence:
                print(f"⚠ Keine Sequenz extrahiert für {chain['unique_chain_id']}")
                continue

            sequence = sequence.replace("O", "X")  # ❌ Fehlerhafte "O"-Zeichen korrigieren

            uniprot_id = run_blast_search(sequence)
            if uniprot_id:
                chain["uniprot_id"] = uniprot_id
                chain["molecule_name"] = f"Matched UniProt: {uniprot_id}"
                chain["molecule_type"] = "Protein"
                print(f"✅ {chain['unique_chain_id']} → {uniprot_id}")
            else:
                chain["molecule_type"] = classify_molecule_type(chain)  # ✅ Alte Klassifikation wiederhergestellt!
                print(
                    f"⚠ Kein UniProt-Match für {chain['unique_chain_id']}, klassifiziert als {chain['molecule_type']}")

    print("🔍 UniProt BLAST-Suche abgeschlossen.")
