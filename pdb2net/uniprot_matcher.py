import os
from Bio import SeqIO
from Bio.Data import IUPACData
from Bio import pairwise2
from config_loader import config
import random

# Lade UniProt FASTA-Datei
UNIPROT_FASTA_PATH = config["uniprot_fasta_path"]

# Liste bekannter Aminosäuren und Nukleotide für die Klassifikation
AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "SEC", "PYL"
}

NUCLEOTIDES = {"A", "T", "C", "G", "U", "DA", "DT", "DC", "DG", "DU"}  # RNA und DNA


def load_uniprot_sequences(fasta_path):
    """Lädt alle UniProt-Sequenzen aus der FASTA-Datei."""
    uniprot_sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        uniprot_id = record.id.split("|")[1]  # UniProt ID aus Header extrahieren
        uniprot_sequences[uniprot_id] = str(record.seq)
    return uniprot_sequences


def extract_sequence_from_parsed_data(chain_data):
    """Konvertiert Residue-Namen (3-Buchstaben-Code) in eine Aminosäuresequenz (1-Buchstaben-Code)."""
    three_to_one = IUPACData.protein_letters_3to1
    sequence = "".join([three_to_one.get(res["residue_name"].capitalize(), "X") for res in chain_data["residues"]])
    return sequence if sequence else None


def classify_molecule_type(chain):
    """Bestimmt, ob eine Kette ein Protein oder eine Nukleinsäure ist, basierend auf Stichproben."""
    residues = chain["residues"]
    if len(residues) == 0:
        return "Unknown"

    # Ziehe zufällig 3 Residuen als Stichprobe oder nimm die ersten, falls weniger als 3 vorhanden sind
    sample_residues = random.sample(residues, min(3, len(residues)))
    sample_names = {res["residue_name"].strip() for res in sample_residues}

    # Falls alle Residuen Aminosäuren sind, ist es ein Protein
    if all(res in AMINO_ACIDS for res in sample_names):
        return "Protein"

    # Falls alle Residuen Nukleotide sind, ist es eine Nukleinsäure
    if all(res in NUCLEOTIDES for res in sample_names):
        return "Nucleic Acid"

    return "Unknown"  # Falls gemischt oder unklar


def match_sequence_to_uniprot(parsed_data):
    """Vergleicht extrahierte Sequenzen mit UniProt, um eine ID zu finden und klassifiziert unbekannte Moleküle."""
    print("🔍 Starte UniProt-Abgleich für Moleküle ohne PDB-ID...")
    uniprot_db = load_uniprot_sequences(UNIPROT_FASTA_PATH)

    for structure in parsed_data:
        if len(structure["pdb_id"]) == 4:
            continue  # Nur für Dateien OHNE PDB-ID!

        for chain in structure["atom_data"]:
            if chain["is_hetatm"]:
                continue  # HETATM-Ketten ignorieren

            extracted_seq = extract_sequence_from_parsed_data(chain)
            if not extracted_seq:
                print(f"⚠ Keine Sequenz extrahiert für {chain['unique_chain_id']}")
                continue

            match_found = False
            for uniprot_id, uniprot_seq in uniprot_db.items():
                alignments = pairwise2.align.globalxx(extracted_seq, uniprot_seq, score_only=True)
                if alignments / len(extracted_seq) > 0.8:  # Mindestens 80% Übereinstimmung
                    chain["uniprot_id"] = uniprot_id
                    chain["molecule_name"] = f"Matched UniProt: {uniprot_id}"
                    chain["molecule_type"] = "Protein"
                    match_found = True
                    break

            if not match_found:
                print(f"⚠ Kein UniProt-Match für {chain['unique_chain_id']} gefunden!")

    print("🔍 UniProt-Abgleich abgeschlossen.")

    # **Neu: Klassifikation unbekannter Moleküle nach der UniProt-Suche**
    print("\n🔎 Klassifiziere verbleibende 'Unknown' Molekültypen...")
    for structure in parsed_data:
        for chain in structure["atom_data"]:
            if chain["molecule_type"] == "Unknown":
                chain["molecule_type"] = classify_molecule_type(chain)
                print(f"  🏷 {chain['unique_chain_id']} klassifiziert als: {chain['molecule_type']}")

    print("✅ Klassifikation abgeschlossen.")
