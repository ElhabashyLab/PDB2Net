import os
from Bio import SeqIO
from config_loader import config


# Lade UniProt FASTA-Datei
UNIPROT_FASTA_PATH = config["uniprot_fasta_path"]


def load_uniprot_sequences(fasta_path):
    """Lädt alle UniProt-Sequenzen aus der FASTA-Datei."""
    uniprot_sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        uniprot_id = record.id.split("|")[1]  # UniProt ID aus Header extrahieren
        uniprot_sequences[uniprot_id] = str(record.seq)
    return uniprot_sequences


def extract_sequence_from_parsed_data(chain_data):
    """Extrahiert die Aminosäuresequenz aus den geparsten Daten einer Kette."""
    sequence = "".join([res["residue_name"] for res in chain_data["residues"]])
    return sequence if sequence else None


def match_sequence_to_uniprot(parsed_data):
    """Vergleicht extrahierte Sequenzen mit UniProt, um eine ID zu finden."""
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
                continue  # Falls keine Sequenz extrahiert werden konnte

            for uniprot_id, uniprot_seq in uniprot_db.items():
                if extracted_seq in uniprot_seq:  # Teilstringsuche
                    chain["uniprot_id"] = uniprot_id
                    chain["molecule_name"] = f"Matched UniProt: {uniprot_id}"
                    chain["molecule_type"] = "Protein"
                    print(f"✅ Match gefunden: {chain['unique_chain_id']} → {uniprot_id}")
                    break  # Ersten Match behalten
    print("🔍 UniProt-Abgleich abgeschlossen.")