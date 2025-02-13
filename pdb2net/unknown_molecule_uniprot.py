import os
from Bio import SeqIO

def load_pdb_fasta(pdb_fasta_path):
    """
    Lädt die PDB FASTA-Datei und speichert die Sequenzen in einem Dictionary.
    Die PDB-IDs werden in Kleinbuchstaben umgewandelt, die Ketten-ID bleibt groß.
    """
    pdb_sequences = {}
    with open(pdb_fasta_path, "r") as f:
        current_key = None
        current_seq = []
        for line in f:
            if line.startswith(">"):
                # Speichere vorherige Sequenz
                if current_key and current_seq:
                    pdb_sequences[current_key]["sequence"] = "".join(current_seq)

                # Extrahiere PDB-ID und Kettenbezeichnung
                parts = line.split()
                fasta_header = parts[0][1:]  # Entferne ">"

                formatted_key = None  # Sicherstellen, dass die Variable existiert
                if "_" in fasta_header:  # Nur wenn "_" vorhanden ist
                    pdb_id_chain = fasta_header.split("_")
                    if len(pdb_id_chain) == 2:
                        pdb_id = pdb_id_chain[0].lower()  # PDB-ID in Kleinbuchstaben
                        chain_id = pdb_id_chain[1].upper()  # Ketten-ID in Großbuchstaben
                        formatted_key = f"{pdb_id}_{chain_id}"
                        pdb_sequences[formatted_key] = {"info": " ".join(parts[1:]), "sequence": ""}

                # Falls kein gültiger Key generiert wurde, Debug-Info ausgeben
                if formatted_key is None:
                    print(f"⚠ Warnung: Ungültiger FASTA-Header erkannt: {fasta_header}")

                # Setze neue aktuelle Sequenz nur, wenn formatted_key gültig ist
                if formatted_key:
                    current_key = formatted_key
                    current_seq = []

            else:
                current_seq.append(line.strip())

        # Speichere letzte Sequenz
        if current_key and current_seq:
            pdb_sequences[current_key]["sequence"] = "".join(current_seq)

    return pdb_sequences


def process_molecule_info(combined_data, pdb_fasta_path):
    """
    Weist jeder Kette einen Namen, Typ und eine Sequenz aus der PDB FASTA-Datei zu.
    """
    print("\n🔍 Bestimme Namen, Typen und Sequenzen aus PDB FASTA...")

    # Lade die PDB FASTA Datei
    pdb_fasta = load_pdb_fasta(pdb_fasta_path)

    # Debugging: Zeige die ersten 10 Keys aus der PDB FASTA
    print(f"✅ Erste 10 Keys in PDB FASTA:\n {list(pdb_fasta.keys())[:10]}")

    # Durchlaufe alle Strukturen
    for structure_data in combined_data:
        pdb_id = os.path.basename(structure_data["file_path"]).split(".")[0].lower()  # PDB-ID in Kleinbuchstaben

        for chain in structure_data["atom_data"]:
            chain_id = chain["chain_id"].upper()  # Ketten-ID in Großbuchstaben
            search_key = f"{pdb_id}_{chain_id}"  # Richtige Formatierung für den Key

            if search_key in pdb_fasta:
                # Hole Daten aus der PDB FASTA Datei
                chain_info = pdb_fasta[search_key]
                fasta_info = chain_info["info"]
                sequence = chain_info["sequence"]

                # Extrahiere Molekülname und Typ
                if "mol:protein" in fasta_info:
                    molecule_type = "Protein"
                elif "mol:na" in fasta_info:
                    if "DNA" in fasta_info:
                        molecule_type = "DNA"
                    elif "RNA" in fasta_info:
                        molecule_type = "RNA"
                    else:
                        molecule_type = "Nucleic Acid"
                else:
                    molecule_type = "Unknown"

                molecule_name = fasta_info.split("length:")[-1].strip()  # Alles nach "length:" als Name nehmen

                # Speichere Informationen in der Datenstruktur
                chain["molecule_name"] = molecule_name
                chain["molecule_type"] = molecule_type
                chain["sequence"] = sequence

                print(f"✅ Treffer für {search_key}: {molecule_name} ({molecule_type})")
            else:
                print(f"⚠ Kein Treffer für: {search_key} (PDB-ID: {pdb_id}, Kette: {chain_id})")


def load_uniprot_fasta(uniprot_fasta_path):
    """
    Lädt die UniProt FASTA-Datei und speichert die Sequenzen mit zugehörigen Namen.
    """
    uniprot_sequences = {}
    for record in SeqIO.parse(uniprot_fasta_path, "fasta"):
        sequence = str(record.seq)
        protein_name = record.description.split(" ", 1)[1]  # Alles nach dem ersten Leerzeichen als Name nehmen
        uniprot_sequences[sequence] = protein_name
    return uniprot_sequences


def find_best_uniprot_match_partial(pdb_sequence, uniprot_fasta):
    """
    Prüft, ob eine PDB FASTA-Sequenz als Teilstring in einer UniProt-Sequenz enthalten ist.
    Gibt den ersten gefundenen Treffer zurück.
    """
    for uniprot_seq, uniprot_name in uniprot_fasta.items():
        if pdb_sequence in uniprot_seq:  # Prüft, ob PDB-Sequenz ein Teil der UniProt-Sequenz ist
            return uniprot_name  # Gibt den Namen des ersten Treffers zurück

    return None  # Kein Treffer gefunden


def update_names_with_uniprot(combined_data, uniprot_fasta_path):
    """
    Aktualisiert die Molekülnamen mit Einträgen aus der UniProt FASTA-Datei, falls eine exakte oder Teilstring-Übereinstimmung gefunden wird.
    """
    print("\n🔍 Überprüfe Namen mit UniProt FASTA...")

    # Lade die UniProt FASTA-Datei
    uniprot_fasta = load_uniprot_fasta(uniprot_fasta_path)

    # Debugging: Zeige die ersten 5 UniProt-Sequenzen
    print(f"✅ Erste 5 UniProt-Sequenzen:\n {list(uniprot_fasta.keys())[:5]}")

    # Durchlaufe alle Strukturen und Ketten
    for structure_data in combined_data:
        for chain in structure_data["atom_data"]:
            pdb_sequence = chain["sequence"]  # PDB-Sequenz der Kette

            # Exaktes Matching zuerst
            if pdb_sequence in uniprot_fasta:
                chain["molecule_name"] = uniprot_fasta[pdb_sequence]
                print(f"🔄 Name (exakt) aktualisiert für {structure_data['pdb_id']}:{chain['chain_id']}")
            else:
                # Teilstring-Suche als zweite Option
                best_match = find_best_uniprot_match_partial(pdb_sequence, uniprot_fasta)
                if best_match:
                    chain["molecule_name"] = best_match
                    print(f"🔄 Name (Teilstring) aktualisiert für {structure_data['pdb_id']}:{chain['chain_id']}")
                else:
                    print(f"⚠ Kein UniProt-Treffer für {structure_data['pdb_id']}:{chain['chain_id']} → Behalte PDB-Name")

    print("\n✅ Namensaktualisierung mit UniProt abgeschlossen.")
