import requests
import time
from Bio.SeqUtils import seq1


def filter_unknown_molecules(combined_data):
    """
    Filtert Ketten aus der kombinierten Datenstruktur, deren Molekülname 'UNKNOWN' ist,
    und fügt die Datei-Information hinzu. Vermeidet doppelte Prüfungen.
    """
    unknown_chains = []
    processed = set()  # Set zum Nachverfolgen bereits geprüfter Datei-Ketten-Kombinationen

    for structure_data in combined_data:
        file_path = structure_data["file_path"]
        for chain in structure_data["atom_data"]:
            if chain["molecule_name"] == "UNKNOWN":
                chain_id = chain["chain_id"]
                unique_key = (file_path, chain_id)
                if unique_key not in processed:
                    sequence = extract_chain_sequence(chain["residues"])
                    print(f"DEBUG: Prüfe Datei {file_path}, Kette {chain_id}, Sequenz: {sequence}")
                    if sequence:
                        unknown_chains.append({
                            "file_path": file_path,
                            "chain_id": chain_id,
                            "sequence": sequence
                        })
                        processed.add(unique_key)  # Markiere als verarbeitet
    return unknown_chains


def extract_chain_sequence(residues):
    """
    Extrahiert die Sequenz einer Kette aus den Residuen.
    """
    sequence = ""
    for residue in residues:
        residue_name = residue["residue_name"]
        try:
            sequence += seq1(residue_name)
        except KeyError:
            pass  # Ignoriere Residuen, die keine Aminosäuren sind
    if not sequence:
        print("WARNUNG: Keine Sequenz extrahiert!")
    return sequence


def fetch_molecule_name_with_blast(sequence):
    """
    Verwendet die aktuelle UniProt BLAST-API, um den Molekülnamen basierend auf einer Sequenz zu finden.
    """
    # UniProt BLAST API-Endpunkt
    url = "https://rest.uniprot.org/blast/run"

    # Parameter für die BLAST-Suche
    headers = {"Content-Type": "application/x-www-form-urlencoded"}
    params = {
        "query": f">Query\n{sequence}",  # FASTA-formatierte Sequenz
        "program": "blastp",
        "database": "uniprotkb"
    }

    try:
        # Schritt 1: Job starten
        response = requests.post(url, headers=headers, data=params)
        response.raise_for_status()
        job_id = response.json().get("jobId")
        if not job_id:
            raise ValueError(f"Kein gültiger Job-ID erhalten für Sequenz: {sequence}")

        # Schritt 2: Ergebnisstatus abfragen
        status_url = f"https://rest.uniprot.org/blast/status/{job_id}"
        while True:
            status_response = requests.get(status_url)
            status_response.raise_for_status()
            status = status_response.json().get("jobStatus")
            if status == "FINISHED":
                break
            elif status == "FAILED":
                raise ValueError(f"BLAST-Job fehlgeschlagen für Sequenz: {sequence}")
            time.sleep(5)  # Wartezeit für Verarbeitung

        # Schritt 3: Ergebnisse abrufen
        result_url = f"https://rest.uniprot.org/blast/results/{job_id}?format=json"
        result_response = requests.get(result_url)
        result_response.raise_for_status()
        results = result_response.json()

        if "hits" in results and results["hits"]:
            # Nimm den ersten Treffer
            top_hit = results["hits"][0]
            protein_name = top_hit["target"]["name"]
            return protein_name

        return "UNKNOWN"  # Keine Treffer gefunden

    except requests.exceptions.RequestException as e:
        print(f"Netzwerkfehler: {e}")
    except ValueError as ve:
        print(f"BLAST-Fehler: {ve}")
    except Exception as e:
        print(f"Fehler: {e}")
    return "UNKNOWN"


def update_combined_data_with_uniprot(combined_data, unknown_chains):
    """
    Aktualisiert die kombinierte Datenstruktur mit den Namen aus UniProt und gibt Kette sowie zugehöriges File aus.
    """
    for unknown_chain in unknown_chains:
        chain_id = unknown_chain["chain_id"]
        sequence = unknown_chain["sequence"]
        file_path = unknown_chain["file_path"]  # File-Information
        print(f"Verarbeite Datei: {file_path}, Kette: {chain_id}, Sequenz: {sequence}")
        molecule_name = fetch_molecule_name_with_blast(sequence)
        print(f"Ergebnis für Datei: {file_path}, Kette: {chain_id}: {molecule_name}")
        # Aktualisiere die Datenstruktur
        for structure_data in combined_data:
            if structure_data["file_path"] == file_path:  # Passendes File finden
                for chain in structure_data["atom_data"]:
                    if chain["chain_id"] == chain_id and chain["molecule_name"] == "UNKNOWN":
                        chain["molecule_name"] = molecule_name


def process_unknown_molecules(combined_data):
    """
    Gesamtprozess zur Verarbeitung von 'UNKNOWN'-Molekülen.
    """
    unknown_chains = filter_unknown_molecules(combined_data)
    if unknown_chains:
        update_combined_data_with_uniprot(combined_data, unknown_chains)
    else:
        print("Alle Ketten haben bekannte Molekülnamen. Keine Aktion erforderlich.")


def list_unknown_molecules(combined_data):
    """
    Extrahiert alle Ketten mit 'UNKNOWN' als Molekülname aus der kombinierten Datenstruktur.
    Gibt eine Liste von Datei-Kette-Kombinationen zurück.
    """
    unknown_chains = []

    for structure_data in combined_data:
        file_path = structure_data["file_path"]
        for chain in structure_data["atom_data"]:
            if chain["molecule_name"] == "UNKNOWN":
                unknown_chains.append({
                    "file_path": file_path,
                    "chain_id": chain["chain_id"]
                })

    return unknown_chains
