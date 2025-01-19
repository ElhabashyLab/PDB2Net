from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import os
import warnings
from Bio import BiopythonWarning
import re

warnings.simplefilter("ignore", BiopythonWarning)


def parse_structure(file_path):
    """
    Parst eine PDB- oder mmCIF-Datei und extrahiert Basisinformationen inklusive ATOM- und HETATM-Daten.

    :param file_path: Pfad zur PDB- oder mmCIF-Datei
    :return: Zwei Listen:
             - ATOM-Daten (Ketteninformationen: Ketten-ID, Residuen, Molekültyp, Molekülname)
             - HETATM-Daten (Residueninformationen für Heteroatome)
    """
    _, ext = os.path.splitext(file_path)

    if ext == ".pdb":
        return parse_pdb_structure(file_path)
    elif ext == ".cif":
        return parse_cif_structure(file_path)
    else:
        raise ValueError(f"Unbekanntes Dateiformat: {ext}")


def parse_pdb_structure(file_path):
    """
    Parst eine PDB-Datei und extrahiert ATOM- und HETATM-Daten.

    :param file_path: Pfad zur PDB-Datei
    :return: Zwei Listen:
             - ATOM-Daten (Ketteninformationen)
             - HETATM-Daten (Residueninformationen für Heteroatome)
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", file_path)

    # COMPND-Block verarbeiten
    chain_to_molecule = extract_pdb_molecule_mapping(file_path)

    # ATOM-Daten extrahieren
    atom_chains = extract_pdb_atom_chains(structure, chain_to_molecule)

    # HETATM-Daten extrahieren
    hetatm_molecules = extract_hetatm_molecules(structure)

    return atom_chains, hetatm_molecules


def extract_pdb_molecule_mapping(file_path):
    """
    Parst den COMPND-Block aus einer PDB-Datei und erstellt ein Mapping von Ketten zu Molekülnamen.

    :param file_path: Pfad zur PDB-Datei
    :return: Dictionary mit Ketten-ID als Schlüssel und (Molekülname, Typ) als Wert
    """
    chain_to_molecule = {}
    current_molecule_name = None

    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("COMPND"):
                # Kombinierte Suche für MOLECULE und CHAIN
                molecule_match = re.search(r"MOLECULE:\s*([^;]+)", line)
                chain_match = re.search(r"CHAIN:\s*([^;]+)", line)

                if molecule_match:
                    current_molecule_name = molecule_match.group(1).strip()
                if chain_match and current_molecule_name:
                    chains = chain_match.group(1).replace(" ", "").split(",")
                    for chain_id in chains:
                        chain_to_molecule[chain_id] = (current_molecule_name, "Protein")

    return chain_to_molecule


def extract_pdb_atom_chains(structure, chain_to_molecule):
    """
    Extrahiert ATOM-Daten (Makromoleküle) inklusive Molekültyp und Molekülname aus einer PDB-Datei.

    :param structure: Geparste Struktur
    :param chain_to_molecule: Mapping von Ketten zu Molekülnamen
    :return: Liste mit Ketteninformationen
    """
    chains_info = []
    for model in structure:
        for chain in model:
            residues = [residue.resname for residue in chain.get_residues() if residue.id[0] == " "]
            chain_id = chain.id.strip()

            # Molekülname und Typ aus dem Mapping holen
            molecule_name, molecule_type = chain_to_molecule.get(chain_id, ("UNKNOWN", "UNKNOWN"))

            chains_info.append({
                "chain_id": chain_id,
                "residues": residues,
                "molecule_type": molecule_type,
                "molecule_name": molecule_name
            })

    return chains_info


def parse_cif_structure(file_path):
    """
    Parst eine mmCIF-Datei und extrahiert ATOM- und HETATM-Daten.

    :param file_path: Pfad zur mmCIF-Datei
    :return: Zwei Listen:
             - ATOM-Daten (Ketteninformationen)
             - HETATM-Daten (Residueninformationen für Heteroatome)
    """
    parser = MMCIFParser()
    structure = parser.get_structure("structure", file_path)

    # ATOM-Daten extrahieren
    atom_chains = extract_atom_chains(structure, file_path)

    # HETATM-Daten extrahieren
    hetatm_molecules = extract_hetatm_molecules(structure)

    return atom_chains, hetatm_molecules


def extract_atom_chains(structure, file_path):
    """
    Extrahiert ATOM-Daten (Makromoleküle) inklusive tatsächlichem Molekültyp und Molekülname aus einer mmCIF-Datei.

    :param structure: Geparste Struktur
    :param file_path: Pfad zur mmCIF-Datei
    :return: Liste mit Ketteninformationen
    """
    chains_info = []

    # Molekültypen und Namen aus der Datei extrahieren
    molecule_data = extract_cif_molecule_data(file_path)

    for model in structure:
        for chain in model:
            # Residuen der Kette extrahieren
            residues = [residue.resname for residue in chain.get_residues() if residue.id[0] == " "]
            chain_id = chain.id

            # Hole Molekültyp und Name aus den extrahierten Daten
            molecule_name, molecule_type = molecule_data.get(chain_id, ("UNKNOWN", "UNKNOWN"))

            chains_info.append({
                "chain_id": chain_id,
                "residues": residues,
                "molecule_type": molecule_type,
                "molecule_name": molecule_name
            })

    return chains_info


def extract_cif_molecule_data(file_path):
    """
    Extrahiert Molekültypen und Namen aus einer mmCIF-Datei.

    :param file_path: Pfad zur mmCIF-Datei
    :return: Dictionary mit Ketten-ID als Schlüssel und (Molekülname, Molekültyp) als Wert
    """
    mmcif_dict = MMCIF2Dict(file_path)

    # Entity-Daten extrahieren
    entity_ids = mmcif_dict["_entity.id"]
    descriptions = mmcif_dict["_entity.pdbx_description"]
    polymer_types = mmcif_dict.get("_entity_poly.type", ["UNKNOWN"] * len(entity_ids))
    chain_mappings = mmcif_dict["_entity_poly.pdbx_strand_id"]

    # Mapping erstellen: Ketten-ID -> (Molekülname, Molekültyp)
    chain_to_molecule = {}
    for entity_id, description, polymer_type, chains in zip(entity_ids, descriptions, polymer_types, chain_mappings):
        for chain in chains.split(","):
            chain_to_molecule[chain.strip()] = (description, polymer_type)

    return chain_to_molecule


def extract_hetatm_molecules(structure):
    """
    Extrahiert einmalige Informationen zu HETATM-Daten (z. B. Liganden, Ionen, Wasser).

    :param structure: Geparste Struktur
    :return: Liste mit HETATM-Informationen (ein Eintrag pro Residuentyp)
    """
    hetatm_info = {}
    for model in structure:
        for chain in model:
            for residue in chain.get_residues():
                if residue.id[0] != " ":  # HETATM-Einträge haben ein anderes ID-Flag
                    residue_name = residue.resname
                    if residue_name not in hetatm_info:
                        hetatm_info[residue_name] = {
                            "chain_id": chain.id,
                            "residues": [residue_name],
                            "molecule_type": "HETATM",
                            "molecule_name": residue_name
                        }

    # Rückgabe als Liste mit einzigartigen HETATM-Daten
    return list(hetatm_info.values())
