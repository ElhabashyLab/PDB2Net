from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1
import os
import re


def generate_atoms_info(residue):
    """
    Generator, der Atominformationen für ein Residuum erzeugt.
    """
    for atom in residue.get_atoms():
        yield {
            "atom_name": atom.get_name(),
            "coordinates": atom.get_coord().tolist(),  # Koordinaten als Liste [x, y, z]
            "element": atom.element
        }


def generate_residues_info(chain):
    """
    Generator, der Residueninformationen für eine Kette erzeugt.
    """
    for residue in chain.get_residues():
        if residue.id[0] == " ":  # Nur ATOM-Residuen
            atoms_info = generate_atoms_info(residue)
            yield {
                "residue_name": residue.resname,
                "residue_id": residue.id[1],  # Residuen-ID
                "atoms": list(atoms_info),   # Enthält Atominformationen inkl. Koordinaten
            }


def generate_atom_data(structure, chain_to_molecule):
    """
    Generator, der ATOM-Daten für Ketten erzeugt.
    """
    for model in structure:
        for chain in model:
            chain_id = chain.id.strip()
            unique_chain_id = f"{structure.id}:{chain_id}"  # Eindeutige Kennzeichnung
            molecule_name, molecule_type = chain_to_molecule.get(chain_id, ("UNKNOWN", "UNKNOWN"))
            residues_info = generate_residues_info(chain)
            yield {
                "chain_id": chain_id,
                "unique_chain_id": unique_chain_id,  # Eindeutige ID
                "residues": list(residues_info),    # Residuen enthalten Atominformationen
                "molecule_type": molecule_type,
                "molecule_name": molecule_name,
            }



def extract_pdb_molecule_mapping(file_path):
    """
    Parst den COMPND-Block aus einer PDB-Datei und erstellt ein Mapping von Ketten zu Molekülnamen.
    """
    chain_to_molecule = {}
    current_molecule_name = None

    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("COMPND"):
                molecule_match = re.search(r"MOLECULE:\s*([^;]+)", line)
                chain_match = re.search(r"CHAIN:\s*([^;]+)", line)
                if molecule_match:
                    current_molecule_name = molecule_match.group(1).strip()
                if chain_match and current_molecule_name:
                    chains = chain_match.group(1).replace(" ", "").split(",")
                    for chain_id in chains:
                        chain_to_molecule[chain_id] = (current_molecule_name, "Protein")
    return chain_to_molecule


def extract_cif_molecule_data(mmcif_dict):
    """
    Extrahiert Molekültypen und Namen aus einer mmCIF-Datei.
    """
    entity_ids = mmcif_dict["_entity.id"]
    descriptions = mmcif_dict["_entity.pdbx_description"]
    polymer_types = mmcif_dict.get("_entity_poly.type", ["UNKNOWN"] * len(entity_ids))
    chain_mappings = mmcif_dict["_entity_poly.pdbx_strand_id"]

    chain_to_molecule = {}
    for entity_id, description, polymer_type, chains in zip(entity_ids, descriptions, polymer_types, chain_mappings):
        for chain in chains.split(","):
            chain_to_molecule[chain.strip()] = (description, polymer_type)
    return chain_to_molecule


def process_structure(structure_data):
    """
    Verarbeitet eine geparste Struktur und erstellt die hierarchische Datenstruktur.
    """
    file_path = structure_data["file_path"]
    structure = structure_data["structure"]
    _, ext = os.path.splitext(file_path)

    # Mapping erstellen (PDB oder mmCIF)
    if ext == ".pdb":
        chain_to_molecule = extract_pdb_molecule_mapping(file_path)
    elif ext in {".cif", ".mmcif"}:
        mmcif_dict = MMCIF2Dict(file_path)
        chain_to_molecule = extract_cif_molecule_data(mmcif_dict)
    else:
        raise ValueError("Unsupported file format")

    # Kombinierte Datenstruktur zurückgeben
    atom_data = list(generate_atom_data(structure, chain_to_molecule))
    return {
        "file_path": file_path,
        "atom_data": atom_data,
    }
