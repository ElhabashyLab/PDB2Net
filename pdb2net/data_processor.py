from structure_extractor import extract_pdb_molecule_mapping, extract_cif_molecule_data
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import os

def process_structure(structure_data):
    """
    Verarbeitet eine geparste Struktur und erstellt die hierarchische Datenstruktur ohne doppelte Berechnungen.
    """
    file_path = structure_data["file_path"]
    structure = structure_data["structure"]
    _, ext = os.path.splitext(file_path)

    # Mapping nur einmal berechnen
    if ext == ".pdb":
        chain_to_molecule = extract_pdb_molecule_mapping(file_path)
    elif ext in {".cif", ".mmcif"}:
        mmcif_dict = MMCIF2Dict(file_path)
        chain_to_molecule = extract_cif_molecule_data(mmcif_dict)
    else:
        raise ValueError("Unsupported file format")

    # Übergabe des fertigen Mappings an die Datenstruktur
    atom_data = list(generate_atom_data(structure, chain_to_molecule))
    return {
        "file_path": file_path,
        "atom_data": atom_data,
    }

def generate_atom_data(structure, chain_to_molecule):
    """
    Generator, der ATOM-Daten für Ketten erzeugt.
    """
    for model in structure:
        for chain in model:
            chain_id = chain.id.strip()
            molecule_name, molecule_type = chain_to_molecule.get(chain_id, ("UNKNOWN", "UNKNOWN"))
            residues_info = generate_residues_info(chain)
            yield {
                "chain_id": chain_id,
                "residues": list(residues_info),
                "molecule_type": molecule_type,
                "molecule_name": molecule_name,
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
                "atoms": list(atoms_info),
            }

def generate_atoms_info(residue):
    """
    Generator, der Atominformationen für ein Residuum erzeugt.
    """
    for atom in residue.get_atoms():
        yield {
            "atom_name": atom.get_name(),
            "coordinates": atom.get_coord().tolist(),
            "element": atom.element,
        }
