import re
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

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
