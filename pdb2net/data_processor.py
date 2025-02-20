from config_loader import config

def process_structure(structure_data):
    """
    Processes a parsed PDB or mmCIF structure and extracts relevant atom and residue information.

    Args:
        structure_data (dict): A dictionary containing the parsed structure with:
            - "pdb_id": The PDB ID of the structure.
            - "structure": The Biopython Structure object.

    Returns:
        dict: A processed data dictionary containing:
            - "file_path": The original file path.
            - "pdb_id": The PDB ID.
            - "atom_data": A list of dictionaries, each representing a chain with:
                - "chain_id": The chain identifier.
                - "unique_chain_id": A unique identifier combining PDB ID and chain ID.
                - "molecule_name": Placeholder (to be updated later).
                - "molecule_type": Placeholder (to be updated later).
                - "sequence": Placeholder for the sequence.
                - "residues": List of residues with atom information.
    """
    pdb_id = structure_data["pdb_id"]
    structure = structure_data["structure"]

    atom_data = []
    for model in structure:
        for chain in model:
            chain_id = chain.id.strip()

            # Extract residues and their atom information
            residues = [
                {
                    "residue_name": residue.resname,
                    "atoms": [
                        {
                            "atom_name": atom.get_name(),
                            "coordinates": atom.get_coord().tolist()
                        }
                        for atom in residue.get_atoms()
                    ]
                }
                for residue in chain.get_residues()
            ]

            # Append processed chain data
            atom_data.append({
                "chain_id": chain_id,
                "unique_chain_id": f"{pdb_id}:{chain_id}",
                "molecule_name": "UNKNOWN",  # To be determined later
                "molecule_type": "UNKNOWN",  # To be determined later
                "sequence": "",  # To be determined later
                "residues": residues
            })

    return {
        "file_path": structure_data["file_path"],
        "pdb_id": pdb_id,
        "atom_data": atom_data
    }
