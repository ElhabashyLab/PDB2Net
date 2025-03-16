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
                - "is_hetatm": True if the chain consists only of HETATM.
                - "residues": List of residues with atom information.
    """
    pdb_id = structure_data["pdb_id"]
    structure = structure_data["structure"]

    raw_atom_data = []  # Stores all ATOM & HETATM entries
    for model in structure:
        for chain in model:
            chain_id = chain.id.strip()

            residues = []
            only_hetatm = True  # Assume the chain consists only of HETATM initially

            for residue in chain.get_residues():
                res_name = residue.resname.strip()

                # Extract residue number and insertion code
                res_number = residue.id[1]  # Residue number from PDB
                insert_code = residue.id[2]  # Optional insertion code

                if insert_code and insert_code != " ":
                    full_res_id = f"{res_number}{insert_code}"
                else:
                    full_res_id = str(res_number)

                # Skip HETATM residues
                if residue.id[0] != " ":
                    continue

                # Check if the residue contains any non-hydrogen atom
                contains_atom = False
                for atom in residue.get_atoms():
                    if atom.element not in ["H", "D"]:  # Heavy atom found
                        contains_atom = True
                        break  # Stop checking further atoms

                if contains_atom:
                    only_hetatm = False

                residue_entry = {
                    "residue_name": res_name,
                    "residue_number": full_res_id,
                    "atoms": [
                        {
                            "atom_name": atom.get_name(),
                            "coordinates": atom.get_coord().tolist()
                        }
                        for atom in residue.get_atoms()
                    ]
                }
                residues.append(residue_entry)

            raw_atom_data.append({
                "chain_id": chain_id,
                "unique_chain_id": f"{pdb_id}:{chain_id}",
                "molecule_name": "UNKNOWN",
                "molecule_type": "UNKNOWN",
                "sequence": "",
                "is_hetatm": only_hetatm,
                "residues": residues
            })

    # Remove chains that contain only HETATM residues
    atom_data = [chain for chain in raw_atom_data if not chain["is_hetatm"]]

    return {
        "file_path": structure_data["file_path"],
        "pdb_id": pdb_id,
        "atom_data": atom_data
    }
