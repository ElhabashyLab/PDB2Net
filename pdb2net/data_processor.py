def process_structure(structure_data):
    """Verarbeitet eine geparste Struktur und erstellt die Datenstruktur."""
    pdb_id = structure_data["pdb_id"]
    structure = structure_data["structure"]

    atom_data = []
    for model in structure:
        for chain in model:
            chain_id = chain.id.strip()
            residues = [
                {"residue_name": residue.resname, "atoms": [{"atom_name": atom.get_name(), "coordinates": atom.get_coord().tolist()} for atom in residue.get_atoms()]}
                for residue in chain.get_residues()
            ]

            atom_data.append({
                "chain_id": chain_id,
                "unique_chain_id": f"{pdb_id}:{chain_id}",
                "molecule_name": "UNKNOWN",
                "molecule_type": "UNKNOWN",
                "sequence": "",
                "residues": residues
            })

    return {"file_path": structure_data["file_path"], "pdb_id": pdb_id, "atom_data": atom_data}
