import gemmi

# Erlaubte Residuen (3-letter-codes) für Proteine und Nukleinsäuren
ALLOWED_RESIDUES = {
    # Proteine
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "SEC", "PYL",
    # RNA
    "A", "U", "G", "C", "I",
    # DNA (Desoxyformen)
    "DA", "DT", "DG", "DC", "DI"
}

def process_structure(structure_data):
    file_path = structure_data["file_path"]
    pdb_id = structure_data["pdb_id"]
    structure = structure_data["structure"]

    atom_data = []

    for model in structure:
        for chain in model:
            chain_id = chain.name.strip()
            residues = []
            is_valid_chain = False

            for res in chain:
                res_name = res.name.upper()
                atoms = []

                for atom in res:
                    if atom.element.name not in ["H", "D"]:
                        atoms.append({
                            "atom_name": atom.name,
                            "coordinates": list(atom.pos)
                        })

                if atoms:
                    residues.append({
                        "residue_name": res_name,
                        "residue_number": res.seqid.num,
                        "atoms": atoms
                    })

                    # Akzeptiere Chain, wenn mind. 1 Residuum aus Protein oder NA besteht
                    if res_name in ALLOWED_RESIDUES:
                        is_valid_chain = True

            if is_valid_chain:
                atom_data.append({
                    "chain_id": chain_id,
                    "unique_chain_id": f"{pdb_id}:{chain_id}",
                    "molecule_name": "UNKNOWN",
                    "molecule_type": "UNKNOWN",
                    "sequence": "",
                    "is_hetatm": False,
                    "residues": residues
                })

    return {
        "file_path": file_path,
        "pdb_id": pdb_id,
        "atom_data": atom_data
    }
