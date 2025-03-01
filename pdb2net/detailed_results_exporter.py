import os
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree

def export_detailed_interactions(structure_data, interactions, run_output_path):
    """
    Exports detailed interaction results, including residues and atoms involved, to a CSV file.
    """
    pdb_id = structure_data["pdb_id"]
    atom_data = structure_data["atom_data"]

    residues_atoms_lookup = {}
    for chain in atom_data:
        chain_id = chain["chain_id"]
        residues_atoms_lookup[chain_id] = []
        for residue in chain["residues"]:
            res_id = residue["residue_name"]
            res_number = residue.get("residue_number", "?")
            for atom in residue["atoms"]:
                residues_atoms_lookup[chain_id].append({
                    "residue": f"{res_number}:{res_id}",
                    "atom_name": atom["atom_name"],
                    "coordinates": atom["coordinates"]
                })

    detailed_interactions = []
    for interaction in interactions:
        chain_a_id = interaction["chain_a"].split(":")[1]
        chain_b_id = interaction["chain_b"].split(":")[1]

        atoms_a = residues_atoms_lookup.get(chain_a_id, [])
        atoms_b = residues_atoms_lookup.get(chain_b_id, [])

        if not atoms_a or not atoms_b:
            continue

        chain_a_data = next((ch for ch in atom_data if ch["chain_id"] == chain_a_id), None)
        chain_b_data = next((ch for ch in atom_data if ch["chain_id"] == chain_b_id), None)

        if not chain_a_data or not chain_b_data:
            continue

        uniprot_a = chain_a_data.get("uniprot_id", "UNKNOWN")
        uniprot_b = chain_b_data.get("uniprot_id", "UNKNOWN")

        interaction_type = interaction["interaction_type"]

        coords_a = np.array([atom["coordinates"] for atom in atoms_a])
        coords_b = np.array([atom["coordinates"] for atom in atoms_b])

        if coords_a.size == 0 or coords_b.size == 0:
            continue

        tree_a = cKDTree(coords_a)
        pairs = tree_a.query_ball_point(coords_b, r=5.0)

        for idx_b, idx_list_a in enumerate(pairs):
            for idx_a in idx_list_a:
                atom_a = atoms_a[idx_a]
                atom_b = atoms_b[idx_b]

                distance = np.linalg.norm(
                    np.array(atom_a["coordinates"]) - np.array(atom_b["coordinates"])
                )

                detailed_interactions.append({
                    "PDB_ID": pdb_id,
                    "Chain_A": chain_a_id,
                    "Residue_A": atom_a["residue"],
                    "Atom_A": atom_a["atom_name"],
                    "Chain_B": chain_b_id,
                    "Residue_B": atom_b["residue"],
                    "Atom_B": atom_b["atom_name"],
                    "Distance": round(distance, 2),
                    "UniProt_A": uniprot_a,
                    "UniProt_B": uniprot_b,
                    "Interaction_Type": interaction_type
                })

    df = pd.DataFrame(detailed_interactions)

    output_dir = os.path.join(run_output_path, pdb_id)
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{pdb_id}_detailed_interactions.csv")
    df.to_csv(output_file, index=False)

    print(f"✅ CSV export successful: {output_file} (showing first 10 rows):")
    print(df.head(10))
