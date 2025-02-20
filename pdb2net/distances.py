from scipy.spatial import cKDTree
import numpy as np
from config_loader import config

def get_nearest_heavy_atom(residue, ca_coord):
    """
    Finds the nearest heavy atom (excluding hydrogen) for a given Cα atom.
    If no heavy atom is found, the Cα coordinates are returned.

    Args:
        residue (dict): A dictionary representing the residue containing atom data.
        ca_coord (list): The coordinates of the Cα atom.

    Returns:
        list: The coordinates of the nearest heavy atom, or the Cα coordinates if none are found.
    """
    heavy_atoms = [
        atom for atom in residue["atoms"]
        if atom["atom_name"] != "CA" and not atom["atom_name"].startswith("H")
    ]

    if not heavy_atoms:
        return ca_coord  # Use Cα if no heavy atom is available

    nearest_atom = min(
        heavy_atoms,
        key=lambda atom: np.linalg.norm(np.array(atom["coordinates"]) - np.array(ca_coord))
    )
    return nearest_atom["coordinates"]

def extract_ca_nn(chain):
    """
    Extracts Cα atoms and their nearest heavy atoms from a protein chain.

    Args:
        chain (dict): Dictionary containing chain information.

    Returns:
        np.ndarray: A NumPy array of shape (N, 3) containing Cα and nearest neighbor (NN) coordinates.
    """
    ca_nn = [
        [ca_atom["coordinates"], get_nearest_heavy_atom(residue, ca_atom["coordinates"])]
        for residue in chain["residues"]
        for ca_atom in residue["atoms"]
        if ca_atom["atom_name"] == "CA"
    ]
    return np.array(ca_nn).reshape(-1, 3) if ca_nn else np.array([])  # Ensure (N,3) shape

def extract_all_atoms(chain):
    """
    Extracts all atomic coordinates from a molecular chain.

    Args:
        chain (dict): Dictionary containing chain information.

    Returns:
        np.ndarray: A NumPy array of shape (N, 3) containing all atom coordinates.
    """
    coords = [
        atom["coordinates"] for residue in chain["residues"] for atom in residue["atoms"]
    ]
    coords_array = np.array(coords)
    if coords_array.size == 0:
        print(f"⚠ WARNING: No atoms found in chain {chain['chain_id']}")
    return coords_array.reshape(-1, 3) if coords_array.size > 0 else np.array([])

def get_close_pairs(tree_a, points_b, radius):
    """
    Computes the number of atoms in `points_b` that are within a given radius of `tree_a`.

    Args:
        tree_a (cKDTree): A KDTree built from the coordinates of chain A.
        points_b (np.ndarray): A NumPy array containing the coordinates of chain B.
        radius (float): The maximum distance threshold for counting close contacts.

    Returns:
        int: The total number of close-contact pairs within the given radius.
    """
    if len(points_b) == 0:
        print(f"⚠ WARNING: Empty point list in distance calculation! Radius: {radius}")
        return 0
    close_points = tree_a.query_ball_point(points_b, r=radius)
    return sum(len(points) for points in close_points)

def calculate_distances_with_ckdtree(combined_data):
    """
    Computes inter-chain atomic distances using cKDTree and filters relevant interactions.

    Args:
        combined_data (list): List of dictionaries containing PDB structures and atom data.

    Returns:
        list: A list of dictionaries containing pairwise interaction results.
    """
    results = []

    for file_data in combined_data:
        file_path = file_data["file_path"]
        atom_data = file_data["atom_data"]

        for i, chain_a in enumerate(atom_data):
            for j, chain_b in enumerate(atom_data):
                if i >= j:
                    continue  # Avoid redundant calculations

                print(f"\n🔍 Computing distances: {chain_a['unique_chain_id']} ↔ {chain_b['unique_chain_id']}")

                # Determine if at least one of the chains is a protein
                protein_interaction = (
                    chain_a["molecule_type"] == "Protein" or chain_b["molecule_type"] == "Protein"
                )

                if protein_interaction:
                    # Protein-Protein or Protein-Nucleic Acid interactions: Use Cα + nearest heavy atom
                    ca_nn_a = extract_ca_nn(chain_a) if chain_a["molecule_type"] == "Protein" else extract_all_atoms(chain_a)
                    ca_nn_b = extract_ca_nn(chain_b) if chain_b["molecule_type"] == "Protein" else extract_all_atoms(chain_b)

                    ca_nn_a_flat = np.array([coord for coord in ca_nn_a]) if protein_interaction else ca_nn_a
                    ca_nn_b_flat = np.array([coord for coord in ca_nn_b]) if protein_interaction else ca_nn_b
                else:
                    # Only nucleic acids: Use all atoms
                    ca_nn_a_flat = extract_all_atoms(chain_a)
                    ca_nn_b_flat = extract_all_atoms(chain_b)

                # Debugging: Print array shapes
                print(f"  ✅ Shape of {chain_a['chain_id']} (A): {ca_nn_a_flat.shape}")
                print(f"  ✅ Shape of {chain_b['chain_id']} (B): {ca_nn_b_flat.shape}")

                # Validate array shapes before proceeding
                if ca_nn_a_flat.ndim != 2 or ca_nn_a_flat.shape[1] != 3:
                    print(f"❌ Error: Invalid shape for {chain_a['chain_id']}: {ca_nn_a_flat.shape}")
                    continue
                if ca_nn_b_flat.ndim != 2 or ca_nn_b_flat.shape[1] != 3:
                    print(f"❌ Error: Invalid shape for {chain_b['chain_id']}: {ca_nn_b_flat.shape}")
                    continue

                # Step 1: Compute distances using cKDTree (radius 15 Å)
                tree_a = cKDTree(ca_nn_a_flat)
                ca_nn_count = get_close_pairs(tree_a, ca_nn_b_flat, radius=15.0)

                if ca_nn_count < 10:
                    print(f"❌ Fewer than 10 Cα/NN pairs found within <15 Å. Skipping.")
                    continue

                # Step 2: Compute all-atom contacts using cKDTree (radius 5 Å)
                all_atoms_a = extract_all_atoms(chain_a)
                all_atoms_b = extract_all_atoms(chain_b)

                if all_atoms_a.size == 0 or all_atoms_b.size == 0:
                    print(f"⚠ WARNING: Empty atom list for {chain_a['chain_id']} or {chain_b['chain_id']}")
                    continue

                tree_atoms_a = cKDTree(all_atoms_a)
                all_atoms_close_count = get_close_pairs(tree_atoms_a, all_atoms_b, radius=5.0)

                if all_atoms_close_count < 10:
                    print(f"❌ Fewer than 10 atom pairs found within <5 Å. Skipping.")
                    continue

                # Store results
                results.append({
                    "file_path": file_path,
                    "chain_a": chain_a["unique_chain_id"],
                    "chain_b": chain_b["unique_chain_id"],
                    "ca_nn_count": ca_nn_count,
                    "all_atoms_close_count": all_atoms_close_count,
                })

                print(f"✅ Distances computed: {ca_nn_count} Cα/NN pairs, {all_atoms_close_count} atom pairs <5 Å")

    return results
