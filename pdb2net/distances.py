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
    return np.array(ca_nn).reshape(-1, 3) if ca_nn else np.array([])

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

                # Verhindere Selbst-Interaktion
                if chain_a["chain_id"] == chain_b["chain_id"]:
                    continue

                print(f"\n🔍 Computing distances: {chain_a['unique_chain_id']} ↔ {chain_b['unique_chain_id']}")

                # Prüfe Molekültyp
                molecule_type_a = chain_a.get("molecule_type", "Unknown")
                molecule_type_b = chain_b.get("molecule_type", "Unknown")

                if molecule_type_a == "Unknown" or molecule_type_b == "Unknown":
                    print(f"⚠ WARNING: Unknown molecule type for {chain_a['chain_id']} or {chain_b['chain_id']}")


                # Klassifiziere Interaktionstyp
                if molecule_type_a == "Protein" and molecule_type_b == "Protein":
                    interaction_type = "Protein-Protein"
                elif (molecule_type_a == "Protein" and molecule_type_b == "Nucleic Acid") or \
                     (molecule_type_a == "Nucleic Acid" and molecule_type_b == "Protein"):
                    interaction_type = "Protein-Nucleic Acid"
                elif molecule_type_a == "Nucleic Acid" and molecule_type_b == "Nucleic Acid":
                    interaction_type = "Nucleic Acid-Nucleic Acid"
                else:
                    interaction_type = "Unknown"

                # Berechnung der Distanzen
                ca_nn_a = extract_ca_nn(chain_a) if molecule_type_a == "Protein" else extract_all_atoms(chain_a)
                ca_nn_b = extract_ca_nn(chain_b) if molecule_type_b == "Protein" else extract_all_atoms(chain_b)

                if ca_nn_a.size == 0 or ca_nn_b.size == 0:
                    print(f"⚠ WARNING: Keine Atome für {chain_a['chain_id']} oder {chain_b['chain_id']}")
                    continue

                tree_a = cKDTree(ca_nn_a)
                ca_nn_count = get_close_pairs(tree_a, ca_nn_b, radius=15.0)

                tree_atoms_a = cKDTree(extract_all_atoms(chain_a))
                all_atoms_close_count = get_close_pairs(tree_atoms_a, extract_all_atoms(chain_b), radius=5.0)

                # ❌ Falls keine Interaktion stattgefunden hat, überspringen
                if ca_nn_count == 0 and all_atoms_close_count == 0:
                    continue  # 🚫 Keine relevante Interaktion → wird nicht gespeichert

                # ✅ Nur Interaktionen mit mindestens einem Kontakt speichern
                results.append({
                    "file_path": file_path,
                    "chain_a": chain_a["unique_chain_id"],
                    "chain_b": chain_b["unique_chain_id"],
                    "ca_nn_count": ca_nn_count,
                    "all_atoms_close_count": all_atoms_close_count,
                    "interaction_type": interaction_type
                })

    return results
