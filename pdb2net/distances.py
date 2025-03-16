from scipy.spatial import cKDTree
import numpy as np
import itertools
from config_loader import config

# Cache for KD-Trees and extracted coordinates
tree_cache, coords_cache = {}, {}

def get_nearest_heavy_atom(residue, ca_coord):
    """
    Finds the nearest heavy atom (excluding hydrogen) to a given Cα atom.

    Args:
        residue (dict): Dictionary containing residue information.
        ca_coord (list): Coordinates of the Cα atom.

    Returns:
        list: Coordinates of the nearest heavy atom, or Cα coordinates if no heavy atom is found.
    """
    heavy_atoms = [
        atom for atom in residue["atoms"] if atom["atom_name"] != "CA" and not atom["atom_name"].startswith("H")
    ]
    return min(heavy_atoms, key=lambda atom: np.linalg.norm(np.array(atom["coordinates"]) - ca_coord),
               default={"coordinates": ca_coord})["coordinates"]

def extract_coordinates(chain, extraction_type):
    """
    Extracts atom coordinates from a chain and caches them.

    Args:
        chain (dict): Chain data containing residue and atom information.
        extraction_type (str): Type of extraction ('ca', 'nn_heavy', or 'all_atoms').

    Returns:
        np.ndarray: Array of extracted coordinates.
    """
    key = (chain["unique_chain_id"], extraction_type)
    if key in coords_cache:
        return coords_cache[key]

    atoms = (atom for residue in chain["residues"] for atom in residue["atoms"])
    if extraction_type == "ca":
        coords = np.array([atom["coordinates"] for atom in atoms if atom["atom_name"] == "CA"])
    elif extraction_type == "nn_heavy":
        coords = np.array(
            [get_nearest_heavy_atom(residue, atom["coordinates"]) for residue in chain["residues"] for atom in
             residue["atoms"] if atom["atom_name"] == "CA"])
    elif extraction_type == "all_atoms":
        coords = np.array([atom["coordinates"] for atom in atoms])
    else:
        return np.array([])

    coords_cache[key] = coords if coords.size else np.array([])
    return coords_cache[key]

def get_or_create_tree(chain, extraction_type):
    """
    Creates or retrieves a stored cKDTree for a specific extraction type.

    Args:
        chain (dict): Chain data containing residue and atom information.
        extraction_type (str): Type of extraction ('ca', 'nn_heavy', 'all_atoms').

    Returns:
        cKDTree or None: The created or cached KD-Tree.
    """
    key = (chain["unique_chain_id"], extraction_type)
    if key not in tree_cache:
        points = extract_coordinates(chain, extraction_type)
        if points.size:
            tree_cache[key] = cKDTree(points)
    return tree_cache.get(key)

def count_nearby_atoms(tree_a, tree_b, radius):
    """
    Counts the number of neighboring atoms within a given radius.

    Args:
        tree_a (cKDTree): KD-Tree of chain A.
        tree_b (cKDTree): KD-Tree of chain B.
        radius (float): Distance threshold for counting neighboring atoms.

    Returns:
        int: Count of nearby atoms.
    """
    return sum(len(p) for p in tree_a.query_ball_point(tree_b.data, r=radius)) if tree_a and tree_b else 0

def determine_interaction_type(mol_type_a, mol_type_b):
    """
    Determines the molecular interaction type based on molecule types.

    Args:
        mol_type_a (str): Molecule type of chain A.
        mol_type_b (str): Molecule type of chain B.

    Returns:
        str or None: Interaction type or None if it cannot be determined.
    """
    if "Unknown" in [mol_type_a, mol_type_b]:
        return None
    if mol_type_a == mol_type_b == "Protein":
        return "Protein-Protein"
    if "Nucleic Acid" in [mol_type_a, mol_type_b]:
        return "Protein-Nucleic Acid" if "Protein" in [mol_type_a, mol_type_b] else "Nucleic Acid-Nucleic Acid"
    return "Unknown"


# Load radius thresholds from config
RADIUS_CA = config["distance_thresholds"]["ca_radius"]
RADIUS_NN_HEAVY = config["distance_thresholds"]["nn_heavy_radius"]
RADIUS_ALL_ATOMS = config["distance_thresholds"]["all_atoms_radius"]


def calculate_distances_with_ckdtree(combined_data):
    """
    Computes inter-chain distances efficiently using cKDTree.

    Args:
        combined_data (list): List of dictionaries containing parsed PDB structures.

    Returns:
        list: A list of dictionaries with interaction results.
    """
    results = []
    for file_data in combined_data:
        file_path, atom_data = file_data["file_path"], file_data["atom_data"]

        # Prepare KD-Trees for each chain
        local_trees = {(c["unique_chain_id"], t): get_or_create_tree(c, t) for c in atom_data for t in
                       ["ca", "nn_heavy", "all_atoms"]}

        for chain_a, chain_b in itertools.combinations(atom_data, 2):
            if chain_a["chain_id"] == chain_b["chain_id"]:
                continue  # Skip self-interactions

            interaction_type = determine_interaction_type(chain_a.get("molecule_type", "Unknown"),
                                                          chain_b.get("molecule_type", "Unknown"))
            if not interaction_type:
                continue

            if "Nucleic Acid" in [chain_a["molecule_type"], chain_b["molecule_type"]]:
                all_atoms_count = count_nearby_atoms(
                    local_trees.get((chain_a["unique_chain_id"], "all_atoms")),
                    local_trees.get((chain_b["unique_chain_id"], "all_atoms")),
                    radius=RADIUS_ALL_ATOMS
                )
                if all_atoms_count:
                    results.append({
                        "file_path": file_path,
                        "chain_a": chain_a["unique_chain_id"],
                        "chain_b": chain_b["unique_chain_id"],
                        "all_atoms_count": all_atoms_count,
                        "interaction_type": interaction_type
                    })
                continue

            # First check Cα distances, then NN-heavy atoms
            if count_nearby_atoms(
                    local_trees.get((chain_a["unique_chain_id"], "ca")),
                    local_trees.get((chain_b["unique_chain_id"], "ca")),
                    radius=RADIUS_CA) >= 10 or \
                    count_nearby_atoms(
                        local_trees.get((chain_a["unique_chain_id"], "nn_heavy")),
                        local_trees.get((chain_b["unique_chain_id"], "nn_heavy")),
                        radius=RADIUS_NN_HEAVY) >= 10:

                all_atoms_count = count_nearby_atoms(
                    local_trees.get((chain_a["unique_chain_id"], "all_atoms")),
                    local_trees.get((chain_b["unique_chain_id"], "all_atoms")),
                    radius=RADIUS_ALL_ATOMS
                )
                if all_atoms_count:
                    results.append({
                        "file_path": file_path,
                        "chain_a": chain_a["unique_chain_id"],
                        "chain_b": chain_b["unique_chain_id"],
                        "all_atoms_count": all_atoms_count,
                        "interaction_type": interaction_type
                    })

    return results