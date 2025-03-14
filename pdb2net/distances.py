from scipy.spatial import cKDTree
import numpy as np
import itertools

# Cache für KD-Trees und extrahierte Koordinaten
tree_cache, coords_cache = {}, {}


def get_nearest_heavy_atom(residue, ca_coord):
    """Findet das nächstgelegene schwere Atom (außer Wasserstoff) zu einem Cα-Atom."""
    heavy_atoms = [
        atom for atom in residue["atoms"] if atom["atom_name"] != "CA" and not atom["atom_name"].startswith("H")
    ]
    return min(heavy_atoms, key=lambda atom: np.linalg.norm(np.array(atom["coordinates"]) - ca_coord),
               default={"coordinates": ca_coord})["coordinates"]


def extract_coordinates(chain, extraction_type):
    """Extrahiert Atomkoordinaten aus einer Kette und speichert sie im Cache."""
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
    """Erstellt oder ruft einen gespeicherten cKDTree für eine bestimmte Extraktion ab."""
    key = (chain["unique_chain_id"], extraction_type)
    if key not in tree_cache:
        points = extract_coordinates(chain, extraction_type)
        if points.size:
            tree_cache[key] = cKDTree(points)
    return tree_cache.get(key)


def count_nearby_atoms(tree_a, tree_b, radius):
    """Zählt die Anzahl der benachbarten Atome innerhalb eines gegebenen Radius."""
    return sum(len(p) for p in tree_a.query_ball_point(tree_b.data, r=radius)) if tree_a and tree_b else 0


def determine_interaction_type(mol_type_a, mol_type_b):
    """Bestimmt den Typ der molekularen Interaktion."""
    if "Unknown" in [mol_type_a, mol_type_b]:
        return None
    if mol_type_a == mol_type_b == "Protein":
        return "Protein-Protein"
    if "Nucleic Acid" in [mol_type_a, mol_type_b]:
        return "Protein-Nucleic Acid" if "Protein" in [mol_type_a, mol_type_b] else "Nucleic Acid-Nucleic Acid"
    return "Unknown"


def calculate_distances_with_ckdtree(combined_data):
    """Berechnet interkettige Distanzen effizient mit cKDTree."""
    results = []
    for file_data in combined_data:
        file_path, atom_data = file_data["file_path"], file_data["atom_data"]

        # KD-Trees für jede Kette vorbereiten
        local_trees = {(c["unique_chain_id"], t): get_or_create_tree(c, t) for c in atom_data for t in
                       ["ca", "nn_heavy", "all_atoms"]}

        for chain_a, chain_b in itertools.combinations(atom_data, 2):
            if chain_a["chain_id"] == chain_b["chain_id"]:
                continue  # Selbst-Interaktion vermeiden

            interaction_type = determine_interaction_type(chain_a.get("molecule_type", "Unknown"),
                                                          chain_b.get("molecule_type", "Unknown"))
            if not interaction_type:
                continue

            if "Nucleic Acid" in [chain_a["molecule_type"], chain_b["molecule_type"]]:
                all_atoms_count = count_nearby_atoms(local_trees.get((chain_a["unique_chain_id"], "all_atoms")),
                                                     local_trees.get((chain_b["unique_chain_id"], "all_atoms")),
                                                     radius=5.0)
                if all_atoms_count:
                    results.append({"file_path": file_path, "chain_a": chain_a["unique_chain_id"],
                                    "chain_b": chain_b["unique_chain_id"], "all_atoms_count": all_atoms_count,
                                    "interaction_type": interaction_type})
                continue

            # Erst Cα-Distanzen prüfen, dann NN-Heavy-Atome
            if count_nearby_atoms(local_trees.get((chain_a["unique_chain_id"], "ca")),
                                  local_trees.get((chain_b["unique_chain_id"], "ca")),
                                  radius=15.0) >= 10 or count_nearby_atoms(
                    local_trees.get((chain_a["unique_chain_id"], "nn_heavy")),
                    local_trees.get((chain_b["unique_chain_id"], "nn_heavy")), radius=15.0) >= 10:
                all_atoms_count = count_nearby_atoms(local_trees.get((chain_a["unique_chain_id"], "all_atoms")),
                                                     local_trees.get((chain_b["unique_chain_id"], "all_atoms")),
                                                     radius=5.0)
                if all_atoms_count:
                    results.append({"file_path": file_path, "chain_a": chain_a["unique_chain_id"],
                                    "chain_b": chain_b["unique_chain_id"], "all_atoms_count": all_atoms_count,
                                    "interaction_type": interaction_type})

    return results
