from scipy.spatial import cKDTree
import numpy as np
import itertools

# Cache für bereits erstellte KD-Trees und extrahierte Koordinaten
tree_cache = {}
coords_cache = {}


def get_nearest_heavy_atom(residue, ca_coord):
    """Findet das nächstgelegene schwere Atom (außer Wasserstoff) für ein gegebenes Cα-Atom."""
    heavy_atoms = [
        atom for atom in residue["atoms"]
        if atom["atom_name"] != "CA" and not atom["atom_name"].startswith("H")
    ]

    if not heavy_atoms:
        return ca_coord  # Falls keine schweren Atome gefunden werden, nehme Cα

    nearest_atom = min(
        heavy_atoms,
        key=lambda atom: np.linalg.norm(np.array(atom["coordinates"]) - np.array(ca_coord))
    )
    return nearest_atom["coordinates"]


def extract_coordinates(chain, extraction_type):
    """Extrahiert Atomkoordinaten aus einer Kette und speichert sie im Cache."""
    key = (chain["unique_chain_id"], extraction_type)

    if key in coords_cache:
        return coords_cache[key]  # Falls schon berechnet, aus Cache abrufen

    if extraction_type == "ca":
        coords = np.array([
            atom["coordinates"]
            for residue in chain["residues"]
            for atom in residue["atoms"]
            if atom["atom_name"] == "CA"
        ])
    elif extraction_type == "nn_heavy":
        coords = np.array([
            get_nearest_heavy_atom(residue, ca_atom["coordinates"])
            for residue in chain["residues"]
            for ca_atom in residue["atoms"]
            if ca_atom["atom_name"] == "CA"
        ])
    elif extraction_type == "all_atoms":
        coords = np.array([
            atom["coordinates"]
            for residue in chain["residues"]
            for atom in residue["atoms"]
        ])
    else:
        return np.array([])

    coords_cache[key] = coords if coords.size > 0 else np.array([])  # Speichern im Cache
    return coords_cache[key]


def get_or_create_tree(chain, extraction_type):
    """Erstellt oder ruft einen gespeicherten cKDTree für eine bestimmte Extraktion ab."""
    key = (chain["unique_chain_id"], extraction_type)

    if key in tree_cache:
        return tree_cache[key]  # Falls bereits vorhanden, aus Cache abrufen

    points = extract_coordinates(chain, extraction_type)
    if points.size > 0:
        tree_cache[key] = cKDTree(points)  # KD-Tree nur einmal erstellen und speichern
        return tree_cache[key]

    return None


def calculate_distances_with_ckdtree(combined_data):
    """Berechnet interkettige Distanzen effizient mit cKDTree."""
    results = []

    for file_data in combined_data:
        file_path = file_data["file_path"]
        atom_data = file_data["atom_data"]

        # Initialisiere KD-Trees für jede Kette nur einmal
        local_trees = {
            (chain["unique_chain_id"], t): get_or_create_tree(chain, t)
            for chain in atom_data
            for t in ["ca", "nn_heavy", "all_atoms"]
        }

        for chain_a, chain_b in itertools.combinations(atom_data, 2):
            if chain_a["chain_id"] == chain_b["chain_id"]:
                continue  # Selbst-Interaktion vermeiden

            molecule_type_a = chain_a.get("molecule_type", "Unknown")
            molecule_type_b = chain_b.get("molecule_type", "Unknown")

            if molecule_type_a == "Unknown" or molecule_type_b == "Unknown":
                continue  # Falls der Typ unbekannt ist, überspringen

            if molecule_type_a == "Protein" and molecule_type_b == "Protein":
                interaction_type = "Protein-Protein"
            elif "Nucleic Acid" in [molecule_type_a, molecule_type_b]:
                interaction_type = "Protein-Nucleic Acid" if "Protein" in [molecule_type_a,
                                                                           molecule_type_b] else "Nucleic Acid-Nucleic Acid"
            else:
                interaction_type = "Unknown"

            # Falls eine Kette eine Nukleinsäure ist, verwende nur `all_atoms`
            if "Nucleic Acid" in [molecule_type_a, molecule_type_b]:
                tree_a = local_trees.get((chain_a["unique_chain_id"], "all_atoms"))
                tree_b = local_trees.get((chain_b["unique_chain_id"], "all_atoms"))

                if not tree_a or not tree_b:
                    continue

                all_atoms_pairs = tree_a.query_ball_point(tree_b.data, r=5.0)
                all_atoms_count = sum(len(p) for p in all_atoms_pairs)

                if all_atoms_count > 0:
                    results.append({
                        "file_path": file_path,
                        "chain_a": chain_a["unique_chain_id"],
                        "chain_b": chain_b["unique_chain_id"],
                        "all_atoms_count": all_atoms_count,
                        "interaction_type": interaction_type
                    })
                continue

            # Cα-zu-Cα-Distanzen berechnen
            tree_ca_a = local_trees.get((chain_a["unique_chain_id"], "ca"))
            tree_ca_b = local_trees.get((chain_b["unique_chain_id"], "ca"))

            if tree_ca_a and tree_ca_b:
                ca_pairs = tree_ca_a.query_ball_point(tree_ca_b.data, r=15.0)
                ca_count = sum(len(p) for p in ca_pairs)

                if ca_count >= 10:
                    tree_atoms_a = local_trees.get((chain_a["unique_chain_id"], "all_atoms"))
                    tree_atoms_b = local_trees.get((chain_b["unique_chain_id"], "all_atoms"))

                    if tree_atoms_a and tree_atoms_b:
                        all_atoms_pairs = tree_atoms_a.query_ball_point(tree_atoms_b.data, r=5.0)
                        all_atoms_count = sum(len(p) for p in all_atoms_pairs)

                        results.append({
                            "file_path": file_path,
                            "chain_a": chain_a["unique_chain_id"],
                            "chain_b": chain_b["unique_chain_id"],
                            "all_atoms_count": all_atoms_count,
                            "interaction_type": interaction_type
                        })
                    continue

            # Falls weniger als 10 Cα-zu-Cα-Distanzen unter 15 Å, prüfe NN-Heavy-Atome
            tree_nn_a = local_trees.get((chain_a["unique_chain_id"], "nn_heavy"))
            tree_nn_b = local_trees.get((chain_b["unique_chain_id"], "nn_heavy"))

            if tree_nn_a and tree_nn_b:
                nn_pairs = tree_nn_a.query_ball_point(tree_nn_b.data, r=15.0)
                nn_count = sum(len(p) for p in nn_pairs)

                if nn_count >= 10:
                    tree_atoms_a = local_trees.get((chain_a["unique_chain_id"], "all_atoms"))
                    tree_atoms_b = local_trees.get((chain_b["unique_chain_id"], "all_atoms"))

                    if tree_atoms_a and tree_atoms_b:
                        all_atoms_pairs = tree_atoms_a.query_ball_point(tree_atoms_b.data, r=5.0)
                        all_atoms_count = sum(len(p) for p in all_atoms_pairs)

                        results.append({
                            "file_path": file_path,
                            "chain_a": chain_a["unique_chain_id"],
                            "chain_b": chain_b["unique_chain_id"],
                            "all_atoms_count": all_atoms_count,
                            "interaction_type": interaction_type
                        })

    return results
