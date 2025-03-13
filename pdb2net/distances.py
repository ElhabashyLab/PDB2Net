from scipy.spatial import cKDTree
import numpy as np
import itertools

# Cache für bereits erstellte KD-Trees
tree_cache = {}


def get_nearest_heavy_atom(residue, ca_coord):
    """
    Findet das nächstgelegene schwere Atom (außer Wasserstoff) für ein gegebenes Cα-Atom.
    Falls kein schweres Atom gefunden wird, werden die Cα-Koordinaten zurückgegeben.
    """
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


def extract_ca_coords(chain):
    """ Extrahiert Cα-Atome aus einer Protein-Kette. """
    ca_coords = [
        atom["coordinates"]
        for residue in chain["residues"]
        for atom in residue["atoms"]
        if atom["atom_name"] == "CA"
    ]
    return np.array(ca_coords) if ca_coords else np.array([])


def extract_nn_heavy_atoms(chain):
    """ Extrahiert die nächstgelegenen schweren Atome zu den Cα-Atomen. """
    nn_heavy_atoms = [
        get_nearest_heavy_atom(residue, ca_atom["coordinates"])
        for residue in chain["residues"]
        for ca_atom in residue["atoms"]
        if ca_atom["atom_name"] == "CA"
    ]
    return np.array(nn_heavy_atoms) if nn_heavy_atoms else np.array([])


def extract_all_atoms(chain):
    """ Extrahiert alle Atomkoordinaten aus einer Kette. """
    coords = [atom["coordinates"] for residue in chain["residues"] for atom in residue["atoms"]]
    return np.array(coords) if coords else np.array([])


def get_or_create_tree(chain, extraction_func, cache_key):
    """
    Erstellt oder ruft einen gespeicherten cKDTree für eine bestimmte Extraktionsfunktion ab.
    """
    key = (chain["unique_chain_id"], cache_key)
    if key in tree_cache:
        return tree_cache[key]

    points = extraction_func(chain)
    if points.size > 0:
        tree = cKDTree(points)
        tree_cache[key] = tree
        return tree
    return None


def calculate_distances_with_ckdtree(combined_data):
    """ Berechnet interkettige Distanzen basierend auf den neuen Kriterien. """
    results = []

    for file_data in combined_data:
        file_path = file_data["file_path"]
        atom_data = file_data["atom_data"]

        tree_cache_local = {
            (chain["unique_chain_id"], "ca"): get_or_create_tree(chain, extract_ca_coords, "ca") for chain in atom_data
        }
        tree_cache_local.update({
            (chain["unique_chain_id"], "nn_heavy"): get_or_create_tree(chain, extract_nn_heavy_atoms, "nn_heavy")
            for chain in atom_data
        })
        tree_cache_local.update({
            (chain["unique_chain_id"], "all_atoms"): get_or_create_tree(chain, extract_all_atoms, "all_atoms")
            for chain in atom_data
        })

        for chain_a, chain_b in itertools.combinations(atom_data, 2):
            if chain_a["chain_id"] == chain_b["chain_id"]:
                continue  # Selbst-Interaktion vermeiden

            # **Interaktionstyp bestimmen**
            molecule_type_a = chain_a.get("molecule_type", "Unknown")
            molecule_type_b = chain_b.get("molecule_type", "Unknown")

            if molecule_type_a == "Unknown" or molecule_type_b == "Unknown":
                continue  # Falls der Typ unbekannt ist, überspringen

            if molecule_type_a == "Protein" and molecule_type_b == "Protein":
                interaction_type = "Protein-Protein"
            elif "Nucleic Acid" in [molecule_type_a, molecule_type_b]:
                interaction_type = "Protein-Nucleic Acid" if "Protein" in [molecule_type_a, molecule_type_b] else "Nucleic Acid-Nucleic Acid"
            else:
                interaction_type = "Unknown"

            # 🔹 Falls eine Kette eine Nukleinsäure ist, verwende nur `extract_all_atoms()`
            if "Nucleic Acid" in [molecule_type_a, molecule_type_b]:
                tree_atoms_a = tree_cache_local.get((chain_a["unique_chain_id"], "all_atoms"))
                tree_atoms_b = tree_cache_local.get((chain_b["unique_chain_id"], "all_atoms"))

                if not tree_atoms_a or not tree_atoms_b:
                    continue  # Falls keine Atome gefunden wurden, ignorieren

                all_atoms_pairs = tree_atoms_a.query_ball_point(tree_atoms_b.data, r=5.0)
                all_atoms_count = sum(len(p) for p in all_atoms_pairs)

                if all_atoms_count > 0:
                    results.append({
                        "file_path": file_path,
                        "chain_a": chain_a["unique_chain_id"],
                        "chain_b": chain_b["unique_chain_id"],
                        "all_atoms_count": all_atoms_count,
                        "interaction_type": interaction_type
                    })
                continue  # Überspringe die Cα- und NN-Prüfung

            # Cα-zu-Cα-Distanzen berechnen
            tree_ca_a = tree_cache_local.get((chain_a["unique_chain_id"], "ca"))
            tree_ca_b = tree_cache_local.get((chain_b["unique_chain_id"], "ca"))

            if not tree_ca_a or not tree_ca_b:
                continue

            ca_pairs = tree_ca_a.query_ball_point(tree_ca_b.data, r=15.0)
            ca_count = sum(len(p) for p in ca_pairs)

            if ca_count >= 10:
                tree_atoms_a = tree_cache_local.get((chain_a["unique_chain_id"], "all_atoms"))
                tree_atoms_b = tree_cache_local.get((chain_b["unique_chain_id"], "all_atoms"))
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
            tree_nn_a = tree_cache_local.get((chain_a["unique_chain_id"], "nn_heavy"))
            tree_nn_b = tree_cache_local.get((chain_b["unique_chain_id"], "nn_heavy"))

            if not tree_nn_a or not tree_nn_b:
                continue

            nn_pairs = tree_nn_a.query_ball_point(tree_nn_b.data, r=15.0)
            nn_count = sum(len(p) for p in nn_pairs)

            if nn_count >= 10:
                tree_atoms_a = tree_cache_local.get((chain_a["unique_chain_id"], "all_atoms"))
                tree_atoms_b = tree_cache_local.get((chain_b["unique_chain_id"], "all_atoms"))
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
