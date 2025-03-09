from scipy.spatial import cKDTree
import numpy as np
from config_loader import config
import itertools

# Cache für bereits erstellte KD-Trees
tree_cache = {}


def get_nearest_heavy_atom(residue, ca_coord):
    """
    Findet das nächstgelegene schwere Atom (außer Wasserstoff) für ein gegebenes Cα-Atom.
    Falls kein schweres Atom gefunden wird, werden die Cα-Koordinaten zurückgegeben.

    Args:
        residue (dict): Ein Dictionary, das Residue-Informationen enthält.
        ca_coord (list): Die Koordinaten des Cα-Atoms.

    Returns:
        list: Die Koordinaten des nächstgelegenen schweren Atoms oder die Cα-Koordinaten.
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


def extract_ca_nn(chain):
    """
    Extrahiert Cα-Atome und ihre nächstgelegenen schweren Atome aus einer Protein-Kette.

    Args:
        chain (dict): Dictionary mit Ketteninformationen.

    Returns:
        np.ndarray: Ein NumPy-Array mit Cα- und NN-Koordinaten.
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
    Extrahiert alle Atomkoordinaten aus einer Kette.

    Args:
        chain (dict): Dictionary mit Ketteninformationen.

    Returns:
        np.ndarray: Ein NumPy-Array mit allen Atom-Koordinaten.
    """
    coords = [atom["coordinates"] for residue in chain["residues"] for atom in residue["atoms"]]
    coords_array = np.array(coords)
    return coords_array.reshape(-1, 3) if coords_array.size > 0 else np.array([])


def get_or_create_tree(chain, use_ca_nn=True):
    """
    Gibt einen zwischengespeicherten cKDTree zurück oder erstellt einen neuen und speichert ihn im Cache.

    Args:
        chain (dict): Ketten-Informationen.
        use_ca_nn (bool): Falls True, wird der Tree für Cα-NN verwendet, sonst für alle Atome.

    Returns:
        cKDTree oder None: KDTree für die gegebenen Koordinaten oder None, falls leer.
    """
    key = (chain["unique_chain_id"], "ca_nn" if use_ca_nn else "all_atoms")
    if key in tree_cache:
        return tree_cache[key]

    points = extract_ca_nn(chain) if use_ca_nn else extract_all_atoms(chain)
    if points.size > 0:
        tree = cKDTree(points)
        tree_cache[key] = tree
        return tree
    return None


def calculate_distances_with_ckdtree(combined_data):
    """
    Berechnet interkettige atomare Distanzen mit cKDTree und speichert detaillierte Interaktionsdaten.

    Args:
        combined_data (list): Liste von Dictionaries mit PDB-Strukturen und Atomdaten.

    Returns:
        list: Eine Liste von Dictionaries mit paarweisen Interaktionsergebnissen.
    """
    results = []

    for file_data in combined_data:
        file_path = file_data["file_path"]
        atom_data = file_data["atom_data"]

        # Lokales Cache für diesen Lauf
        tree_cache_local = {
            (chain["unique_chain_id"], "ca_nn"): get_or_create_tree(chain, use_ca_nn=True) for chain in atom_data
        }
        tree_cache_local.update({
            (chain["unique_chain_id"], "all_atoms"): get_or_create_tree(chain, use_ca_nn=False) for chain in atom_data
        })

        # Nur einzigartige Kombinationen iterieren
        for chain_a, chain_b in itertools.combinations(atom_data, 2):
            if chain_a["chain_id"] == chain_b["chain_id"]:
                continue  # Selbst-Interaktion vermeiden

            molecule_type_a = chain_a.get("molecule_type", "Unknown")
            molecule_type_b = chain_b.get("molecule_type", "Unknown")

            # 🚨 Logging für Unknown-Ketten
            if molecule_type_a == "Unknown":
                print(f"⚠️ SKIPPING: Unknown molecule in {file_path}, Chain: {chain_a['chain_id']}")
            if molecule_type_b == "Unknown":
                print(f"⚠️ SKIPPING: Unknown molecule in {file_path}, Chain: {chain_b['chain_id']}")

            if molecule_type_a == "Unknown" or molecule_type_b == "Unknown":
                continue

            # Bestimme Interaktionstyp
            if molecule_type_a == "Protein" and molecule_type_b == "Protein":
                interaction_type = "Protein-Protein"
            elif (molecule_type_a == "Protein" and molecule_type_b == "Nucleic Acid") or \
                    (molecule_type_a == "Nucleic Acid" and molecule_type_b == "Protein"):
                interaction_type = "Protein-Nucleic Acid"
            elif molecule_type_a == "Nucleic Acid" and molecule_type_b == "Nucleic Acid":
                interaction_type = "Nucleic Acid-Nucleic Acid"
            else:
                interaction_type = "Unknown"

            # Hole gespeicherte KD-Trees
            tree_a = tree_cache_local.get((chain_a["unique_chain_id"], "ca_nn"))
            tree_b = tree_cache_local.get((chain_b["unique_chain_id"], "ca_nn"))

            if not tree_a or not tree_b:
                continue  # Falls eine Kette leer ist, überspringen

            # Berechnung mit `query_ball_point()` für detaillierte Paar-Informationen
            ca_nn_pairs = tree_a.query_ball_point(tree_b.data, r=15.0)

            # Berechnung der Anzahl aller Atome im Nahbereich
            tree_atoms_a = tree_cache_local[(chain_a["unique_chain_id"], "all_atoms")]
            tree_atoms_b = tree_cache_local[(chain_b["unique_chain_id"], "all_atoms")]
            all_atoms_close_pairs = tree_atoms_a.query_ball_point(tree_atoms_b.data, r=5.0)

            ca_nn_count = sum(len(p) for p in ca_nn_pairs)
            all_atoms_close_count = sum(len(p) for p in all_atoms_close_pairs)

            if ca_nn_count == 0 and all_atoms_close_count == 0:
                continue  # Falls wirklich keine Kontakte → Skip

            results.append({
                "file_path": file_path,
                "chain_a": chain_a["unique_chain_id"],
                "chain_b": chain_b["unique_chain_id"],
                "ca_nn_count": ca_nn_count,
                "ca_nn_pairs": ca_nn_pairs,  # 🔥 Enthält detaillierte Interaktionspaare
                "all_atoms_close_count": all_atoms_close_count,
                "interaction_type": interaction_type
            })

    return results
