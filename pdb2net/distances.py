from scipy.spatial import cKDTree
import numpy as np


def get_nearest_heavy_atom(residue, ca_coord):
    """
    Findet das nächstgelegene schwere Atom (außer Wasserstoff) für ein gegebenes Cα.
    Falls kein weiteres Atom gefunden wird, wird nur das Cα zurückgegeben.
    """
    heavy_atoms = [
        atom for atom in residue["atoms"]
        if atom["atom_name"] != "CA" and not atom["atom_name"].startswith("H")
    ]

    if not heavy_atoms:
        return ca_coord  # Falls kein schweres Atom da ist, verwende nur Cα

    nearest_atom = min(
        heavy_atoms,
        key=lambda atom: np.linalg.norm(np.array(atom["coordinates"]) - np.array(ca_coord))
    )
    return nearest_atom["coordinates"]


def extract_ca_nn(chain):
    """Extrahiert Cα-Atome und deren nächstgelegene schwere Atome aus einer Kette."""
    ca_nn = [
        [ca_atom["coordinates"], get_nearest_heavy_atom(residue, ca_atom["coordinates"])]
        for residue in chain["residues"]
        for ca_atom in residue["atoms"]
        if ca_atom["atom_name"] == "CA"
    ]
    return np.array(ca_nn).reshape(-1, 3) if ca_nn else np.array([])  # Sicherstellen, dass es (N,3) ist


def extract_all_atoms(chain):
    """Extrahiert alle Atome einer Kette als NumPy-Array (N,3)."""
    coords = [
        atom["coordinates"] for residue in chain["residues"] for atom in residue["atoms"]
    ]
    coords_array = np.array(coords)
    if coords_array.size == 0:
        print(f"⚠ WARNUNG: Keine Atome in Kette {chain['chain_id']}")
    return coords_array.reshape(-1, 3) if coords_array.size > 0 else np.array([])


def get_close_pairs(tree_a, points_b, radius):
    """Berechnet die Anzahl der Punkte in B, die innerhalb des Radius um A liegen."""
    if len(points_b) == 0:
        print(f"⚠ WARNUNG: Leere Punkte-Liste bei Distanzberechnung! Radius: {radius}")
        return 0
    close_points = tree_a.query_ball_point(points_b, r=radius)
    return sum(len(points) for points in close_points)


def calculate_distances_with_ckdtree(combined_data):
    """Berechnet Atomentfernungen zwischen Ketten und speichert die Ergebnisse."""
    results = []

    for file_data in combined_data:
        file_path = file_data["file_path"]
        atom_data = file_data["atom_data"]

        for i, chain_a in enumerate(atom_data):
            for j, chain_b in enumerate(atom_data):
                if i >= j:
                    continue

                print(f"\n🔍 Berechne Distanzen: {chain_a['unique_chain_id']} ↔ {chain_b['unique_chain_id']}")

                # Prüfe, ob mindestens eine der Ketten ein Protein ist
                protein_interaction = (
                    chain_a["molecule_type"] == "Protein" or chain_b["molecule_type"] == "Protein"
                )

                if protein_interaction:
                    # Protein-Protein oder Protein-DNA/RNA: Nutze Cα + NN Heavy Atom
                    ca_nn_a = extract_ca_nn(chain_a) if chain_a["molecule_type"] == "Protein" else extract_all_atoms(chain_a)
                    ca_nn_b = extract_ca_nn(chain_b) if chain_b["molecule_type"] == "Protein" else extract_all_atoms(chain_b)

                    ca_nn_a_flat = np.array([coord for coord in ca_nn_a]) if protein_interaction else ca_nn_a
                    ca_nn_b_flat = np.array([coord for coord in ca_nn_b]) if protein_interaction else ca_nn_b
                else:
                    # Nur Nukleinsäuren: Alle Atome verwenden
                    ca_nn_a_flat = extract_all_atoms(chain_a)
                    ca_nn_b_flat = extract_all_atoms(chain_b)

                # Debugging: Shape überprüfen
                print(f"  ✅ Shape von {chain_a['chain_id']} (A): {ca_nn_a_flat.shape}")
                print(f"  ✅ Shape von {chain_b['chain_id']} (B): {ca_nn_b_flat.shape}")

                # Shape-Validierung
                if ca_nn_a_flat.ndim != 2 or ca_nn_a_flat.shape[1] != 3:
                    print(f"❌ Fehler: Falsches Shape für {chain_a['chain_id']}: {ca_nn_a_flat.shape}")
                    continue
                if ca_nn_b_flat.ndim != 2 or ca_nn_b_flat.shape[1] != 3:
                    print(f"❌ Fehler: Falsches Shape für {chain_b['chain_id']}: {ca_nn_b_flat.shape}")
                    continue

                # Schritt 1: Berechne Distanzen mit cKDTree (Radius 15 Å)
                tree_a = cKDTree(ca_nn_a_flat)
                ca_nn_count = get_close_pairs(tree_a, ca_nn_b_flat, radius=15.0)

                if ca_nn_count < 10:
                    print(f"❌ Weniger als 10 Cα/NN-Paare mit Abstand <15 Å gefunden. Überspringe.")
                    continue

                # Schritt 2: Berechne alle Atome mit cKDTree (Radius 5 Å)
                all_atoms_a = extract_all_atoms(chain_a)
                all_atoms_b = extract_all_atoms(chain_b)

                if all_atoms_a.size == 0 or all_atoms_b.size == 0:
                    print(f"⚠ WARNUNG: Leere Atomliste für {chain_a['chain_id']} oder {chain_b['chain_id']}")
                    continue

                tree_atoms_a = cKDTree(all_atoms_a)
                all_atoms_close_count = get_close_pairs(tree_atoms_a, all_atoms_b, radius=5.0)

                if all_atoms_close_count < 10:
                    print(f"❌ Weniger als 10 Atome mit Abstand <5 Å gefunden. Überspringe.")
                    continue

                results.append({
                    "file_path": file_path,
                    "chain_a": chain_a["unique_chain_id"],
                    "chain_b": chain_b["unique_chain_id"],
                    "ca_nn_count": ca_nn_count,
                    "all_atoms_close_count": all_atoms_close_count,
                })

                print(f"✅ Distanzen berechnet: {ca_nn_count} Cα/NN-Paare, {all_atoms_close_count} Atompunkte <5 Å")

    return results
