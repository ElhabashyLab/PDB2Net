from scipy.spatial import cKDTree
import numpy as np


def calculate_distances_with_ckdtree(combined_data):
    """
    Geht die kombinierte Datenstruktur durch, berechnet Atomentfernungen mithilfe von cKDTree
    und filtert Kettenpaare, die bestimmte Distanzkriterien erfüllen.
    """
    results = []

    for file_data in combined_data:
        file_path = file_data["file_path"]
        atom_data = file_data["atom_data"]

        # Liste aller Ketten durchgehen
        for i, chain_a in enumerate(atom_data):
            for j, chain_b in enumerate(atom_data):
                if i >= j:
                    continue  # Nur einmal jedes Paar prüfen

                # Extrahiere die Cα- und Cβ-Koordinaten der beiden Ketten
                ca_cb_a = [
                    atom["coordinates"]
                    for residue in chain_a["residues"]
                    for atom in residue["atoms"]
                    if atom["atom_name"] in {"CA", "CB"}
                ]
                ca_cb_b = [
                    atom["coordinates"]
                    for residue in chain_b["residues"]
                    for atom in residue["atoms"]
                    if atom["atom_name"] in {"CA", "CB"}
                ]

                # Validierung: Listen dürfen nicht leer sein
                if not ca_cb_a or not ca_cb_b:
                    print(f"WARNUNG: Leere Koordinatenliste in Datei {file_path}, Ketten {chain_a['chain_id']} oder {chain_b['chain_id']}")
                    continue

                # Schritt 1: Berechne Distanzen von Cα und Cβ mit cKDTree
                tree_a = cKDTree(ca_cb_a)
                close_points = tree_a.query_ball_point(ca_cb_b, r=15.0)  # Radius: 15 Å
                ca_cb_count = sum(len(points) for points in close_points)

                # Schritt 2: Prüfe auf mindestens 10 Cα/Cβ mit Abstand <15 Å
                if ca_cb_count >= 10:
                    # Alle Atome extrahieren
                    all_atoms_a = np.array([
                        atom["coordinates"]
                        for residue in chain_a["residues"]
                        for atom in residue["atoms"]
                    ])
                    all_atoms_b = np.array([
                        atom["coordinates"]
                        for residue in chain_b["residues"]
                        for atom in residue["atoms"]
                    ])

                    # Validierung: Listen dürfen nicht leer sein
                    if all_atoms_a.size == 0 or all_atoms_b.size == 0:
                        print(f"WARNUNG: Keine Atome in Ketten {chain_a['chain_id']} oder {chain_b['chain_id']}")
                        continue

                    # Schritt 3: Berechne Distanzen aller Atome mit cKDTree
                    tree_atoms_a = cKDTree(all_atoms_a)
                    close_points_atoms = tree_atoms_a.query_ball_point(all_atoms_b, r=5.0)  # Radius: 5 Å
                    all_atoms_close_count = sum(len(points) for points in close_points_atoms)

                    # Schritt 4: Prüfe auf mindestens 10 Atome mit Abstand <5 Å
                    if all_atoms_close_count >= 10:
                        results.append({
                            "file_path": file_path,
                            "chain_a": chain_a["unique_chain_id"],
                            "chain_b": chain_b["unique_chain_id"],
                            "ca_cb_count": ca_cb_count,
                            "all_atoms_close_count": all_atoms_close_count,
                        })

    return results
