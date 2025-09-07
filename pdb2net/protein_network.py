from cytoscape import create_cytoscape_network


def create_protein_network(results, combined_data, run_output_path, network_config):
    """
    Erzeugt Protein-Level Netzwerke (UniProt-Knoten) – pro PDB und/oder kombiniert.
    NEU: Auch wenn es keine Protein-Protein-Kanten gibt, wird ein Nodes-only Netzwerk gebaut
    (z. B. wenn eine PDB nur ein einziges UniProt-Protein enthält oder nur Homomere vorliegen).
    """
    if not network_config.get("protein_per_pdb", False) and not network_config.get("combined_protein_network", False):
        print("Protein network creation is disabled.")
        return

    # --------- Sammle UniProt-Informationen ----------
    chain_to_uniprot = {}       # "PDB:Chain" -> UniProt
    uniprot_to_pdb_ids = {}     # UniProt -> {PDBs}
    uniprot_to_name = {}        # UniProt -> schöner Name
    pdb_to_uniprots = {}        # PDB -> {UniProt-IDs in dieser PDB}

    for structure in combined_data:
        pdb_id = structure["pdb_id"]
        for chain in structure["atom_data"]:
            unique_chain_id = f"{pdb_id}:{chain['chain_id']}"
            uniprot_id = chain.get("uniprot_id")
            if not uniprot_id:
                continue
            chain_to_uniprot[unique_chain_id] = uniprot_id
            uniprot_to_name[uniprot_id] = chain.get("molecule_name", "Unknown")
            uniprot_to_pdb_ids.setdefault(uniprot_id, set()).add(pdb_id)
            pdb_to_uniprots.setdefault(pdb_id, set()).add(uniprot_id)

    # --------- Aggregiere Interaktionen auf UniProt-Ebene ----------
    protein_interactions = set()
    interaction_data = {}

    for entry in results:
        if entry.get("all_atoms_count", 0) <= 0:
            continue

        up_a = chain_to_uniprot.get(entry["chain_a"])
        up_b = chain_to_uniprot.get(entry["chain_b"])
        if not up_a or not up_b:
            continue

        # Homomere (up_a == up_b) werden als Kanten weggelassen; für Einzel-Protein-PDBs wollen wir Nodes-only
        if up_a == up_b:
            continue

        if network_config.get("protein_per_pdb", False):
            pdb_id = entry["chain_a"].split(":")[0]
            inter_key = (pdb_id, tuple(sorted([up_a, up_b])))
        else:
            inter_key = tuple(sorted([up_a, up_b]))

        if inter_key not in protein_interactions:
            protein_interactions.add(inter_key)
            interaction_data[inter_key] = {
                "uniprot_a": up_a,
                "uniprot_b": up_b,
                "all_atoms_count": entry["all_atoms_count"]
            }

    # --------- Hilfsfunktionen für Nodes ----------
    def get_color_group(uniprot_id: str) -> str:
        """Ein PDB-Code, wenn das UniProt nur in genau 1 PDB vorkommt; sonst 'Multi'."""
        pdbs = uniprot_to_pdb_ids.get(uniprot_id, set())
        return "Multi" if len(pdbs) != 1 else next(iter(pdbs))

    def nodes_from_uniprots(uniprot_ids):
        return [{
            "id": up,
            "color_group": get_color_group(up),
            "molecule_name": uniprot_to_name.get(up, up)
        } for up in sorted(uniprot_ids)]

    # --------- Netzwerke pro PDB ----------
    if network_config.get("protein_per_pdb", False):
        # Edges je PDB
        results_by_pdb = {}
        for (pdb_id, (a, b)), inter in interaction_data.items():
            results_by_pdb.setdefault(pdb_id, []).append({
                "chain_a": inter["uniprot_a"],
                "chain_b": inter["uniprot_b"],
                "all_atoms_count": inter["all_atoms_count"]
            })

        # 1) PDBs MIT Kanten: nutze alle in der PDB vorkommenden UniProts als Nodes
        for pdb_id, edges in results_by_pdb.items():
            node_set = pdb_to_uniprots.get(pdb_id, set())
            nodes = nodes_from_uniprots(node_set)
            network_title = f"Protein_Network_{pdb_id}"
            create_cytoscape_network(edges, network_title, run_output_path, nodes_data=nodes)

        # 2) PDBs OHNE Kanten (z. B. nur 1 UniProt oder nur Homomere) -> Nodes-only
        for pdb_id, node_set in pdb_to_uniprots.items():
            if pdb_id not in results_by_pdb:
                if not node_set:
                    continue
                nodes = nodes_from_uniprots(node_set)
                network_title = f"Protein_Network_{pdb_id}"
                create_cytoscape_network([], network_title, run_output_path, nodes_data=nodes)

    # --------- Kombiniertes Protein-Netzwerk ----------
    if network_config.get("combined_protein_network", False):
        combined_edges = []
        for inter in interaction_data.values():
            combined_edges.append({
                "chain_a": inter["uniprot_a"],
                "chain_b": inter["uniprot_b"],
                "all_atoms_count": inter["all_atoms_count"]
            })

        # Alle UniProts, die überhaupt im Datensatz vorkommen – auch wenn sie keine Kante haben
        all_uniprots = set(uniprot_to_name.keys())
        nodes = nodes_from_uniprots(all_uniprots)
        create_cytoscape_network(combined_edges, "Combined_Protein_Network", run_output_path, nodes_data=nodes)
