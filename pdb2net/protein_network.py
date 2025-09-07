from cytoscape import create_cytoscape_network

def create_protein_network(results, combined_data, run_output_path, network_config):
    """
    Build protein-level networks (UniProt nodes) per PDB and/or as a combined network.

    - Per-PDB protein networks: force UniProt nodes to color_group="Protein" (blue).
    - Include interacting nucleic-acid chains (DNA/RNA/DNA-RNA/Nucleic Acid) as extra nodes,
      colored by their molecule_type (same palette as chain networks).
    - Avoid duplicate exports (no overwrite dialog).
    """
    make_per_pdb = network_config.get("protein_per_pdb", False)
    make_combined = network_config.get("combined_protein_network", False)

    if not make_per_pdb and not make_combined:
        print("Protein network creation is disabled.")
        return

    # --------- Collect UniProt and chain metadata ----------
    chain_to_uniprot = {}       # "PDB:Chain" -> UniProt
    uniprot_to_pdb_ids = {}     # UniProt -> {PDBs}
    uniprot_to_name = {}        # UniProt -> display name
    pdb_to_uniprots = {}        # PDB -> {UniProt-IDs present}
    chain_info = {}             # "PDB:Chain" -> {"molecule_type": str, "molecule_name": str}

    for structure in combined_data:
        pdb_id = structure["pdb_id"]
        for chain in structure["atom_data"]:
            uid = f"{pdb_id}:{chain['chain_id']}"
            chain_info[uid] = {
                "molecule_type": chain.get("molecule_type", "Unknown"),
                "molecule_name": chain.get("molecule_name", "Unknown"),
            }
            up = chain.get("uniprot_id")
            if not up:
                continue
            chain_to_uniprot[uid] = up
            uniprot_to_name[up] = chain.get("molecule_name", up)
            uniprot_to_pdb_ids.setdefault(up, set()).add(pdb_id)
            pdb_to_uniprots.setdefault(pdb_id, set()).add(up)

    # Helpers
    def get_color_group_for_combined(up_id: str) -> str:
        """For the combined network keep PDB/Multi coloring."""
        pdbs = uniprot_to_pdb_ids.get(up_id, set())
        return "Multi" if len(pdbs) != 1 else next(iter(pdbs))

    def nodes_from_uniprots(uniprot_ids, force_protein_color=False):
        nodes = []
        for up in sorted(uniprot_ids):
            cg = "Protein" if force_protein_color else get_color_group_for_combined(up)
            nodes.append({
                "id": up,
                "color_group": cg,
                "molecule_name": uniprot_to_name.get(up, up)
            })
        return nodes

    # --------- Aggregate edges ----------
    # Protein–Protein edges per PDB and for combined
    pp_edges_by_pdb = {}        # pdb_id -> {(up_a, up_b): weight}
    pp_edges_combined = {}      # (up_a, up_b) -> weight

    # Protein–Nucleic Acid edges per PDB + NA nodes per PDB
    na_types = {"DNA", "RNA", "DNA/RNA", "Nucleic Acid"}
    na_edges_by_pdb = {}        # pdb_id -> {(uniprot, na_uid): weight}
    na_nodes_by_pdb = {}        # pdb_id -> {na_uid: node_dict}

    for entry in results:
        cnt = entry.get("all_atoms_count", 0)
        if cnt <= 0:
            continue

        a, b = entry["chain_a"], entry["chain_b"]
        # PDB-ID aus A ableiten (Ergebnisse stammen i.d.R. aus einer Struktur)
        pdb_id = a.split(":")[0] if ":" in a else b.split(":")[0]

        ia = chain_info.get(a, {})
        ib = chain_info.get(b, {})
        ta = (ia.get("molecule_type") or "Unknown").strip()
        tb = (ib.get("molecule_type") or "Unknown").strip()

        up_a = chain_to_uniprot.get(a)
        up_b = chain_to_uniprot.get(b)

        # --- Protein–Protein ---
        if up_a and up_b and up_a != up_b:
            key = tuple(sorted((up_a, up_b)))
            pp_edges_combined[key] = pp_edges_combined.get(key, 0) + cnt
            pp_edges_by_pdb.setdefault(pdb_id, {})
            pp_edges_by_pdb[pdb_id][key] = pp_edges_by_pdb[pdb_id].get(key, 0) + cnt

        # --- Protein–NA (genau ein Protein + eine NA) ---
        if ta == "Protein" and tb in na_types:
            up, na_uid, na_type, na_name = up_a, b, tb, ib.get("molecule_name", "Unknown")
        elif tb == "Protein" and ta in na_types:
            up, na_uid, na_type, na_name = up_b, a, ta, ia.get("molecule_name", "Unknown")
        else:
            up = None

        if up:
            # Falls Proteinchain keine UniProt-ID hat, können wir nicht ins Protein-Netz integrieren
            if not up:
                pass
            else:
                na_nodes_by_pdb.setdefault(pdb_id, {})
                if na_uid not in na_nodes_by_pdb[pdb_id]:
                    na_nodes_by_pdb[pdb_id][na_uid] = {
                        "id": na_uid,
                        "color_group": na_type,
                        "molecule_name": na_name
                    }
                na_edges_by_pdb.setdefault(pdb_id, {})
                k = (up, na_uid)
                na_edges_by_pdb[pdb_id][k] = na_edges_by_pdb[pdb_id].get(k, 0) + cnt

    # --------- PER-PDB protein networks (with NA neighbors) ----------
    if make_per_pdb:
        # PDBs, in denen es überhaupt Kanten gibt (PP oder P–NA)
        pdbs_with_edges = set(pp_edges_by_pdb.keys()) | set(na_edges_by_pdb.keys())

        # 1) PDBs MIT Kanten
        for pdb_id in sorted(pdbs_with_edges):
            prot_set = pdb_to_uniprots.get(pdb_id, set())
            if not prot_set:
                continue

            # Protein–Protein-Kanten
            edges_pp = []
            for (up1, up2), w in pp_edges_by_pdb.get(pdb_id, {}).items():
                edges_pp.append({"chain_a": up1, "chain_b": up2, "all_atoms_count": w})

            # Protein–NA-Kanten
            edges_pna = []
            for (up, na_uid), w in na_edges_by_pdb.get(pdb_id, {}).items():
                edges_pna.append({"chain_a": up, "chain_b": na_uid, "all_atoms_count": w})

            edges = edges_pp + edges_pna
            if not edges:
                continue  # Sicherheitsgurt

            # Nodes: Proteine (immer blau) + NA-Nodes (farben nach Typ)
            nodes = nodes_from_uniprots(prot_set, force_protein_color=True)
            nodes.extend(list(na_nodes_by_pdb.get(pdb_id, {}).values()))

            create_cytoscape_network(edges, f"Protein_Network_{pdb_id}", run_output_path, nodes_data=nodes)

        # 2) PDBs OHNE Kanten -> nur Protein-Knoten (blau)
        for pdb_id in sorted(set(pdb_to_uniprots.keys()) - pdbs_with_edges):
            prot_set = pdb_to_uniprots.get(pdb_id, set())
            if not prot_set:
                continue
            nodes = nodes_from_uniprots(prot_set, force_protein_color=True)
            create_cytoscape_network([], f"Protein_Network_{pdb_id}", run_output_path, nodes_data=nodes)

    # --------- Combined protein network (unchanged coloring logic) ----------
    if make_combined:
        combined_edges = [
            {"chain_a": up1, "chain_b": up2, "all_atoms_count": w}
            for (up1, up2), w in pp_edges_combined.items()
        ]
        all_uniprots = set(uniprot_to_name.keys())
        nodes = nodes_from_uniprots(all_uniprots, force_protein_color=False)
        create_cytoscape_network(combined_edges, "Combined_Protein_Network", run_output_path, nodes_data=nodes)
