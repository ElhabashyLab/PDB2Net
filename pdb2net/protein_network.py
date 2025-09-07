from cytoscape import create_cytoscape_network

def create_protein_network(results, combined_data, run_output_path, network_config):
    """
    Protein-level networks (UniProt-Knoten). 
    - Per-PDB: Protein-Knoten blau (color_group="Protein"); zusätzlich NA-Ketten (DNA/RNA/...) als Nachbarn.
    - Combined: unverändert (PDB/Multi-Farb-Logik).
    - NEU: Tooltips für Protein- und NA-Knoten (Labels bleiben unverändert = Standard aus create_cytoscape_network).
    - Kein doppelter Export.
    """
    make_per_pdb = network_config.get("protein_per_pdb", False)
    make_combined = network_config.get("combined_protein_network", False)
    if not make_per_pdb and not make_combined:
        print("Protein network creation is disabled.")
        return

    # --- Residue-Klassen für Längenangaben in Tooltips ---
    dna_set = {"DA", "DT", "DG", "DC", "DI"}
    rna_set = {"A", "U", "G", "C", "I"}
    aa_set = {
        "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS",
        "ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP",
        "TYR","VAL","SEC","PYL"
    }

    def count_lengths(res_list):
        aa = nt = 0
        for res in (res_list or []):
            rn = (res.get("residue_name") or "").upper()
            if rn in aa_set:
                aa += 1
            elif rn in dna_set or rn in rna_set:
                nt += 1
        return aa, nt

    # --- Metadaten sammeln ---
    chain_to_uniprot = {}              # "PDB:Chain" -> UniProt
    uniprot_to_pdb_ids = {}            # UniProt -> {PDBs}
    uniprot_to_name = {}               # UniProt -> display name (voll)
    pdb_to_uniprots = {}               # PDB -> {UniProt-IDs}
    chain_info = {}                    # "PDB:Chain" -> dict(type, name, aa_len, nt_len, pdb_id, chain_id)
    uniprot_pdb_to_chains = {}         # (UniProt, PDB) -> {ChainIDs}

    for structure in combined_data:
        pdb_id = structure["pdb_id"]
        for chain in structure["atom_data"]:
            chain_id = chain["chain_id"]
            uid = f"{pdb_id}:{chain_id}"
            aa_len, nt_len = count_lengths(chain.get("residues"))
            chain_info[uid] = {
                "molecule_type": (chain.get("molecule_type") or "Unknown").strip(),
                "molecule_name": chain.get("molecule_name", "Unknown"),
                "aa_len": aa_len,
                "nt_len": nt_len,
                "pdb_id": pdb_id,
                "chain_id": chain_id,
            }
            up = chain.get("uniprot_id")
            if up:
                chain_to_uniprot[uid] = up
                uniprot_to_name[up] = chain.get("molecule_name", up)
                uniprot_to_pdb_ids.setdefault(up, set()).add(pdb_id)
                pdb_to_uniprots.setdefault(pdb_id, set()).add(up)
                uniprot_pdb_to_chains.setdefault((up, pdb_id), set()).add(chain_id)

    # --- Helper für Farb-Logik im Combined-Netz ---
    def get_color_group_for_combined(up_id: str) -> str:
        pdbs = uniprot_to_pdb_ids.get(up_id, set())
        return "Multi" if len(pdbs) != 1 else next(iter(pdbs))

    # --- Tooltips bauen ---
    def protein_tooltip(up_id: str, pdb_scope: str | None) -> str:
        """
        pdb_scope=None => Combined-Netz: zeige Anzahl/Beispiele PDBs
        pdb_scope='1TUP' => per-PDB: zeige PDB und Ketten, Länge (aa) summiert über zugeordnete Ketten
        """
        lines = [str(uniprot_to_name.get(up_id, up_id))]
        lines.append("Type: Protein")
        lines.append(f"UniProt: {up_id}")

        pdbs = sorted(uniprot_to_pdb_ids.get(up_id, []))
        if pdb_scope:
            chains = sorted(uniprot_pdb_to_chains.get((up_id, pdb_scope), []))
            # Länge als Summe AA-Längen über alle Ketten dieses Proteins in der PDB
            aa_total = 0
            for ch in chains:
                info = chain_info.get(f"{pdb_scope}:{ch}", {})
                aa_total += int(info.get("aa_len", 0))
            lines.append(f"PDB: {pdb_scope}")
            if chains:
                lines.append(f"Chains: {', '.join(chains)}")
            if aa_total:
                lines.append(f"Length: {aa_total} aa")
        else:
            lines.append(f"PDBs: {len(pdbs)}")
            if pdbs:
                # nur wenige Beispiele zeigen, um Tooltips kurz zu halten
                sample = ", ".join(pdbs[:5])
                more = "" if len(pdbs) <= 5 else f" (+{len(pdbs)-5} more)"
                lines.append(f"Examples: {sample}{more}")

        return "\n".join(lines)

    def na_tooltip(na_uid: str) -> str:
        info = chain_info.get(na_uid, {})
        name = info.get("molecule_name", "Unknown")
        mtype = info.get("molecule_type", "Unknown")
        aa_len = info.get("aa_len", 0)
        nt_len = info.get("nt_len", 0)
        pdb_id = info.get("pdb_id", "")
        chain_id = info.get("chain_id", "")
        lines = [str(name)]
        lines.append(f"Type: {mtype}")
        if aa_len:
            lines.append(f"Length: {aa_len} aa")
        if nt_len:
            lines.append(f"Length: {nt_len} nt")
        if pdb_id:
            lines.append(f"PDB: {pdb_id}:{chain_id}")
        return "\n".join(lines)

    # --- Kanten aggregieren ---
    pp_edges_by_pdb = {}     # pdb_id -> {(up_a, up_b): weight}
    pp_edges_combined = {}   # (up_a, up_b) -> weight
    na_types = {"DNA", "RNA", "DNA/RNA", "Nucleic Acid"}
    na_edges_by_pdb = {}     # pdb_id -> {(uniprot, na_uid): weight}
    na_nodes_by_pdb = {}     # pdb_id -> {na_uid: node_dict}

    for entry in results:
        cnt = entry.get("all_atoms_count", 0)
        if cnt <= 0:
            continue

        a, b = entry["chain_a"], entry["chain_b"]
        pdb_id = a.split(":")[0] if ":" in a else b.split(":")[0]

        ia = chain_info.get(a, {})
        ib = chain_info.get(b, {})
        ta = (ia.get("molecule_type") or "Unknown").strip()
        tb = (ib.get("molecule_type") or "Unknown").strip()

        up_a = chain_to_uniprot.get(a)
        up_b = chain_to_uniprot.get(b)

        # Protein–Protein
        if up_a and up_b and up_a != up_b:
            key = tuple(sorted((up_a, up_b)))
            pp_edges_combined[key] = pp_edges_combined.get(key, 0) + cnt
            pp_edges_by_pdb.setdefault(pdb_id, {})
            pp_edges_by_pdb[pdb_id][key] = pp_edges_by_pdb[pdb_id].get(key, 0) + cnt

        # Protein–NA (genau ein Protein + eine NA)
        up = None
        na_uid = None
        na_type = None
        if ta == "Protein" and tb in na_types:
            up, na_uid, na_type = up_a, b, tb
        elif tb == "Protein" and ta in na_types:
            up, na_uid, na_type = up_b, a, ta

        if up and na_uid:
            na_nodes_by_pdb.setdefault(pdb_id, {})
            if na_uid not in na_nodes_by_pdb[pdb_id]:
                na_nodes_by_pdb[pdb_id][na_uid] = {
                    "id": na_uid,
                    "color_group": na_type,
                    "molecule_name": chain_info.get(na_uid, {}).get("molecule_name", "Unknown"),
                    "tooltip": na_tooltip(na_uid)   # NEU: informativer Tooltip
                }
            na_edges_by_pdb.setdefault(pdb_id, {})
            k = (up, na_uid)
            na_edges_by_pdb[pdb_id][k] = na_edges_by_pdb[pdb_id].get(k, 0) + cnt

    # --- Node-Helfer ---
    def nodes_from_uniprots(uniprot_ids, force_protein_color=False, pdb_scope=None):
        nodes = []
        for up in sorted(uniprot_ids):
            cg = "Protein" if force_protein_color else get_color_group_for_combined(up)
            nodes.append({
                "id": up,
                "color_group": cg,
                "molecule_name": uniprot_to_name.get(up, up),
                "tooltip": protein_tooltip(up, pdb_scope=pdb_scope)  # NEU: informativer Tooltip
                # kein "name": create_cytoscape_network setzt Label=ID (alte Darstellung)
            })
        return nodes

    # --- PER-PDB: mit NA-Nachbarn, Proteine blau ---
    if make_per_pdb:
        pdbs_with_edges = set(pp_edges_by_pdb.keys()) | set(na_edges_by_pdb.keys())

        # 1) PDBs MIT Kanten
        for pdb_id in sorted(pdbs_with_edges):
            prot_set = pdb_to_uniprots.get(pdb_id, set())
            if not prot_set:
                continue

            # Edges sammeln
            edges_pp = [
                {"chain_a": up1, "chain_b": up2, "all_atoms_count": w}
                for (up1, up2), w in pp_edges_by_pdb.get(pdb_id, {}).items()
            ]
            edges_pna = [
                {"chain_a": up, "chain_b": na_uid, "all_atoms_count": w}
                for (up, na_uid), w in na_edges_by_pdb.get(pdb_id, {}).items()
            ]
            edges = edges_pp + edges_pna
            if not edges:
                continue

            # Nodes: Proteine (blau) + NA-Nodes
            nodes = nodes_from_uniprots(prot_set, force_protein_color=True, pdb_scope=pdb_id)
            nodes.extend(list(na_nodes_by_pdb.get(pdb_id, {}).values()))

            create_cytoscape_network(edges, f"Protein_Network_{pdb_id}", run_output_path, nodes_data=nodes)

        # 2) PDBs OHNE Kanten -> nur Protein-Knoten (blau)
        for pdb_id in sorted(set(pdb_to_uniprots.keys()) - pdbs_with_edges):
            prot_set = pdb_to_uniprots.get(pdb_id, set())
            if not prot_set:
                continue
            nodes = nodes_from_uniprots(prot_set, force_protein_color=True, pdb_scope=pdb_id)
            create_cytoscape_network([], f"Protein_Network_{pdb_id}", run_output_path, nodes_data=nodes)

    # --- COMBINED: unveränderte Farb-Logik, aber Tooltip enthalten ---
    if make_combined:
        combined_edges = [
            {"chain_a": up1, "chain_b": up2, "all_atoms_count": w}
            for (up1, up2), w in pp_edges_combined.items()
        ]
        all_uniprots = set(uniprot_to_name.keys())
        nodes = nodes_from_uniprots(all_uniprots, force_protein_color=False, pdb_scope=None)
        create_cytoscape_network(combined_edges, "Combined_Protein_Network", run_output_path, nodes_data=nodes)
