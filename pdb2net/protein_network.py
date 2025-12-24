"""Protein-level network construction.

This module builds protein-centric networks (UniProt nodes) from chain-level
interaction results and parsed structure metadata.

Modes
-----
- Per-PDB networks:
  * Protein nodes are colored uniformly as "Protein" (blue).
  * Nucleic acid (NA) chains (DNA/RNA/...) from the same PDB are added as
    neighbor nodes with their own molecule-type colors.
  * Tooltips are informative for both protein and NA nodes.
- Combined protein network:
  * Keeps the existing color logic (per-PDB vs Multi) and aggregates all
    protein–protein edges across PDBs.
- No duplicate exports are produced.

Notes
-----
- Labels remain unchanged (matching `create_cytoscape_network` behavior).
- The function expects `results` (chain-level interactions) and `combined_data`
  (per-PDB parsed chains with molecule typing and optional UniProt IDs).
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from cytoscape import create_cytoscape_network


def create_protein_network(
    results: List[Dict[str, Any]],
    combined_data: List[Dict[str, Any]],
    run_output_path: str,
    network_config: Dict[str, Any],
) -> None:
    """Create protein-level networks (per PDB and/or combined).

    Parameters
    ----------
    results : list[dict]
        Chain-level interaction records with at least:
        - "chain_a": "PDB:Chain"
        - "chain_b": "PDB:Chain"
        - "all_atoms_count": int
    combined_data : list[dict]
        Per-structure payloads as emitted by the parsing/processing stage.
        Each entry includes "pdb_id" and "atom_data" (list of chain dicts with
        "molecule_type", "molecule_name", "residues", and optional "uniprot_id").
    run_output_path : str
        Output directory passed through to Cytoscape export.
    network_config : dict
        Contains booleans:
        - "protein_per_pdb": build per-PDB protein networks
        - "combined_protein_network": build one combined protein network
    """
    make_per_pdb = network_config.get("protein_per_pdb", False)
    make_combined = network_config.get("combined_protein_network", False)
    if not make_per_pdb and not make_combined:
        print("Protein network creation is disabled.")
        return

    # --- Residue classes for quick length annotations in tooltips ---
    dna_set = {"DA", "DT", "DG", "DC", "DI"}
    rna_set = {"A", "U", "G", "C", "I"}
    aa_set = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
        "TYR", "VAL", "SEC", "PYL",
    }

    def count_lengths(res_list: Optional[List[Dict[str, Any]]]) -> Tuple[int, int]:
        """Return (aa_len, nt_len) for a residue list."""
        aa = nt = 0
        for res in (res_list or []):
            rn = (res.get("residue_name") or "").upper()
            if rn in aa_set:
                aa += 1
            elif rn in dna_set or rn in rna_set:
                nt += 1
        return aa, nt

    # --- Collect chain/protein metadata ---
    chain_to_uniprot: Dict[str, str] = {}            # "PDB:Chain" -> UniProt
    uniprot_to_pdb_ids: Dict[str, Set[str]] = {}     # UniProt -> {PDBs}
    uniprot_to_name: Dict[str, str] = {}             # UniProt -> display name
    pdb_to_uniprots: Dict[str, Set[str]] = {}        # PDB -> {UniProt-IDs}
    chain_info: Dict[str, Dict[str, Any]] = {}       # "PDB:Chain" -> info dict
    uniprot_pdb_to_chains: Dict[Tuple[str, str], Set[str]] = {}  # (UniProt, PDB) -> {ChainIDs}

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

    # --- Color-group logic for the combined network ---
    def get_color_group_for_combined(up_id: str) -> str:
        pdbs = uniprot_to_pdb_ids.get(up_id, set())
        return "Multi" if len(pdbs) != 1 else next(iter(pdbs))

    # --- Tooltips ---
    def protein_tooltip(up_id: str, pdb_scope: Optional[str]) -> str:
        """Format a tooltip for a protein node.

        pdb_scope = None
            Combined network: show number/examples of PDBs.
        pdb_scope = "1TUP"
            Per-PDB: show PDB, chains, and summed AA length across chains.
        """
        lines = [str(uniprot_to_name.get(up_id, up_id))]
        lines.append("Type: Protein")
        lines.append(f"UniProt: {up_id}")

        pdbs = sorted(uniprot_to_pdb_ids.get(up_id, []))
        if pdb_scope:
            chains = sorted(uniprot_pdb_to_chains.get((up_id, pdb_scope), []))
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
                sample = ", ".join(pdbs[:5])
                more = "" if len(pdbs) <= 5 else f" (+{len(pdbs) - 5} more)"
                lines.append(f"Examples: {sample}{more}")

        return "\n".join(lines)

    def na_tooltip(na_uid: str) -> str:
        """Format a tooltip for a nucleic-acid chain node."""
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

    # --- Aggregate edges ---
    pp_edges_by_pdb: Dict[str, Dict[Tuple[str, str], int]] = {}  # pdb_id -> {(up_a, up_b): weight}
    pp_edges_combined: Dict[Tuple[str, str], int] = {}           # (up_a, up_b) -> weight

    na_types = {"DNA", "RNA", "DNA/RNA", "Nucleic Acid"}
    na_edges_by_pdb: Dict[str, Dict[Tuple[str, str], int]] = {}  # pdb_id -> {(uniprot, na_uid): weight}
    na_nodes_by_pdb: Dict[str, Dict[str, Dict[str, Any]]] = {}   # pdb_id -> {na_uid: node_dict}

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

        # Protein–NA (exactly one protein + one NA)
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
                    "tooltip": na_tooltip(na_uid),
                }
            na_edges_by_pdb.setdefault(pdb_id, {})
            k = (up, na_uid)
            na_edges_by_pdb[pdb_id][k] = na_edges_by_pdb[pdb_id].get(k, 0) + cnt

    # --- Node builders ---
    def nodes_from_uniprots(
        uniprot_ids: Iterable[str],
        force_protein_color: bool = False,
        pdb_scope: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        nodes: List[Dict[str, Any]] = []
        for up in sorted(uniprot_ids):
            cg = "Protein" if force_protein_color else get_color_group_for_combined(up)
            nodes.append({
                "id": up,
                "color_group": cg,
                "molecule_name": uniprot_to_name.get(up, up),
                "tooltip": protein_tooltip(up, pdb_scope=pdb_scope),
                # No explicit "name": create_cytoscape_network uses 'id' as label (legacy style)
            })
        return nodes

    # --- Per-PDB networks: proteins in blue + NA neighbors ---
    if make_per_pdb:
        pdbs_with_edges = set(pp_edges_by_pdb.keys()) | set(na_edges_by_pdb.keys())

        # 1) PDBs WITH edges
        for pdb_id in sorted(pdbs_with_edges):
            prot_set = pdb_to_uniprots.get(pdb_id, set())
            if not prot_set:
                continue

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

            nodes = nodes_from_uniprots(prot_set, force_protein_color=True, pdb_scope=pdb_id)
            nodes.extend(list(na_nodes_by_pdb.get(pdb_id, {}).values()))
            create_cytoscape_network(edges, f"Protein_Network_{pdb_id}", run_output_path, nodes_data=nodes)

        # 2) PDBs WITHOUT edges → nodes-only protein graphs
        for pdb_id in sorted(set(pdb_to_uniprots.keys()) - pdbs_with_edges):
            prot_set = pdb_to_uniprots.get(pdb_id, set())
            if not prot_set:
                continue
            nodes = nodes_from_uniprots(prot_set, force_protein_color=True, pdb_scope=pdb_id)
            create_cytoscape_network([], f"Protein_Network_{pdb_id}", run_output_path, nodes_data=nodes)

    # --- Combined protein network ---
    if make_combined:
        combined_edges = [
            {"chain_a": up1, "chain_b": up2, "all_atoms_count": w}
            for (up1, up2), w in pp_edges_combined.items()
        ]
        all_uniprots = set(uniprot_to_name.keys())
        nodes = nodes_from_uniprots(all_uniprots, force_protein_color=False, pdb_scope=None)
        create_cytoscape_network(combined_edges, "Combined_Protein_Network", run_output_path, nodes_data=nodes)
