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

from typing import Any, Dict, Iterable, List, Mapping, Optional, Set, Tuple

from .components import build_identity_edges, find_linked_components, make_component_title
from .cytoscape import create_cytoscape_network
from .graph_limits import combined_graph_skip, normalize_combined_graph_limits
from .network_annotations import annotation_node_metadata
from .residue_types import count_polymer_lengths

PROTEIN_NETWORK_TITLE = "Protein_Network"
COMBINED_PROTEIN_NETWORK_TITLE = "Combined_Protein_Network"


def create_protein_network(
    results: List[Dict[str, Any]],
    combined_data: List[Dict[str, Any]],
    run_output_path: str,
    network_config: Dict[str, Any],
    *,
    per_pdb_output_path: Optional[str] = None,
    combined_output_path: Optional[str] = None,
    combined_graph_limits: Mapping[str, Any] | None = None,
) -> List[Dict[str, Any]]:
    """Create protein-level networks (per PDB and/or combined).

    per_pdb_output_path:
        Output folder for Protein_Network_<pdb>.cx2
    combined_output_path:
        Output folder for Combined_Protein_Network.cx2
    """
    skipped_outputs: List[Dict[str, Any]] = []
    graph_limits = normalize_combined_graph_limits(combined_graph_limits)
    make_per_pdb = network_config.get("protein_per_pdb", False)
    make_combined = network_config.get("combined_protein_network", False)
    if not make_per_pdb and not make_combined:
        print("Protein network creation is disabled.")
        return skipped_outputs

    # Default fallbacks (backward compatible)
    per_pdb_out = per_pdb_output_path or run_output_path
    combined_out = combined_output_path or run_output_path

    def count_lengths(chain: Dict[str, Any]) -> Tuple[int, int]:
        if "aa_len" in chain or "nt_len" in chain:
            return int(chain.get("aa_len", 0)), int(chain.get("nt_len", 0))
        return count_polymer_lengths(
            (res.get("residue_name") for res in chain.get("residues") or [])
        )

    # --- Collect chain/protein metadata ---
    chain_to_uniprot: Dict[str, str] = {}
    uniprot_to_pdb_ids: Dict[str, Set[str]] = {}
    uniprot_to_file_labels: Dict[str, Set[str]] = {}
    pdb_to_uniprots: Dict[str, Set[str]] = {}
    chain_info: Dict[str, Dict[str, Any]] = {}
    uniprot_pdb_to_chains: Dict[Tuple[str, str], Set[str]] = {}
    uniprot_pdb_to_uids: Dict[Tuple[str, str], Set[str]] = {}
    chain_to_pdb_id: Dict[str, str] = {}

    for structure in combined_data:
        pdb_id = structure["pdb_id"]
        file_path = str(structure.get("file_path") or "")
        file_label = file_path.rsplit("/", 1)[-1] if file_path else pdb_id
        for chain in structure["atom_data"]:
            chain_id = chain["chain_id"]
            uid = str(chain.get("unique_chain_id") or f"{pdb_id}:{chain_id}")
            chain_to_pdb_id[uid] = pdb_id
            aa_len, nt_len = count_lengths(chain)
            chain_info[uid] = {
                "molecule_type": (chain.get("molecule_type") or "Unknown").strip(),
                "molecule_name": chain.get("molecule_name", "Unknown"),
                "aa_len": aa_len,
                "nt_len": nt_len,
                "pdb_id": pdb_id,
                "chain_id": chain_id,
                "file_label": file_label,
                "embedded_annotations": chain.get("embedded_annotations", {}),
                "embedded_annotation_source": chain.get("embedded_annotation_source"),
            }
            up = chain.get("uniprot_id")
            if up:
                chain_to_uniprot[uid] = up
                uniprot_to_pdb_ids.setdefault(up, set()).add(pdb_id)
                uniprot_to_file_labels.setdefault(up, set()).add(file_label)
                pdb_to_uniprots.setdefault(pdb_id, set()).add(up)
                uniprot_pdb_to_chains.setdefault((up, pdb_id), set()).add(chain_id)
                uniprot_pdb_to_uids.setdefault((up, pdb_id), set()).add(uid)

    def get_color_group_for_combined(up_id: str) -> str:
        pdbs = uniprot_to_pdb_ids.get(up_id, set())
        return "Multi" if len(pdbs) != 1 else next(iter(pdbs))

    def source_chain_uids_for_uniprot(up_id: str, pdb_scope: Optional[str]) -> List[str]:
        """Return internal chain keys without reconstructing them from labels."""
        pdbs = [pdb_scope] if pdb_scope else sorted(uniprot_to_pdb_ids.get(up_id, []))
        return sorted(
            {
                uid
                for pdb_id in pdbs
                for uid in uniprot_pdb_to_uids.get((up_id, pdb_id), set())
            }
        )

    def source_chains_for_uniprot(up_id: str, pdb_scope: Optional[str]) -> List[str]:
        labels = {
            f"{chain_info[uid]['pdb_id']}:{chain_info[uid]['chain_id']}"
            for uid in source_chain_uids_for_uniprot(up_id, pdb_scope)
            if uid in chain_info
        }
        return sorted(labels)

    def protein_name(up_id: str, pdb_scope: Optional[str]) -> str:
        for uid in source_chain_uids_for_uniprot(up_id, pdb_scope):
            name = chain_info[uid].get("molecule_name")
            if name and name != "Unknown":
                return str(name)
        return up_id

    def protein_tooltip(up_id: str, pdb_scope: Optional[str]) -> str:
        lines = [protein_name(up_id, pdb_scope)]
        lines.append("Type: Protein")
        lines.append(f"UniProt: {up_id}")

        pdbs = sorted(uniprot_to_pdb_ids.get(up_id, []))
        if pdb_scope:
            chains = sorted(uniprot_pdb_to_chains.get((up_id, pdb_scope), []))
            aa_total = sum(
                int(chain_info.get(uid, {}).get("aa_len", 0))
                for uid in source_chain_uids_for_uniprot(up_id, pdb_scope)
            )
            lines.append(f"PDB: {pdb_scope}")
            if chains:
                lines.append(f"Chains: {', '.join(chains)}")
            if aa_total:
                lines.append(f"Length: {aa_total} aa")
        else:
            lines.append(f"PDBs: {len(pdbs)}")
            if pdbs:
                lines.append(f"PDB IDs: {', '.join(pdbs)}")
            source_chains = source_chains_for_uniprot(up_id, pdb_scope=None)
            if source_chains:
                lines.append(f"Source chains: {', '.join(source_chains)}")
            source_files = sorted(uniprot_to_file_labels.get(up_id, []))
            if source_files:
                lines.append(f"Source files: {', '.join(source_files)}")

        annotation_metadata = protein_annotation_metadata(up_id, pdb_scope)
        lines.extend(annotation_metadata.get("tooltip_lines", []))

        return "\n".join(lines)

    def protein_annotation_metadata(up_id: str, pdb_scope: Optional[str]) -> Dict[str, Any]:
        records = [
            chain_info[chain_uid]
            for chain_uid in source_chain_uids_for_uniprot(up_id, pdb_scope=pdb_scope)
            if chain_uid in chain_info
        ]
        return annotation_node_metadata(records, uniprot_accession=up_id)

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
        lines.extend(annotation_node_metadata([info]).get("tooltip_lines", []))
        return "\n".join(lines)

    # --- Aggregate edges ---
    pp_edges_by_pdb: Dict[str, Dict[Tuple[str, str], int]] = {}

    na_types = {"DNA", "RNA", "DNA/RNA", "Nucleic Acid"}
    na_edges_by_pdb: Dict[str, Dict[Tuple[str, str], int]] = {}
    na_nodes_by_pdb: Dict[str, Dict[str, Dict[str, Any]]] = {}

    for entry in results:
        cnt = entry.get("all_atoms_count", 0)
        if cnt <= 0:
            continue

        a, b = entry["chain_a"], entry["chain_b"]
        pdb_id = chain_to_pdb_id.get(a) or chain_to_pdb_id.get(b) or str(entry.get("pdb_id") or "")

        ia = chain_info.get(a, {})
        ib = chain_info.get(b, {})
        ta = (ia.get("molecule_type") or "Unknown").strip()
        tb = (ib.get("molecule_type") or "Unknown").strip()

        up_a = chain_to_uniprot.get(a)
        up_b = chain_to_uniprot.get(b)

        if up_a and up_b and up_a != up_b:
            key = tuple(sorted((up_a, up_b)))
            pp_edges_by_pdb.setdefault(pdb_id, {})
            pp_edges_by_pdb[pdb_id][key] = pp_edges_by_pdb[pdb_id].get(key, 0) + cnt

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
                    "name": na_uid,
                    "color_group": na_type,
                    "pdb_id": chain_info.get(na_uid, {}).get("pdb_id", ""),
                    "chain_id": chain_info.get(na_uid, {}).get("chain_id", ""),
                    "source_file": chain_info.get(na_uid, {}).get("file_label", ""),
                    "molecule_type": chain_info.get(na_uid, {}).get("molecule_type", "Unknown"),
                    "molecule_name": chain_info.get(na_uid, {}).get("molecule_name", "Unknown"),
                    "aa_len": chain_info.get(na_uid, {}).get("aa_len", 0),
                    "nt_len": chain_info.get(na_uid, {}).get("nt_len", 0),
                    "uniprot_id": "",
                    "node_kind": "nucleic_acid",
                    "tooltip": na_tooltip(na_uid),
                }
                annotation_metadata = annotation_node_metadata([chain_info.get(na_uid, {})])
                annotation_metadata.pop("tooltip_lines", None)
                for key, value in annotation_metadata.items():
                    if key.startswith("annotation_"):
                        na_nodes_by_pdb[pdb_id][na_uid][key] = value
            na_edges_by_pdb.setdefault(pdb_id, {})
            k = (up, na_uid)
            na_edges_by_pdb[pdb_id][k] = na_edges_by_pdb[pdb_id].get(k, 0) + cnt

    def nodes_from_uniprots(
        uniprot_ids: Iterable[str],
        force_protein_color: bool = False,
        pdb_scope: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        nodes: List[Dict[str, Any]] = []
        for up in sorted(uniprot_ids):
            cg = "Protein" if force_protein_color else get_color_group_for_combined(up)
            pdbs = [pdb_scope] if pdb_scope is not None else sorted(uniprot_to_pdb_ids.get(up, []))
            source_chains = source_chains_for_uniprot(up, pdb_scope=pdb_scope)
            source_files = sorted({
                chain_info[uid]["file_label"]
                for uid in source_chain_uids_for_uniprot(up, pdb_scope)
            })
            display_name = up
            if pdb_scope is None and len(pdbs) > 1:
                display_name = f"{up}\n{len(pdbs)} PDBs"
            nodes.append({
                "id": up,
                "name": display_name,
                "color_group": cg,
                "uniprot_id": up,
                "node_kind": "protein",
                "pdb_count": len(pdbs),
                "pdb_ids": ", ".join(pdbs),
                "source_chains": ", ".join(source_chains),
                "source_files": ", ".join(source_files),
                "molecule_name": protein_name(up, pdb_scope),
                "tooltip": protein_tooltip(up, pdb_scope=pdb_scope),
            })
            annotation_metadata = protein_annotation_metadata(up, pdb_scope)
            annotation_metadata.pop("tooltip_lines", None)
            for key, value in annotation_metadata.items():
                if key.startswith("annotation_"):
                    nodes[-1][key] = value
        return nodes

    # --- Per-PDB networks: proteins in blue + NA neighbors ---
    if make_per_pdb:
        pdbs_with_edges = set(pp_edges_by_pdb.keys()) | set(na_edges_by_pdb.keys())

        for pdb_id in sorted(pdbs_with_edges):
            prot_set = pdb_to_uniprots.get(pdb_id, set())
            if not prot_set:
                continue

            edges_pp = [{"chain_a": up1, "chain_b": up2, "all_atoms_count": w}
                        for (up1, up2), w in pp_edges_by_pdb.get(pdb_id, {}).items()]
            edges_pna = [{"chain_a": up, "chain_b": na_uid, "all_atoms_count": w}
                         for (up, na_uid), w in na_edges_by_pdb.get(pdb_id, {}).items()]
            edges = edges_pp + edges_pna
            if not edges:
                continue

            nodes = nodes_from_uniprots(prot_set, force_protein_color=True, pdb_scope=pdb_id)
            nodes.extend(list(na_nodes_by_pdb.get(pdb_id, {}).values()))
            # NEW: protein per-PDB folder
            create_cytoscape_network(
                edges,
                f"{PROTEIN_NETWORK_TITLE}_{pdb_id}",
                per_pdb_out,
                nodes_data=nodes,
            )

        for pdb_id in sorted(set(pdb_to_uniprots.keys()) - pdbs_with_edges):
            prot_set = pdb_to_uniprots.get(pdb_id, set())
            if not prot_set:
                continue
            nodes = nodes_from_uniprots(prot_set, force_protein_color=True, pdb_scope=pdb_id)
            # NEW: protein per-PDB folder
            create_cytoscape_network(
                [],
                f"{PROTEIN_NETWORK_TITLE}_{pdb_id}",
                per_pdb_out,
                nodes_data=nodes,
            )

    def make_combined_component_title(component_uniprots: Set[str]) -> str:
        return make_component_title(
            COMBINED_PROTEIN_NETWORK_TITLE, component_uniprots
        )

    def find_interlinked_chain_components() -> List[Set[str]]:
        """Return chain components connected by interactions plus cross-PDB UniProt identity."""
        uniprot_to_chain_ids: Dict[str, Set[str]] = {}
        chain_to_pdb: Dict[str, str] = {}
        valid_nodes = set(chain_info.keys())

        for up_id, pdb_ids in uniprot_to_pdb_ids.items():
            for pdb_id in pdb_ids:
                for uid in uniprot_pdb_to_uids.get((up_id, pdb_id), set()):
                    if uid in chain_info:
                        uniprot_to_chain_ids.setdefault(up_id, set()).add(uid)
                        chain_to_pdb[uid] = pdb_id

        identity_edges = build_identity_edges(uniprot_to_chain_ids, chain_to_pdb)
        return find_linked_components(results, identity_edges, valid_nodes=valid_nodes)

    # --- Combined protein networks: same interlinked components as combined chain networks ---
    if make_combined:
        component_node_sets = find_interlinked_chain_components()

        for component_nodes in component_node_sets:
            component_uniprots = {
                chain_to_uniprot[node]
                for node in component_nodes
                if node in chain_to_uniprot
            }
            if not component_uniprots:
                continue

            component_edges_by_pair: Dict[Tuple[str, str], int] = {}
            for entry in results:
                a, b = entry["chain_a"], entry["chain_b"]
                if a not in component_nodes or b not in component_nodes:
                    continue

                cnt = entry.get("all_atoms_count", 0)
                if cnt <= 0:
                    continue

                up_a = chain_to_uniprot.get(a)
                up_b = chain_to_uniprot.get(b)
                if not up_a or not up_b or up_a == up_b:
                    continue

                key = tuple(sorted((up_a, up_b)))
                component_edges_by_pair[key] = component_edges_by_pair.get(key, 0) + cnt

            combined_edges = [
                {"chain_a": up1, "chain_b": up2, "all_atoms_count": weight}
                for (up1, up2), weight in sorted(component_edges_by_pair.items())
            ]
            nodes = nodes_from_uniprots(component_uniprots, force_protein_color=False, pdb_scope=None)
            network_title = make_combined_component_title(component_uniprots)
            skipped = combined_graph_skip(
                network_kind="combined_protein_network",
                name=network_title,
                node_count=len(nodes),
                edge_count=len(combined_edges),
                limits=graph_limits,
            )
            if skipped:
                skipped_outputs.append(skipped)
                continue

            create_cytoscape_network(
                combined_edges,
                network_title,
                combined_out,
                nodes_data=nodes,
            )

    return skipped_outputs
