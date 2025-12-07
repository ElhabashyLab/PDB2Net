"""Cytoscape integration helpers for PDB2Net.

This module handles:
- Launching/connecting to Cytoscape (via py4cytoscape).
- Computing deterministic node coordinates for headless exports.
- Building a portable CX2 file without requiring a running Cytoscape UI.
- Creating and styling networks in Cytoscape when a UI session is available.
- Generating node data dictionaries from parsed atom data.

UPDATES:
- Restores ORIGINAL style names and defaults for Protein/Chain networks.
- Applies new 'Linked Identity' style ONLY to the 'Combined_Network'.
- Fixes Collection grouping.
"""

import os
import json
import time
import subprocess
import re

import pandas as pd
import py4cytoscape as p4c
from matplotlib import cm
from matplotlib.colors import to_hex

from config_loader import config


def ensure_cytoscape_running():
    """Ensure a Cytoscape instance is reachable via CyREST."""
    try:
        p4c.cytoscape_ping()
        print("Cytoscape is already running.")
        return
    except Exception:
        pass

    cyto_path = config.get("cytoscape_path")
    if not cyto_path:
        print("Error: Cytoscape is not running and 'cytoscape_path' is not configured.")
        raise SystemExit(1)

    print(f"Starting Cytoscape using: {cyto_path!r}")
    try:
        subprocess.Popen([cyto_path])
    except Exception as e:
        print(f"Error: Failed to start Cytoscape: {e}")
        raise SystemExit(1)

    wait_total = float(config.get("cytoscape_wait_seconds", 60))
    deadline = time.time() + wait_total
    while time.time() < deadline:
        try:
            p4c.cytoscape_ping()
            print("Cytoscape started successfully.")
            return
        except Exception:
            time.sleep(5.0)

    print("Error: Cytoscape did not respond within timeout.")
    raise SystemExit(1)


def compute_preset_positions_spring(nodes_df, edges_df, network_title, scale=1000.0):
    """Compute fast, deterministic force-directed positions."""
    import math
    import zlib
    
    node_ids = list(nodes_df["id"]) if len(nodes_df) else []
    N = len(node_ids)
    if N == 0: return {}
    if N == 1: return {node_ids[0]: {"x": 0.0, "y": 0.0}}

    try:
        import networkx as nx
        G = nx.Graph()
        G.add_nodes_from(node_ids)
        if edges_df is not None and len(edges_df) > 0:
            for _, e in edges_df.iterrows():
                s, t = e["source"], e["target"]
                if s != t:
                    w = 1.0 + math.log10(max(1.0, float(e.get("all_atoms_count", 1.0))))
                    G.add_edge(s, t, weight=w)
        
        seed = zlib.adler32(str(network_title).encode("utf-8")) & 0xFFFFFFFF
        pos = nx.spring_layout(G, seed=seed, weight="weight", scale=scale * 0.8, k=1.0/math.sqrt(max(N, 1)))
        return {n: {"x": float(x), "y": float(y)} for n, (x, y) in pos.items()}
    except Exception:
        return {nid: {"x": math.cos(2*math.pi*i/N)*scale, "y": math.sin(2*math.pi*i/N)*scale} for i, nid in enumerate(node_ids)}


def _export_cx2_headless(network_title, run_output_path, nodes_df, edges_df_for_export, color_map, positions):
    """Write a portable CX2 file without requiring a running Cytoscape session."""
    os.makedirs(run_output_path, exist_ok=True)
    out_path = os.path.join(run_output_path, f"{network_title}.cx2")

    # Only "Combined_Network" gets the special Linked Identity style in Headless too
    is_linked_identity_network = (network_title == "Combined_Network")

    # --- Node Mapping ---
    node_ids = list(nodes_df["id"].astype(str))
    nid_map = {nid: i + 1 for i, nid in enumerate(node_ids)}

    nodes_aspect = []
    for _, row in nodes_df.iterrows():
        nid = str(row["id"])
        pos = positions.get(nid, {"x": 0.0, "y": 0.0})
        
        # Border Color Logic
        if is_linked_identity_network:
            border_color = str(row.get("uniprot_border_color", "#555555"))
        else:
            border_color = "#000000"

        nodes_aspect.append({
            "id": nid_map[nid],
            "x": float(pos["x"]),
            "y": float(pos["y"]),
            "v": {
                "name": str(row.get("name", nid)),
                "tooltip": str(row.get("tooltip", "")),
                "color_group": str(row.get("color_group", "Unknown")),
                "uniprot_border_color": border_color
            },
        })

    # --- Edge Mapping ---
    edges_aspect = []
    max_w = 1
    for _, r in edges_df_for_export.iterrows():
        try: max_w = max(max_w, int(r.get("all_atoms_count", 1)))
        except: pass

    for i, r in edges_df_for_export.iterrows():
        s = nid_map.get(str(r["source"]))
        t = nid_map.get(str(r["target"]))
        if s is None or t is None or s == t: continue
        edges_aspect.append({
            "id": i + 1, "s": s, "t": t,
            "v": {
                "interaction": str(r.get("interaction", "interacts_with")),
                "all_atoms_count": int(r.get("all_atoms_count", 1)),
            },
        })

    # --- Visual Properties ---
    discrete_map = [{"v": k, "vp": v} for k, v in color_map.items()]

    if is_linked_identity_network:
        # Linked Identity Style
        border_width = 5.0
        border_paint_mapping = {"type": "PASSTHROUGH", "definition": {"attribute": "uniprot_border_color", "type": "string"}}
        edge_style_mapping = {"type": "DISCRETE", "definition": {"attribute": "interaction", "type": "string", "map": [{"v": "identity", "vp": "DOT"}, {"v": "interacts_with", "vp": "SOLID"}]}}
    else:
        # Original Headless Defaults (black borders 2.0 were standard here)
        border_width = 2.0
        border_paint_mapping = None
        edge_style_mapping = None 

    visual_props = {
        "visualProperties": [{
            "default": {
                "network": {"NETWORK_BACKGROUND_COLOR": "#FFFFFF"},
                "node": {
                    "NODE_SHAPE": "ellipse",
                    "NODE_WIDTH": 40.0, "NODE_HEIGHT": 40.0,
                    "NODE_BORDER_WIDTH": border_width,
                    "NODE_BORDER_OPACITY": 1.0,
                    "NODE_BORDER_COLOR": "#000000" if not is_linked_identity_network else "#555555"
                },
                "edge": {
                    "EDGE_LINE_COLOR": "#888888", "EDGE_OPACITY": 0.6, "EDGE_CURVED": True, "EDGE_LINE_STYLE": "SOLID"
                },
            },
            "nodeMapping": {
                "NODE_LABEL": { "type": "PASSTHROUGH", "definition": {"attribute": "name", "type": "string"} },
                "NODE_TOOLTIP": { "type": "PASSTHROUGH", "definition": {"attribute": "tooltip", "type": "string"} },
                "NODE_BACKGROUND_COLOR": { "type": "DISCRETE", "definition": { "attribute": "color_group", "type": "string", "map": discrete_map } }
            },
            "edgeMapping": {
                "EDGE_WIDTH": {
                    "type": "CONTINUOUS",
                    "definition": {
                        "attribute": "all_atoms_count", "type": "integer",
                        "map": [{"min": 1.0, "includeMin": True, "max": float(max_w), "includeMax": True, "minVPValue": 1.0, "maxVPValue": 6.0}],
                    },
                },
            },
        }]
    }
    
    if border_paint_mapping:
        visual_props["visualProperties"][0]["nodeMapping"]["NODE_BORDER_COLOR"] = border_paint_mapping
    if edge_style_mapping:
        visual_props["visualProperties"][0]["edgeMapping"]["EDGE_LINE_STYLE"] = edge_style_mapping

    attr_decls = {
        "attributeDeclarations": [{
            "networkAttributes": {"name": {"d": "string"}},
            "nodes": {"name": {"d": "string"}, "tooltip": {"d": "string"}, "color_group": {"d": "string"}, "uniprot_border_color": {"d": "string"}},
            "edges": {"interaction": {"d": "string"}, "all_atoms_count": {"d": "integer"}},
        }]
    }

    cx = [{"CXVersion": "2.0", "hasFragments": False}, attr_decls, {"networkAttributes": [{"name": network_title}]}, {"nodes": nodes_aspect}, {"edges": edges_aspect}, visual_props, {"status": [{"error": "", "success": True}]}]
    
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(cx, f, ensure_ascii=False, indent=2)


def create_cytoscape_network(results, network_title="Protein_Interaction_Network", run_output_path=".", nodes_data=None):
    """Create a network either headlessly (CX2 only) or inside Cytoscape, then export CX2."""
    
    # --- Prepare Data ---
    unique_nodes = set()
    edges = []
    for entry in results:
        if entry.get("all_atoms_count", 0) > 0:
            a, b = entry["chain_a"], entry["chain_b"]
            unique_nodes.add(a); unique_nodes.add(b)
            edges.append({
                "chain_a": a, "chain_b": b, 
                "all_atoms_count": entry["all_atoms_count"],
                "interaction": entry.get("interaction_type", "interacts_with")
            })

    if nodes_data:
        nodes_df = pd.DataFrame(nodes_data).copy()
        if "name" not in nodes_df.columns: nodes_df["name"] = nodes_df["id"]
        if "tooltip" not in nodes_df.columns: nodes_df["tooltip"] = nodes_df.get("molecule_name", nodes_df["name"])
    else:
        nodes_df = pd.DataFrame({"id": list(unique_nodes)})
        nodes_df["name"] = nodes_df["id"]
        nodes_df["tooltip"] = nodes_df["id"]

    if "uniprot_border_color" not in nodes_df.columns:
        nodes_df["uniprot_border_color"] = "#000000"

    edges_df = pd.DataFrame(edges).rename(columns={"chain_a": "source", "chain_b": "target"}) if len(edges) > 0 else None
    if edges_df is not None and "interaction" not in edges_df.columns: 
        edges_df["interaction"] = "interacts_with"

    # --- Setup Colors ---
    fixed_colors = {"Protein": "#1f77b4", "DNA": "#ff7f0e", "RNA": "#2ca02c", "DNA/RNA": "#a2a200", "Nucleic Acid": "#9467bd", "Unknown": "#7f7f7f"}
    color_groups = sorted(nodes_df["color_group"].dropna().unique()) if "color_group" in nodes_df.columns else []
    base_color_groups = [g for g in color_groups if g not in fixed_colors and g != "Multi"]
    cmap = cm.get_cmap("tab20", max(1, len(base_color_groups)))
    auto_colors = {group: to_hex(cmap(i)) for i, group in enumerate(base_color_groups)}
    if "Multi" in color_groups and "Multi" not in fixed_colors: auto_colors["Multi"] = "#FF0000"
    color_map = {**fixed_colors, **auto_colors}

    # === PATH 1: Headless ===
    if not config.get("open_in_cytoscape", True):
        positions = compute_preset_positions_spring(nodes_df, edges_df, network_title, scale=1000.0)
        _export_cx2_headless(network_title, run_output_path, nodes_df, edges_df, color_map, positions)
        return

    # === PATH 2: Live Cytoscape ===
    existing_networks = p4c.get_network_list()
    while len(existing_networks) > config["keep_last_n_networks"]:
        p4c.delete_network(existing_networks.pop(0))

    # --- Collection Naming Logic ---
    title_lower = str(network_title).lower()
    if "combined" in title_lower:
        collection_name = "PDB2Net — Combined"
    elif "protein" in title_lower:
        collection_name = "PDB2Net — Protein"
    else:
        collection_name = "PDB2Net — Chain"

    try:
        p4c.create_network_from_data_frames(nodes=nodes_df, edges=edges_df, title=network_title, collection=collection_name)
    except Exception as e: print(f"Error creating network: {e}")

    try:
        cols = ["id", "name", "color_group", "tooltip", "uniprot_border_color", "uniprot_id"]
        load_cols = [c for c in cols if c in nodes_df.columns]
        p4c.load_table_data(nodes_df[load_cols], data_key_column="id", table="node", table_key_column="id")
    except Exception as e: print(f"Error loading data: {e}")

    # --- Style Naming Logic (RESTORED ORIGINAL LOGIC) ---
    is_combined_chain_network = (network_title == "Combined_Network")
    is_protein_network = "Protein" in network_title
    is_combined_protein = is_protein_network and "combined" in title_lower
    is_chain_network = "Chain" in network_title

    if is_combined_chain_network:
        style_name = "PDB2Net_Linked_Identity_Style"
    elif is_combined_protein:
        style_name = "PDB2Net_Protein_Combined_Style"
    elif is_protein_network:
        style_name = "PDB2Net_Protein_Style"
    else:
        style_name = "PDB2Net_Chain_Style"

    # Define defaults and mappings
    if style_name not in p4c.get_visual_style_names():
        
        # Base Defaults (Original)
        defaults = {
            "NODE_SHAPE": "ELLIPSE",
            "NODE_SIZE": 40 if not is_combined_chain_network else 45,
            "NODE_LABEL_POSITION": "C,C,c,0.00,0.00",
            "EDGE_TRANSPARENCY": 120,
        }

        # Base Mappings (Original)
        mappings = [
            {"mappingType": "passthrough", "mappingColumn": "name", "mappingColumnType": "String", "visualProperty": "NODE_LABEL"},
            {"mappingType": "passthrough", "mappingColumn": "tooltip", "mappingColumnType": "String", "visualProperty": "NODE_TOOLTIP"},
            {"mappingType": "discrete", "mappingColumn": "color_group", "mappingColumnType": "String", "visualProperty": "NODE_FILL_COLOR", "map": [{"key": k, "value": v} for k, v in color_map.items()]},
        ]

        # ONLY for the new Linked Identity Network add the border/dash customization
        if is_combined_chain_network:
            defaults["NODE_BORDER_WIDTH"] = 5.0
            defaults["NODE_BORDER_PAINT"] = "#555555"
            mappings.append({"mappingType": "passthrough", "mappingColumn": "uniprot_border_color", "mappingColumnType": "String", "visualProperty": "NODE_BORDER_PAINT"})
            mappings.append({"mappingType": "discrete", "mappingColumn": "interaction", "mappingColumnType": "String", "visualProperty": "EDGE_LINE_TYPE", "map": [{"key": "identity", "value": "LONG_DASH"}, {"key": "interacts_with", "value": "SOLID"}]})
        
        p4c.create_visual_style(style_name, mappings=mappings, defaults=defaults)

    p4c.set_visual_style(style_name)
    p4c.layout_network(layout_name="force-directed")

    try:
        os.makedirs(run_output_path, exist_ok=True)
        p4c.export_network(os.path.join(run_output_path, f"{network_title}.cx2"), type="cx2")
    except Exception as e: print(f"Error exporting: {e}")


def generate_nodes_from_atom_data(atom_data, pdb_id=None):
    """Create Cytoscape chain nodes from parsed atom/chain data."""
    dna_set = {"DA", "DT", "DG", "DC", "DI"}
    rna_set = {"A", "U", "G", "C", "I"}
    protein_set = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
        "TYR", "VAL", "SEC", "PYL",
    }

    def count_lengths(res_list):
        aa = nt = 0
        for res in (res_list or []):
            rn = (res.get("residue_name") or "").upper()
            if rn in protein_set: aa += 1
            elif rn in dna_set or rn in rna_set: nt += 1
        return aa, nt

    nodes = []
    for chain in atom_data:
        uid = chain.get("unique_chain_id") or chain.get("id")
        mol_type = (chain.get("molecule_type") or "Unknown").strip()
        mol_name_full = chain.get("molecule_name") or "Unknown"
        up_id = chain.get("uniprot_id")
        aa_len, nt_len = count_lengths(chain.get("residues"))

        details = [str(mol_name_full)]
        details.append(f"Type: {mol_type}")
        if aa_len: details.append(f"Length: {aa_len} aa")
        if nt_len: details.append(f"Length: {nt_len} nt")
        if uid: details.append(f"PDB: {uid}")
        if up_id: details.append(f"UniProt: {up_id}")
        tooltip = "\n".join(details)

        nodes.append({
            "id": uid,
            "name": uid,
            "tooltip": tooltip,
            "color_group": mol_type or "Unknown",
            "molecule_name": mol_name_full,
        })

    return nodes