"""Cytoscape integration helpers for PDB2Net.

This module handles:
- Launching/connecting to Cytoscape (via py4cytoscape).
- Computing deterministic node coordinates for headless exports.
- Building a portable CX2 file without requiring a running Cytoscape UI.
- Creating and styling networks in Cytoscape when a UI session is available.

Notes
-----
- Behavior is preserved; only documentation/comments and message strings are
  standardized to English.
- Headless path writes CX2 directly (no *.cyjs export).
"""

import os
import json
import time
import subprocess
import re  # kept because it may be used by callers and to avoid accidental removal

import pandas as pd
import py4cytoscape as p4c
from matplotlib import cm
from matplotlib.colors import to_hex

from config_loader import config


def ensure_cytoscape_running():
    """Ensure a Cytoscape instance is reachable.

    Tries to ping the Cytoscape REST API. If not available, launches Cytoscape
    using the path configured in `config["cytoscape_path"]`, waits a short
    period, and pings again.

    Exits the program if Cytoscape could not be started.
    """
    try:
        p4c.cytoscape_ping()
        print("Cytoscape is already running.")
    except Exception:
        print("Starting Cytoscape...")
        subprocess.Popen(config["cytoscape_path"])
        time.sleep(40)
        try:
            p4c.cytoscape_ping()
            print("Cytoscape started successfully.")
        except Exception:
            print("Error: Cytoscape could not be started. Check the path in config.json.")
            exit(1)


def compute_preset_positions_spring(nodes_df, edges_df, network_title, scale=1000.0):
    """Compute fast, deterministic force-directed positions (Fruchterman–Reingold style).

    The layout intentionally avoids post-hoc min/max normalization to preserve
    natural spacing. Returns a mapping suitable for CX2 export:

        {node_id: {"x": float, "y": float}}

    Parameters
    ----------
    nodes_df : pandas.DataFrame
        Must contain a column 'id' with node identifiers (strings).
    edges_df : pandas.DataFrame | None
        Optional DataFrame with 'source' and 'target' and optionally
        'all_atoms_count' used as an edge-weight proxy.
    network_title : str
        Used to derive a stable seed so layouts are reproducible across runs.
    scale : float, default 1000.0
        Overall coordinate scale factor.

    Returns
    -------
    dict[str, dict[str, float]]
        Mapping of node id to x/y coordinates.
    """
    import math
    import zlib

    node_ids = list(nodes_df["id"]) if len(nodes_df) else []
    N = len(node_ids)
    if N == 0:
        return {}
    if N == 1:
        return {node_ids[0]: {"x": 0.0, "y": 0.0}}
    if N == 2:
        d = scale * 0.35
        return {node_ids[0]: {"x": -d, "y": 0.0}, node_ids[1]: {"x": d, "y": 0.0}}

    # Prefer NetworkX if available; otherwise fall back to a circle layout.
    try:
        import networkx as nx
    except Exception:
        r = scale * 0.35
        out = {}
        for i, nid in enumerate(node_ids):
            a = 2.0 * math.pi * i / N
            out[nid] = {"x": float(r * math.cos(a)), "y": float(r * math.sin(a))}
        return out

    # Build graph with gently compressed weights to avoid extreme distortions.
    G = nx.Graph()
    G.add_nodes_from(node_ids)

    if edges_df is not None and len(edges_df) > 0:
        for _, e in edges_df.iterrows():
            s, t = e["source"], e["target"]
            if s in G and t in G and s != t:
                w = float(e.get("all_atoms_count", 1.0))
                w = 1.0 + math.log10(max(1.0, w))  # compress large counts
                if G.has_edge(s, t):
                    if w > G[s][t].get("weight", 1.0):
                        G[s][t]["weight"] = w
                else:
                    G.add_edge(s, t, weight=w)

    # Stable seed derived from network title.
    seed = zlib.adler32(str(network_title).encode("utf-8")) & 0xFFFFFFFF

    # Iterations: slightly higher for smaller networks to improve stability.
    iters = max(100, min(250, 8 * N))

    # Let NetworkX scale; center around (0, 0); use half of the provided scale.
    pos = nx.spring_layout(
        G,
        seed=seed,
        weight="weight",
        iterations=iters,
        dim=2,
        center=(0.0, 0.0),
        scale=scale * 0.5,
        k=None,  # let NetworkX choose a good default, works well for N < ~50
    )

    return {n: {"x": float(x), "y": float(y)} for n, (x, y) in pos.items()}


def _export_cx2_headless(network_title, run_output_path, nodes_df, edges_df_for_export, color_map, positions):
    """Write a portable CX2 file without requiring a running Cytoscape session.

    The CX2 will include:
      - Nodes/edges with attributes (name, tooltip, color_group, interaction, all_atoms_count)
      - Node positions (x, y) from `positions`
      - Visual style encoded via visualProperties:

        NODE_LABEL           (Passthrough: name)
        NODE_TOOLTIP         (Passthrough: tooltip)
        NODE_BACKGROUND_COLOR(Discrete: color_group → color)
        EDGE_WIDTH           (Continuous: all_atoms_count → 1..6)

      - Reasonable defaults (node size, edge opacity, curvature, etc.)

    Parameters
    ----------
    network_title : str
        Title to store in networkAttributes.
    run_output_path : str
        Output directory for the CX2 file.
    nodes_df : pandas.DataFrame
        DataFrame with at least 'id', optional 'name', 'tooltip', 'color_group'.
    edges_df_for_export : pandas.DataFrame
        DataFrame with 'source', 'target', 'interaction', 'all_atoms_count'.
    color_map : dict[str, str]
        Mapping from color_group to hex color string.
    positions : dict[str, dict[str, float]]
        Mapping from node id to {'x': float, 'y': float}.
    """
    os.makedirs(run_output_path, exist_ok=True)
    out_path = os.path.join(run_output_path, f"{network_title}.cx2")

    # --- Map string node ids to compact numeric ids required by CX2 ---
    node_ids = list(nodes_df["id"].astype(str))
    nid_map = {nid: i + 1 for i, nid in enumerate(node_ids)}

    # --- Nodes aspect (with x, y and attribute bag 'v') ---
    nodes_aspect = []
    for _, row in nodes_df.iterrows():
        nid = str(row["id"])
        pos = positions.get(nid, {"x": 0.0, "y": 0.0})
        nodes_aspect.append({
            "id": nid_map[nid],
            "x": float(pos["x"]),
            "y": float(pos["y"]),
            "v": {
                "name": str(row.get("name", nid)),
                "tooltip": str(row.get("tooltip", "")),
                "color_group": str(row.get("color_group", "Unknown")),
            },
        })

    # --- Edges aspect (numeric s/t and attribute bag 'v') ---
    edges_aspect = []
    max_w = 1
    for _, r in edges_df_for_export.iterrows():
        try:
            max_w = max(max_w, int(r.get("all_atoms_count", 1)))
        except Exception:
            pass

    for i, r in edges_df_for_export.iterrows():
        s = nid_map.get(str(r["source"]))
        t = nid_map.get(str(r["target"]))
        if s is None or t is None or s == t:
            continue
        edges_aspect.append({
            "id": i + 1,
            "s": s,
            "t": t,
            "v": {
                "interaction": str(r.get("interaction", "interacts_with")),
                "all_atoms_count": int(r.get("all_atoms_count", 1)),
            },
        })

    # --- Attribute declarations (types for attributes) ---
    attr_decls = {
        "attributeDeclarations": [{
            "networkAttributes": {"name": {"d": "string"}},
            "nodes": {
                "name": {"d": "string"},
                "tooltip": {"d": "string"},
                "color_group": {"d": "string"},
            },
            "edges": {
                "interaction": {"d": "string"},
                "all_atoms_count": {"d": "integer"},
            },
        }]
    }

    # --- Visual properties: discrete colors for nodes; continuous width for edges ---
    discrete_map = [{"v": k, "vp": v} for k, v in color_map.items()]

    visual_props = {
        "visualProperties": [{
            "default": {
                "network": {"NETWORK_BACKGROUND_COLOR": "#FFFFFF"},
                "node": {
                    "NODE_SHAPE": "ellipse",
                    "NODE_WIDTH": 40.0,
                    "NODE_HEIGHT": 40.0,
                    "NODE_BORDER_COLOR": "#555555",
                },
                "edge": {
                    "EDGE_LINE_COLOR": "#888888",
                    "EDGE_OPACITY": 0.6,
                    "EDGE_CURVED": True,
                },
            },
            "nodeMapping": {
                "NODE_LABEL": {
                    "type": "PASSTHROUGH",
                    "definition": {"attribute": "name", "type": "string"},
                },
                "NODE_TOOLTIP": {
                    "type": "PASSTHROUGH",
                    "definition": {"attribute": "tooltip", "type": "string"},
                },
                "NODE_BACKGROUND_COLOR": {
                    "type": "DISCRETE",
                    "definition": {
                        "attribute": "color_group",
                        "type": "string",
                        "map": discrete_map,
                    },
                },
            },
            "edgeMapping": {
                "EDGE_WIDTH": {
                    "type": "CONTINUOUS",
                    "definition": {
                        "attribute": "all_atoms_count",
                        "type": "integer",
                        "map": [{
                            "min": 1.0, "includeMin": True,
                            "max": float(max_w), "includeMax": True,
                            "minVPValue": 1.0, "maxVPValue": 6.0,
                        }],
                    },
                },
            },
        }]
    }

    # --- Network attributes (store the network title) ---
    network_attrs = {"networkAttributes": [{"name": network_title}]}

    # --- CX2 container ---
    cx = [
        {"CXVersion": "2.0", "hasFragments": False},
        attr_decls,
        network_attrs,
        {"nodes": nodes_aspect},
        {"edges": edges_aspect},
        visual_props,
        {"status": [{"error": "", "success": True}]},
    ]

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(cx, f, ensure_ascii=False, indent=2)


def create_cytoscape_network(results, network_title="Protein_Interaction_Network", run_output_path=".", nodes_data=None):
    """Create a network either headlessly (CX2 only) or inside Cytoscape, then export CX2.

    Parameters
    ----------
    results : list[dict]
        Edge-like dictionaries including 'chain_a', 'chain_b', and 'all_atoms_count'.
    network_title : str, default "Protein_Interaction_Network"
        Title to assign to the created network / exported CX2 file.
    run_output_path : str, default "."
        Output directory used for CX2 export.
    nodes_data : list[dict] | None
        Optional precomputed nodes (each dict should include 'id', optionally
        'name', 'tooltip', 'color_group', 'molecule_name').

    Notes
    -----
    - When `config["open_in_cytoscape"]` is falsy, no Cytoscape connection is
      required; a CX2 file is written directly.
    - When Cytoscape is open, a visual style is created/applied and the network
      is exported as CX2 via the REST API.
    """
    # --- Collect edges and unique node ids from `results` ---
    unique_nodes = set()
    edges = []
    for entry in results:
        if entry.get("all_atoms_count", 0) > 0:
            a = entry["chain_a"]
            b = entry["chain_b"]
            unique_nodes.add(a)
            unique_nodes.add(b)
            edges.append({"chain_a": a, "chain_b": b, "all_atoms_count": entry["all_atoms_count"]})

    # --- Build node DataFrame (prefer provided nodes_data) ---
    if nodes_data:
        nodes_df = pd.DataFrame(nodes_data).copy()
        if "name" not in nodes_df.columns:
            nodes_df["name"] = nodes_df["id"]
        if "tooltip" not in nodes_df.columns:
            # fallback tooltip to molecule_name or name if not provided
            nodes_df["tooltip"] = nodes_df.get("molecule_name", nodes_df["name"])
    else:
        nodes_df = pd.DataFrame({"id": list(unique_nodes)})
        nodes_df["name"] = nodes_df["id"]
        nodes_df["tooltip"] = nodes_df["id"]

    # --- Build edge DataFrames (one for export, one for Cytoscape UI) ---
    if len(edges) > 0:
        edges_df_for_export = pd.DataFrame(edges).rename(columns={"chain_a": "source", "chain_b": "target"})
        edges_df_for_export["interaction"] = "interacts_with"
        edges_df_for_cyto = edges_df_for_export
    else:
        edges_df_for_export = pd.DataFrame(columns=["source", "target", "interaction", "all_atoms_count"])
        edges_df_for_cyto = None

    # === PATH 1: Headless → write CX2 only (no web.cyjs) ===
    if not config.get("open_in_cytoscape", True):
        os.makedirs(run_output_path, exist_ok=True)

        # Create color mapping: fixed palette + automatic colors for any extras.
        fixed_colors = {
            "Protein": "#1f77b4",
            "DNA": "#ff7f0e",
            "RNA": "#2ca02c",
            "DNA/RNA": "#a2a200",
            "Nucleic Acid": "#9467bd",
            "Unknown": "#7f7f7f",
        }
        color_groups = sorted(nodes_df["color_group"].dropna().unique()) if "color_group" in nodes_df.columns else []
        base_color_groups = [g for g in color_groups if g not in fixed_colors and g != "Multi"]
        cmap = cm.get_cmap("tab20", max(1, len(base_color_groups)))
        auto_colors = {group: to_hex(cmap(i)) for i, group in enumerate(base_color_groups)}
        if "Multi" in color_groups and "Multi" not in fixed_colors:
            auto_colors["Multi"] = "#FF0000"
        color_map = {**fixed_colors, **auto_colors}

        # Deterministic node positions for portable CX2 exports.
        positions = compute_preset_positions_spring(
            nodes_df, edges_df_for_export, network_title, scale=1000.0
        )

        _export_cx2_headless(
            network_title=network_title,
            run_output_path=run_output_path,
            nodes_df=nodes_df,
            edges_df_for_export=edges_df_for_export,
            color_map=color_map,
            positions=positions,
        )
        return

    # === PATH 2: Cytoscape UI available → create network and export CX2 ===
    existing_networks = p4c.get_network_list()
    while len(existing_networks) > config["keep_last_n_networks"]:
        oldest_network = existing_networks.pop(0)
        p4c.delete_network(oldest_network)

    # Choose a collection name by network flavor (Chain / Protein / Combined).
    title_lower = str(network_title).lower()
    if "combined" in title_lower:
        collection_name = "PDB2Net — Combined"
    elif "protein" in title_lower:
        collection_name = "PDB2Net — Protein"
    else:
        collection_name = "PDB2Net — Chain"

    # Create the network from DataFrames.
    try:
        if edges_df_for_cyto is not None and len(edges_df_for_cyto) > 0:
            p4c.create_network_from_data_frames(
                nodes=nodes_df,
                edges=edges_df_for_cyto,
                title=network_title,
                collection=collection_name,
            )
        else:
            p4c.create_network_from_data_frames(
                nodes=nodes_df,
                edges=None,
                title=network_title,
                collection=collection_name,
            )
    except Exception as e:
        print(f"Error while creating network: {e}")

    # Push node attributes into Cytoscape tables.
    try:
        if "name" in nodes_df.columns:
            p4c.load_table_data(
                data=nodes_df[["id", "name"]],
                data_key_column="id",
                table="node",
                table_key_column="id",
            )
        if "color_group" in nodes_df.columns:
            p4c.load_table_data(
                data=nodes_df[["id", "color_group"]],
                data_key_column="id",
                table="node",
                table_key_column="id",
            )
        if "tooltip" in nodes_df.columns:
            p4c.load_table_data(
                data=nodes_df[["id", "tooltip"]],
                data_key_column="id",
                table="node",
                table_key_column="id",
            )
    except Exception as e:
        print(f"Error while loading node data: {e}")

    # Create/apply visual style and run a quick layout.
    try:
        color_groups = sorted(nodes_df["color_group"].dropna().unique()) if "color_group" in nodes_df.columns else []
        is_chain_network = "Chain" in network_title or network_title == "Combined_Network"
        is_protein_network = "Protein" in network_title
        is_combined_protein = is_protein_network and "combined" in network_title.lower()

        style_name = (
            "PDB2Net_Chain_Style" if is_chain_network
            else "PDB2Net_Protein_Combined_Style" if is_combined_protein
            else "PDB2Net_Protein_Style"
        )

        fixed_colors = {
            "Protein": "#1f77b4",
            "DNA": "#ff7f0e",
            "RNA": "#2ca02c",
            "DNA/RNA": "#a2a200",
            "Nucleic Acid": "#9467bd",
        }
        base_color_groups = [g for g in color_groups if g not in fixed_colors and g != "Multi"]
        cmap = cm.get_cmap("tab20", max(1, len(base_color_groups)))
        auto_colors = {group: to_hex(cmap(i)) for i, group in enumerate(base_color_groups)}
        if "Multi" in color_groups and "Multi" not in fixed_colors:
            auto_colors["Multi"] = "#FF0000"
        color_map = {**fixed_colors, **auto_colors}

        if style_name not in p4c.get_visual_style_names():
            defaults = {
                "NODE_SHAPE": "ELLIPSE",
                "NODE_SIZE": 40,
                "NODE_LABEL_POSITION": "C,C,c,0.00,0.00",
                "EDGE_TRANSPARENCY": 120,
            }
            mappings = [
                {
                    "mappingType": "passthrough",
                    "mappingColumn": "name",
                    "mappingColumnType": "String",
                    "visualProperty": "NODE_LABEL",
                },
                {
                    "mappingType": "passthrough",
                    "mappingColumn": "tooltip",
                    "mappingColumnType": "String",
                    "visualProperty": "NODE_TOOLTIP",
                },
                {
                    "mappingType": "discrete",
                    "mappingColumn": "color_group",
                    "mappingColumnType": "String",
                    "visualProperty": "NODE_FILL_COLOR",
                    "map": [{"key": k, "value": v} for k, v in color_map.items()],
                },
            ]
            p4c.create_visual_style(style_name, mappings=mappings, defaults=defaults)

        p4c.set_current_network(network_title)
        p4c.set_visual_style(style_name)
        p4c.layout_network(layout_name="force-directed")
    except Exception as e:
        print(f"Error while applying style: {e}")

    # Export CX2 (fallback to CX if the server lacks CX2 support).
    try:
        os.makedirs(run_output_path, exist_ok=True)
        cx2_file = os.path.join(run_output_path, f"{network_title}.cx2")
        try:
            p4c.export_network(cx2_file, type="cx2")
        except Exception:
            try:
                p4c.export_network(cx2_file, type="cx")
            except Exception:
                p4c.export_network(cx2_file, type="CX")
    except Exception as e:
        print(f"Error while exporting files: {e}")


def generate_nodes_from_atom_data(atom_data, pdb_id=None):
    """Create Cytoscape chain nodes from parsed atom/chain data.

    Produces:
      - name: classic label equal to 'id' ('PDBID:Chain')
      - tooltip: multi-line info (full name, type, length, PDB, UniProt)
      - color_group: molecule type (Protein, DNA, RNA, DNA/RNA, Nucleic Acid, Unknown)

    Parameters
    ----------
    atom_data : list[dict]
        Each chain dict should include keys like:
        'unique_chain_id', 'id', 'molecule_type', 'molecule_name',
        'uniprot_id', 'residues' (list with 'residue_name').
    pdb_id : str | None
        Unused here, kept for API compatibility with callers.

    Returns
    -------
    list[dict]
        Node attribute dicts compatible with create_cytoscape_network().
    """
    # Residue-class sets used for quick length classification.
    dna_set = {"DA", "DT", "DG", "DC", "DI"}
    rna_set = {"A", "U", "G", "C", "I"}
    protein_set = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
        "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
        "TYR", "VAL", "SEC", "PYL",
    }

    def count_lengths(res_list):
        """Count protein (aa) and nucleic-acid (nt) residues in a residue list."""
        aa = nt = 0
        for res in (res_list or []):
            rn = (res.get("residue_name") or "").upper()
            if rn in protein_set:
                aa += 1
            elif rn in dna_set or rn in rna_set:
                nt += 1
        return aa, nt

    nodes = []
    for chain in atom_data:
        # Unique ID used across the pipeline (usually 'PDBID:Chain').
        uid = chain.get("unique_chain_id") or chain.get("id")

        mol_type = (chain.get("molecule_type") or "Unknown").strip()
        mol_name_full = chain.get("molecule_name") or "Unknown"
        up_id = chain.get("uniprot_id")  # optional UniProt accession
        aa_len, nt_len = count_lengths(chain.get("residues"))

        # Build a rich multi-line tooltip.
        details = [str(mol_name_full)]
        details.append(f"Type: {mol_type}")
        if aa_len:
            details.append(f"Length: {aa_len} aa")
        if nt_len:
            details.append(f"Length: {nt_len} nt")
        if uid:
            details.append(f"PDB: {uid}")
        if up_id:
            details.append(f"UniProt: {up_id}")
        tooltip = "\n".join(details)

        nodes.append({
            "id": uid,
            "name": uid,            # keep legacy label equal to the id
            "tooltip": tooltip,     # rich, human-readable tooltip
            "color_group": mol_type or "Unknown",
            "molecule_name": mol_name_full,
        })

    return nodes
