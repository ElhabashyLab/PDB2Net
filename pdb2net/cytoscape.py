"""Cytoscape integration helpers for PDB2Net.

This module handles:
- Launching/connecting to Cytoscape (via py4cytoscape).
- Computing deterministic node coordinates for headless exports.
- Building a portable CX2 file without requiring a running Cytoscape UI.
- Creating and styling networks in Cytoscape when a UI session is available.
- Generating node data dictionaries from parsed atom data.

UPDATES:
- "Grid Packing" algorithm implemented for Combined Networks:
  Separates disconnected islands, layouts them individually, and packs them tightly.
- Standard Networks (Protein/Chain per PDB) use the EXACT original scaling logic.
"""

import os
import json
import time
import subprocess
import re
import math
import zlib

import pandas as pd
import py4cytoscape as p4c
from matplotlib import cm
from matplotlib.colors import to_hex

from config_loader import config


# Shared style constants for the Linked Identity Combined Network
LINKED_IDENTITY_BORDER_WIDTH = 12.0

# Unique tag per Python run so styles are always fresh across runs
STYLE_RUN_TAG = f"{int(time.time() * 1000)}_{os.getpid()}"


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
    """
    Compute positions using specific strategies for different network types.
    """
    import networkx as nx

    node_ids = list(nodes_df["id"]) if len(nodes_df) else []
    N = len(node_ids)

    if N == 0:
        return {}
    if N == 1:
        return {node_ids[0]: {"x": 0.0, "y": 0.0}}

    # 1. Build Graph
    G = nx.Graph()
    G.add_nodes_from(node_ids)
    if edges_df is not None and len(edges_df) > 0:
        for _, e in edges_df.iterrows():
            s, t = e["source"], e["target"]
            if s != t:
                # Logarithmic weight
                w = float(e.get("all_atoms_count", 1.0))
                w = 1.0 + math.log10(max(1.0, w))
                G.add_edge(s, t, weight=w)

    # 2. Determine Strategy
    is_combined = "combined" in str(network_title).lower()

    # Deterministic Seed
    seed = zlib.adler32(str(network_title).encode("utf-8")) & 0xFFFFFFFF

    # === STRATEGY A: GRID PACKING (For Combined Networks) ===
    # Solves the "Exploded Graph" problem by packing components tightly.
    if is_combined:
        components = list(nx.connected_components(G))
        # Sort components: largest first to anchor the grid
        components.sort(key=len, reverse=True)

        # Grid Setup
        num_comps = len(components)
        grid_cols = math.ceil(math.sqrt(num_comps))

        final_pos = {}

        # Spacing parameters
        # We assume a base size per node to calculate cell offsets
        base_padding = 150.0

        for i, comp in enumerate(components):
            sub_G = G.subgraph(comp)
            sub_N = len(comp)

            # Local layout for this island
            # Scale depends on island size to keep edge lengths consistent
            # Small islands get small space, big ones get more
            sub_scale = 100.0 + (math.sqrt(sub_N) * 60.0)

            sub_seed = (seed + i) & 0xFFFFFFFF
            sub_pos = nx.spring_layout(
                sub_G,
                seed=sub_seed,
                weight="weight",
                scale=sub_scale,
                center=(0, 0),  # Center at 0,0 relative to cell
                k=None,  # Let nx decide optimal k locally
            )

            # Calculate Grid Cell Position
            row = i // grid_cols
            col = i % grid_cols

            # We use a fixed stride for simplicity, assuming the largest component defines the stride.
            # (A refined version would do bin-packing, but grid is robust enough)
            # We assume a "max cell size" based on the overall N roughly
            stride = 600.0 + (math.sqrt(N / grid_cols) * 50.0)

            offset_x = col * stride
            offset_y = row * stride

            # Apply offset
            for node, coords in sub_pos.items():
                final_pos[node] = {
                    "x": float(coords[0] + offset_x),
                    "y": float(coords[1] + offset_y),
                }

        return final_pos

    # === STRATEGY B: ORIGINAL LOGIC (For Standard Networks) ===
    # EXACT copy of original logic for backward compatibility
    else:
        if N == 2:
            d = scale * 0.08
            return {
                node_ids[0]: {"x": -d, "y": 0.0},
                node_ids[1]: {"x": d, "y": 0.0},
            }

        iters = max(100, min(250, 8 * N))

        scale_used = scale * 0.8
        if N < 40:
            scale_used *= max(0.12, (N / 40.0))

        k = 1.0 / math.sqrt(max(N, 100))
        if N < 40:
            k *= 0.75

        pos = nx.spring_layout(
            G,
            seed=seed,
            weight="weight",
            iterations=iters,
            dim=2,
            center=(0.0, 0.0),
            scale=scale_used,
            k=k,
        )
        return {n: {"x": float(x), "y": float(y)} for n, (x, y) in pos.items()}


def _export_cx2_headless(
    network_title, run_output_path, nodes_df, edges_df_for_export, color_map, positions
):
    """Write a portable CX2 file without requiring a running Cytoscape session."""
    os.makedirs(run_output_path, exist_ok=True)
    out_path = os.path.join(run_output_path, f"{network_title}.cx2")

    # Only "Combined_Network" gets the special Linked Identity style
    is_linked_identity_network = network_title == "Combined_Network"

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

        nodes_aspect.append(
            {
                "id": nid_map[nid],
                "x": float(pos["x"]),
                "y": float(pos["y"]),
                "v": {
                    "name": str(row.get("name", nid)),
                    "tooltip": str(row.get("tooltip", "")),
                    "color_group": str(row.get("color_group", "Unknown")),
                    "uniprot_border_color": border_color,
                },
            }
        )

    # --- Edge Mapping ---
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
        edges_aspect.append(
            {
                "id": i + 1,
                "s": s,
                "t": t,
                "v": {
                    "interaction": str(r.get("interaction", "interacts_with")),
                    "all_atoms_count": int(r.get("all_atoms_count", 1)),
                },
            }
        )

    # --- Visual Properties ---
    discrete_map = [{"v": k, "vp": v} for k, v in color_map.items()]

    if is_linked_identity_network:
        # Style: Linked Identity
        border_width = LINKED_IDENTITY_BORDER_WIDTH
        border_paint_mapping = {
            "type": "PASSTHROUGH",
            "definition": {"attribute": "uniprot_border_color", "type": "string"},
        }
        edge_style_mapping = {
            "type": "DISCRETE",
            "definition": {
                "attribute": "interaction",
                "type": "string",
                "map": [
                    {"v": "identity", "vp": "DOT"},
                    {"v": "interacts_with", "vp": "SOLID"},
                ],
            },
        }
    else:
        # Style: Standard
        border_width = 2.0
        border_paint_mapping = None
        edge_style_mapping = None

    visual_props = {
        "visualProperties": [
            {
                "default": {
                    "network": {"NETWORK_BACKGROUND_COLOR": "#FFFFFF"},
                    "node": {
                        "NODE_SHAPE": "ellipse",
                        "NODE_WIDTH": 40.0,
                        "NODE_HEIGHT": 40.0,
                        "NODE_BORDER_WIDTH": border_width,
                        "NODE_BORDER_OPACITY": 1.0,
                        "NODE_BORDER_COLOR": "#000000"
                        if not is_linked_identity_network
                        else "#555555",
                    },
                    "edge": {
                        "EDGE_LINE_COLOR": "#888888",
                        "EDGE_OPACITY": 0.6,
                        "EDGE_CURVED": True,
                        "EDGE_LINE_STYLE": "SOLID",
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
                            "map": [
                                {
                                    "min": 1.0,
                                    "includeMin": True,
                                    "max": float(max_w),
                                    "includeMax": True,
                                    "minVPValue": 1.0,
                                    "maxVPValue": 6.0,
                                }
                            ],
                        },
                    },
                },
            }
        ]
    }

    if border_paint_mapping:
        visual_props["visualProperties"][0]["nodeMapping"][
            "NODE_BORDER_COLOR"
        ] = border_paint_mapping
    if edge_style_mapping:
        visual_props["visualProperties"][0]["edgeMapping"][
            "EDGE_LINE_STYLE"
        ] = edge_style_mapping

    attr_decls = {
        "attributeDeclarations": [
            {
                "networkAttributes": {"name": {"d": "string"}},
                "nodes": {
                    "name": {"d": "string"},
                    "tooltip": {"d": "string"},
                    "color_group": {"d": "string"},
                    "uniprot_border_color": {"d": "string"},
                },
                "edges": {
                    "interaction": {"d": "string"},
                    "all_atoms_count": {"d": "integer"},
                },
            }
        ]
    }

    # Some viewers require this even though it is optional in the spec
    meta = {
        "metaData": [
            {"name": "attributeDeclarations", "elementCount": 1},
            {"name": "networkAttributes", "elementCount": 1},
            {"name": "nodes", "elementCount": len(nodes_aspect)},
            {"name": "edges", "elementCount": len(edges_aspect)},
            {
                "name": "visualProperties",
                "elementCount": len(visual_props.get("visualProperties", [])),
            },
        ]
    }

    cx = [
        {"CXVersion": "2.0", "hasFragments": False},
        meta,
        attr_decls,
        {"networkAttributes": [{"name": network_title}]},
        {"nodes": nodes_aspect},
        {"edges": edges_aspect},
        visual_props,
        {"status": [{"error": "", "success": True}]},
    ]

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(cx, f, ensure_ascii=False, indent=2)


def create_cytoscape_network(
    results, network_title="Protein_Interaction_Network", run_output_path=".", nodes_data=None
):
    """Create a network either headlessly (CX2 only) or inside Cytoscape, then export CX2."""

    def _verbose_enabled() -> bool:
        return os.environ.get("PDB2NET_VERBOSE", "").strip().lower() in {
            "1",
            "true",
            "yes",
            "on",
        }

    def _compute_edge_name(source: str, target: str, interaction: str) -> str:
        return f"{source} ({interaction}) {target}"

    def _style_signature(
        color_map: dict, extra_tag: str = "", border_width: float | None = None
    ) -> str:
        """
        Deterministic short signature for the current visual mapping/config.
        """
        import zlib

        keys = sorted([str(k) for k in color_map.keys()])
        bw = "" if border_width is None else f"|bw={border_width}"
        payload = (extra_tag + bw + "|" + "|".join(keys)).encode("utf-8")
        return f"{zlib.adler32(payload) & 0xFFFFFFFF:08x}"

    def _ensure_style(style_name: str, mappings: list, defaults: dict) -> None:
        """
        Create the style if it does not already exist in the current Cytoscape session.
        Because STYLE_RUN_TAG is unique per Python run, styles are fresh across runs,
        while still avoiding duplicate-creation errors within the same run.
        """
        try:
            existing = set(p4c.get_visual_style_names())
        except Exception:
            existing = set()

        if style_name not in existing:
            p4c.create_visual_style(style_name, mappings=mappings, defaults=defaults)

    def _load_edge_attributes_by_name(network_suid: int, edges_full_df: pd.DataFrame) -> None:
        """Load edge attributes keyed by edge 'name' (stable), avoiding SUID join issues."""
        if edges_full_df is None or edges_full_df.empty:
            return

        base_cols = {"source", "target", "interaction"}
        extra_cols = [c for c in edges_full_df.columns if c not in base_cols]
        if not extra_cols:
            return

        df = edges_full_df.copy()
        df["name"] = [
            _compute_edge_name(str(s), str(t), str(i))
            for s, t, i in zip(df["source"], df["target"], df["interaction"])
        ]
        load_df = df[["name"] + extra_cols].copy()

        last_exc = None
        for _ in range(8):
            try:
                p4c.load_table_data(
                    load_df,
                    data_key_column="name",
                    table="edge",
                    table_key_column="name",
                    network=network_suid,
                )
                return
            except Exception as e:
                last_exc = e
                time.sleep(0.25)

        if _verbose_enabled() and last_exc is not None:
            print(f"[cytoscape] Edge attribute import failed (non-fatal): {last_exc}")

    # --- Prepare Data ---
    unique_nodes = set()
    edges = []
    for entry in results:
        if entry.get("all_atoms_count", 0) > 0:
            a, b = entry["chain_a"], entry["chain_b"]
            unique_nodes.add(a)
            unique_nodes.add(b)
            edges.append(
                {
                    "chain_a": a,
                    "chain_b": b,
                    "all_atoms_count": entry["all_atoms_count"],
                    "interaction": entry.get("interaction_type", "interacts_with"),
                }
            )

    if nodes_data:
        nodes_df = pd.DataFrame(nodes_data).copy()
        if "name" not in nodes_df.columns:
            nodes_df["name"] = nodes_df["id"]
        if "tooltip" not in nodes_df.columns:
            nodes_df["tooltip"] = nodes_df.get("molecule_name", nodes_df["name"])
    else:
        nodes_df = pd.DataFrame({"id": list(unique_nodes)})
        nodes_df["name"] = nodes_df["id"]
        nodes_df["tooltip"] = nodes_df["id"]

    if "uniprot_border_color" not in nodes_df.columns:
        nodes_df["uniprot_border_color"] = "#000000"

    if len(edges) > 0:
        edges_df = pd.DataFrame(edges).rename(columns={"chain_a": "source", "chain_b": "target"})
        if "interaction" not in edges_df.columns:
            edges_df["interaction"] = "interacts_with"
    else:
        edges_df = pd.DataFrame(
            columns=["source", "target", "interaction", "all_atoms_count"]
        )

    # --- Setup Colors ---
    fixed_colors = {
        "Protein": "#1f77b4",
        "DNA": "#ff7f0e",
        "RNA": "#2ca02c",
        "DNA/RNA": "#a2a200",
        "Nucleic Acid": "#9467bd",
        "Unknown": "#7f7f7f",
    }
    color_groups = (
        sorted(nodes_df["color_group"].dropna().unique())
        if "color_group" in nodes_df.columns
        else []
    )
    base_color_groups = [g for g in color_groups if g not in fixed_colors and g != "Multi"]
    cmap = cm.get_cmap("tab20", max(1, len(base_color_groups)))
    auto_colors = {group: to_hex(cmap(i)) for i, group in enumerate(base_color_groups)}
    if "Multi" in color_groups and "Multi" not in fixed_colors:
        auto_colors["Multi"] = "#FF0000"

    color_map = {**fixed_colors, **auto_colors}

    # === PATH 1: Headless ===
    if not config.get("open_in_cytoscape", True):
        positions = compute_preset_positions_spring(nodes_df, edges_df, network_title)
        _export_cx2_headless(network_title, run_output_path, nodes_df, edges_df, color_map, positions)
        return

    # === PATH 2: Live Cytoscape ===
    existing_networks = p4c.get_network_list()
    while len(existing_networks) > config["keep_last_n_networks"]:
        p4c.delete_network(existing_networks.pop(0))

    # --- Collection Naming ---
    title_lower = str(network_title).lower()
    if "combined" in title_lower:
        collection_name = "PDB2Net — Combined"
    elif "protein" in title_lower:
        collection_name = "PDB2Net — Protein"
    else:
        collection_name = "PDB2Net — Chain"

    # Create network WITHOUT extra edge attribute columns to avoid py4c SUID-join noise
    edges_df_for_create = None
    if edges_df is not None and not edges_df.empty:
        edges_df_for_create = edges_df[["source", "target", "interaction"]].copy()

    try:
        network_suid = p4c.create_network_from_data_frames(
            nodes=nodes_df,
            edges=edges_df_for_create,
            title=network_title,
            collection=collection_name,
        )
    except Exception as e:
        if _verbose_enabled():
            print(f"[cytoscape] Error creating network: {e}")
        return

    # Load edge attributes after creation (stable)
    _load_edge_attributes_by_name(network_suid, edges_df)

    # --- Style Naming Strategy ---
    is_linked_identity_network = network_title == "Combined_Network"
    is_combined_protein_network = network_title == "Combined_Protein_Network"

    if is_linked_identity_network:
        style_name = (
            f"PDB2Net_Linked_Identity_{STYLE_RUN_TAG}_"
            f"{_style_signature(color_map, extra_tag='linked', border_width=LINKED_IDENTITY_BORDER_WIDTH)}"
        )
    elif is_combined_protein_network:
        style_name = (
            f"PDB2Net_Combined_Protein_{STYLE_RUN_TAG}_"
            f"{_style_signature(color_map, extra_tag='combined_protein')}"
        )
    else:
        style_name = (
            f"PDB2Net_Standard_{STYLE_RUN_TAG}_"
            f"{_style_signature(color_map, extra_tag='standard')}"
        )

    # Defaults
    defaults = {
        "NODE_SHAPE": "ELLIPSE",
        "NODE_SIZE": 45,
        "NODE_LABEL_POSITION": "C,C,c,0.00,0.00",
        "EDGE_TRANSPARENCY": 120,
        "NODE_BORDER_WIDTH": LINKED_IDENTITY_BORDER_WIDTH if is_linked_identity_network else 2.0,
        "NODE_BORDER_PAINT": "#555555" if is_linked_identity_network else "#000000",
    }

    if is_linked_identity_network or is_combined_protein_network:
        # Neutral default for Combined so unmapped values do not appear "all red"
        defaults["NODE_FILL_COLOR"] = "#BDBDBD"

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

    if is_linked_identity_network:
        mappings.append(
            {
                "mappingType": "passthrough",
                "mappingColumn": "uniprot_border_color",
                "mappingColumnType": "String",
                "visualProperty": "NODE_BORDER_PAINT",
            }
        )
        mappings.append(
            {
                "mappingType": "discrete",
                "mappingColumn": "interaction",
                "mappingColumnType": "String",
                "visualProperty": "EDGE_LINE_TYPE",
                "map": [
                    {"key": "identity", "value": "LONG_DASH"},
                    {"key": "interacts_with", "value": "SOLID"},
                ],
            }
        )

    _ensure_style(style_name, mappings=mappings, defaults=defaults)
    p4c.set_visual_style(style_name)

    # Layout
    p4c.layout_network(layout_name="force-directed")

    # Export
    try:
        os.makedirs(run_output_path, exist_ok=True)
        p4c.export_network(os.path.join(run_output_path, f"{network_title}.cx2"), type="cx2")
    except Exception as e:
        if _verbose_enabled():
            print(f"[cytoscape] Error exporting: {e}")


def generate_nodes_from_atom_data(atom_data, pdb_id=None):
    """Create Cytoscape chain nodes from parsed atom/chain data."""
    dna_set = {"DA", "DT", "DG", "DC", "DI"}
    rna_set = {"A", "U", "G", "C", "I"}
    protein_set = {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "SEC",
        "PYL",
    }

    def count_lengths(res_list):
        aa = nt = 0
        for res in res_list or []:
            rn = (res.get("residue_name") or "").upper()
            if rn in protein_set:
                aa += 1
            elif rn in dna_set or rn in rna_set:
                nt += 1
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
        if aa_len:
            details.append(f"Length: {aa_len} aa")
        if nt_len:
            details.append(f"Length: {nt_len} nt")
        if uid:
            details.append(f"PDB: {uid}")
        if up_id:
            details.append(f"UniProt: {up_id}")
        tooltip = "\n".join(details)

        nodes.append(
            {
                "id": uid,
                "name": uid,
                "tooltip": tooltip,
                "color_group": mol_type or "Unknown",
                "molecule_name": mol_name_full,
            }
        )

    return nodes