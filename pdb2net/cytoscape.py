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
import math
import zlib
import colorsys

import pandas as pd
import py4cytoscape as p4c
from matplotlib import cm
from matplotlib.colors import to_hex

from config_loader import config


# ---------------------------------------------------------------------------
# Visual tuning block
# Change these values when you want to test different colors / widths / styles.
# The rest of the file reads from these constants so headless + Cytoscape stay
# aligned as closely as possible.
# ---------------------------------------------------------------------------

NODE_COLOR_FIXED = {
    "Protein": "#1f77b4",
    "DNA": "#ff7f0e",
    "RNA": "#2ca02c",
    "DNA/RNA": "#a2a200",
    "Nucleic Acid": "#9467bd",
    "Unknown": "#7f7f7f",
}
MULTI_NODE_COLOR = "#FF0000"

COMBINED_FILE_PALETTE = [
    "#FF8DA1",
    "#a8e953",
    "#50caec",
    "#bec0c0",
    "#EE4E6A",  
    "#6ae771",  
    "#6593ec",
    "#a4b1a0",
    "#E359c2",
    "#cded64",
]

VISUAL_TUNING = {
    "network_background_color": "#FFFFFF",
    "combined_default_node_fill": "#BDBDBD",
    "node_shape": "ellipse",
    "node_width": 40.0,
    "node_height": 40.0,
    "node_size_live": 45,
    "node_border_opacity": 1.0,
    "node_label_position_live": "C,C,c,0.00,0.00",
    "edge_opacity_headless": 0.6,
    "edge_transparency_live": 120,
    "standard": {
        "node_border_width": 2.0,
        "node_border_color": "#000000",
        "edge_width": 2.5,
        "edge_color": "#888888",
        "edge_style": "SOLID",
    },
    "linked_identity": {
        "node_border_width": 13.0,
        "node_border_color": "#555555",
        "node_border_transparency_headless": 110,
        "node_border_transparency_identity_headless": 245,
        "node_border_transparency_live": 110,
        "node_border_transparency_identity_live": 245,
        "edge_width_standard": 2.5,
        "edge_width_identity": 6.0,
        "edge_color_standard": "#888888",
        "edge_color_identity": "#FF0000",
        "edge_style_standard": "SOLID",
        "edge_style_standard_headless": "solid",
        "edge_style_identity_headless": "dotted",
        "edge_style_identity_live": "DOT",
        "edge_opacity_identity_headless": 1.0,
        "edge_opacity_identity_live": 255,
    },
}

# Backward-compatible aliases used throughout the file
LINKED_IDENTITY_BORDER_WIDTH = VISUAL_TUNING["linked_identity"]["node_border_width"]
LINKED_IDENTITY_EDGE_WIDTH = VISUAL_TUNING["linked_identity"]["edge_width_identity"]
STANDARD_EDGE_WIDTH = VISUAL_TUNING["standard"]["edge_width"]
LINKED_IDENTITY_EDGE_STYLE = VISUAL_TUNING["linked_identity"]["edge_style_identity_live"]
STANDARD_EDGE_STYLE = VISUAL_TUNING["standard"]["edge_style"]

# Unique tag per Python run so styles are always fresh across runs
STYLE_RUN_TAG = f"{int(time.time() * 1000)}_{os.getpid()}"


def _is_linked_identity_network(network_title: str) -> bool:
    return str(network_title).startswith("Combined_Network")


def _is_combined_protein_network(network_title: str) -> bool:
    return network_title == "Combined_Protein_Network"


def _get_network_visual_profile(network_title: str) -> dict:
    linked = _is_linked_identity_network(network_title)
    combined_protein = _is_combined_protein_network(network_title)

    standard = VISUAL_TUNING["standard"]
    linked_identity = VISUAL_TUNING["linked_identity"]

    if linked:
        return {
            "is_linked_identity_network": True,
            "is_combined_protein_network": False,
            "is_combined_network": True,
            "node_border_width": linked_identity["node_border_width"],
            "node_border_color": linked_identity["node_border_color"],
            "node_border_transparency_headless": linked_identity["node_border_transparency_headless"],
            "node_border_transparency_live": linked_identity["node_border_transparency_live"],
            "default_edge_width": linked_identity["edge_width_standard"],
            "default_edge_color": linked_identity["edge_color_standard"],
            "default_edge_style": linked_identity["edge_style_standard"],
            "default_edge_style_headless": linked_identity["edge_style_standard_headless"],
            "identity_edge_width": linked_identity["edge_width_identity"],
            "identity_edge_color": linked_identity["edge_color_identity"],
            "identity_edge_style": linked_identity["edge_style_identity_live"],
            "identity_edge_style_headless": linked_identity["edge_style_identity_headless"],
            "default_edge_opacity_headless": VISUAL_TUNING["edge_opacity_headless"],
            "identity_edge_opacity_headless": linked_identity["edge_opacity_identity_headless"],
            "combined_default_node_fill": VISUAL_TUNING["combined_default_node_fill"],
        }

    return {
        "is_linked_identity_network": False,
        "is_combined_protein_network": combined_protein,
        "is_combined_network": combined_protein,
        "node_border_width": standard["node_border_width"],
        "node_border_color": standard["node_border_color"],
        "node_border_transparency_headless": 255,
        "node_border_transparency_live": 255,
        "default_edge_width": standard["edge_width"],
        "default_edge_color": standard["edge_color"],
        "default_edge_style": standard["edge_style"],
        "default_edge_style_headless": standard["edge_style"].lower(),
        "identity_edge_width": standard["edge_width"],
        "identity_edge_color": standard["edge_color"],
        "identity_edge_style": standard["edge_style"],
        "identity_edge_style_headless": standard["edge_style"].lower(),
        "default_edge_opacity_headless": VISUAL_TUNING["edge_opacity_headless"],
        "identity_edge_opacity_headless": VISUAL_TUNING["edge_opacity_headless"],
        "combined_default_node_fill": VISUAL_TUNING["combined_default_node_fill"],
    }


def _build_color_map(color_groups, network_title: str = "") -> dict:
    if _is_linked_identity_network(network_title):
        ordered_groups = [str(g) for g in color_groups]
        return {
            group: COMBINED_FILE_PALETTE[i % len(COMBINED_FILE_PALETTE)]
            for i, group in enumerate(ordered_groups)
        }

    base_color_groups = [g for g in color_groups if g not in NODE_COLOR_FIXED and g != "Multi"]
    cmap = cm.get_cmap("tab20", max(1, len(base_color_groups)))
    auto_colors = {group: to_hex(cmap(i)) for i, group in enumerate(base_color_groups)}
    if "Multi" in color_groups:
        auto_colors["Multi"] = MULTI_NODE_COLOR
    return {**NODE_COLOR_FIXED, **auto_colors}


def _normalize_uniprot_id(value) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return None
    return text



def _hex_from_hsv(h: float, s: float, v: float) -> str:
    r, g, b = colorsys.hsv_to_rgb(h, s, v)
    return f"#{int(round(r * 255)):02x}{int(round(g * 255)):02x}{int(round(b * 255)):02x}"



def _build_linked_identity_uniprot_color_map(uniprot_ids) -> dict:
    ordered_ids = sorted({_normalize_uniprot_id(uid) for uid in uniprot_ids if _normalize_uniprot_id(uid)})
    if not ordered_ids:
        return {}

    if len(ordered_ids) == 1:
        return {ordered_ids[0]: _hex_from_hsv(0.60, 0.72, 0.82)}

    hue_start = 0.10
    hue_end = 0.82
    span = hue_end - hue_start

    color_map = {}
    for i, up_id in enumerate(ordered_ids):
        hue = hue_start + (span * i / (len(ordered_ids) - 1))
        color_map[up_id] = _hex_from_hsv(hue, 0.72, 0.82)
    return color_map



def _annotate_linked_identity_node_borders(nodes_df: pd.DataFrame, edges_df: pd.DataFrame) -> pd.DataFrame:
    """Color all linked-relevant UniProt groups red and assign distinct colors to the remaining UniProt groups."""
    if nodes_df is None or nodes_df.empty:
        return nodes_df

    df = nodes_df.copy()

    if "uniprot_border_color" not in df.columns:
        df["uniprot_border_color"] = "#555555"

    if "uniprot_id" in df.columns:
        df["uniprot_id"] = df["uniprot_id"].apply(_normalize_uniprot_id)
    else:
        df["uniprot_id"] = None

    regular_uniprot_color_map = _build_linked_identity_uniprot_color_map(df["uniprot_id"].dropna().tolist())

    df["linked_identity_border_color"] = df["uniprot_id"].map(regular_uniprot_color_map).fillna("#555555")
    df["linked_identity_border_transparency"] = VISUAL_TUNING["linked_identity"]["node_border_transparency_headless"]

    linked_uniprot_ids = set()
    if edges_df is not None and not edges_df.empty and "interaction" in edges_df.columns:
        identity_rows = edges_df[edges_df["interaction"].astype(str) == "identity"]
        if not identity_rows.empty:
            id_to_uniprot = {
                str(node_id): up_id
                for node_id, up_id in zip(df["id"].astype(str), df["uniprot_id"])
            }
            for node_id in pd.concat([identity_rows["source"], identity_rows["target"]]).astype(str):
                up_id = id_to_uniprot.get(node_id)
                if up_id:
                    linked_uniprot_ids.add(up_id)

    if linked_uniprot_ids:
        linked_mask = df["uniprot_id"].isin(linked_uniprot_ids)
        df.loc[linked_mask, "linked_identity_border_color"] = VISUAL_TUNING["linked_identity"]["edge_color_identity"]
        df.loc[linked_mask, "linked_identity_border_transparency"] = VISUAL_TUNING["linked_identity"]["node_border_transparency_identity_headless"]

    return df


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

        for i, comp in enumerate(components):
            sub_G = G.subgraph(comp)
            sub_N = len(comp)

            # Local layout for this island
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

            # Approximate stride
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
    profile = _get_network_visual_profile(network_title)

    # --- Node Mapping ---
    node_ids = list(nodes_df["id"].astype(str))
    nid_map = {nid: i + 1 for i, nid in enumerate(node_ids)}

    nodes_aspect = []
    for _, row in nodes_df.iterrows():
        nid = str(row["id"])
        pos = positions.get(nid, {"x": 0.0, "y": 0.0})
        border_color = (
            str(row.get("linked_identity_border_color", row.get("uniprot_border_color", profile["node_border_color"])))
            if profile["is_linked_identity_network"]
            else profile["node_border_color"]
        )
        border_transparency = (
            int(row.get("linked_identity_border_transparency", profile["node_border_transparency_headless"]))
            if profile["is_linked_identity_network"]
            else profile["node_border_transparency_headless"]
        )

        nodes_aspect.append(
            {
                "id": nid_map[nid],
                "x": float(pos["x"]),
                "y": float(pos["y"]),
                "v": {
                    "name": str(row.get("name", nid)),
                    "tooltip": str(row.get("tooltip", "")),
                    "color_group": str(row.get("color_group", "Unknown")),
                    "uniprot_border_color": str(row.get("uniprot_border_color", border_color)),
                    "linked_identity_border_color": border_color,
                    "linked_identity_border_transparency": border_transparency,
                },
            }
        )

    # --- Edge Mapping ---
    edges_aspect = []
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

    discrete_map = [{"v": k, "vp": v} for k, v in color_map.items()]

    visual_props = {
        "visualProperties": [
            {
                "default": {
                    "network": {
                        "NETWORK_BACKGROUND_COLOR": VISUAL_TUNING["network_background_color"]
                    },
                    "node": {
                        "NODE_SHAPE": VISUAL_TUNING["node_shape"],
                        "NODE_WIDTH": VISUAL_TUNING["node_width"],
                        "NODE_HEIGHT": VISUAL_TUNING["node_height"],
                        "NODE_BORDER_WIDTH": profile["node_border_width"],
                        "NODE_BORDER_TRANSPARENCY": profile["node_border_transparency_headless"],
                        "NODE_BORDER_COLOR": profile["node_border_color"],
                    },
                    "edge": {
                        "EDGE_LINE_COLOR": profile["default_edge_color"],
                        "EDGE_OPACITY": profile["default_edge_opacity_headless"],
                        "EDGE_CURVED": True,
                        "EDGE_LINE_STYLE": profile["default_edge_style_headless"],
                        "EDGE_WIDTH": profile["default_edge_width"],
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
                "edgeMapping": {},
            }
        ]
    }

    if profile["is_linked_identity_network"]:
        visual_props["visualProperties"][0]["nodeMapping"]["NODE_BORDER_COLOR"] = {
            "type": "PASSTHROUGH",
            "definition": {"attribute": "linked_identity_border_color", "type": "string"},
        }
        visual_props["visualProperties"][0]["nodeMapping"]["NODE_BORDER_TRANSPARENCY"] = {
            "type": "PASSTHROUGH",
            "definition": {"attribute": "linked_identity_border_transparency", "type": "integer"},
        }
        visual_props["visualProperties"][0]["edgeMapping"]["EDGE_LINE_STYLE"] = {
            "type": "DISCRETE",
            "definition": {
                "attribute": "interaction",
                "type": "string",
                "map": [
                    {"v": "identity", "vp": profile["identity_edge_style_headless"]},
                    {"v": "interacts_with", "vp": profile["default_edge_style_headless"]},
                ],
            },
        }
        visual_props["visualProperties"][0]["edgeMapping"]["EDGE_LINE_COLOR"] = {
            "type": "DISCRETE",
            "definition": {
                "attribute": "interaction",
                "type": "string",
                "map": [
                    {"v": "identity", "vp": profile["identity_edge_color"]},
                    {"v": "interacts_with", "vp": profile["default_edge_color"]},
                ],
            },
        }
        visual_props["visualProperties"][0]["edgeMapping"]["EDGE_WIDTH"] = {
            "type": "DISCRETE",
            "definition": {
                "attribute": "interaction",
                "type": "string",
                "map": [
                    {"v": "identity", "vp": profile["identity_edge_width"]},
                    {"v": "interacts_with", "vp": profile["default_edge_width"]},
                ],
            },
        }
        visual_props["visualProperties"][0]["edgeMapping"]["EDGE_OPACITY"] = {
            "type": "DISCRETE",
            "definition": {
                "attribute": "interaction",
                "type": "string",
                "map": [
                    {"v": "identity", "vp": profile["identity_edge_opacity_headless"]},
                    {"v": "interacts_with", "vp": profile["default_edge_opacity_headless"]},
                ],
            },
        }

    attr_decls = {
        "attributeDeclarations": [
            {
                "networkAttributes": {"name": {"d": "string"}},
                "nodes": {
                    "name": {"d": "string"},
                    "tooltip": {"d": "string"},
                    "color_group": {"d": "string"},
                    "uniprot_border_color": {"d": "string"},
                    "linked_identity_border_color": {"d": "string"},
                    "linked_identity_border_transparency": {"d": "integer"},
                },
                "edges": {
                    "interaction": {"d": "string"},
                    "all_atoms_count": {"d": "integer"},
                },
            }
        ]
    }

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

    def _apply_identity_edge_live_bypasses(network_suid: int, edges_full_df: pd.DataFrame) -> None:
        """Make identity edges thick, dotted, and fully opaque in Cytoscape live mode."""
        if edges_full_df is None or edges_full_df.empty:
            return

        try:
            identity_rows = edges_full_df[edges_full_df["interaction"].astype(str) == "identity"]
        except Exception:
            return

        if identity_rows.empty:
            return

        edge_names = [
            _compute_edge_name(str(s), str(t), "identity")
            for s, t in zip(identity_rows["source"], identity_rows["target"])
        ]
        widths = [LINKED_IDENTITY_EDGE_WIDTH] * len(edge_names)
        styles = [LINKED_IDENTITY_EDGE_STYLE] * len(edge_names)
        opacities = [VISUAL_TUNING["linked_identity"]["edge_opacity_identity_live"]] * len(edge_names)

        try:
            p4c.set_edge_line_width_bypass(edge_names, widths, network=network_suid)
        except Exception:
            try:
                p4c.set_edge_line_width_bypass(edge_names, widths)
            except Exception:
                if _verbose_enabled():
                    print("[cytoscape] Identity width bypass could not be applied; continuing without it.")

        try:
            p4c.set_edge_line_style_bypass(edge_names, styles, network=network_suid)
        except Exception:
            try:
                p4c.set_edge_line_style_bypass(edge_names, styles)
            except Exception:
                if _verbose_enabled():
                    print("[cytoscape] Identity line-style bypass could not be applied; continuing without it.")

        try:
            p4c.set_edge_opacity_bypass(edge_names, opacities, network=network_suid)
        except Exception:
            try:
                p4c.set_edge_opacity_bypass(edge_names, opacities)
            except Exception:
                if _verbose_enabled():
                    print("[cytoscape] Identity opacity bypass could not be applied; continuing without it.")

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

    if _is_linked_identity_network(network_title):
        nodes_df = _annotate_linked_identity_node_borders(nodes_df, edges_df)

    color_groups = (
        sorted(nodes_df["color_group"].dropna().unique())
        if "color_group" in nodes_df.columns
        else []
    )
    color_map = _build_color_map(color_groups, network_title)
    profile = _get_network_visual_profile(network_title)

    positions = compute_preset_positions_spring(nodes_df, edges_df, network_title)

    if not config.get("open_in_cytoscape", True):
        _export_cx2_headless(network_title, run_output_path, nodes_df, edges_df, color_map, positions)
        return

    is_linked_identity_network = profile["is_linked_identity_network"]
    is_combined_protein_network = profile["is_combined_protein_network"]
    is_combined_network = profile["is_combined_network"]

    existing_networks = p4c.get_network_list()
    while len(existing_networks) > config["keep_last_n_networks"]:
        p4c.delete_network(existing_networks.pop(0))

    title_lower = str(network_title).lower()
    if "combined" in title_lower:
        collection_name = "PDB2Net — Combined"
    elif "protein" in title_lower:
        collection_name = "PDB2Net — Protein"
    else:
        collection_name = "PDB2Net — Chain"

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

    if is_combined_network:
        _load_edge_attributes_by_name(network_suid, edges_df)

    if is_linked_identity_network:
        style_name = (
            f"PDB2Net_Linked_Identity_{STYLE_RUN_TAG}_"
            f"{_style_signature(color_map, extra_tag='linked', border_width=profile['node_border_width'])}"
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

    defaults = {
        "NODE_SHAPE": VISUAL_TUNING["node_shape"].upper(),
        "NODE_SIZE": VISUAL_TUNING["node_size_live"],
        "NODE_LABEL_POSITION": VISUAL_TUNING["node_label_position_live"],
        "EDGE_TRANSPARENCY": VISUAL_TUNING["edge_transparency_live"],
        "EDGE_STROKE_UNSELECTED_PAINT": profile["default_edge_color"],
        "EDGE_WIDTH": profile["default_edge_width"],
        "NODE_BORDER_WIDTH": profile["node_border_width"],
        "NODE_BORDER_PAINT": profile["node_border_color"],
        "NODE_BORDER_TRANSPARENCY": profile["node_border_transparency_live"],
    }

    if is_linked_identity_network or is_combined_protein_network:
        defaults["NODE_FILL_COLOR"] = profile["combined_default_node_fill"]

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
                "mappingColumn": "linked_identity_border_color",
                "mappingColumnType": "String",
                "visualProperty": "NODE_BORDER_PAINT",
            }
        )
        mappings.append(
            {
                "mappingType": "passthrough",
                "mappingColumn": "linked_identity_border_transparency",
                "mappingColumnType": "Integer",
                "visualProperty": "NODE_BORDER_TRANSPARENCY",
            }
        )
        mappings.append(
            {
                "mappingType": "discrete",
                "mappingColumn": "interaction",
                "mappingColumnType": "String",
                "visualProperty": "EDGE_STROKE_UNSELECTED_PAINT",
                "map": [
                    {"key": "identity", "value": profile["identity_edge_color"]},
                    {"key": "interacts_with", "value": profile["default_edge_color"]},
                ],
            }
        )

    _ensure_style(style_name, mappings=mappings, defaults=defaults)
    p4c.set_visual_style(style_name)

    if is_linked_identity_network:
        _apply_identity_edge_live_bypasses(network_suid, edges_df)

    if is_combined_network:
        try:
            p4c.layout_network(layout_name="force-directed")
        except Exception as e:
            if _verbose_enabled():
                print(f"[cytoscape] Layout failed (non-fatal): {e}")

    try:
        _export_cx2_headless(network_title, run_output_path, nodes_df, edges_df, color_map, positions)
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