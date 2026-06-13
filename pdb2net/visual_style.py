"""Visual style helpers shared by headless CX2 and live Cytoscape export."""

from __future__ import annotations

import colorsys

import pandas as pd
from matplotlib import cm, colormaps
from matplotlib.colors import to_hex


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

LINKED_IDENTITY_EDGE_WIDTH = VISUAL_TUNING["linked_identity"]["edge_width_identity"]
LINKED_IDENTITY_EDGE_STYLE = VISUAL_TUNING["linked_identity"]["edge_style_identity_live"]


def is_linked_identity_network(network_title: str) -> bool:
    return str(network_title).startswith("Combined_Network")


def is_combined_protein_network(network_title: str) -> bool:
    return network_title == "Combined_Protein_Network"


def get_network_visual_profile(network_title: str) -> dict:
    linked = is_linked_identity_network(network_title)
    combined_protein = is_combined_protein_network(network_title)

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


def build_color_map(color_groups, network_title: str = "") -> dict:
    if is_linked_identity_network(network_title):
        ordered_groups = [str(group) for group in color_groups]
        return {
            group: COMBINED_FILE_PALETTE[index % len(COMBINED_FILE_PALETTE)]
            for index, group in enumerate(ordered_groups)
        }

    base_color_groups = [group for group in color_groups if group not in NODE_COLOR_FIXED and group != "Multi"]
    color_count = max(1, len(base_color_groups))
    try:
        cmap = colormaps.get_cmap("tab20").resampled(color_count)
    except AttributeError:
        cmap = cm.get_cmap("tab20", color_count)
    auto_colors = {group: to_hex(cmap(index)) for index, group in enumerate(base_color_groups)}
    if "Multi" in color_groups:
        auto_colors["Multi"] = MULTI_NODE_COLOR
    return {**NODE_COLOR_FIXED, **auto_colors}


def normalize_uniprot_id(value) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return None
    return text


def _hex_from_hsv(h: float, s: float, v: float) -> str:
    r, g, b = colorsys.hsv_to_rgb(h, s, v)
    return f"#{int(round(r * 255)):02x}{int(round(g * 255)):02x}{int(round(b * 255)):02x}"


def build_linked_identity_uniprot_color_map(uniprot_ids) -> dict:
    ordered_ids = sorted({normalize_uniprot_id(uid) for uid in uniprot_ids if normalize_uniprot_id(uid)})
    if not ordered_ids:
        return {}

    if len(ordered_ids) == 1:
        return {ordered_ids[0]: _hex_from_hsv(0.60, 0.72, 0.82)}

    hue_start = 0.10
    hue_end = 0.82
    span = hue_end - hue_start

    color_map = {}
    for index, up_id in enumerate(ordered_ids):
        hue = hue_start + (span * index / (len(ordered_ids) - 1))
        color_map[up_id] = _hex_from_hsv(hue, 0.72, 0.82)
    return color_map


def annotate_linked_identity_node_borders(nodes_df: pd.DataFrame, edges_df: pd.DataFrame) -> pd.DataFrame:
    """Color identity-linked UniProt groups red and color remaining UniProt groups distinctly."""
    if nodes_df is None or nodes_df.empty:
        return nodes_df

    df = nodes_df.copy()

    if "uniprot_border_color" not in df.columns:
        df["uniprot_border_color"] = "#555555"

    if "uniprot_id" in df.columns:
        df["uniprot_id"] = df["uniprot_id"].apply(normalize_uniprot_id)
    else:
        df["uniprot_id"] = None

    regular_uniprot_color_map = build_linked_identity_uniprot_color_map(df["uniprot_id"].dropna().tolist())

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
