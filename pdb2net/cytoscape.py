"""Cytoscape integration helpers for PDB2Net.

This module handles:
- Launching/connecting to Cytoscape (via py4cytoscape).
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
import zlib
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "pdb2net-matplotlib"))

import pandas as pd

from .config_loader import config
from .cx2_export import export_cx2_headless
from .layout_engine import calculate_positions
from .network_annotations import annotation_node_metadata
from .residue_types import count_polymer_lengths
from .visual_style import (
    LINKED_IDENTITY_EDGE_STYLE,
    LINKED_IDENTITY_EDGE_WIDTH,
    VISUAL_TUNING,
    annotate_linked_identity_node_borders,
    build_color_map,
    get_network_visual_profile,
    is_linked_identity_network,
)

# Unique tag per Python run so styles are always fresh across runs
STYLE_RUN_TAG = f"{int(time.time() * 1000)}_{os.getpid()}"


def _get_py4cytoscape():
    import py4cytoscape as p4c

    return p4c


def ensure_cytoscape_running():
    """Ensure a Cytoscape instance is reachable via CyREST."""
    p4c = _get_py4cytoscape()
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
        # Preserve the distinction between an absent optional field and an
        # explicitly supplied non-finite number.  Pandas otherwise turns both
        # into NaN before the CX2 serializer can enforce finite numerics.
        node_fields = sorted(
            {str(field) for node in nodes_data for field in node}
        )
        normalized_nodes = [
            {field: node[field] if field in node else "" for field in node_fields}
            for node in nodes_data
        ]
        nodes_df = pd.DataFrame(normalized_nodes).copy()
        if "name" not in nodes_df.columns:
            nodes_df["name"] = nodes_df["id"]
        else:
            nodes_df["name"] = nodes_df["name"].where(nodes_df["name"].notna(), nodes_df["id"])
        if "tooltip" not in nodes_df.columns:
            nodes_df["tooltip"] = nodes_df.get("molecule_name", nodes_df["name"])
        else:
            nodes_df["tooltip"] = nodes_df["tooltip"].where(nodes_df["tooltip"].notna(), nodes_df["name"])
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

    if is_linked_identity_network(network_title):
        nodes_df = annotate_linked_identity_node_borders(nodes_df, edges_df)

    color_groups = (
        sorted(nodes_df["color_group"].dropna().unique())
        if "color_group" in nodes_df.columns
        else []
    )
    color_map = build_color_map(color_groups, network_title)
    profile = get_network_visual_profile(network_title)

    positions = calculate_positions(
        nodes_df=nodes_df,
        edges_df=edges_df,
        network_title=network_title,
        layout_mode=config.get("layout_mode", "python_fast"),
    )

    if not config.get("open_in_cytoscape", True):
        export_cx2_headless(network_title, run_output_path, nodes_df, edges_df, color_map, positions)
        return

    p4c = _get_py4cytoscape()
    linked_identity = profile["is_linked_identity_network"]
    combined_protein = profile["is_combined_protein_network"]
    combined_network = profile["is_combined_network"]

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
        # The portable CX2 artifact is the output contract even when the
        # optional live Cytoscape integration is unavailable.
        export_cx2_headless(network_title, run_output_path, nodes_df, edges_df, color_map, positions)
        return

    if combined_network:
        _load_edge_attributes_by_name(network_suid, edges_df)

    if linked_identity:
        style_name = (
            f"PDB2Net_Linked_Identity_{STYLE_RUN_TAG}_"
            f"{_style_signature(color_map, extra_tag='linked', border_width=profile['node_border_width'])}"
        )
    elif combined_protein:
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

    if linked_identity or combined_protein:
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

    if linked_identity:
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

    if linked_identity:
        _apply_identity_edge_live_bypasses(network_suid, edges_df)

    if combined_network:
        try:
            p4c.layout_network(layout_name="force-directed")
        except Exception as e:
            if _verbose_enabled():
                print(f"[cytoscape] Layout failed (non-fatal): {e}")

    # Do not turn serialization/disk failures into a successful run with
    # silently missing artifacts.
    export_cx2_headless(network_title, run_output_path, nodes_df, edges_df, color_map, positions)


def generate_nodes_from_atom_data(atom_data, pdb_id=None):
    """Create Cytoscape chain nodes from parsed atom/chain data."""
    def count_lengths(chain):
        if "aa_len" in chain or "nt_len" in chain:
            return int(chain.get("aa_len", 0)), int(chain.get("nt_len", 0))
        return count_polymer_lengths(
            (res.get("residue_name") for res in chain.get("residues") or [])
        )

    nodes = []
    for chain in atom_data:
        uid = chain.get("unique_chain_id") or chain.get("id")
        mol_type = (chain.get("molecule_type") or "Unknown").strip()
        mol_name_full = chain.get("molecule_name") or "Unknown"
        up_id = chain.get("uniprot_id")
        annotation_source = str(chain.get("annotation_source") or "")
        matched_database = str(chain.get("matched_database") or "")
        matched_id = str(chain.get("matched_id") or "")
        representative_accession = str(chain.get("representative_accession") or "")
        annotation_confidence = str(chain.get("annotation_confidence") or "")
        aa_len, nt_len = count_lengths(chain)
        chain_identity = chain.get("chain_identity", {})
        identity_display = (
            str(chain_identity.get("structure_display_id") or "")
            if isinstance(chain_identity, dict)
            else ""
        )
        node_pdb_id = str(pdb_id or chain.get("_parent_pdb_id") or identity_display or "")
        node_chain_id = str(chain.get("chain_id") or "")
        source_file = str(
            chain.get("_parent_file_label")
            or Path(str(chain.get("_parent_file_path") or "")).name
            or ""
        )

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
        if annotation_source:
            annotation_label = annotation_source
            if matched_database:
                annotation_label = f"{annotation_label} / {matched_database}"
            if annotation_confidence:
                annotation_label = f"{annotation_label} ({annotation_confidence})"
            details.append(f"Annotation: {annotation_label}")
        if matched_id:
            details.append(f"Matched ID: {matched_id}")
        if representative_accession and representative_accession != up_id:
            details.append(f"Representative accession: {representative_accession}")
        embedded_metadata = annotation_node_metadata([chain])
        details.extend(embedded_metadata.pop("tooltip_lines", []))
        tooltip = "\n".join(details)

        node = {
                "id": uid,
                "name": uid,
                "tooltip": tooltip,
                "color_group": mol_type or "Unknown",
                "molecule_name": mol_name_full,
                "pdb_id": node_pdb_id,
                "chain_id": node_chain_id,
                "source_file": source_file,
                "molecule_type": mol_type or "Unknown",
                "uniprot_id": up_id or "",
                "annotation_source": annotation_source,
                "matched_database": matched_database,
                "matched_id": matched_id,
                "representative_accession": representative_accession,
                "annotation_confidence": annotation_confidence,
                "aa_len": aa_len,
                "nt_len": nt_len,
                "node_kind": "chain",
            }
        if chain.get("embedded_annotation_source"):
            node["embedded_annotation_source"] = str(chain["embedded_annotation_source"])
        if chain.get("embedded_uniprot_status"):
            node["embedded_uniprot_status"] = str(chain["embedded_uniprot_status"])
        if chain.get("embedded_uniprot_accessions"):
            node["embedded_uniprot_accessions"] = json.dumps(
                sorted(str(value) for value in chain["embedded_uniprot_accessions"]),
                separators=(",", ":"),
            )
        for key, value in embedded_metadata.items():
            if key.startswith("annotation_"):
                node[key] = value
        nodes.append(node)

    return nodes
