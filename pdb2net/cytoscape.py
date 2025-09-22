import os
import pandas as pd
from config_loader import config
from matplotlib import cm
from matplotlib.colors import to_hex
import json
import time
import py4cytoscape as p4c
import subprocess
from config_loader import config
import re

def ensure_cytoscape_running():
    try:
        p4c.cytoscape_ping()
        print("Cytoscape is already running.")
    except:
        print("Starting Cytoscape...")
        subprocess.Popen(config["cytoscape_path"])
        time.sleep(40)
        try:
            p4c.cytoscape_ping()
            print("Cytoscape started successfully.")
        except:
            print("Error: Cytoscape could not be started. Check the path in config.json.")
            exit(1)

def compute_preset_positions_spring(nodes_df, edges_df, network_title, scale=1000.0):
    """
    Schnelles, deterministisches Force-Directed-Layout (Fruchterman-Reingold),
    ohne nachträgliche Min/Max-Normalisierung (dadurch natürlichere Abstände).
    Gibt {node_id: {"x":..., "y":...}} zurück.
    """
    import math, zlib

    node_ids = list(nodes_df["id"]) if len(nodes_df) else []
    N = len(node_ids)
    if N == 0:
        return {}
    if N == 1:
        return {node_ids[0]: {"x": 0.0, "y": 0.0}}
    if N == 2:
        d = scale * 0.35
        return {node_ids[0]: {"x": -d, "y": 0.0}, node_ids[1]: {"x": d, "y": 0.0}}

    try:
        import networkx as nx
    except Exception:
        # Fallback: Kreislayout
        r = scale * 0.35
        out = {}
        for i, nid in enumerate(node_ids):
            a = 2.0 * math.pi * i / N
            out[nid] = {"x": float(r * math.cos(a)), "y": float(r * math.sin(a))}
        return out

    # Graph bauen
    G = nx.Graph()
    G.add_nodes_from(node_ids)

    # Kanten mit "sanfter" Gewichtung (log-Kompression vermeidet extreme Zerrungen)
    if edges_df is not None and len(edges_df) > 0:
        for _, e in edges_df.iterrows():
            s, t = e["source"], e["target"]
            if s in G and t in G and s != t:
                w = float(e.get("all_atoms_count", 1.0))
                w = 1.0 + math.log10(max(1.0, w))  # 1.. ~
                if G.has_edge(s, t):
                    if w > G[s][t].get("weight", 1.0):
                        G[s][t]["weight"] = w
                else:
                    G.add_edge(s, t, weight=w)

    # Reproduzierbarer Seed
    seed = zlib.adler32(str(network_title).encode("utf-8")) & 0xFFFFFFFF

    # Iterationen: etwas höher für stabilere Konstellation bei kleinen Netzen
    iters = max(100, min(250, 8 * N))

    # WICHTIG: NetworkX skaliert selbst; wir geben scale/2 vor und zentrieren auf (0,0)
    # -> natürliche Distanzen bleiben erhalten, keine Nach-Normalisierung nötig
    pos = nx.spring_layout(
        G,
        seed=seed,
        weight="weight",
        iterations=iters,
        dim=2,
        center=(0.0, 0.0),
        scale=scale * 0.5,
        k=None  # Auto-K; für N<50 idR gut
    )

    return {n: {"x": float(x), "y": float(y)} for n, (x, y) in pos.items()}


def export_custom_cyjs(nodes_df, edges_df, network_title, run_output_path):
    """
    Saves a Cytoscape-compatible .cyjs file without using Cytoscape.
    Includes molecule_name, tooltip, color_group etc. in the node attributes.
    """
    import json
    os.makedirs(run_output_path, exist_ok=True)
    network_file = os.path.join(run_output_path, f"{network_title}.cyjs")

    nodes = []
    for _, row in nodes_df.iterrows():
        node = {
            "data": {
                "id": row["id"],
                "name": row["name"],
                "molecule_name": row.get("molecule_name", "Unknown")
            }
        }
        if "color_group" in row:
            node["data"]["color_group"] = row["color_group"]
        if "tooltip" in row:
            node["data"]["tooltip"] = row["tooltip"]
        nodes.append(node)

    edges = []
    for idx, row in edges_df.iterrows():
        edge = {
            "data": {
                "id": f"edge_{idx}",
                "source": row["source"],
                "target": row["target"],
                "interaction": row.get("interaction", "interacts_with"),
                "all_atoms_count": row.get("all_atoms_count", 0)
            }
        }
        edges.append(edge)

    cyjs_data = {
        "data": {"name": network_title},
        "elements": {"nodes": nodes, "edges": edges}
    }

    with open(network_file, "w") as f:
        json.dump(cyjs_data, f, indent=2)

def create_cytoscape_network(results, network_title="Protein_Interaction_Network", run_output_path=".", nodes_data=None):
    import py4cytoscape as p4c
    import pandas as pd
    import json, os
    from matplotlib import cm
    from matplotlib.colors import to_hex

    # --- Edges aus results sammeln ---
    unique_nodes = set()
    edges = []
    for entry in results:
        if entry.get("all_atoms_count", 0) > 0:
            a = entry["chain_a"]; b = entry["chain_b"]
            unique_nodes.add(a); unique_nodes.add(b)
            edges.append({"chain_a": a, "chain_b": b, "all_atoms_count": entry["all_atoms_count"]})

    # --- Nodes-DF ---
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

    # --- Edges-DF (Export & Cytoscape) ---
    if edges:
        edges_df_for_export = pd.DataFrame(edges).rename(columns={"chain_a": "source", "chain_b": "target"})
        edges_df_for_export["interaction"] = "interacts_with"
        edges_df_for_cyto = edges_df_for_export
    else:
        edges_df_for_export = pd.DataFrame(columns=["source", "target", "interaction", "all_atoms_count"])
        edges_df_for_cyto = None

    # === PFAD 1: Headless → schreibe .cyjs mit JS-Style + preset ===
    if not config.get("open_in_cytoscape", True):
        os.makedirs(run_output_path, exist_ok=True)
        network_file = os.path.join(run_output_path, f"{network_title}.cyjs")

        fixed_colors = {
            "Protein": "#1f77b4",
            "DNA": "#ff7f0e",
            "RNA": "#2ca02c",
            "DNA/RNA": "#a2a200",
            "Nucleic Acid": "#9467bd",
            "Unknown": "#7f7f7f"
        }
        color_groups = sorted(nodes_df["color_group"].dropna().unique()) if "color_group" in nodes_df.columns else []
        base_color_groups = [g for g in color_groups if g not in fixed_colors and g != "Multi"]
        cmap = cm.get_cmap('tab20', max(1, len(base_color_groups)))
        auto_colors = {group: to_hex(cmap(i)) for i, group in enumerate(base_color_groups)}
        if "Multi" in color_groups and "Multi" not in fixed_colors:
            auto_colors["Multi"] = "#FF0000"
        color_map = {**fixed_colors, **auto_colors}

        positions = compute_preset_positions_spring(nodes_df, edges_df_for_export, network_title, scale=1000.0)

        nodes_elems = []
        for _, row in nodes_df.iterrows():
            nid = row["id"]
            cg = row.get("color_group", "Unknown")
            node_color = color_map.get(cg, "#7f7f7f")
            nodes_elems.append({
                "data": {
                    "id": nid,
                    "name": row.get("name", nid),
                    "tooltip": row.get("tooltip", ""),
                    "color_group": cg,
                    "node_color": node_color
                },
                "position": positions.get(nid, {"x": 0.0, "y": 0.0})
            })

        edges_elems = []
        max_w = 1
        for _, r in edges_df_for_export.iterrows():
            try:
                max_w = max(max_w, int(r.get("all_atoms_count", 1)))
            except Exception:
                pass
        for idx, row in edges_df_for_export.iterrows():
            edges_elems.append({
                "data": {
                    "id": f"edge_{idx}",
                    "source": row["source"],
                    "target": row["target"],
                    "interaction": row.get("interaction", "interacts_with"),
                    "all_atoms_count": int(row.get("all_atoms_count", 1))
                }
            })

        style = [
            {"selector": "node", "style": {
                "shape": "ellipse",
                "width": 40, "height": 40,
                "background-color": "data(node_color)",
                "label": "data(name)",
                "font-size": 10,
                "text-valign": "center", "text-halign": "center",
                "border-width": 1, "border-color": "#555"
            }},
            {"selector": "edge", "style": {
                "curve-style": "bezier",
                "opacity": 0.6,
                "width": f"mapData(all_atoms_count, 1, {max_w}, 1, 6)",
                "line-color": "#888"
            }}
        ]

        cyjs_data = {
            "data": {"name": network_title},
            "elements": {"nodes": nodes_elems, "edges": edges_elems},
            "style": style,
            "layout": {"name": "preset"}
        }

        with open(network_file, "w", encoding="utf-8") as f:
            json.dump(cyjs_data, f, indent=2)
        return

    # === PFAD 2: Cytoscape offen → Desktop-Workflow + Datei-Exporte ===
    existing_networks = p4c.get_network_list()
    while len(existing_networks) > config["keep_last_n_networks"]:
        oldest_network = existing_networks.pop(0)
        p4c.delete_network(oldest_network)

    p4c.create_network_from_data_frames(nodes_df, edges_df_for_cyto, title=network_title)

    # Node-Tabellenspalten nachladen (Label/Tooltip/Farbgruppe)
    try:
        if "color_group" in nodes_df.columns:
            p4c.load_table_data(
                data=nodes_df[["id", "color_group"]],
                data_key_column="id",
                table="node", table_key_column="id"
            )
        p4c.load_table_data(
            data=nodes_df[["id", "name"]],
            data_key_column="id",
            table="node", table_key_column="id"
        )
        if "tooltip" in nodes_df.columns:
            p4c.load_table_data(
                data=nodes_df[["id", "tooltip"]],
                data_key_column="id",
                table="node", table_key_column="id"
            )
    except Exception as e:
        print(f"Error while loading node data: {e}")

    # Visual Style erstellen/anwenden + Layout
    try:
        color_groups = sorted(nodes_df["color_group"].dropna().unique()) if "color_group" in nodes_df.columns else []
        is_chain_network = "Chain" in network_title or network_title == "Combined_Network"
        is_protein_network = "Protein" in network_title
        is_combined_protein = is_protein_network and "combined" in network_title.lower()

        style_name = (
            "PDB2Net_Chain_Style" if is_chain_network else
            "PDB2Net_Protein_Combined_Style" if is_combined_protein else
            "PDB2Net_Protein_Style"
        )

        fixed_colors = {
            "Protein": "#1f77b4",
            "DNA": "#ff7f0e",
            "RNA": "#2ca02c",
            "DNA/RNA": "#a2a200",
            "Nucleic Acid": "#9467bd"
        }
        base_color_groups = [g for g in color_groups if g not in fixed_colors and g != "Multi"]
        cmap = cm.get_cmap('tab20', len(base_color_groups))
        auto_colors = {group: to_hex(cmap(i)) for i, group in enumerate(base_color_groups)}
        if is_combined_protein and "Multi" in color_groups:
            auto_colors["Multi"] = "#FF0000"
        color_map = {**auto_colors, **fixed_colors}

        if style_name not in p4c.get_visual_style_names():
            defaults = {
                "NODE_SHAPE": "ELLIPSE",
                "NODE_SIZE": 40,
                "NODE_LABEL_POSITION": "C,C,c,0.00,0.00",
                "EDGE_TRANSPARENCY": 120
            }
            mappings = [
                {"mappingType": "passthrough", "mappingColumn": "name",    "mappingColumnType": "String", "visualProperty": "NODE_LABEL"},
                {"mappingType": "passthrough", "mappingColumn": "tooltip", "mappingColumnType": "String", "visualProperty": "NODE_TOOLTIP"},
                {
                    "mappingType": "discrete",
                    "mappingColumn": "color_group",
                    "mappingColumnType": "String",
                    "visualProperty": "NODE_FILL_COLOR",
                    "map": [{"key": k, "value": v} for k, v in color_map.items()]
                }
            ]
            p4c.create_visual_style(style_name, mappings=mappings, defaults=defaults)

        p4c.set_current_network(network_title)
        p4c.set_visual_style(style_name)
        p4c.layout_network(layout_name="force-directed")
    except Exception as e:
        print(f"Error while applying style: {e}")

    # === NEU: Immer exportieren (Desktop-kompatibel) ===
    try:
        os.makedirs(run_output_path, exist_ok=True)
        cyjs_desktop = os.path.join(run_output_path, f"{network_title}.desktop.cyjs")
        cx2_file     = os.path.join(run_output_path, f"{network_title}.cx2")

        # CyJS mit angewendetem Desktop-Style
        p4c.export_network(cyjs_desktop, type="cyjs")

        # Robust: CX2 (falls 'cx2' nicht akzeptiert wird, auf 'cx'/'CX' fallbacken)
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
    """
    Chain-Nodes für Cytoscape erzeugen.
    - Label (name): wieder klassisch = id ('PDBID:Chain')
    - Tooltip: reichhaltig (voller Name, Typ, Länge, PDB, UniProt)
    - color_group: Molekültyp (Protein, DNA, RNA, DNA/RNA, Nucleic Acid, Unknown)
    """
    dna_set = {"DA", "DT", "DG", "DC", "DI"}
    rna_set = {"A", "U", "G", "C", "I"}
    protein_set = {
        "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS",
        "ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP",
        "TYR","VAL","SEC","PYL"
    }

    def count_lengths(res_list):
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
        uid = chain.get("unique_chain_id") or chain.get("id")  # 'PDBID:Chain'
        mol_type = (chain.get("molecule_type") or "Unknown").strip()
        mol_name_full = chain.get("molecule_name") or "Unknown"
        up_id = chain.get("uniprot_id")  # optional
        aa_len, nt_len = count_lengths(chain.get("residues"))

        # Tooltip (reichhaltig)
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
            "name": uid,                   # << altes Label beibehalten
            "tooltip": tooltip,            # << neue, reichhaltige Tooltips
            "color_group": mol_type or "Unknown",
            "molecule_name": mol_name_full
        })

    return nodes
