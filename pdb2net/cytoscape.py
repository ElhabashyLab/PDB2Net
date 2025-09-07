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
    from matplotlib import cm
    from matplotlib.colors import to_hex
    import json
    import os

    # --- Edges sammeln (nur mit all_atoms_count > 0) ---
    unique_nodes = set()
    edges = []
    for entry in results:
        if entry.get("all_atoms_count", 0) > 0:
            chain_a = entry["chain_a"]
            chain_b = entry["chain_b"]
            unique_nodes.add(chain_a)
            unique_nodes.add(chain_b)
            edges.append({
                "chain_a": chain_a,
                "chain_b": chain_b,
                "all_atoms_count": entry["all_atoms_count"]
            })

    # --- Nodes-DataFrame herstellen ---
    if nodes_data:
        nodes_df = pd.DataFrame(nodes_data)
        # WICHTIG: vorhandene Felder NICHT überschreiben
        if "name" not in nodes_df.columns:
            nodes_df["name"] = nodes_df["id"]
        if "tooltip" not in nodes_df.columns:
            nodes_df["tooltip"] = nodes_df.get("molecule_name", "Unknown")
    else:
        nodes_df = pd.DataFrame({"id": list(unique_nodes)})
        nodes_df["name"] = nodes_df["id"]
        nodes_df["tooltip"] = "Unknown"

    # --- Edges-DFs für Export vs. Cytoscape ---
    # Für .cyjs-Export darf ein leerer DF existieren; für Cytoscape muss dann None übergeben werden
    if edges:
        edges_df_for_export = pd.DataFrame(edges).rename(columns={"chain_a": "source", "chain_b": "target"})
        edges_df_for_export["interaction"] = "interacts_with"
        edges_df_for_cyto = edges_df_for_export
    else:
        edges_df_for_export = pd.DataFrame(columns=["source", "target", "interaction", "all_atoms_count"])
        edges_df_for_cyto = None  # verhindert CyError bei Kanten=0

    # --- Pfad 1: kein Cytoscape – direkt .cyjs schreiben und fertig ---
    if not config.get("open_in_cytoscape", True):
        os.makedirs(run_output_path, exist_ok=True)
        network_file = os.path.join(run_output_path, f"{network_title}.cyjs")

        # Minimales Cytoscape-JSON (Nodes/Edges)
        nodes_elems = []
        for _, row in nodes_df.iterrows():
            nodes_elems.append({
                "data": {
                    "id": row["id"],
                    "name": row.get("name", row["id"]),
                    "tooltip": row.get("tooltip", ""),
                    "color_group": row.get("color_group", "")
                }
            })

        edges_elems = []
        for idx, row in edges_df_for_export.iterrows():
            edges_elems.append({
                "data": {
                    "id": f"edge_{idx}",
                    "source": row["source"],
                    "target": row["target"],
                    "interaction": row.get("interaction", "interacts_with"),
                    "all_atoms_count": row.get("all_atoms_count", 0)
                }
            })

        cyjs_data = {
            "data": {"name": network_title},
            "elements": {"nodes": nodes_elems, "edges": edges_elems}
        }
        with open(network_file, "w", encoding="utf-8") as f:
            json.dump(cyjs_data, f, indent=2)
        return  # kein Cytoscape nötig

    # --- Pfad 2: Cytoscape offen – Netzwerk bauen, Style anwenden, exportieren ---
    # ggf. alte Netze einkürzen
    existing_networks = p4c.get_network_list()
    while len(existing_networks) > config["keep_last_n_networks"]:
        oldest_network = existing_networks.pop(0)
        p4c.delete_network(oldest_network)

    # edges_df_for_cyto ist None, wenn keine Kanten vorhanden sind (ist ok)
    p4c.create_network_from_data_frames(nodes_df, edges_df_for_cyto, title=network_title)

    # --- Node-Tabellen in Cytoscape updaten (immer über id mappen) ---
    try:
        if "color_group" in nodes_df.columns:
            p4c.load_table_data(
                data=nodes_df[["id", "color_group"]],
                data_key_column="id",
                table="node",
                table_key_column="id"     # <— vorher "name"; jetzt robust über "id"
            )

        p4c.load_table_data(
            data=nodes_df[["id", "name"]],
            data_key_column="id",
            table="node",
            table_key_column="id"
        )

        if "tooltip" in nodes_df.columns:
            p4c.load_table_data(
                data=nodes_df[["id", "tooltip"]],
                data_key_column="id",
                table="node",
                table_key_column="id"     # <— vorher "name"; jetzt "id"
            )
    except Exception as e:
        print(f"Error while loading node data: {e}")

    # --- Style bestimmen & anwenden ---
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

        # Feste Farben (wie bisher)
        fixed_colors = {
            "Protein": "#1f77b4",
            "DNA": "#ff7f0e",
            "RNA": "#2ca02c",
            "DNA/RNA": "#d62728",
            "Nucleic Acid": "#9467bd"
        }

        # Automatische Palette für sonstige color_groups (z. B. PDB-IDs)
        base_color_groups = [g for g in color_groups if g not in fixed_colors and g != "Multi"]
        cmap = cm.get_cmap('tab20', len(base_color_groups))
        auto_colors = {group: to_hex(cmap(i)) for i, group in enumerate(base_color_groups)}
        if is_combined_protein and "Multi" in color_groups:
            auto_colors["Multi"] = "#8c564b"
        color_map = {**auto_colors, **fixed_colors}

        # Style einmalig anlegen (falls nicht vorhanden)
        if style_name not in p4c.get_visual_style_names():
            defaults = {
                "NODE_SHAPE": "ELLIPSE",
                "NODE_SIZE": 40,
                "NODE_LABEL_POSITION": "C,C,c,0.00,0.00",
                "EDGE_TRANSPARENCY": 120
            }
            mappings = [
                # Label und Tooltip als Passthrough
                {"mappingType": "passthrough", "mappingColumn": "name",    "mappingColumnType": "String", "visualProperty": "NODE_LABEL"},
                {"mappingType": "passthrough", "mappingColumn": "tooltip", "mappingColumnType": "String", "visualProperty": "NODE_TOOLTIP"},
                # Farbe diskret aus color_group
                {
                    "mappingType": "discrete",
                    "mappingColumn": "color_group",
                    "mappingColumnType": "String",
                    "visualProperty": "NODE_FILL_COLOR",
                    "map": [{"key": k, "value": v} for k, v in color_map.items()]
                }
            ]
            p4c.create_visual_style(style_name, mappings=mappings, defaults=defaults)

        # Style (ggf. neu) anwenden
        p4c.set_current_network(network_title)
        p4c.set_visual_style(style_name)
        p4c.layout_network(layout_name="force-directed")
    except Exception as e:
        print(f"Error while applying style: {e}")

    # --- Export in Datei (.cyjs) ---
    os.makedirs(run_output_path, exist_ok=True)
    network_file = os.path.join(run_output_path, f"{network_title}.cyjs")
    p4c.export_network(network_file, type="cyjs")


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
