"""Reproduce the approved Q02 network semantics through actual CX2 export.

Run from the repository root:
  .venv/bin/python docs/examples/q02_network_examples.py

Tiny synthetic, already annotated chains isolate aggregation from annotation.
Weights are explicit inputs, not newly measured physical contacts. Exports use
the real headless CX2 path. Assembly probes use real PDB/mmCIF parsing.
Only a temporary directory is written; no reference data or Cytoscape is needed.
"""

import json
from pathlib import Path
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import gemmi
import pandas as pd

from pdb2net import config_loader, distances, file_parser, pipeline, protein_network
from pdb2net.cytoscape import create_cytoscape_network
from pdb2net.cx2_export import export_cx2_headless
from pdb2net.outputs import _cx2_counts
from pdb2net.structure_identity import identity_from_official_id


def chain(pdb, label, accession=None, kind="Protein", offset=0.):
    return {
        "chain_id": label, "unique_chain_id": f"{pdb}:{label}", "model_index": 1,
        "molecule_type": kind, "molecule_name": kind, "uniprot_id": accession,
        "residues": [{"residue_name": "ALA" if kind in {"Protein", "Unknown"} else "DA",
                      "residue_number": i, "is_polymer": True,
                      "atoms": [{"atom_name": "CA", "coordinates": [30.*i, offset, 0.]}]}
                     for i in range(10)],
    }


def structure(pdb, *chains):
    return {"pdb_id": pdb, "file_path": f"{pdb}.pdb", "atom_data": list(chains),
            "structure_identity": identity_from_official_id(pdb).as_dict()}


def edge(a, b, count, kind="Protein-Protein"):
    return {"chain_a": a, "chain_b": b, "all_atoms_count": count, "interaction_type": kind}


def read_networks(root):
    result = {}
    for path in sorted(root.rglob("*.cx2")):
        _cx2_counts(path)  # Actual Core artifact validation, not only JSON decoding.
        aspects = {key: value for block in json.loads(path.read_text()) for key, value in block.items()}
        labels = {n["id"]: n["v"]["name"].split("\n")[0] for n in aspects["nodes"]}
        result[path.stem] = {
            "nodes": sorted(labels.values()),
            "edges": sorted([*sorted((labels[e["s"]], labels[e["t"]])),
                             e["v"]["all_atoms_count"], e["v"]["interaction"]]
                            for e in aspects["edges"]),
        }
    return result


def exports(root, structures, edges):
    pipeline._export_chain_networks(structures, edges, str(root / "chain"))
    pipeline._create_linked_identity_network(edges, structures, str(root / "combined"))
    protein_network.create_protein_network(
        edges, structures, str(root / "protein"),
        {"protein_per_pdb": True, "combined_protein_network": True},
    )
    return read_networks(root)


def assembly_pdb():
    lines = ["HEADER".ljust(62) + "1ABC\n", "REMARK 350 BIOMOLECULE: 1\n",
             "REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: DIMERIC\n",
             "REMARK 350 APPLY THE FOLLOWING TO CHAINS: A\n"]
    for operation, y in [(1, 0.), (2, 3.)]:
        for row in range(3):
            matrix = [float(i == row) for i in range(3)]
            lines.append(f"REMARK 350   BIOMT{row+1} {operation:3d}  "
                         + " ".join(f"{v:9.6f}" for v in matrix)
                         + f"     {y if row == 1 else 0.:10.5f}\n")
    for i in range(10):
        lines.append(f"ATOM  {i+1:5d}  CA  ALA A{i+1:4d}    "
                     f"{i*30.:8.3f}{0.:8.3f}{0.:8.3f}  1.00 20.00           C  \n")
    return "".join(lines) + "TER\nEND\n"


def main():
    cfg = json.loads((Path(__file__).resolve().parents[2] / "pdb2net/configs/config.base.json").read_text())
    cfg["open_in_cytoscape"] = False
    config_loader.activate_config(cfg)
    observations = {}
    with tempfile.TemporaryDirectory(prefix="pdb2net-q02-examples-") as temporary:
        root = Path(temporary)
        for case, second in [("heteromer", "Q99999"), ("homomer", "P12345"), ("unmapped_protein", None)]:
            data = [structure("1ABC", chain("1ABC", "A", "P12345"), chain("1ABC", "B", second))]
            out = exports(root / case, data, [edge("1ABC:A", "1ABC:B", 12)])
            assert out["Chain_Interaction_Network_1ABC"]["edges"] == [["1ABC:A", "1ABC:B", 12, "Protein-Protein"]]
            expected = [["P12345", "Q99999", 12, "interacts_with"]] if second == "Q99999" else []
            assert out["Protein_Network_1ABC"]["edges"] == expected
            assert out["Protein_Network_1ABC"]["nodes"] == (["P12345", "Q99999"] if second == "Q99999" else ["P12345"])
            observations[case] = out

        data = [structure("1ABC", chain("1ABC", "A", "P12345"), chain("1ABC", "B", "P12345"),
                          chain("1ABC", "C", "Q99999"))]
        rows = [edge("1ABC:A", "1ABC:C", 12), edge("1ABC:B", "1ABC:C", 7), edge("1ABC:A", "1ABC:B", 5)]
        out = exports(root / "summation", data, rows)
        assert out["Protein_Network_1ABC"]["edges"] == [["P12345", "Q99999", 19, "interacts_with"]]
        observations["summation"] = out

        data = [structure("1ABC", chain("1ABC", "A", "P12345"), chain("1ABC", "B", "P12345"),
                          chain("1ABC", "D", kind="DNA"), chain("1ABC", "E", kind="RNA")),
                structure("2XYZ", chain("2XYZ", "A", "P12345"), chain("2XYZ", "B", "Q99999"))]
        rows = [edge("1ABC:A", "1ABC:D", 3, "Protein-DNA"), edge("1ABC:B", "1ABC:D", 4, "Protein-DNA"),
                edge("1ABC:D", "1ABC:E", 6, "DNA-RNA"), edge("2XYZ:A", "2XYZ:B", 8)]
        out = exports(root / "combined", data, rows)
        assert out["Protein_Network_1ABC"]["edges"] == [["1ABC:D", "P12345", 7, "interacts_with"]]
        assert out["Combined_Protein_Network_P12345_Q99999"]["edges"] == [["P12345", "Q99999", 8, "interacts_with"]]
        linked, = [value for name, value in out.items() if name.startswith("Combined_Network_")]
        identity = [row for row in linked["edges"] if row[3] == "identity"]
        assert len(identity) == 2 and all(row[2] == 0 for row in identity)
        assert len([row for row in linked["edges"] if row[3] != "identity"]) == 4
        observations["combined_and_nucleic_acids"] = out

        # Identity zero is valid at both the serializer and the public exporter.
        export_cx2_headless(
            "Combined_Network_identity_zero", str(root / "zero-candidate"),
            pd.DataFrame([{"id": n, "name": n, "tooltip": n, "color_group": "Protein"} for n in ("A", "B")]),
            pd.DataFrame([{"source": "A", "target": "B", "interaction": "identity", "all_atoms_count": 0}]),
            {}, {"A": {"x": 0., "y": 0.}, "B": {"x": 1., "y": 0.}},
        )
        candidate = read_networks(root / "zero-candidate")["Combined_Network_identity_zero"]
        assert candidate["edges"] == [["A", "B", 0, "identity"]]
        create_cytoscape_network([edge("A", "B", 0, "identity")], "identity_zero_current",
                                str(root / "zero-current"), nodes_data=[{"id": n, "name": n} for n in ("A", "B")])
        assert read_networks(root / "zero-current")["identity_zero_current"]["edges"] == [["A", "B", 0, "identity"]]
        observations["identity_zero_export"] = {"low_level_export_valid": True, "high_level_retains_edge": True}

        unknown = structure("1ABC", chain("1ABC", "A", "P12345"), chain("1ABC", "B", kind="Unknown", offset=3))
        assert distances.calculate_distances_with_ckdtree([unknown]) == []
        unknown["atom_data"][1]["molecule_type"] = "Protein"
        assert distances.calculate_distances_with_ckdtree([unknown])[0]["all_atoms_count"] == 10
        observations["unknown_type_vs_missing_accession"] = {"unknown_type_contacts": 0, "protein_without_accession_contacts": 10}

        pdb = assembly_pdb()
        # mmCIF is a Gemmi conversion of this probe, not an independent format oracle.
        document = gemmi.read_pdb_string(pdb).make_mmcif_document()
        document[0].name = "1abc"
        cif = document.as_string()
        for extension, text in [("pdb", pdb), ("cif", cif)]:
            source = root / f"1abc.{extension}"
            source.write_text(text)
            parsed = file_parser.parse_structure_input(str(source))
            st = parsed["structure"]
            assert len(st.assemblies) == 1
            assert len(st.assemblies[0].generators[0].operators) == 2
            result = pipeline.process_single_file(str(source))
            assert len(result["atom_data"]) == 1
            result["atom_data"][0]["molecule_type"] = "Protein"
            assert distances.calculate_distances_with_ckdtree([result]) == []
            out = exports(root / f"assembly-{extension}", [result], [])
            assert out["Chain_Interaction_Network_1ABC"] == {"nodes": ["1ABC:A"], "edges": []}
            observations[f"assembly_{extension}"] = {"recognized_assemblies": 1, "operators": 2, "networks": out}

        try:
            create_cytoscape_network([edge("P12345", "P12345", 12)], "self_loop_probe", str(root / "selfloop"),
                                    nodes_data=[{"id": "P12345", "name": "P12345", "color_group": "Protein"}])
        except ValueError as exc:
            assert "self-loop" in str(exc)
            observations["self_loop_export"] = str(exc)
        else:
            raise AssertionError("Self-loop unexpectedly accepted")

        # Artifact validation independently rejects a self-loop in otherwise valid CX2.
        source = root / "heteromer/chain/Chain_Interaction_Network_1ABC.cx2"
        document = json.loads(source.read_text())
        for block in document:
            for row in block.get("edges", []):
                row["t"] = row["s"]
        altered = root / "self_loop_artifact.cx2"
        altered.write_text(json.dumps(document))
        try:
            _cx2_counts(altered)
        except ValueError as exc:
            observations["self_loop_validation"] = str(exc)
        else:
            raise AssertionError("Self-loop artifact unexpectedly accepted")
    print(json.dumps(observations, indent=2))


if __name__ == "__main__":
    main()
