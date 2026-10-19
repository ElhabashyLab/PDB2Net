"""Exercise an installed wheel outside the checkout using frozen scientific facts.

Run with the isolated environment's Python and -I. Only the standard library
bootstraps this checker; both the installed console command and Python API
exercise runtime packages.
"""

import argparse
from collections import Counter
import csv
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


def graph(path):
    aspects = {key: value for block in json.loads(path.read_text()) for key, value in block.items()}
    nodes = {node["id"]: node["v"] for node in aspects["nodes"]}
    assert len(nodes) == len(aspects["nodes"])
    edges = {}
    for edge in aspects["edges"]:
        a, b = nodes[edge["s"]]["name"], nodes[edge["t"]]["name"]
        pair = tuple(sorted((a, b)))
        assert pair not in edges
        edges[pair] = (edge["v"]["interaction"], edge["v"]["all_atoms_count"])
    return nodes, edges


def verify_outputs(web, facts, contacts):
    summary = json.loads((web / "summary.json").read_text())
    assert summary["status"] == "success" and summary["errors"] == []
    assert summary["counts"]["structures"] == 2 and summary["counts"]["chains"] == 6
    expected_files = set()
    for filename in ("4hhb.pdb", "1bna.cif"):
        structure = facts[filename]
        pdb_id = structure["display_id"]
        name = f"Chain_Interaction_Network_{pdb_id}.cx2"
        expected_files.add(name)
        nodes, edges = graph(web / "networks" / name)
        assert {n["chain_id"] for n in nodes.values()} == set(structure["chains"])
        for node in nodes.values():
            expected = structure["chains"][node["chain_id"]]
            assert node["molecule_type"] == expected["molecule_type"]
            assert (node.get("uniprot_id") or None) == expected["uniprot_id"]
            length = "aa_len" if expected["molecule_type"] == "Protein" else "nt_len"
            assert node[length] == expected["polymer_residues"]
        assert edges == {
            tuple(sorted(f"{pdb_id}:{c}" for c in row["chains"])): (row["interaction"], row["all_atoms_count"])
            for row in contacts[filename] if row["accepted"]
        }
        proteins = {c["uniprot_id"] for c in structure["chains"].values() if c["uniprot_id"]}
        if proteins:
            name = f"Protein_Network_{pdb_id}.cx2"
            expected_files.add(name)
            nodes, edges = graph(web / "networks" / name)
            assert {n["uniprot_id"] for n in nodes.values()} == proteins
            totals = Counter()
            for row in contacts[filename]:
                a, b = (structure["chains"][c]["uniprot_id"] for c in row["chains"])
                if row["accepted"] and a and b and a != b:
                    totals[tuple(sorted((a, b)))] += row["all_atoms_count"]
            assert edges == {pair: ("interacts_with", count) for pair, count in totals.items()}
    assert {p.name for p in (web / "networks").glob("*.cx2")} == expected_files
    assert {Path(p).name for p in summary["networks"]} == expected_files


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixtures", required=True, type=Path)
    parser.add_argument("--work-dir", required=True, type=Path, help="A new directory outside the checkout")
    args = parser.parse_args()
    fixtures = args.fixtures.resolve()
    checkout = Path(__file__).resolve().parents[1]
    work = args.work_dir.resolve()
    assert not Path.cwd().resolve().is_relative_to(checkout), "Run outside the checkout"
    assert not work.is_relative_to(checkout), "Work must remain outside the checkout"
    assert sys.prefix != sys.base_prefix, "Use a separate virtual environment"
    import pdb2net
    package = Path(pdb2net.__file__).resolve().parent
    assert package.is_relative_to(Path(sys.prefix).resolve()) and not package.is_relative_to(checkout)
    assert importlib.util.find_spec("pdb2net.precomputed") is None
    assert importlib.util.find_spec("py4cytoscape") is None, "Check the base package without the Cytoscape extra"
    for name in ("base", "linux", "darwin", "windows"):
        assert (package / f"configs/config.{name}.json").is_file()
    assert not (package / "configs/config.local.json").exists()
    work.mkdir(parents=True, exist_ok=False)
    inputs, references = work / "inputs", work / "references"
    inputs.mkdir()
    shutil.copytree(fixtures / "references", references)
    for filename in ("4hhb.pdb", "1bna.cif"):
        shutil.copyfile(fixtures / "structures" / filename, inputs / filename)
    facts = json.loads((fixtures / "expected/structures.json").read_text())
    contacts = json.loads((fixtures / "expected/networks.json").read_text())
    cfg = json.loads((fixtures / "config.json").read_text())
    cfg.update(
        pdb_fasta_path=str(references / "pdb_seqres.txt"), sifts_tsv_path=str(references / "sifts.tsv"),
        uniprot_fasta_path=str(references / "swissprot.fasta"), reference_manifest_id="installed-wheel-frozen-fixtures",
        blastp_executable=str(work / "absent-blastp"), blast_db_path=str(work / "absent-db"),
        blast_cache_path=str(work / "cache.sqlite3"), open_in_cytoscape=False,
        workers={"parsing": 1, "blast_threads": 1},
        networks={"chain_per_pdb": True, "protein_per_pdb": True,
                  "combined_chain_network": False, "combined_protein_network": False},
    )
    config = work / "config.json"
    config.write_text(json.dumps(cfg))
    # Remove machine overrides and Python import-path injection for every CLI.
    env = {k: v for k, v in os.environ.items() if not k.startswith(("PDB2NET_", "PYTHON"))}
    env.update(PYTHONNOUSERSITE="1", MPLCONFIGDIR=str(work / "matplotlib"))
    cli = Path(sys.executable).parent / ("pdb2net.exe" if os.name == "nt" else "pdb2net")
    assert cli.is_file()

    def command(name, arguments):
        result = subprocess.run([str(cli), *arguments, "--config", str(config), "--headless"],
                                cwd=work, env=env, capture_output=True, text=True, timeout=90)
        (work / f"{name}.log").write_text(result.stdout + result.stderr)
        assert result.returncode == 0, result.stdout + result.stderr

    command("raw", ["run", "--input-dir", str(inputs), "--output-dir", str(work / "raw"),
                    "--web-output-dir", str(work / "raw-web"), "--export-detailed-interactions"])
    verify_outputs(work / "raw-web", facts, contacts)
    for filename in ("4hhb.pdb", "1bna.cif"):
        path = work / "raw-web/interactions" / f"{facts[filename]['display_id']}_detailed_interactions.csv"
        with path.open() as stream:
            count = sum(1 for _ in csv.DictReader(stream))
        assert count == sum(row["all_atoms_count"] for row in contacts[filename] if row["accepted"])
    from pdb2net.config_loader import activate_config
    from pdb2net.structure_analysis import prepare_structure, analyze_prepared
    from pdb2net.structure_facts import canonical_json
    cfg.update(output_path=str(work / "prepared"), export_detailed_interactions=True)
    activate_config(cfg)
    prepared, streams = [], {}
    for source in sorted(inputs.iterdir()):
        rows = []
        item = prepare_structure(source, contact_sink=rows.append)
        prepared.append(item)
        streams[item["structure_identity"]["canonical_id"]] = lambda rows=rows: iter(rows)
    before = canonical_json(prepared)
    shutil.rmtree(inputs)
    analyze_prepared(prepared, cfg, contact_sources=streams, web_output_dir=str(work / "prepared-web"))
    verify_outputs(work / "prepared-web", facts, contacts)
    assert canonical_json(prepared) == before
    for source in (work / "raw-web/interactions").glob("*.csv"):
        assert source.read_bytes() == (work / "prepared-web/interactions" / source.name).read_bytes()
    print(json.dumps({"status": "passed", "package": str(package), "structures": 2,
                      "networks_per_run": 3, "immutable_facts": True, "fresh_annotation": True, "cytoscape_extra": False}))


if __name__ == "__main__":
    main()
