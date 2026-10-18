"""Exercise an installed wheel outside the checkout using frozen scientific facts.

Run with the isolated environment's Python and -I. Only the standard library
is used by this checker; the installed console command loads runtime packages.
"""

import argparse
from collections import Counter
import csv
import hashlib
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
    store = work / "store"
    command("precompute", ["precompute", "--input-dir", str(inputs), "--store", str(store)])
    shutil.rmtree(inputs)
    shutil.rmtree(references)

    def hashes():
        return {str(p.relative_to(store)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in store.rglob("*") if p.is_file()}

    before = hashes()
    if os.name == "posix":
        for path in store.rglob("*"):
            path.chmod(0o555 if path.is_dir() else 0o444)
        store.chmod(0o555)
    try:
        command("assemble", ["assemble", "--store", str(store), "--pdb-id", "4hhb", "--pdb-id", "1bna",
                             "--output-dir", str(work / "assembled"), "--web-output-dir", str(work / "assembled-web")])
        verify_outputs(work / "assembled-web", facts, contacts)
        assert hashes() == before, "Assembly changed the store"
    finally:
        if os.name == "posix":
            store.chmod(0o755)
            for path in store.rglob("*"):
                path.chmod(0o755 if path.is_dir() else 0o644)
    print(json.dumps({"status": "passed", "package": str(package), "structures": 2,
                      "networks_per_run": 3, "read_only_store": True, "cytoscape_extra": False}))


if __name__ == "__main__":
    main()
