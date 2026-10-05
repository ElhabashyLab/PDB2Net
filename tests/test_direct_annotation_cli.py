"""Exercise filename versus content annotation through the real headless CLI."""

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


@pytest.mark.parametrize(
    ("kind", "residue", "atom", "element", "sifts", "expected_id"),
    [
        ("Protein", "ALA", "CA", "C", "1abc\tA\tP12345\n1abc\tA\tQ99999\n", None),
        ("RNA", "U", "P", "P", "# PDB\tCHAIN\tSP_PRIMARY\n", None),
        ("Protein", "ALA", "CA", "C", "1abc\tA\tQ99999\n", "Q99999"),
    ],
)
def test_af_filename_cannot_override_real_content_annotation(
    tmp_path, kind, residue, atom, element, sifts, expected_id
):
    root = Path(__file__).resolve().parents[1]
    (tmp_path / "sifts.tsv").write_text(sifts)
    (tmp_path / "uniprot.fasta").write_text(
        ">sp|P12345|A Protein A\nAG\n>sp|Q99999|B Protein B\nAG\n"
    )
    molecule = "mol:na length:1 RNA" if kind == "RNA" else "mol:protein length:1 Protein"
    (tmp_path / "pdb.fasta").write_text(f">1abc_A {molecule}\nA\n")
    config = {
        "pdb_fasta_path": str(tmp_path / "pdb.fasta"),
        "uniprot_fasta_path": str(tmp_path / "uniprot.fasta"),
        "sifts_tsv_path": str(tmp_path / "sifts.tsv"),
        "blast_db_path": str(tmp_path / "no-blast-database"),
        "blast_cache_path": str(tmp_path / "cache.sqlite3"),
        "reference_manifest_id": "direct-annotation-cli-fixture",
        "open_in_cytoscape": False,
        "workers": {"parsing": 1, "blast_threads": 1},
        "diamond": {"enabled": False},
        "networks": {
            "chain_per_pdb": True, "protein_per_pdb": False,
            "combined_chain_network": False, "combined_protein_network": False,
        },
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))
    pdb = "HEADER".ljust(62) + "1ABC\n" + (
        f"ATOM  {1:5d} {atom:^4s} {residue:3s} A{1:4d}    "
        f"{0.:8.3f}{0.:8.3f}{0.:8.3f}  1.00 20.00          {element:>2s}  \nEND\n"
    )
    env = {key: value for key, value in os.environ.items() if not key.startswith("PDB2NET_")}
    env["MPLCONFIGDIR"] = str(tmp_path / "matplotlib")
    defaults = tmp_path / "defaults"
    defaults.mkdir()
    shutil.copyfile(root / "pdb2net/configs/config.base.json", defaults / "config.base.json")
    env["PDB2NET_CONFIG_DIR"] = str(defaults)
    for name in ["1abc.pdb", "AF-P69905-F1-model_v6.pdb"]:
        case = tmp_path / name.removesuffix(".pdb")
        inputs = case / "inputs"
        inputs.mkdir(parents=True)
        (inputs / name).write_text(pdb)
        result = subprocess.run(
            [sys.executable, "-m", "pdb2net", "run", "--input-dir", str(inputs),
             "--output-dir", str(case / "work"), "--web-output-dir", str(case / "web"),
             "--config", str(config_path), "--headless"],
            cwd=root, env=env, capture_output=True, text=True, timeout=30,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        summary = json.loads((case / "web/summary.json").read_text())
        assert summary["status"] == "success" and summary["errors"] == []
        assert summary["identities"][0]["canonical_id"] == "pdb_00001abc"
        (network_path,) = summary["networks"]
        document = json.loads((case / "web" / network_path).read_text())
        (node,) = [node for block in document for node in block.get("nodes", [])]
        assert node["v"]["molecule_type"] == kind
        assert node["v"].get("uniprot_id") == expected_id
