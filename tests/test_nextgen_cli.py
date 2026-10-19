"""Run the checked-in real NextGen excerpt through CLI and coordinate-free analysis."""

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/nextgen/pdb_00008aly_xyz-enrich.cif"


def test_checked_in_nextgen_cli_preserves_independent_annotation_expectations(tmp_path):
    assert hashlib.sha256(FIXTURE.read_text(encoding="utf-8").encode("utf-8")).hexdigest() == (
        "c53dbbbc5629cc0899afcc16ae5770e0632e400bad625ebf63be930d51030661"
    )
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    shutil.copyfile(FIXTURE, inputs / FIXTURE.name)
    # Synthetic lookup files for this excerpt, not a full reference generation.
    (tmp_path / "pdb.fasta").write_text(
        ">8aly_E mol:protein length:2 Mini-fixture tankyrase-2\nDF\n"
    )
    (tmp_path / "sifts.tsv").write_text("# PDB\tCHAIN\tSP_PRIMARY\n")
    (tmp_path / "uniprot.fasta").write_text(
        ">sp|Q9H2K2|TNKS2 Mini-fixture tankyrase-2\nDF\n"
    )
    config = {
        "pdb_fasta_path": str(tmp_path / "pdb.fasta"),
        "sifts_tsv_path": str(tmp_path / "sifts.tsv"),
        "uniprot_fasta_path": str(tmp_path / "uniprot.fasta"),
        "blast_db_path": str(tmp_path / "no-blast-database"),
        "blast_cache_path": str(tmp_path / "cache.sqlite3"),
        "reference_manifest_id": "checked-in-nextgen-excerpt-v1",
        "open_in_cytoscape": False,
        "workers": {"parsing": 1, "blast_threads": 1},
        "diamond": {"enabled": False},
        "networks": {
            "chain_per_pdb": True, "protein_per_pdb": False,
            "combined_chain_network": False, "combined_protein_network": False,
        },
        "network_annotations": {"tooltip_fields": ["uniprot", "pfam"]},
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config))
    env = {key: value for key, value in os.environ.items() if not key.startswith("PDB2NET_")}
    env["MPLCONFIGDIR"] = str(tmp_path / "matplotlib")
    defaults = tmp_path / "defaults"
    defaults.mkdir()
    shutil.copyfile(ROOT / "pdb2net/configs/config.base.json", defaults / "config.base.json")
    env["PDB2NET_CONFIG_DIR"] = str(defaults)

    def cli(command):
        result = subprocess.run(
            [sys.executable, "-m", "pdb2net", *command, "--config", str(config_path), "--headless"],
            cwd=ROOT, env=env, capture_output=True, text=True, timeout=30,
        )
        assert result.returncode == 0, result.stdout + result.stderr

    cli(["run", "--input-dir", str(inputs), "--output-dir", str(tmp_path / "raw"),
         "--web-output-dir", str(tmp_path / "raw-web")])
    script = r"""
import json,sys
from pathlib import Path
from pdb2net.config_loader import load_config,activate_config
from pdb2net.structure_analysis import prepare_structure,analyze_prepared
root,source=map(Path,sys.argv[1:])
cfg=load_config(explicit_file=root/'config.json');cfg['output_path']=str(root/'assembled');activate_config(cfg)
facts=prepare_structure(source);source.unlink()
analyze_prepared([facts],cfg,web_output_dir=str(root/'assembled-web'))
"""
    result=subprocess.run([sys.executable,"-c",script,str(tmp_path),str(inputs/FIXTURE.name)],
                          cwd=ROOT,env=env,capture_output=True,text=True,timeout=30)
    assert result.returncode==0,result.stdout+result.stderr

    # Both paths are checked against the original fixture facts, not each other.
    for directory in [tmp_path / "raw-web", tmp_path / "assembled-web"]:
        summary = json.loads((directory / "summary.json").read_text())
        assert summary["status"] == "success" and summary["errors"] == []
        assert summary["counts"]["structures"] == 1
        assert summary["counts"]["chains"] == 1
        assert summary["annotations"]["by_source"] == {"embedded_sifts": 1}
        (name,) = summary["networks"]
        aspects = {key: value for block in json.loads((directory / name).read_text()) for key, value in block.items()}
        assert aspects["edges"] == []
        (node,) = aspects["nodes"]
        values = node["v"]
        assert (values["name"], values["chain_id"], values["uniprot_id"], values["aa_len"]) == (
            "8ALY:E", "E", "Q9H2K2", 2,
        )
        assert values["molecule_type"] == "Protein"
        assert [
            (s["accession"], s["chain_id"], s["label_asym_id"], s["pdb_start"],
             s["pdb_end"], s["database_start"], s["database_end"])
            for s in json.loads(values["annotation_uniprot"])
        ] == [("Q9H2K2", "E", "A", 4, 299, 867, 1162)]
        assert {
            (s["accession"], s["pdb_start"], s["pdb_end"])
            for s in json.loads(values["annotation_pfam"])
        } == {("PF00644", 99, 270), ("PF07647", 8, 71)}
