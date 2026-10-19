from pathlib import Path
import json
import os
import subprocess
import sys

from pdb2net.reference_data import (
    load_pdb_fasta_headers,
    load_sifts_mapping,
    load_sifts_segments,
)


def test_pdb_header_lookup_does_not_retain_reference_sequences(tmp_path: Path) -> None:
    fasta = tmp_path / "pdb_seqres.txt"
    fasta.write_text(
        ">1abc_A mol:protein length:6 Example protein\n"
        "AAAAAA\n"
        ">2xyz_B mol:na length:4 Example DNA\n"
        "ATGC\n"
        ">3def_C mol:protein length:3 Unrequested protein\n"
        "GGG\n",
        encoding="utf-8",
    )
    load_pdb_fasta_headers.cache_clear()

    headers = load_pdb_fasta_headers(str(fasta), ("1ABC", "2xyz"))

    assert headers == {
        "1abc_A": {
            "info": "mol:protein length:6 Example protein",
        },
        "2xyz_B": {
            "info": "mol:na length:4 Example DNA",
        },
    }


def test_real_1ao7_chimera_segments_are_retained_and_compatibility_mapping_is_ambiguous(
    tmp_path: Path,
) -> None:
    """Chain D reflects the two accessions in the official PDBe-SIFTS record."""
    mapping = tmp_path / "uniprot_segments_observed.tsv"
    mapping.write_text(
        "# PDB\tCHAIN\tSP_PRIMARY\tRES_BEG\tRES_END\tPDB_BEG\tPDB_END\tSP_BEG\tSP_END\n"
        "1ao7\tD\tA0A075B6T6\t1\t92\t1\t92\t23\t112\n"
        "1ao7\tD\tP01848\t118\t210\t118\t210\t1\t93\n",
        encoding="utf-8",
    )
    load_sifts_segments.cache_clear()
    load_sifts_mapping.cache_clear()

    segments = load_sifts_segments(str(mapping))
    compatibility = load_sifts_mapping(str(mapping))

    assert [segment["accession"] for segment in segments["1ao7_D"]] == [
        "A0A075B6T6",
        "P01848",
    ]
    assert "1ao7_D" not in compatibility


def test_repeated_external_sifts_segments_for_one_accession_remain_unique(
    tmp_path: Path,
) -> None:
    mapping = tmp_path / "uniprot_segments_observed.tsv"
    mapping.write_text(
        "# PDB\tCHAIN\tSP_PRIMARY\tRES_BEG\tRES_END\tPDB_BEG\tPDB_END\tSP_BEG\tSP_END\n"
        "8aly\tE\tQ9H2K2\t4\t71\t4\t71\t867\t934\n"
        "8aly\tE\tQ9H2K2\t99\t270\t99\t270\t962\t1133\n",
        encoding="utf-8",
    )
    load_sifts_segments.cache_clear()
    load_sifts_mapping.cache_clear()

    assert len(load_sifts_segments(str(mapping))["8aly_E"]) == 2
    assert load_sifts_mapping(str(mapping))["8aly_E"] == "Q9H2K2"


def test_new_process_uses_its_reference_generation_and_search_cache_namespace(tmp_path):
    """The supported lifecycle, with real annotation and no sequence searches."""
    repo = Path(__file__).resolve().parents[1]
    configs = []
    for number, accession, name in [(1, "P12345", "First name"), (2, "Q99999", "Second name")]:
        generation = tmp_path / f"generation-{number}"
        generation.mkdir()
        (generation / "sifts.tsv").write_text(f"1abc\tA\t{accession}\n")
        (generation / "uniprot.fasta").write_text(f">sp|{accession}|EXAMPLE {name}\nAAA\n")
        (generation / "seqres.fasta").write_text(">1abc_A mol:protein length:3 Example\nAAA\n")
        cfg = json.loads((repo / "pdb2net/configs/config.base.json").read_text())
        cfg.update(sifts_tsv_path=str(generation / "sifts.tsv"),
                   pdb_fasta_path=str(generation / "seqres.fasta"),
                   uniprot_fasta_path=str(generation / "uniprot.fasta"),
                   blast_db_path=str(generation), blast_cache_path=str(tmp_path / "search.sqlite3"),
                   reference_manifest_id=f"generation-{number}")
        config_file = generation / "config.json"
        config_file.write_text(json.dumps(cfg))
        configs.append(config_file)
    code = '''
import json, sys
from pathlib import Path
from pdb2net import config_loader, uniprot_matcher as search
from pdb2net.unknown_molecule_uniprot import process_molecule_info
from pdb2net.structure_facts import profile_id
config_loader.activate_config(json.loads(Path(sys.argv[1]).read_text()))
search._CACHE_ENABLED = True
annotations = []
for _ in range(2):
    data = [{"pdb_id":"1ABC", "file_path":"1abc.pdb", "atom_data":[
        {"chain_id":"A", "unique_chain_id":"1ABC:A", "uniprot_id":None,
         "residues":[{"residue_name":"ALA", "is_polymer":True}]}]}]
    process_molecule_info(data)
    chain = data[0]["atom_data"][0]
    annotations.append([chain["uniprot_id"], chain["molecule_name"]])
cached = search._cache_get("generation-probe")[0]
search._cache_put("generation-probe", None)
print(json.dumps({"annotations":annotations, "cached":cached, "profile":profile_id()}))
'''
    env = {k: v for k, v in os.environ.items() if not k.startswith("PDB2NET_")}
    observed = []
    for config_file in [configs[0], configs[0], configs[1]]:
        result = subprocess.run([sys.executable, "-c", code, str(config_file)], cwd=repo,
                                env=env, capture_output=True, text=True, check=True, timeout=15)
        observed.append(json.loads(result.stdout))
    assert observed[0]["annotations"] == observed[1]["annotations"] == [["P12345", "First name"]] * 2
    assert observed[2]["annotations"] == [["Q99999", "Second name"]] * 2
    assert [result["cached"] for result in observed] == [False, True, False]
    assert observed[0]["profile"] == observed[1]["profile"] == observed[2]["profile"]
