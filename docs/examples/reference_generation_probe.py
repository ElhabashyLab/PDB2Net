"""Reproduce D02 with tiny temporary references; do not change cache policy.

Run: .venv/bin/python docs/examples/reference_generation_probe.py
The probe changes only its own reference files/config in an isolated process.
No structures, searches, downloads or real reference generations are touched.
"""

import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from pdb2net import config_loader, reference_data, unknown_molecule_uniprot, uniprot_matcher

REPO = Path(__file__).resolve().parents[2]


def annotate():
    data = [{"pdb_id": "1ABC", "file_path": "1abc.pdb", "atom_data": [
        {"chain_id": "A", "unique_chain_id": "1ABC:A", "molecule_type": "Unknown",
         "molecule_name": "Unknown", "uniprot_id": None,
         "residues": [{"residue_name": "ALA", "is_polymer": True}]}]}]
    unknown_molecule_uniprot.process_molecule_info(data)
    c = data[0]["atom_data"][0]
    return {"uniprot_id": c["uniprot_id"], "name": c["molecule_name"]}


def main():
    cfg = json.loads((REPO/'pdb2net/configs/config.base.json').read_text())
    with tempfile.TemporaryDirectory(prefix="pdb2net-d02-probe-") as temporary:
        root = Path(temporary)
        first, second = root/'generation-one', root/'generation-two'
        first.mkdir()
        sifts, names, fasta = first/'sifts.tsv', first/'names.fasta', first/'seqres.fasta'
        sifts.write_text("1abc\tA\tP12345\n")
        names.write_text(">sp|P12345|FIRST Name old\nAAA\n>sp|Q99999|SECOND Name new\nAAA\n")
        fasta.write_text(">1abc_A mol:protein length:3 Header old\nAAA\n")
        cfg.update(sifts_tsv_path=str(sifts), uniprot_fasta_path=str(names), pdb_fasta_path=str(fasta),
                   reference_manifest_id="one", blast_db_path=str(first), blast_cache_path=str(root/'search.sqlite3'))
        config_loader.activate_config(cfg)
        initial = annotate()
        header_before = reference_data.load_pdb_fasta_headers(str(fasta), ("1abc",))["1abc_A"]["info"]
        signature_before = uniprot_matcher._db_signature()
        uniprot_matcher._CACHE_ENABLED = True
        uniprot_matcher._cache_put("probe-key", None)
        assert uniprot_matcher._cache_get("probe-key")[0]

        sifts.write_text("1abc\tA\tQ99999\n")
        names.write_text(">sp|P12345|FIRST Name replaced\nAAA\n>sp|Q99999|SECOND Name new\nAAA\n")
        fasta.write_text(">1abc_A mol:protein length:3 Header new\nAAA\n")
        cfg["reference_manifest_id"] = "two"
        config_loader.activate_config(cfg)
        same_process = annotate()
        assert same_process == initial == {"uniprot_id": "P12345", "name": "Name old"}
        assert reference_data.load_pdb_fasta_headers(str(fasta), ("1abc",))["1abc_A"]["info"] == header_before
        assert reference_data.load_uniprot_names(str(names))["P12345"] == "Name old"
        assert uniprot_matcher._db_signature() != signature_before
        assert not uniprot_matcher._cache_get("probe-key")[0]

        config_file = root/'config.json'
        config_file.write_text(json.dumps(cfg))
        child = subprocess.run(
            [sys.executable, "-c", "import json,sys; from pathlib import Path; "
             "sys.path.insert(0,sys.argv[1]); import reference_generation_probe as probe; "
             "probe.config_loader.activate_config(json.loads(Path(sys.argv[2]).read_text())); "
             "print(json.dumps(probe.annotate()))", str(Path(__file__).parent), str(config_file)],
            capture_output=True, text=True, check=True, timeout=10,
        )
        fresh = json.loads(child.stdout)
        assert fresh == {"uniprot_id": "Q99999", "name": "Name new"}
        shutil.copytree(first, second)
        cfg.update(sifts_tsv_path=str(second/'sifts.tsv'), uniprot_fasta_path=str(second/'names.fasta'),
                   pdb_fasta_path=str(second/'seqres.fasta'))
        config_loader.activate_config(cfg)
        changed_paths = annotate()
        assert changed_paths == fresh
        observations = {
            "initial": initial, "same_process_same_paths_after_replacement": same_process,
            "fresh_process_same_paths": fresh, "same_process_new_generation_paths": changed_paths,
            "old_seqres_header_cached": True, "search_cache_signature_changed": True,
            "old_search_cache_row_not_reused": True,
            "sifts_mapping_cached_paths": reference_data.load_sifts_mapping.cache_info().currsize,
        }
        # Release the probe's SQLite handle before deleting its temporary files.
        if uniprot_matcher._CACHE_CONN is not None:
            uniprot_matcher._CACHE_CONN.close()
            uniprot_matcher._CACHE_CONN = None
    print(json.dumps(observations, indent=2))


if __name__ == "__main__":
    main()
