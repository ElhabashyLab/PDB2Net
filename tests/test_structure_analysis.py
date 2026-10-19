"""Prepared inputs must reproduce independent scientific fixtures without sources."""
from __future__ import annotations

import copy
import gzip
import json
from pathlib import Path

import pytest

from pdb2net import pipeline
from pdb2net.config_loader import activate_config, config
from pdb2net.input_contract import InputValidationError
from pdb2net.structure_analysis import analyze_prepared, prepare_structure
from pdb2net.structure_facts import canonical_json, geometry_profile, validate_contacts, validate_facts
from scripts.compare_cx2_outputs import compare_directories

DATA = Path(__file__).parent / "fixtures/core_reference"


@pytest.mark.parametrize("special,expected", [("SEC", "AAAU"), ("PYL", "AAAO")])
def test_prepared_sequences_preserve_rare_amino_acids_without_sources(tmp_path, monkeypatch, special, expected):
    from pdb2net.structure_facts import materialize_facts
    from pdb2net.uniprot_matcher import extract_sequence_from_parsed_data
    monkeypatch.setitem(config, "structure_model_policy", "first")
    source = tmp_path / "1abc.pdb"
    source.write_text("HEADER".ljust(62) + "1ABC\n" + "".join(
        f"ATOM  {i:5d}  CA  {name} A{i:4d}    {i * 3.:8.3f}{0.:8.3f}{0.:8.3f}  1.00 20.00           C  \n"
        for i, name in enumerate(("ALA", "ALA", "ALA", special), 1)) + "TER\nEND\n")
    saved = canonical_json(prepare_structure(source))
    source.unlink()
    restored = materialize_facts(json.loads(saved))["atom_data"][0]
    assert restored["structure_sequence"] == expected
    assert extract_sequence_from_parsed_data(restored) == expected
    assert restored["aa_len"] == 4


def test_prepared_exact_boundary_contacts_survive_serialization(tmp_path, monkeypatch):
    from pdb2net.structure_facts import profile_id
    monkeypatch.setitem(config, "structure_model_policy", "first")
    source = tmp_path / "1abc.pdb"
    source.write_text("HEADER".ljust(62) + "1ABC\n" + "".join(
        f"ATOM  {i:5d}  CA  ALA {chain}{1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  \nTER\n"
        for i, (chain, x, y, z) in enumerate([
            ("A", -4.533, 66.745, -21.694), ("B", -1.533, 70.745, -21.694)], 1)) + "END\n")
    rows = []
    facts = json.loads(canonical_json(prepare_structure(source, contact_sink=rows.append)))
    source.unlink()
    assert facts["pairs"] == [{"chain_a": "1ABC:A", "chain_b": "1ABC:B",
                               "all_atoms_count": 1, "ca_neighbors": 1}]
    assert len(list(validate_contacts(facts, json.loads(json.dumps(rows))))) == 1
    assert rows[0][-1] == 5.
    # Even internally consistent old profiles must fail rather than silently
    # reusing facts with lost sequences or the old numerical contact policy.
    for field, old in [("contacts", "pdb2net-distance-classification-v2"),
                       ("sequence_type_policy", "polymer-residues-v1")]:
        stale = copy.deepcopy(facts)
        stale["profile"][field] = old
        stale["profile_id"] = profile_id(stale["profile"])
        with pytest.raises(InputValidationError, match="profile is incompatible"):
            validate_facts(stale)


@pytest.fixture
def settings(tmp_path, integration_config):
    previous = copy.deepcopy(dict(config))
    cfg = copy.deepcopy(integration_config)
    cfg.update(output_path=str(tmp_path / "results"), export_detailed_interactions=True)
    cfg["workers"]["parsing"] = 1
    activate_config(cfg)
    yield cfg
    activate_config(previous)


@pytest.mark.integration
@pytest.mark.parametrize("names", [("1bna.cif",), ("4hhb.pdb", "1bna.cif"), ("6m17.cif",), ("1a8o.cif",)])
def test_prepared_real_outputs_match_raw_without_sources(tmp_path, monkeypatch, settings, names):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    payloads, sources = [], {}
    for name in names:
        path = inputs / name
        path.write_bytes((DATA / "structures" / name).read_bytes())
        rows = []
        facts = prepare_structure(path, contact_sink=rows.append)
        payloads.append(facts)
        sources[facts["structure_identity"]["canonical_id"]] = lambda rows=rows: iter(rows)
        assert list(validate_contacts(facts, rows)) == rows
    before = canonical_json(payloads)
    raw = pipeline.run_pipeline(str(inputs))
    for path in inputs.iterdir():
        path.unlink()
    monkeypatch.setattr(pipeline, "process_single_file", lambda *a, **k: pytest.fail("Prepared input reparsed a structure"))
    prepared = analyze_prepared(payloads, settings, contact_sources=sources)
    comparison = compare_directories(Path(prepared.run_output_path), Path(raw.run_output_path), 1.)
    assert comparison["status"] == "PASS", comparison
    assert comparison["files"] and not comparison["missing_expected_files"] and not comparison["extra_actual_files"]
    assert {p.name: p.read_bytes() for p in Path(raw.distances_dir).glob("*.csv")} == {
        p.name: p.read_bytes() for p in Path(prepared.distances_dir).glob("*.csv")}
    assert canonical_json(payloads) == before
    fixed = json.loads((DATA / "expected/networks.json").read_text())
    for facts, name in zip(payloads, names):
        measured = {tuple(sorted([p["chain_a"].split(":")[-1], p["chain_b"].split(":")[-1]])): p["all_atoms_count"] for p in facts["pairs"]}
        for pair in fixed[name]:
            if pair["accepted"]:
                assert measured[tuple(pair["chains"])] == pair["all_atoms_count"]


@pytest.mark.integration
def test_compression_and_references_do_not_define_geometry(tmp_path, settings):
    source = DATA / "structures/1bna.cif"
    facts = prepare_structure(source)
    zipped = tmp_path / "1bna.cif.gz"
    zipped.write_bytes(gzip.compress(source.read_bytes()))
    config["reference_manifest_id"] = "another-generation"
    config["sifts_tsv_path"] = str(tmp_path / "absent")
    config["uniprot_fasta_path"] = str(tmp_path / "absent-fasta")
    config["blast"] = {"max_target_seqs_default": 1}
    assert prepare_structure(zipped) == facts
    assert facts["source_sha256"] == json.loads((DATA / "provenance.json").read_text())["structures"]["1bna.cif"]["sha256"]
    assert "reference_manifest_id" not in canonical_json(geometry_profile()).decode()


@pytest.mark.integration
@pytest.mark.parametrize("corruption", ["name", "sequence", "identity", "pair", "model", "nested_model", "evidence", "profile", "source"])
def test_corrupt_facts_rejected_before_annotation(settings, corruption):
    facts = prepare_structure(DATA / "structures/1bna.cif")
    if corruption == "name":
        facts["chains"][0]["molecule_name"] = "old annotation"
    elif corruption == "sequence":
        facts["chains"][0]["structure_sequence"] = "a*"
    elif corruption == "identity":
        facts["chains"][0]["unique_chain_id"] = "1BNA:missing"
    elif corruption == "pair":
        facts["pairs"].append(copy.deepcopy(facts["pairs"][0]))
    elif corruption == "model":
        facts["chains"][0]["model_index"] = 2
    elif corruption == "nested_model":
        facts["chains"][0]["chain_identity"]["model_index"] = True
    elif corruption == "evidence":
        facts["chains"][0]["embedded_annotations"] = []
    elif corruption == "profile":
        facts["profile"]["all_atoms_radius"] = 4.
    else:
        facts["source_sha256"] = "not-a-hash"
    with pytest.raises(InputValidationError):
        validate_facts(facts)


@pytest.mark.integration
@pytest.mark.parametrize("corruption", ["missing", "duplicate", "endpoint", "distance", "insertion"])
def test_corrupt_contacts_rejected(settings, corruption):
    rows = []
    facts = prepare_structure(DATA / "structures/1bna.cif", contact_sink=rows.append)
    assert len(rows) == 652  # independent frozen 1BNA oracle
    if corruption == "missing":
        rows.pop()
    elif corruption == "duplicate":
        rows.insert(1, rows[0])
    elif corruption == "endpoint":
        rows[0][0] = "missing"
    elif corruption == "distance":
        rows[0][-1] = float("nan")
    else:
        rows[0][3] = "wrong"
    with pytest.raises(InputValidationError):
        list(validate_contacts(facts, rows))


@pytest.mark.integration
def test_prepared_detail_failure_is_not_published(tmp_path, settings):
    rows = []
    facts = prepare_structure(DATA / "structures/1bna.cif", contact_sink=rows.append)
    with pytest.raises(InputValidationError):
        analyze_prepared([facts], settings, contact_sources={"pdb_00001bna": lambda: iter(rows[:-1])})
    assert not list((tmp_path / "results").rglob("*.csv"))
    summary, = (tmp_path / "results").rglob("run_summary.json")
    assert json.loads(summary.read_text())["status"] == "failed"


@pytest.mark.integration
def test_unsupported_radius_fails_before_export(settings):
    facts = prepare_structure(DATA / "structures/1bna.cif")
    settings["distance_thresholds"]["ca_radius"] = 10.
    with pytest.raises(InputValidationError) as caught:
        analyze_prepared([facts], settings)
    assert caught.value.code == "STRUCTURE_PROFILE_UNSUPPORTED"


@pytest.mark.integration
def test_new_search_database_changes_annotation_without_repreparation(tmp_path, settings):
    import os
    import subprocess
    import sys

    facts = prepare_structure(DATA / "structures/1a8o.cif")
    facts_path = tmp_path / "facts.json"
    facts_path.write_bytes(canonical_json(facts))
    original = facts_path.read_bytes()
    observations = []
    for stage, accessions in (("before", {"P69905"}), ("after", {"P69905", "P12497"})):
        root = tmp_path / stage
        root.mkdir()
        (root / "sifts.tsv").write_text("# PDB\tCHAIN\tSP_PRIMARY\n")
        (root / "pdb_seqres.txt").write_text(">9xyz_A mol:protein length:1 Other\nA\n")
        records = (DATA / "references/swissprot.fasta").read_text().split(">")
        fasta = root / "swissprot.fasta"
        fasta.write_text("".join(">" + record for record in records if record and record.split("|", 2)[1] in accessions))
        database = root / "db"
        database.mkdir()
        subprocess.run(["makeblastdb", "-in", str(fasta), "-dbtype", "prot", "-out", str(database / "uniprot_db")],
                       check=True, capture_output=True, timeout=30)
        cfg = copy.deepcopy(settings)
        cfg.update(output_path=str(root / "out"), sifts_tsv_path=str(root / "sifts.tsv"),
                   pdb_fasta_path=str(root / "pdb_seqres.txt"), uniprot_fasta_path=str(fasta),
                   blast_db_path=str(database), reference_manifest_id=stage, export_detailed_interactions=False)
        cfg_path = root / "config.json"
        cfg_path.write_text(json.dumps(cfg))
        script = '''
import json,sys
from pathlib import Path
from pdb2net.structure_analysis import analyze_prepared
from pdb2net import pipeline
pipeline.process_single_file = lambda *a, **k: (_ for _ in ()).throw(AssertionError("No parsing allowed"))
paths=analyze_prepared([json.loads(Path(sys.argv[1]).read_text())],json.loads(Path(sys.argv[2]).read_text()))
Path(sys.argv[3]).write_text(paths.chain_dir)
'''
        result = subprocess.run([sys.executable, "-c", script, str(facts_path), str(cfg_path), str(root / "path")],
                                env={**{k:v for k,v in os.environ.items() if not k.startswith("PDB2NET_")}, "PDB2NET_BLAST_CACHE":"0"},
                                capture_output=True, text=True, timeout=30)
        assert result.returncode == 0, result.stdout + result.stderr
        network, = Path((root / "path").read_text()).glob("*.cx2")
        nodes = next(block["nodes"] for block in json.loads(network.read_text()) if "nodes" in block)
        observations.append(nodes[0]["v"].get("uniprot_id"))
    assert observations == [None, "P12497"]
    assert facts_path.read_bytes() == original


@pytest.mark.integration
def test_embedded_evidence_is_preserved_without_final_assignments(tmp_path, settings):
    source = Path(__file__).parent / "fixtures/nextgen/pdb_00008aly_xyz-enrich.cif"
    facts = prepare_structure(source)
    assert all("uniprot_id" not in c and "molecule_name" not in c for c in facts["chains"])
    assert facts["chains"][0]["embedded_annotations"]["uniprot"]
    config["network_annotations"]["use_embedded_sifts"] = False
    assert prepare_structure(source) == facts
    paths = analyze_prepared([facts], settings, contact_sources={"pdb_00008aly": lambda: iter(())})
    network, = Path(paths.chain_dir).glob("*.cx2")
    nodes = next(block["nodes"] for block in json.loads(network.read_text()) if "nodes" in block)
    assert nodes[0]["v"]["uniprot_id"] == "Q9H2K2"
    assert nodes[0]["v"]["annotation_source"] == "embedded_sifts"


@pytest.mark.integration
def test_changed_sifts_refilters_existing_geometry(tmp_path, settings):
    import os
    import subprocess
    import sys

    facts = prepare_structure(DATA / "structures/1bna.cif")
    assert facts["pairs"][0]["all_atoms_count"] == 652
    assert facts["pairs"][0]["ca_neighbors"] == 0
    payload = tmp_path / "facts.json"
    payload.write_bytes(canonical_json(facts))
    counts = []
    for stage in ("before", "after"):
        root = tmp_path / stage
        root.mkdir()
        sifts = root / "sifts.tsv"
        sifts.write_text("# PDB\tCHAIN\tSP_PRIMARY\n" + ("1bna\tA\tP69905\n1bna\tB\tP68871\n" if stage == "after" else ""))
        cfg = copy.deepcopy(settings)
        cfg.update(sifts_tsv_path=str(sifts), output_path=str(root / "out"),
                   export_detailed_interactions=False, reference_manifest_id=stage)
        cfg_path = root / "config.json"
        cfg_path.write_text(json.dumps(cfg))
        script = '''
import json,sys
from pathlib import Path
from pdb2net.structure_analysis import analyze_prepared
from pdb2net.config_loader import activate_config
from pdb2net import pipeline
cfg=json.loads(Path(sys.argv[2]).read_text()); activate_config(cfg)
raw=pipeline.run_pipeline([sys.argv[3]])
pipeline.process_single_file = lambda *a, **k: (_ for _ in ()).throw(AssertionError("No parsing allowed"))
prepared=analyze_prepared([json.loads(Path(sys.argv[1]).read_text())],cfg)
a=json.loads(Path(raw.summary_file).read_text());b=json.loads(Path(prepared.summary_file).read_text())
assert a['counts']['interactions']==b['counts']['interactions']
Path(sys.argv[4]).write_text(str(b['counts']['interactions']))
'''
        result = subprocess.run([sys.executable, "-c", script, str(payload), str(cfg_path),
                                 str((DATA / "structures/1bna.cif").resolve()), str(root / "count")],
                                env={**{k:v for k,v in os.environ.items() if not k.startswith("PDB2NET_")}, "PDB2NET_BLAST_CACHE":"0"},
                                capture_output=True, text=True, timeout=30)
        assert result.returncode == 0, result.stdout + result.stderr
        counts.append(int((root / "count").read_text()))
    assert counts == [1, 0]  # new protein typing needs CA contacts, which DNA coordinates lack
    assert payload.read_bytes() == canonical_json(facts)


@pytest.mark.integration
def test_serialized_contacts_preserve_minimum_witness_insertion_and_zero_occupancy(tmp_path, settings):
    import csv
    lines = ["HEADER    CONTACT WITNESSES                     01-JAN-00   1ABC\n"]
    serial = 1
    for chain, alternatives in (("A", (("A", 0., .8), ("B", 4., .2))), ("B", (("A", 7., .8), ("B", 6., .2)))):
        for n in range(1, 12):
            insertion = "A" if n == 1 else " "
            for alt, offset, occupancy in alternatives:
                x = (n-1)*30 + offset if n != 11 else offset
                occupancy = 0. if n == 11 else occupancy
                lines.append(f"ATOM  {serial:5d}  CA {alt}ALA {chain}{n:4d}{insertion}   {x:8.3f}{0.:8.3f}{0.:8.3f}{occupancy:6.2f} 20.00           C  \n")
                serial += 1
        lines.append("TER\n")
    source = tmp_path / "1abc.pdb"
    source.write_text("".join(lines) + "END\n")
    rows = []
    facts = prepare_structure(source, contact_sink=rows.append)
    assert facts["pairs"] == [{"chain_a":"1ABC:A", "chain_b":"1ABC:B", "all_atoms_count":10, "ca_neighbors":10}]
    assert len(list(validate_contacts(facts, rows))) == 10
    assert {row[-1] for row in rows} == {2.}
    inserted, = [row for row in rows if row[2] == "1"]
    assert inserted[2:6] == ["1", "A", "1A:ALA", "CA"]
    assert inserted[6:10] == ["1", "A", "1A:ALA", "CA"]
    assert all(row[2] != "11" and row[6] != "11" for row in rows)
    sifts = tmp_path / "sifts.tsv"
    sifts.write_text("1abc\tA\tP69905\n1abc\tB\tP68871\n")
    settings["sifts_tsv_path"] = str(sifts)
    out = analyze_prepared([facts], settings, contact_sources={"pdb_00001abc":lambda:iter(rows)})
    table, = Path(out.distances_dir).glob("*.csv")
    with table.open() as stream:
        exported = list(csv.DictReader(stream))
    assert len(exported) == 10 and {r["Distance"] for r in exported} == {"2.0"}
    assert {r["UniProt_A"] for r in exported} == {"P69905"}
    assert next(r for r in exported if r["Residue_A"] == "1A:ALA")["Residue_B"] == "1A:ALA"


@pytest.mark.integration
def test_prepared_public_provenance_matches_facts_and_rejects_forged_source(tmp_path, settings):
    from pdb2net.structure_facts import prepared_provenance, facts_id
    facts = prepare_structure(DATA / 'structures/1bna.cif')
    settings['export_detailed_interactions'] = False
    provenance = prepared_provenance([facts])
    provenance.update(catalog_id='a'*64, checked_at='2026-10-19T00:00:00+00:00')
    web = tmp_path / 'web'
    analyze_prepared([facts], settings, provenance=provenance, web_output_dir=str(web))
    summary = json.loads((web/'summary.json').read_text())
    assert summary['output_contract_version'] == '3.0'
    assert summary['prepared_inputs'] == provenance
    forged = copy.deepcopy(provenance)
    forged['entries'][0]['source_sha256'] = '0'*64
    forged['entries'][0]['revision_id'] = facts_id({**facts, 'source_sha256':'0'*64})
    with pytest.raises(InputValidationError, match='differs from analyzed inputs'):
        analyze_prepared([facts], settings, provenance=forged)


@pytest.mark.integration
def test_prepared_missing_output_path_keeps_original_error(tmp_path, settings, monkeypatch):
    settings['output_path'] = ''
    monkeypatch.chdir(tmp_path)
    with pytest.raises(InputValidationError) as error:
        analyze_prepared([], settings)
    assert error.value.code == 'OUTPUT_PATH_MISSING'
    assert list(tmp_path.iterdir()) == []


@pytest.mark.integration
def test_prepared_report_failure_does_not_mask_interrupt(settings, monkeypatch):
    from pdb2net import structure_analysis
    facts = prepare_structure(DATA / 'structures/1bna.cif')
    def interrupted(*args, **kwargs):
        raise KeyboardInterrupt
    def cannot_report(*args, **kwargs):
        raise OSError('Report destination unavailable')
    monkeypatch.setattr(pipeline, '_run_blast_annotation', interrupted)
    monkeypatch.setattr(structure_analysis, 'write_failed_run_manifest', cannot_report)
    with pytest.raises(KeyboardInterrupt):
        analyze_prepared([facts], settings)
