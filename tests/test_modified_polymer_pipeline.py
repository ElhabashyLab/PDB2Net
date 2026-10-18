"""Declared polymers must survive selection even without canonical residue names."""

import copy
import csv
import hashlib
import json
from pathlib import Path

import pytest

from pdb2net import config_loader, pipeline, precomputed
from pdb2net.uniprot_matcher import classify_molecule_type


def modified_cif(names, kind="polyribonucleotide", canonical="A"):
    header = (
        "data_1ABC\n_entry.id 1ABC\n_entity.id 1\n_entity.type polymer\n"
        f"_entity_poly.entity_id 1\n_entity_poly.type '{kind}'\n"
        "loop_\n_struct_asym.id\n_struct_asym.entity_id\nA 1\nB 1\nloop_\n"
    )
    fields = (
        "group_PDB id type_symbol label_atom_id label_alt_id label_comp_id label_asym_id "
        "label_entity_id label_seq_id pdbx_PDB_ins_code Cartn_x Cartn_y Cartn_z occupancy "
        "B_iso_or_equiv auth_seq_id auth_comp_id auth_asym_id auth_atom_id pdbx_PDB_model_num"
    )
    text = header + "".join(f"_atom_site.{field}\n" for field in fields.split())
    serial = 0
    for chain, residues, y in (("A", names, 0), ("B", [canonical], 1)):
        for number, residue in enumerate(residues, 1):
            serial += 1
            text += (
                f"ATOM {serial} P P . {residue} {chain} 1 {number} ? "
                f"{number * 3} {y} 0 1 20 {number} {residue} {chain} P 1\n"
            )
    return text


@pytest.fixture
def configured(tmp_path):
    original = copy.deepcopy(config_loader.get_config())
    cfg = json.loads((Path(__file__).parents[1] / "pdb2net/configs/config.base.json").read_text())
    references = {
        "pdb_fasta_path": ("pdb.fasta", ">9xyz_A mol:protein length:1 Unrelated\nA\n"),
        "sifts_tsv_path": ("sifts.tsv", "# PDB\tCHAIN\tSP_PRIMARY\n"),
        "uniprot_fasta_path": ("uniprot.fasta", ">sp|P12345|OTHER Unrelated\nA\n"),
    }
    for key, (name, contents) in references.items():
        path = tmp_path / name
        path.write_text(contents)
        cfg[key] = str(path)
    cfg.update(
        open_in_cytoscape=False, output_path=str(tmp_path / "raw"),
        blastp_executable=str(tmp_path / "absent-blastp"), blast_db_path=str(tmp_path / "absent-db"),
        blast_cache_path=str(tmp_path / "cache.sqlite3"), reference_manifest_id="modified-polymer-test",
        workers={"parsing": 1, "blast_threads": 1}, export_detailed_interactions=True,
        networks={"chain_per_pdb": True, "combined_chain_network": False,
                  "protein_per_pdb": False, "combined_protein_network": False},
    )
    config_loader.activate_config(cfg)
    try:
        yield tmp_path
    finally:
        config_loader.activate_config(original)


@pytest.mark.parametrize(("names", "kind", "canonical", "molecule"), [
    (["PSU"] * 3, "polyribonucleotide", "A", "RNA"),
    (["A", "PSU", "PSU"], "polyribonucleotide", "A", "RNA"),
    (["5MC"] * 3, "polydeoxyribonucleotide", "DA", "DNA"),
])
def test_modified_polymer_contacts_and_lengths_survive_raw_and_store(configured, names, kind, canonical, molecule):
    root = configured
    inputs = root / "inputs"
    inputs.mkdir()
    source = inputs / "1abc.cif"
    source.write_text(modified_cif(names, kind, canonical))
    raw = pipeline.run_pipeline(str(inputs))
    manifest = json.loads(Path(raw.manifest_file).read_text())
    assert manifest["status"] == "success" and manifest["warnings"] == []
    assert (manifest["counts"]["chains"], manifest["counts"]["interactions"]) == (2, 1)
    csv_path, = Path(raw.distances_dir).glob("*.csv")
    with csv_path.open() as stream:
        rows = list(csv.DictReader(stream))
    assert {(row["Residue_A"], row["Residue_B"], float(row["Distance"])) for row in rows} == {
        (f"1:{names[0]}", f"1:{canonical}", 1.0),
        (f"2:{names[1]}", f"1:{canonical}", 3.16),
    }
    assert len(rows) == 2

    cfg = copy.deepcopy(config_loader.get_config())
    cfg.update(export_detailed_interactions=False, output_path=str(root / "assembled"))
    config_loader.activate_config(cfg)
    store = root / "store"
    report = precomputed.precompute_sources(store, [source])
    assert report["failed"] == 0 and report["written"] == 1
    source.unlink()
    before = {p.relative_to(store): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in store.rglob("*") if p.is_file()}
    assembled = precomputed.run_assemble_pipeline(store, ["1abc"])
    assert before == {p.relative_to(store): hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in store.rglob("*") if p.is_file()}
    for run in (raw, assembled):
        cx2, = Path(run.chain_dir).glob("*.cx2")
        aspects = {key: value for block in json.loads(cx2.read_text()) for key, value in block.items()}
        assert {n["v"]["chain_id"]: (n["v"]["molecule_type"], n["v"]["nt_len"])
                for n in aspects["nodes"]} == {"A": (molecule, 3), "B": (molecule, 1)}
        edge, = aspects["edges"]
        assert edge["v"]["all_atoms_count"] == 2
        assert edge["v"]["interaction"] == f"{molecule}-{molecule}"


@pytest.mark.parametrize(("records", "expected"), [
    ([("UNK", "PeptideL")], "Protein"),
    ([("PSU", "Rna"), ("5MC", "Dna")], "DNA/RNA"),
    ([("PSU", "Rna"), ("UNK", "PeptideL")], "Unknown"),
    ([("PSU", "Rna"), ("UNK", "Unknown")], "Unknown"),
    ([("PSU", "Other")], "Unknown"),
    ([("ALA", "Rna")], "Protein"),  # Existing determinate name inference stays intact.
])
def test_polymer_type_only_fills_unambiguous_unknown_inference(records, expected):
    chain = {"residues": [{"residue_name": name, "is_polymer": True, "polymer_type": kind}
                          for name, kind in records]}
    before = copy.deepcopy(chain)
    assert classify_molecule_type(chain) == expected
    assert chain == before


def test_raw_run_with_only_quality_excluded_queries_needs_no_blast(configured):
    root = configured
    source = root / "1abc.pdb"
    source.write_text("HEADER".ljust(62) + "1ABC\n" + "".join(
        f"ATOM  {i:5d}  CA  {name} A{i:4d}    {i * 3.:8.3f}{0.:8.3f}{0.:8.3f}  1.00 20.00           C  \n"
        for i, name in enumerate(("ALA", "UNK"), 1)
    ) + "END\n")
    run = pipeline.run_pipeline([str(source)])
    manifest = json.loads(Path(run.manifest_file).read_text())
    assert manifest["status"] == "success" and manifest["counts"]["chains"] == 1
    assert not (root / "cache.sqlite3").exists()
