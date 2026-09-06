"""Opt-in integration checks using the documented external real-file dataset."""

from collections import Counter
import copy
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

import gemmi
import numpy as np
import pytest

from pdb2net import config_loader, distances, file_parser, pipeline

REPO = Path(__file__).resolve().parents[2]
NEXTGEN_NAME = "pdb_00008aly_xyz-enrich.cif.gz"
AF = "AF-P69905-F1-model_v6"
pytestmark = pytest.mark.integration


@pytest.fixture(autouse=True)
def isolated_config(integration_config, tmp_path):
    original = copy.deepcopy(config_loader.get_config())
    cfg = copy.deepcopy(integration_config)
    cfg["blast_cache_path"] = str(tmp_path / "cache.sqlite3")
    config_loader.activate_config(cfg)
    distances.coords_cache.clear()
    distances.tree_cache.clear()
    try:
        yield
    finally:
        config_loader.activate_config(original)
        distances.coords_cache.clear()
        distances.tree_cache.clear()


def command(case, arguments, *, env_extra=None):
    env = {k: v for k, v in os.environ.items() if not k.startswith("PDB2NET_")}
    env["MPLCONFIGDIR"] = str(case / "matplotlib")
    env["PYTHONIOENCODING"] = "utf-8"
    env.update(env_extra or {})
    case.mkdir(exist_ok=True, parents=True)
    result = subprocess.run(
        [sys.executable, *arguments],
        cwd=REPO,
        env=env,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=90,
    )
    (case / "command.json").write_text(
        json.dumps(
            {
                "argv": [sys.executable, *arguments],
                "returncode": result.returncode,
                "env_overrides": env_extra or {},
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    (case / "stdout.log").write_text(result.stdout, encoding="utf-8")
    (case / "stderr.log").write_text(result.stderr, encoding="utf-8")
    return result


def run(case, files, overrides=None, *, expect=0, legacy=False):
    inputs = case / "inputs"
    inputs.mkdir(parents=True)
    for item in files:
        source, name = item if isinstance(item, tuple) else (item, item.name)
        shutil.copyfile(source, inputs / name)
    cfg = copy.deepcopy(config_loader.get_config())
    cfg.update(overrides or {})
    cfg.update(
        input_folder_path=str(inputs),
        output_path=str(case / "work"),
        blast_cache_path=str(case / "blast-cache.sqlite3"),
    )
    config_path = case / "config.json"
    config_path.write_text(json.dumps(cfg, indent=2), encoding="utf-8")
    if legacy:
        result = command(
            case,
            ["-m", "pdb2net.main"],
            env_extra={"PDB2NET_CONFIG_FILE": str(config_path)},
        )
        summaries = list((case / "work").glob("*/run_summary.json"))
        summary = (
            json.loads(summaries[0].read_text(encoding="utf-8")) if summaries else None
        )
    else:
        result = command(
            case,
            [
                "-m",
                "pdb2net",
                "run",
                "--input-dir",
                str(inputs),
                "--output-dir",
                str(case / "work"),
                "--web-output-dir",
                str(case / "web"),
                "--config",
                str(config_path),
                "--headless",
            ],
        )
        path = case / "web/summary.json"
        summary = (
            json.loads(path.read_text(encoding="utf-8")) if path.exists() else None
        )
    assert result.returncode == expect, result.stdout[-4000:] + result.stderr
    assert summary is not None
    assert summary["status"] == ("success" if expect == 0 else "failed")
    if not legacy and expect == 0:
        assert summary["errors"] == []
        assert summary["output_contract_version"] == "2.0"
        assert summary["counts"]["networks"] == len(
            list((case / "web/networks").glob("*.cx2"))
        )
        for artifact in (
            summary["artifacts"]["networks"] + summary["artifacts"]["interactions"]
        ):
            path = case / "web" / artifact["path"]
            assert hashlib.sha256(path.read_bytes()).hexdigest() == artifact["sha256"]
        for path in (case / "web/networks").glob("*.cx2"):
            network(path)
    return summary


def network(path):
    blocks = json.loads(
        path.read_text(encoding="utf-8"),
        parse_constant=lambda value: (_ for _ in ()).throw(AssertionError(value)),
    )
    assert blocks[0] == {"CXVersion": "2.0", "hasFragments": False}
    aspects = {key: value for block in blocks[1:] for key, value in block.items()}
    ids = {node["id"] for node in aspects["nodes"]}
    assert len(ids) == len(aspects["nodes"])
    names = {node["id"]: node["v"]["name"] for node in aspects["nodes"]}
    for node in aspects["nodes"]:
        assert math.isfinite(node["x"]) and math.isfinite(node["y"])
    for edge in aspects["edges"]:
        assert edge["s"] in ids and edge["t"] in ids
    fields = (
        "name",
        "molecule_type",
        "uniprot_id",
        "aa_len",
        "nt_len",
        "node_kind",
        "pdb_count",
        "source_chains",
        "color_group",
    )
    nodes = sorted(
        ({key: node["v"].get(key) for key in fields} for node in aspects["nodes"]),
        key=lambda row: str(row["name"]),
    )
    edges = sorted(
        tuple(sorted((names[edge["s"]], names[edge["t"]])))
        + (json.dumps(edge["v"], sort_keys=True),)
        for edge in aspects["edges"]
    )
    return nodes, edges


def network_set(case):
    return {path.name: network(path) for path in (case / "web/networks").glob("*.cx2")}


REAL_CASES = (
    "4hhb.pdb",
    "9rat.pdb",
    "1TUP.pdb",
    "1tup.cif",
    "1a34.cif",
    "124d.cif",
    "170d.cif",
    "5rlu.pdb",
    "6w41.cif",
    "6m17.cif",
    "1a1t.pdb",
    "1a1t.cif",
    "1a8o.pdb",
    "1a8o.cif",
    "8aly.cif",
    NEXTGEN_NAME,
    AF + ".pdb",
    AF + ".cif",
)


@pytest.mark.parametrize("filename", REAL_CASES)
def test_real_file_through_public_cli(filename, tmp_path, data_files):
    source = data_files[filename]
    reference = gemmi.read_structure(str(source))
    actual = file_parser.parse_structure_input(str(source))["structure"]
    assert len(actual) == len(reference)
    for model, expected in zip(actual, reference):
        assert [chain.name for chain in model] == [chain.name for chain in expected]

        def atoms(chains):
            return Counter(
                (
                    c.name,
                    r.name,
                    str(r.seqid),
                    a.name,
                    a.altloc,
                    a.element.name,
                    tuple(a.pos),
                    a.occ,
                )
                for c in chains
                for r in c
                for a in r
            )

        assert atoms(model) == atoms(expected)
    processed = pipeline.process_single_file(str(source))
    retained = {
        chain["chain_id"]: sum(len(res["atoms"]) for res in chain["residues"])
        for chain in processed["atom_data"]
    }
    assert len(retained) == len(processed["atom_data"])
    expected_atoms = {
        chain.name: sum(
            atom.element.name not in ("H", "D") for res in chain for atom in res
        )
        for chain in reference[0]
    }
    for chain_id, count in retained.items():
        assert count == expected_atoms[chain_id]
    summary = run(tmp_path / "run", [source])
    assert summary["counts"]["structures"] == 1
    assert summary["counts"]["chains"] > 0
    if source == data_files[NEXTGEN_NAME]:
        assert summary["structure_inputs"][0]["kind"] == "nextgen_enriched_mmcif"
        assert (
            summary["annotations"]["embedded_sifts"]["segments_by_database"]["uniprot"]
            > 0
        )
    if source.name.startswith(AF):
        assert summary["annotations"]["by_source"].get("direct_uniprot") == 1


def test_1tup_pdb_mmcif_network_equivalence(tmp_path, data_files):
    for label, source in [
        ("pdb", data_files["1TUP.pdb"]),
        ("cif", data_files["1tup.cif"]),
    ]:
        run(tmp_path / label, [source])
    assert network_set(tmp_path / "pdb") == network_set(tmp_path / "cif")


def test_only_renamed_official_file_keeps_identity_and_graph(tmp_path, data_files):
    run(tmp_path / "original", [data_files["1tup.cif"]])
    summary = run(
        tmp_path / "renamed", [(data_files["1tup.cif"], "freier_name_123456.cif")]
    )
    assert summary["identities"][0]["canonical_id"] == "pdb_00001tup"
    assert any(
        row["code"] == "STRUCTURE_FILENAME_ID_MISMATCH" for row in summary["warnings"]
    )
    assert network_set(tmp_path / "original") == network_set(tmp_path / "renamed")


@pytest.mark.parametrize("extension", ["pdb", "cif"])
def test_renamed_alphafold_is_annotated(extension, tmp_path, data_files):
    summary = run(
        tmp_path / "renamed",
        [(data_files[f"{AF}.{extension}"], f"prediction.{extension}")],
    )
    assert summary["counts"]["chains"] == 1
    assert summary["annotations"]["chains_annotated"] == 1
    # CIF retains its AF entry identity; PDB needs sequence-search fallback.
    expected = "direct_uniprot" if extension == "cif" else "blastp_swissprot"
    assert summary["annotations"]["by_source"].get(expected) == 1


def test_conflict_and_duplicate_inputs_are_diagnostic_only(tmp_path, data_files):
    for label, files, code in [
        ("conflict", [data_files["123456.cif"]], "CONFLICTING_STRUCTURE_IDENTITY"),
        (
            "duplicate",
            [data_files["1TUP.pdb"], data_files["1tup.cif"]],
            "DUPLICATE_STRUCTURE_IDENTITY",
        ),
    ]:
        summary = run(tmp_path / label, files, expect=1)
        assert summary["errors"][0]["code"] == code
        assert summary["networks"] == [] and summary["interactions"] == []


def test_full_nextgen_geometry_segments_and_gzip(tmp_path, data_files):
    uncompressed = tmp_path / "pdb_00008aly_xyz-enrich.cif"
    payload = gzip.decompress(data_files[NEXTGEN_NAME].read_bytes())
    uncompressed.write_bytes(payload)
    parsed = file_parser.parse_structure_input(str(data_files[NEXTGEN_NAME]))
    block = gemmi.cif.read_string(payload.decode())[0]
    raw = block.get_mmcif_category("_pdbx_sifts_unp_segments.")
    observed = [
        segment
        for chain in parsed["embedded_annotations"]["by_chain"].values()
        for segment in chain.get("uniprot", [])
    ]
    assert len(observed) == len(raw["unp_acc"])
    assert Counter(segment["accession"] for segment in observed) == Counter(
        raw["unp_acc"]
    )

    # Compare the actual files, without relying on a previous review inventory.
    def atom_records(structure):
        return Counter(
            (
                i,
                c.name,
                r.name,
                str(r.seqid),
                a.name,
                a.altloc,
                a.element.name,
                tuple(a.pos),
                a.occ,
            )
            for i, model in enumerate(structure)
            for c in model
            for r in c
            for a in r
        )

    regular = gemmi.read_structure(str(data_files["8aly.cif"]))
    enriched = gemmi.read_structure(str(data_files[NEXTGEN_NAME]))
    assert atom_records(regular) == atom_records(enriched)
    for label, source in [("gzip", data_files[NEXTGEN_NAME]), ("plain", uncompressed)]:
        run(tmp_path / label, [source])
    assert network_set(tmp_path / "gzip") == network_set(tmp_path / "plain")


def test_nextgen_tooltip_selection_preserves_networks(tmp_path, data_files):
    run(tmp_path / "uniprot", [data_files[NEXTGEN_NAME]])
    summary = run(
        tmp_path / "all_fields",
        [data_files[NEXTGEN_NAME]],
        {
            "network_annotations": {
                "use_embedded_sifts": True,
                "tooltip_fields": ["uniprot", "pfam", "cath", "scop2"],
                "max_tooltip_segments_per_database": 20,
            }
        },
    )
    assert network_set(tmp_path / "uniprot") == network_set(tmp_path / "all_fields")
    assert summary["annotations"]["embedded_sifts"]["segments_by_database"]["pfam"] > 0
    assert any(
        "Pfam" in path.read_text(encoding="utf-8")
        for path in (tmp_path / "all_fields/web/networks").glob("*.cx2")
    )


def test_real_nmr_all_models_do_not_interact_with_each_other(tmp_path, data_files):
    summary = run(
        tmp_path / "all",
        [data_files["1a1t.cif"]],
        {
            "structure_model_policy": "all",
            "networks": {
                "chain_per_pdb": True,
                "combined_chain_network": True,
                "protein_per_pdb": False,
                "combined_protein_network": False,
            },
        },
    )
    assert summary["counts"]["chains"] == 50
    config_loader.config["structure_model_policy"] = "all"
    parsed = pipeline.process_single_file(str(data_files["1a1t.cif"]))
    from pdb2net.uniprot_matcher import classify_molecule_type

    for chain in parsed["atom_data"]:
        chain["molecule_type"] = classify_molecule_type(chain)
    interactions = distances.calculate_distances_with_ckdtree([parsed])
    chain_models = {
        chain["unique_chain_id"]: chain["model_index"] for chain in parsed["atom_data"]
    }
    assert {edge["model_index"] for edge in interactions} == set(range(1, 26))
    assert all(
        chain_models[edge["chain_a"]] == chain_models[edge["chain_b"]]
        for edge in interactions
    )


def test_real_mse_normalization_and_hydrogen_filter(data_files):
    parsed = pipeline.process_single_file(str(data_files["1a8o.cif"]))
    residues = [res for chain in parsed["atom_data"] for res in chain["residues"]]
    modified = [res for res in residues if res.get("original_residue_name") == "MSE"]
    assert len(modified) == 4 and all(res["residue_name"] == "MET" for res in modified)
    source = data_files["124d.cif"]
    structure = gemmi.read_structure(str(source))
    atoms = [atom for chain in structure[0] for residue in chain for atom in residue]
    assert any(atom.element.name == "H" for atom in atoms)
    expected = sum(atom.element.name not in ("H", "D") for atom in atoms)
    parsed = pipeline.process_single_file(str(source))
    assert (
        sum(
            len(res["atoms"])
            for chain in parsed["atom_data"]
            for res in chain["residues"]
        )
        == expected
    )


def test_real_nmr_csv_normalizes_model_policy(tmp_path, data_files):
    summary = run(
        tmp_path / "all_csv",
        [data_files["1a1t.cif"]],
        {
            "structure_model_policy": " ALL ",
            "export_detailed_interactions": True,
            "networks": {
                "chain_per_pdb": True,
                "combined_chain_network": True,
                "protein_per_pdb": False,
                "combined_protein_network": False,
            },
        },
    )
    assert summary["counts"]["chains"] == 50
    (name,) = summary["interactions"]
    with (tmp_path / "all_csv/web" / name).open(encoding="utf-8") as stream:
        pairs = {
            tuple(sorted((row["Chain_A"], row["Chain_B"])))
            for row in csv.DictReader(stream)
        }
    assert pairs == {(f"1A1T:model{i}:A", f"1A1T:model{i}:B") for i in range(1, 26)}


def test_independent_contact_count_oracle_and_rigid_transform():
    rng = np.random.default_rng(728)
    a, b = rng.normal(size=(23, 3)), rng.normal(size=(19, 3))
    for radius in [0.0, 0.5, 1.0, 5.0]:
        expected = sum(float(np.linalg.norm(x - y)) <= radius for x in a for y in b)
        for left, right in [(a, b), (a[::-1], b[::-1]), (a + [5, 8, 2], b + [5, 8, 2])]:
            assert (
                distances.count_nearby_atoms(
                    distances.cKDTree(left), distances.cKDTree(right), radius
                )
                == expected
            )
    for coordinate, expected in [(1 - 1e-9, 1), (1.0, 1), (1 + 1e-9, 0)]:
        assert (
            distances.count_nearby_atoms(
                distances.cKDTree([[0, 0, 0]]),
                distances.cKDTree([[coordinate, 0, 0]]),
                1.0,
            )
            == expected
        )


def test_legacy_cli_and_real_detailed_csv(tmp_path, data_files):
    summary = run(tmp_path / "legacy", [data_files["1TUP.pdb"]], legacy=True)
    assert summary["counts"]["structures"] == 1
    summary = run(
        tmp_path / "csv",
        [data_files["170d.cif"]],
        {"export_detailed_interactions": True},
    )
    assert summary["interactions"]
    for name in summary["interactions"]:
        with (tmp_path / "csv/web" / name).open(encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        assert rows


def test_worker_and_batch_scheduling_preserves_networks(tmp_path, data_files):
    sources = [data_files["1TUP.pdb"], data_files["4hhb.pdb"]]
    for label, workers, batch in [("one", 1, 550_000), ("two", 2, 2_000_000)]:
        limits = copy.deepcopy(config_loader.config["resource_limits"])
        limits["max_processing_batch_bytes"] = batch
        run(
            tmp_path / label,
            sources if workers == 1 else sources[::-1],
            {
                "workers": {"parsing": workers, "blast_threads": 1},
                "resource_limits": limits,
            },
        )
    assert network_set(tmp_path / "one") == network_set(tmp_path / "two")


def test_real_precompute_assemble_equivalence_and_readonly_store(tmp_path, data_files):
    raw = tmp_path / "raw"
    run(raw, [data_files["6m17.cif"], data_files["6w41.cif"], data_files[NEXTGEN_NAME]])
    store = tmp_path / "store"
    result = command(
        tmp_path / "precompute",
        [
            "-m",
            "pdb2net",
            "precompute",
            "--input-dir",
            str(raw / "inputs"),
            "--store",
            str(store),
            "--config",
            str(raw / "config.json"),
            "--headless",
        ],
    )
    assert result.returncode == 0, result.stdout[-4000:] + result.stderr

    def snapshot():
        return {
            str(p.relative_to(store)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in store.rglob("*")
            if p.is_file()
        }

    before = snapshot()
    for path in store.rglob("*"):
        path.chmod(0o555 if path.is_dir() else 0o444)
    store.chmod(0o555)
    assembled = tmp_path / "assembled"
    args = [
        "-m",
        "pdb2net",
        "assemble",
        "--store",
        str(store),
        "--pdb-id",
        "6m17",
        "--pdb-id",
        "6w41",
        "--pdb-id",
        "8aly",
        "--output-dir",
        str(assembled / "work"),
        "--web-output-dir",
        str(assembled / "web"),
        "--config",
        str(raw / "config.json"),
        "--headless",
    ]
    result = command(assembled, args)
    assert result.returncode == 0, result.stdout[-4000:] + result.stderr
    assert before == snapshot()
    assert network_set(raw) == network_set(assembled)
    assert any("Protein" in name for name in network_set(raw))
    assert any("Combined" in name for name in network_set(raw))
    comparison = command(
        tmp_path / "cx2-comparison",
        [
            "scripts/compare_cx2_outputs.py",
            "compare",
            "--actual",
            str(raw / "web/networks"),
            "--expected",
            str(assembled / "web/networks"),
            "--out",
            str(tmp_path / "cx2-comparison.md"),
            "--json-out",
            str(tmp_path / "cx2-comparison.json"),
        ],
    )
    assert comparison.returncode == 0, (tmp_path / "cx2-comparison.md").read_text(
        encoding="utf-8"
    )
    assert (
        json.loads((tmp_path / "cx2-comparison.json").read_text(encoding="utf-8"))[
            "status"
        ]
        == "PASS"
    )
    missing = command(
        tmp_path / "missing",
        [
            "-m",
            "pdb2net",
            "assemble",
            "--store",
            str(store),
            "--pdb-id",
            "9rat",
            "--output-dir",
            str(tmp_path / "missing/work"),
            "--config",
            str(raw / "config.json"),
            "--headless",
        ],
    )
    assert missing.returncode != 0
    assert "PRECOMPUTED_ENTRY_MISSING" in missing.stdout + missing.stderr
    assert before == snapshot()
