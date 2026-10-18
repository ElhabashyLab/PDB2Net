"""Real-file integration checks using the frozen, reviewed repository dataset."""

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
import pytest

from pdb2net import config_loader, distances, file_parser, pipeline
from pdb2net.contact_atoms import ContactView
from pdb2net.uniprot_matcher import extract_sequence_from_parsed_data
from pdb2net.unknown_molecule_uniprot import process_molecule_info
from reference_oracle import DATASET, STANDARD, contact_pairs, fasta_records, read_model

REPO = Path(__file__).resolve().parents[2]
NEXTGEN_NAME = "pdb_00008aly_xyz-enrich.cif.gz"
AF = "AF-P69905-F1-model_v6"
pytestmark = pytest.mark.integration


@pytest.fixture(scope="module")
def reference_expected():
    return {name: json.loads((DATASET / "expected" / f"{name}.json").read_text())
            for name in ("structures", "networks", "search")}


def assert_extracted_facts(processed, expected):
    chains = {chain["chain_id"]: chain for chain in processed["atom_data"]}
    assert chains.keys() == expected["chains"].keys()
    for label, actual in chains.items():
        facts = expected["chains"][label]
        polymer = [r for r in actual["residues"] if r["is_polymer"]]
        assert sum(len(r["atoms"]) for r in actual["residues"]) == facts["raw_heavy_atoms"]
        assert len(polymer) == facts["polymer_residues"]
        view = ContactView(actual)
        assert len(view.positions) == facts["contact_positions"]
        assert len(view.groups) == facts["contact_sites"]
        assert sum(a["occupancy"] == 0 for r in polymer for a in r["atoms"]) == facts["zero_occupancy_polymer_atoms"]
        def residue_label(r):
            return f"{r['residue_number']}{r.get('insertion_code', '')}"
        assert [residue_label(r) for r in polymer if r.get("insertion_code")] == facts["insertion_codes"]
        assert [[residue_label(r), r.get("original_residue_name", r["residue_name"]), r["residue_name"]]
                for r in polymer if r.get("original_residue_name", r["residue_name"]) not in STANDARD] == facts["modified_residues"]
        if "protein_sequence" in facts:
            assert extract_sequence_from_parsed_data(actual) == facts["protein_sequence"]


def assert_contacts(actual, expected, *, model=1):
    observed = {tuple(sorted([r["chain_a"].split(":")[-1], r["chain_b"].split(":")[-1]])): r
                for r in actual if r.get("model_index", 1) == model}
    accepted = {tuple(row["chains"]): row for row in expected if row["accepted"]}
    assert observed.keys() == accepted.keys()
    for pair, row in accepted.items():
        assert observed[pair]["all_atoms_count"] == row["all_atoms_count"]
        assert observed[pair]["interaction_type"] == row["interaction"]
        if row["interaction"] == "Protein-Protein":
            assert observed[pair]["ca_neighbors"] == row["ca_neighbors"]


def chain_export(case, summary):
    names = [name for name in summary["networks"] if Path(name).name.startswith("Chain_Interaction_Network_")]
    assert len(names) == 1
    return network(case / "web" / names[0])


def protein_projection(facts, contacts):
    """Apply the documented projection to independent chain expectations."""
    chains = facts["chains"]
    nodes = {c["uniprot_id"] for c in chains.values() if c["uniprot_id"]}
    edges = Counter()
    for row in contacts:
        if not row["accepted"]:
            continue
        a, b = row["chains"]
        labels = []
        for chain_id in (a, b):
            c = chains[chain_id]
            labels.append(c["uniprot_id"] if c["molecule_type"] == "Protein"
                          else f"{facts['display_id']}:{chain_id}")
        if "Protein" not in (chains[a]["molecule_type"], chains[b]["molecule_type"]):
            continue
        if all(labels) and labels[0] != labels[1]:
            nodes.update(labels)
            edges[tuple(sorted(labels))] += row["all_atoms_count"]
    return nodes, dict(edges)


def exported_weights(graph, *, decorated=False):
    nodes, edges = graph
    label = (lambda name: name.split("\n", 1)[0]) if decorated else (lambda name: name)
    pairs = {tuple(sorted([label(a), label(b)])): json.loads(v)["all_atoms_count"] for a, b, v in edges}
    assert len(pairs) == len(edges)
    return {label(n["name"]) for n in nodes}, pairs


def assert_protein_export(case, summary, facts, contacts):
    expected_nodes, expected_edges = protein_projection(facts, contacts)
    names = [n for n in summary["networks"] if Path(n).name.startswith("Protein_Network_")]
    assert len(names) == int(bool(expected_nodes))
    if names:
        assert exported_weights(network(case / "web" / names[0])) == (expected_nodes, expected_edges)


def assert_combined_exports(case, filenames, expected):
    """Independent components and sums from fixed chain facts, including zero bridges."""
    nodes, edges = {}, {}
    for filename in filenames:
        facts = expected["structures"][filename]
        for label, chain in facts["chains"].items():
            nodes[f"{facts['display_id']}:{label}"] = (filename, chain["uniprot_id"])
        for row in expected["networks"][filename]:
            if row["accepted"]:
                pair = tuple(sorted(f"{facts['display_id']}:{c}" for c in row["chains"]))
                edges[pair] = row["all_atoms_count"]
    identity = set()
    labels = sorted(nodes)
    for i, a in enumerate(labels):
        for b in labels[i + 1:]:
            if nodes[a][0] != nodes[b][0] and nodes[a][1] and nodes[a][1] == nodes[b][1]:
                identity.add((a, b))
                edges[a, b] = 0
    adjacency = {n: set() for n in nodes}
    for a, b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    components = []
    unseen = set(nodes)
    while unseen:
        pending = [min(unseen)]
        component = set()
        while pending:
            n = pending.pop()
            if n not in component:
                component.add(n)
                pending.extend(adjacency[n] - component)
        unseen -= component
        if len({nodes[n][0] for n in component}) > 1:
            components.append(component)
    actual = network_set(case)
    chain_graphs = []
    actual_identity = set()
    for name, graph in actual.items():
        if not name.startswith("Combined_Network_"):
            continue
        exported_nodes, weights = exported_weights(graph, decorated=True)
        bridges = {tuple(sorted((a.split("\n", 1)[0], b.split("\n", 1)[0])))
                   for a, b, value in graph[1] if json.loads(value)["interaction"] == "identity"}
        assert all(weights[pair] == 0 for pair in bridges)
        assert not actual_identity.intersection(bridges)
        actual_identity.update(bridges)
        chain_graphs.append((exported_nodes, {p: w for p, w in weights.items() if p not in bridges}))
    protein_graphs = [exported_weights(v, decorated=True) for k, v in actual.items() if k.startswith("Combined_Protein_Network_")]
    assert len(chain_graphs) == len(protein_graphs) == len(components)
    for component in components:
        contact_edges = {pair: count for pair, count in edges.items()
                         if set(pair) <= component and pair not in identity}
        assert (component, contact_edges) in chain_graphs
        proteins = {nodes[n][1] for n in component if nodes[n][1]}
        totals = Counter()
        for (a, b), count in contact_edges.items():
            x, y = nodes[a][1], nodes[b][1]
            if x and y and x != y:
                totals[tuple(sorted([x, y]))] += count
        assert (proteins, dict(totals)) in protein_graphs
    # Complete cross-file pairs above define connectivity independently. The
    # export must instead contain a spanning tree for each shared accession.
    assert actual_identity <= identity
    for accession in {up for _filename, up in nodes.values() if up}:
        group = {node for node, (_filename, up) in nodes.items() if up == accession}
        if len({nodes[node][0] for node in group}) < 2:
            continue
        bridges = {pair for pair in actual_identity if set(pair) <= group}
        assert len(bridges) == len(group) - 1
        reached = {min(group)}
        while True:
            expanded = reached | {node for pair in bridges if reached.intersection(pair) for node in pair}
            if expanded == reached:
                break
            reached = expanded
        assert reached == group


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
    # Isolate machine-local configs as well as PDB2NET_* environment overrides.
    defaults = case / "defaults"
    defaults.mkdir(exist_ok=True)
    shutil.copyfile(DATASET / "config.json", defaults / "config.base.json")
    env["PDB2NET_CONFIG_DIR"] = str(defaults)
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
    "1bna.cif",
    "1ehz.cif",
    "2erk.cif",
    "1ao7.cif",
)


@pytest.mark.parametrize("filename", REAL_CASES)
def test_real_file_through_public_cli(filename, tmp_path, data_files, reference_expected):
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
    expected = reference_expected["structures"][filename]
    assert len(reference) == expected["models"]
    assert_extracted_facts(processed, expected)
    process_molecule_info([processed])
    pipeline._run_blast_annotation([processed])
    for chain in processed["atom_data"]:
        facts = expected["chains"][chain["chain_id"]]
        assert chain["molecule_type"] == facts["molecule_type"]
        assert chain.get("uniprot_id") == facts["uniprot_id"]
        if facts["annotation_source"] and facts["annotation_source"] != "direct_uniprot":
            assert chain["annotation_source"] == facts["annotation_source"]
        if len(facts["sifts_accessions"]) > 1:
            assert sorted(chain["external_sifts_accessions"]) == facts["sifts_accessions"]
    assert_contacts(distances.calculate_distances_with_ckdtree([processed]), reference_expected["networks"][filename])
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
    nodes, edges = chain_export(tmp_path / "run", summary)
    assert {n["name"].split(":")[-1] for n in nodes} == set(expected["chains"])
    for node in nodes:
        facts = expected["chains"][node["name"].split(":")[-1]]
        assert node["molecule_type"] == facts["molecule_type"]
        assert node["uniprot_id"] == facts["uniprot_id"]
        field = "aa_len" if facts["molecule_type"] == "Protein" else "nt_len"
        assert node[field] == facts["polymer_residues"]
    assert sorted((a.split(":")[-1], b.split(":")[-1], json.loads(v)["all_atoms_count"])
                  for a, b, v in edges) == sorted((*row["chains"], row["all_atoms_count"])
                                                 for row in reference_expected["networks"][filename] if row["accepted"])
    assert_protein_export(tmp_path / "run", summary, expected, reference_expected["networks"][filename])
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
        if label == "conflict":
            assert summary["errors"][0]["code"] == "NO_PARSEABLE_STRUCTURES"
            assert summary["warnings"][0]["code"] == "INPUT_FILE_SKIPPED"
            assert code in summary["warnings"][0]["message"]
            assert "123456.cif" in summary["warnings"][0]["message"]
        else:
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


def test_real_nmr_all_models_do_not_interact_with_each_other(tmp_path, data_files, reference_expected):
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
    for model, expected in reference_expected["networks"]["1a1t.cif#all"].items():
        assert_contacts(interactions, expected, model=int(model))
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


@pytest.mark.parametrize("extension", ["pdb", "cif"])
def test_real_hydrated_protein_fallback_uses_polymer_sequence(tmp_path, data_files, extension, reference_expected):
    from pdb2net.uniprot_matcher import classify_molecule_type, extract_sequence_from_parsed_data

    source = data_files[f"1a8o.{extension}"]
    (chain,) = pipeline.process_single_file(str(source))["atom_data"]
    assert len(chain["residues"]) == 158
    assert sum(res["is_polymer"] for res in chain["residues"]) == 70
    assert sum(len(res["atoms"]) for res in chain["residues"]) == 644
    # 556 polymer heavy atoms minus 3 explicitly unoccupied positions; 88 waters
    # remain in raw data but cannot contribute contact coordinates.
    assert len(distances.extract_coordinates(chain, "all_atoms")) == 553
    modified = [res for res in chain["residues"] if res.get("original_residue_name") == "MSE"]
    assert len(modified) == 4 and all(res["is_polymer"] for res in modified)
    # The 70-position expectation comes from the fixture's independent SEQRES
    # reference, not from a previous PDB2Net run; its four MSE residues map to M.
    expected_sequence = "MDIRQGPKEPFRDYVDRFYKTLRAEQASQEVKNWMTETLLVQNANPDCKTILKALGPGATLEEMMTACQG"
    assert extract_sequence_from_parsed_data(chain) == expected_sequence
    assert classify_molecule_type(chain) == "Protein"
    waterless = {
        **chain,
        "residues": [res for res in chain["residues"] if res["residue_name"] != "HOH"],
    }
    assert extract_sequence_from_parsed_data(waterless) == expected_sequence

    # Force actual BLAST fallback instead of letting known SIFTS/FASTA labels
    # hide the bug. Swiss-Prot and the tiny real BLAST database remain unchanged.
    sifts = tmp_path / "no-mapping.tsv"
    sifts.write_text("# PDB\tCHAIN\tSP_PRIMARY\n", encoding="utf-8")
    fasta = tmp_path / "no-mapping.fasta"
    fasta.write_text(">1abc_A mol:protein length:1 Unrelated\nA\n", encoding="utf-8")
    case = tmp_path / "fallback"
    summary = run(case, [source], {
        "sifts_tsv_path": str(sifts), "pdb_fasta_path": str(fasta),
    })
    assert summary["annotations"]["by_source"] == {"blastp_swissprot": 1}
    nodes, edges = network(case / "web/networks/Chain_Interaction_Network_1A8O.cx2")
    (node,) = nodes
    assert (node["molecule_type"], node["uniprot_id"], node["aa_len"]) == (
        "Protein", reference_expected["search"]["1a8o_fallback"]["A"], 70,
    )
    assert edges == []


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


def test_pure_dna_csv_matches_every_independent_atom_pair(tmp_path, data_files):
    source = data_files["1bna.cif"]
    _, chains = read_model(source)
    assert set(chains) == {"A", "B"}
    expected = contact_pairs(chains["A"]["positions"], chains["B"]["positions"], 5.0)
    summary = run(tmp_path / "dna-csv", [source], {"export_detailed_interactions": True})
    assert len(summary["interactions"]) == 1
    with (tmp_path / "dna-csv/web" / summary["interactions"][0]).open() as stream:
        rows = list(csv.DictReader(stream))
    labels_a = {p["site"]: p["residue"] for p in chains["A"]["positions"]}
    labels_b = {p["site"]: p["residue"] for p in chains["B"]["positions"]}
    observed = {}
    for row in rows:
        assert (row["PDB_ID"], row["Chain_A"], row["Chain_B"]) == ("1BNA", "A", "B")
        assert row["Interaction_Type"] == "DNA-DNA"
        a = (int(row["Residue_A"].split(":")[0]), "", row["Atom_A"])
        b = (int(row["Residue_B"].split(":")[0]), "", row["Atom_B"])
        assert row["Residue_A"] == labels_a[a] and row["Residue_B"] == labels_b[b]
        assert (a, b) not in observed
        observed[a, b] = float(row["Distance"])
    assert observed.keys() == expected.keys()
    assert observed == {pair: round(distance, 2) for pair, distance in expected.items()}


@pytest.mark.parametrize(("filename", "lengths"), [("1ehz.cif", [76]), ("170d.cif", [12, 12])])
def test_modified_nucleotide_lengths_include_all_observed_polymer_residues(filename, lengths, tmp_path, data_files):
    summary = run(tmp_path / "length", [data_files[filename]])
    nodes, _ = chain_export(tmp_path / "length", summary)
    actual = sorted(n["nt_len"] for n in nodes)
    assert actual == lengths


def test_nonpolymer_amino_acids_do_not_inflate_protein_length(tmp_path, data_files):
    summary = run(tmp_path / "length", [data_files["6m17.cif"]])
    nodes, _ = chain_export(tmp_path / "length", summary)
    actual = {n["name"].split(":")[-1]: n["aa_len"] for n in nodes}
    assert actual.keys() == {"A", "B", "C", "D", "E", "F"}
    assert (actual["A"], actual["C"]) == (605, 605)


def test_polymer_lengths_survive_precompute_and_readonly_assembly(tmp_path, data_files):
    from pdb2net import precomputed

    expected = {
        "1ehz": {"A": (0, 76)},
        "170d": {"A": (0, 12), "B": (0, 12)},
        "6m17": {"A": (605, 0), "B": (748, 0), "C": (605, 0),
                 "D": (748, 0), "E": (183, 0), "F": (183, 0)},
    }
    raw = tmp_path / "raw-lengths"
    run(raw, [data_files[f"{pdbid}.cif"] for pdbid in expected])
    store = tmp_path / "length-store"
    result = command(tmp_path / "build-lengths", [
        "-m", "pdb2net", "precompute", "--input-dir", str(raw / "inputs"),
        "--store", str(store), "--config", str(raw / "config.json"), "--headless"])
    assert result.returncode == 0, result.stdout + result.stderr
    for pdbid, lengths in expected.items():
        chains = precomputed.load_entry(store, pdbid)["geometry"]["structure"]["atom_data"]
        assert {c["chain_id"]: (c["aa_len"], c["nt_len"]) for c in chains} == lengths
    before = {str(p.relative_to(store)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in store.rglob("*") if p.is_file()}
    shutil.rmtree(raw / "inputs")
    if os.name == "posix":
        for path in store.rglob("*"):
            path.chmod(0o555 if path.is_dir() else 0o444)
        store.chmod(0o555)
    assembled = tmp_path / "assembled-lengths"
    try:
        args = ["-m", "pdb2net", "assemble", "--store", str(store)]
        for pdbid in expected:
            args.extend(["--pdb-id", pdbid])
        args.extend(["--output-dir", str(assembled / "work"), "--web-output-dir", str(assembled / "web"),
                     "--config", str(raw / "config.json"), "--headless"])
        result = command(assembled, args)
        assert result.returncode == 0, result.stdout + result.stderr
        assert before == {str(p.relative_to(store)): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in store.rglob("*") if p.is_file()}
    finally:
        if os.name == "posix":
            store.chmod(0o755)
            for path in store.rglob("*"):
                path.chmod(0o755 if path.is_dir() else 0o644)
    assert network_set(raw) == network_set(assembled)
    for case in (raw, assembled):
        for pdbid, lengths in expected.items():
            nodes, _ = network(case / f"web/networks/Chain_Interaction_Network_{pdbid.upper()}.cx2")
            assert {n["name"].split(":")[-1]: (n["aa_len"], n["nt_len"]) for n in nodes} == lengths


@pytest.mark.parametrize(("filename", "case"), [("4hhb.pdb", "4hhb_competitor"), ("9rat.pdb", "9rat_without_target")])
def test_controlled_search_competitor_and_miss(filename, case, tmp_path, data_files, reference_expected):
    pdbid = filename.split(".")[0]
    pdb_fasta = tmp_path / "pdb.fasta"
    pdb_fasta.write_text("".join(f">{header}\n{sequence}\n" for header, sequence in fasta_records(DATASET / "references/pdb_seqres.txt")
                               if not header.lower().startswith(pdbid + "_")))
    sifts = tmp_path / "sifts.tsv"
    sifts.write_text("".join(line for line in (DATASET / "references/sifts.tsv").read_text().splitlines(True)
                            if not line.lower().startswith(pdbid + "\t")))
    overrides = {"pdb_fasta_path": str(pdb_fasta), "sifts_tsv_path": str(sifts),
                 "reference_manifest_id": f"core-reference-v1-{case}"}
    if case == "9rat_without_target":
        fasta = tmp_path / "without-target.fasta"
        fasta.write_text("".join(f">{header}\n{sequence}\n" for header, sequence in fasta_records(DATASET / "references/swissprot.fasta")
                                if header.split("|")[1] != "P61823"))
        db = tmp_path / "database"
        db.mkdir()
        subprocess.run([shutil.which("makeblastdb"), "-in", str(fasta), "-dbtype", "prot", "-out", str(db / "uniprot_db")],
                       capture_output=True, text=True, check=True, timeout=30)
        overrides.update(uniprot_fasta_path=str(fasta), blast_db_path=str(db))
    else:
        query = tmp_path / "alpha.fasta"
        sequence = reference_expected["structures"][filename]["chains"]["A"]["protein_sequence"]
        query.write_text(f">alpha\n{sequence}\n")
        result = subprocess.run([shutil.which("blastp"), "-query", str(query),
                                 "-db", str(Path(config_loader.config["blast_db_path"]) / "uniprot_db"),
                                 "-outfmt", "6 sseqid pident qcovs bitscore evalue", "-evalue", "1e-10"],
                                capture_output=True, text=True, check=True, timeout=30)
        candidates = {fields[0].split("|")[1] for line in result.stdout.splitlines()
                      if (fields := line.split()) and float(fields[1]) >= 30 and float(fields[2]) >= 70
                      and float(fields[3]) >= 70 and float(fields[4]) <= 1e-10}
        assert {"P69905", "P01942"} <= candidates
    summary = run(tmp_path / "search", [data_files[filename]], overrides)
    nodes, _ = chain_export(tmp_path / "search", summary)
    assert {n["name"].split(":")[-1]: n["uniprot_id"] for n in nodes} == reference_expected["search"][case]
    if case == "4hhb_competitor":
        assert summary["annotations"]["by_source"] == {"blastp_swissprot": 4}


def test_real_precompute_assemble_equivalence_and_readonly_store(tmp_path, data_files, reference_expected):
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
    for case in (raw, assembled):
        assert_combined_exports(case, ["6m17.cif", "6w41.cif", NEXTGEN_NAME], reference_expected)
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
