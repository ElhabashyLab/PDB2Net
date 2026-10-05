"""Check the independent comparison tool using explicitly specified graphs."""

import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts/compare_cx2_outputs.py"
SPEC = importlib.util.spec_from_file_location("cx2_comparison", SCRIPT)
comparison = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(comparison)


@pytest.fixture
def network():
    # The expected graph is written independently of the production exporter.
    return [
        {"CXVersion": "2.0", "hasFragments": False},
        {"nodes": [
            {"id": 0, "x": 0.0, "y": 0.0, "v": {
                "name": "1ABC:A", "molecule_type": "Protein",
                "molecule_name": "Protein A", "uniprot_id": "P12345",
                "tooltip": "Type: Protein\nUniProt: P12345",
                "aa_len": 10, "nt_len": 0, "node_kind": "chain",
                "source_chains": "1ABC:A", "source_file": "1abc.pdb",
                "annotation_source": "embedded_sifts",
                "annotation_confidence": "high", "matched_id": "P12345",
                "representative_accession": "P12345",
                "embedded_uniprot_status": "unique_best_mapping",
                "annotation_pfam": '[{"accession":"PF00001","pdb_start":1,"pdb_end":10}]',
            }},
            {"id": 1, "x": 10.0, "y": 0.0, "v": {
                "name": "1ABC:B", "molecule_type": "Protein", "uniprot_id": "Q99999",
            }},
        ]},
        {"edges": [
            {"id": 0, "s": 0, "t": 1, "v": {
                "interaction": "Protein-Protein", "all_atoms_count": 4,
            }},
        ]},
    ]


def _aspects(document):
    return {key: value for block in document for key, value in block.items()}


def _compare(actual, expected):
    return comparison.compare_networks(
        comparison.extract_canonical_network(actual),
        comparison.extract_canonical_network(expected),
    )


@pytest.mark.parametrize(("field", "value"), [
    ("aa_len", 100000), ("nt_len", 7), ("node_kind", "protein"),
    ("uniprot_id", "Q99999"), ("molecule_name", "Wrong protein"),
    ("source_chains", "2XYZ:A"), ("annotation_source", "diamond_uniref90"),
    ("annotation_confidence", "low"), ("matched_id", "UniRef90_Q99999"),
    ("representative_accession", "Q99999"),
    ("embedded_uniprot_status", "ambiguous_multi_mapping"),
    ("annotation_pfam", "WRONG"), ("annotation_pfam_total", 17),
])
def test_changed_scientific_node_attribute_fails(network, field, value):
    actual = copy.deepcopy(network)
    _aspects(actual)["nodes"][0]["v"][field] = value
    result = _compare(actual, network)
    assert result["status"] == "FAIL"
    assert any(item["field"] == field for item in result["node_differences"])


@pytest.mark.parametrize("field", ["aa_len", "uniprot_id", "molecule_type", "annotation_pfam"])
def test_missing_scientific_attribute_fails_even_with_matching_tooltip(network, field):
    actual = copy.deepcopy(network)
    del _aspects(actual)["nodes"][0]["v"][field]
    assert _compare(actual, network)["status"] == "FAIL"


@pytest.mark.parametrize("kind", [
    "node_id", "node_name", "edge_id", "edge_pair", "dangling_endpoint",
    "missing_endpoint", "missing_node_id",
])
def test_invalid_graph_cannot_be_hidden_by_normalization(network, kind):
    actual = copy.deepcopy(network)
    aspects = _aspects(actual)
    if kind.startswith("node_"):
        node = copy.deepcopy(aspects["nodes"][0])
        if kind == "node_name":
            node["id"] = 2
        aspects["nodes"].append(node)
    elif kind.startswith("edge_"):
        edge = copy.deepcopy(aspects["edges"][0])
        if kind == "edge_pair":
            edge["id"] = 2
            edge["s"], edge["t"] = edge["t"], edge["s"]
        else:
            edge["v"]["interaction"] = "identity"
        aspects["edges"].append(edge)
    elif kind == "dangling_endpoint":
        aspects["edges"][0]["t"] = 99
    elif kind == "missing_endpoint":
        del aspects["edges"][0]["t"]
    else:
        del aspects["nodes"][0]["id"]
    # Check both sides: an invalid reference is not a valid goldstandard.
    assert _compare(actual, network)["status"] == "FAIL"
    assert _compare(network, actual)["status"] == "FAIL"


@pytest.mark.parametrize("change", ["count", "missing_count", "missing_edge"])
def test_changed_or_missing_contact_data_fails(network, change):
    actual = copy.deepcopy(network)
    edges = _aspects(actual)["edges"]
    if change == "count":
        edges[0]["v"]["all_atoms_count"] += 1
    elif change == "missing_count":
        del edges[0]["v"]["all_atoms_count"]
    else:
        edges.clear()
    assert _compare(actual, network)["status"] == "FAIL"


def test_export_ids_order_and_source_file_provenance_do_not_change_science(network):
    actual = copy.deepcopy(network)
    aspects = _aspects(actual)
    for node in aspects["nodes"]:
        node["id"] += 20
        node["v"]["source_file"] = "pdb:pdb_00001abc"
        node["v"]["source_files"] = "store-generation"
    aspects["nodes"].reverse()
    edge = aspects["edges"][0]
    edge.update(id=37, s=21, t=20)
    assert _compare(actual[::-1], network)["status"] == "PASS"


def test_layout_only_change_remains_a_warning(network):
    actual = copy.deepcopy(network)
    _aspects(actual)["nodes"][0]["x"] += 1.0
    result = _compare(actual, network)
    assert result["status"] == "WARN"
    assert result["failures"] == []


def test_comparison_cli_writes_failure_report_for_invalid_graph(network, tmp_path):
    expected = tmp_path / "expected"
    actual = tmp_path / "actual"
    expected.mkdir()
    actual.mkdir()
    (expected / "network.cx2").write_text(json.dumps(network))
    invalid = copy.deepcopy(network)
    _aspects(invalid)["nodes"].append(copy.deepcopy(_aspects(invalid)["nodes"][0]))
    (actual / "network.cx2").write_text(json.dumps(invalid))
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "compare", "--actual", str(actual),
         "--expected", str(expected), "--out", str(tmp_path / "report.md"),
         "--json-out", str(tmp_path / "report.json")],
        capture_output=True, text=True, timeout=10,
    )
    assert result.returncode == 1, result.stderr
    assert "Traceback" not in result.stderr
    report = json.loads((tmp_path / "report.json").read_text())
    assert report["status"] == "FAIL"
    assert report["files"][0]["comparison"]["failures"]


@pytest.mark.parametrize(("left", "right"), [
    ("file", "file"), ("missing", "missing"), ("empty", "empty"),
    ("file", "valid"), ("valid", "file"),
    ("missing", "valid"), ("valid", "missing"),
    ("empty", "valid"), ("valid", "empty"),
    ("non_cx2", "non_cx2"), ("cx2_directory", "cx2_directory"),
])
def test_comparison_rejects_invalid_roots_in_api_and_cli(network, tmp_path, left, right):
    paths = []
    for side, kind in [("actual", left), ("expected", right)]:
        path = tmp_path / side
        if kind == "file":
            path.write_text(json.dumps(network))
        elif kind != "missing":
            path.mkdir()
            if kind == "valid":
                (path / "network.cx2").write_text(json.dumps(network))
            elif kind == "non_cx2":
                (path / "notes.txt").write_text("not a network")
            elif kind == "cx2_directory":
                (path / "not-a-file.cx2").mkdir()
        paths.append(path)
    with pytest.raises(ValueError, match="directory|CX2 files"):
        comparison.compare_directories(*paths, overlap_threshold=20.)
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "compare", "--actual", str(paths[0]),
         "--expected", str(paths[1]), "--out", str(tmp_path / "report.md"),
         "--json-out", str(tmp_path / "report.json")],
        capture_output=True, text=True, timeout=10,
    )
    assert result.returncode == 2
    assert "error:" in result.stderr and "Traceback" not in result.stderr
    assert not (tmp_path / "report.md").exists()
    assert not (tmp_path / "report.json").exists()


def test_comparison_recurses_and_compares_real_empty_graph_files(network, tmp_path):
    actual, expected = tmp_path / "actual", tmp_path / "expected"
    for root in [actual, expected]:
        (root / "nested").mkdir(parents=True)
        (root / "nested/network.cx2").write_text(json.dumps(network))
        (root / "empty.cx2").write_text(json.dumps([{"nodes": []}, {"edges": []}]))
    report = comparison.compare_directories(actual, expected, 20.)
    assert report["status"] == "PASS"
    assert len(report["files"]) == 2
    (actual / "nested/network.cx2").rename(actual / "different.cx2")
    report = comparison.compare_directories(actual, expected, 20.)
    assert report["status"] == "FAIL"
    assert len(report["files"]) == 1
    assert report["missing_expected_files"] == ["nested/network.cx2"]


def test_goldstandard_wrapper_cannot_report_success_for_empty_comparison(tmp_path):
    actual, expected = tmp_path / "actual", tmp_path / "expected"
    actual.mkdir()
    expected.mkdir()
    result = subprocess.run(
        [sys.executable, str(SCRIPT.with_name("run_goldstandard_check.py")),
         "--actual", str(actual), "--expected", str(expected), "--report-dir", str(tmp_path / "reports")],
        capture_output=True, text=True, timeout=10,
    )
    assert result.returncode != 0
    assert "Goldstandard status: PASS" not in result.stdout
