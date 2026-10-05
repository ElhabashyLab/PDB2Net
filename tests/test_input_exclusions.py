"""Known bad input files must not silently disappear or invalidate good inputs."""

import copy
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import sys

import gemmi
import pytest

from pdb2net import config_loader, pipeline
from pdb2net.input_contract import InputValidationError


def _pdb(identity="1ABC"):
    return "HEADER".ljust(62) + identity + "\n" + "".join(
        f"ATOM  {i:5d}  P     U {chain}{1:4d}    "
        f"{x:8.3f}{0.:8.3f}{0.:8.3f}  1.00 20.00           P  \n"
        for i, chain, x in [(1, "A", 0.), (2, "B", 3.)]
    ) + "END\n"


def _missing_chain_input(extension, identity="3ABC"):
    # Keep one usable named chain: losing the whole file must be explicit.
    text = _pdb(identity).replace(" U B", " U  ")
    if extension == "cif":
        document = gemmi.read_pdb_string(text).make_mmcif_document()
        document[0].name = identity.lower()
        text = document.as_string()
    return text


@pytest.fixture
def setup_run(tmp_path, monkeypatch):
    original = copy.deepcopy(config_loader.get_config())
    cfg = json.loads((Path(__file__).parents[1] / "pdb2net/configs/config.base.json").read_text())
    (tmp_path / "pdb.fasta").write_text(">1abc_A mol:na length:1 RNA\nU\n>1abc_B mol:na length:1 RNA\nU\n")
    (tmp_path / "sifts.tsv").write_text("# PDB\tCHAIN\tSP_PRIMARY\n")
    (tmp_path / "uniprot.fasta").write_text(">sp|P12345|UNRELATED Other protein\nAG\n")
    cfg.update(
        pdb_fasta_path=str(tmp_path / "pdb.fasta"), sifts_tsv_path=str(tmp_path / "sifts.tsv"),
        uniprot_fasta_path=str(tmp_path / "uniprot.fasta"),
        blast_cache_path=str(tmp_path / "cache.sqlite3"),
        output_path=str(tmp_path / "work"), open_in_cytoscape=False,
        reference_manifest_id="input-exclusions-fixture", workers={"parsing": 1, "blast_threads": 1},
    )
    cfg["networks"].update(protein_per_pdb=False, combined_protein_network=False)
    config_loader.activate_config(cfg)
    # Keep the real bounded scheduler and parser; threads permit call counting.
    monkeypatch.setattr(pipeline, "ProcessPoolExecutor", lambda **_kw: ThreadPoolExecutor(max_workers=1))
    try:
        yield tmp_path
    finally:
        config_loader.activate_config(original)


def _run(root, name, files):
    inputs = root / name
    inputs.mkdir()
    for label, contents in files.items():
        (inputs / label).write_bytes(contents if isinstance(contents, bytes) else contents.encode())
    web = root / (name + "-web")
    paths = pipeline.run_pipeline(str(inputs), web_output_dir=str(web))
    return json.loads((web / "summary.json").read_text()), web, paths


@pytest.mark.parametrize(("name", "contents", "code"), [
    ("empty.pdb", "", "EMPTY_INPUT_FILE"),
    ("no-atoms.pdb", "END\n", "NO_VALID_CHAINS"),
    ("no-atoms.cif", "data_local\n", "NO_VALID_CHAINS"),
    ("water.pdb", _pdb("3ABC").replace("  U ", "HOH "), "NO_VALID_CHAINS"),
    ("hydrogen.pdb", _pdb("3ABC").replace("           P  ", "           H  "), "NO_VALID_CHAINS"),
    ("zero.pdb", _pdb("3ABC").replace("  1.00 20.00", "  0.00 20.00"), "NO_VALID_CHAINS"),
    ("occupancy.pdb", _pdb("3ABC").replace("  1.00 20.00", "   bad 20.00"), "INVALID_ATOM_OCCUPANCY"),
    ("duplicate.pdb", _pdb("3ABC").replace(" U B", " U A"), "AMBIGUOUS_ATOM_IDENTITY"),
    ("missing-chain.pdb", _missing_chain_input("pdb"), "MISSING_CHAIN_ID"),
    ("missing-chain.cif", _missing_chain_input("cif"), "MISSING_CHAIN_ID"),
    ("encoding.cif", b"\xff", "INVALID_STRUCTURE_TEXT_ENCODING"),
    ("broken.cif", "data_broken\nloop_\n_atom_site.id\n_atom_site.type_symbol\n1\n", "INVALID_MMCIF_INPUT"),
    ("blocks.cif", "data_one\n#\ndata_two\n", "INVALID_MMCIF_DATA_BLOCK_COUNT"),
    ("conflict.cif", "data_1abc\n_entry.id 2xyz\n", "CONFLICTING_STRUCTURE_IDENTITY"),
])
def test_bad_input_is_reported_and_good_network_matches_control(setup_run, name, contents, code):
    root = setup_run
    control, control_web, _ = _run(root, "control", {"good.pdb": _pdb()})
    mixed, web, paths = _run(root, "mixed", {name: contents, "good.pdb": _pdb()})
    assert mixed["status"] == "success" and mixed["errors"] == []
    assert mixed["counts"]["structures"] == 1
    assert mixed["resources"]["input"]["files"] == 2
    assert mixed["input_files"] == ["good.pdb"]
    skipped = [w for w in mixed["warnings"] if w["code"] == "INPUT_FILE_SKIPPED"]
    assert len(skipped) == 1
    assert name in skipped[0]["message"] and code in skipped[0]["message"]
    assert str(root) not in skipped[0]["message"]
    assert control["networks"]
    for network in mixed["networks"]:
        assert json.loads((web / network).read_text()) == json.loads((control_web / network).read_text())
    assert mixed["networks"] == control["networks"]
    for report in [paths.manifest_file, paths.summary_file]:
        assert json.loads(Path(report).read_text())["warnings"][0] == skipped[0]


def test_partial_chain_exclusion_reaches_reports_without_becoming_file_exclusion(setup_run):
    text = _pdb().replace("  1.00 20.00", "  0.00 20.00", 1)
    summary, _, paths = _run(setup_run, "partial-chain", {"1abc.pdb": text})
    assert summary["status"] == "success" and summary["counts"]["structures"] == 1
    warnings = summary["warnings"]
    assert len(warnings) == 1 and warnings[0]["code"] == "CHAIN_WITHOUT_CONTACT_ATOMS"
    assert "1abc.pdb" in warnings[0]["message"] and "chain 'A'" in warnings[0]["message"]
    for report in [paths.manifest_file, paths.summary_file]:
        assert json.loads(Path(report).read_text())["warnings"] == warnings


def test_precompute_reports_empty_chain_and_rejects_empty_structure(setup_run):
    from pdb2net import precomputed

    partial = setup_run / "1abc.pdb"
    partial.write_text(_pdb().replace("  1.00 20.00", "  0.00 20.00", 1))
    report = precomputed.precompute_sources(setup_run / "partial-store", [partial])
    assert report["failed"] == 0 and report["written"] == 1
    assert report["warnings"][0]["code"] == "CHAIN_WITHOUT_CONTACT_ATOMS"
    entry = precomputed.load_entry(setup_run / "partial-store", "1abc")
    assert entry["counts"]["chains"] == 1
    assert "input_warnings" not in entry["geometry"]["structure"]  # Schema remains compact.

    empty = setup_run / "3abc.pdb"
    empty.write_text(_pdb("3ABC").replace("  1.00 20.00", "  0.00 20.00"))
    report = precomputed.precompute_sources(setup_run / "empty-store", [empty])
    assert report["written"] == 0 and report["failed"] == 1
    assert report["errors"][0]["code"] == "NO_VALID_CHAINS"
    assert not (setup_run / "empty-store/manifest.json").exists()


@pytest.mark.parametrize("extension", ["pdb", "cif"])
def test_precompute_reports_missing_chain_id_without_publishing_store(setup_run, extension):
    from pdb2net import precomputed

    good = setup_run / "1abc.pdb"
    good.write_text(_pdb())
    bad = setup_run / f"3abc.{extension}"
    bad.write_text(_missing_chain_input(extension))
    store = setup_run / "missing-chain-store"
    report = precomputed.precompute_sources(store, [good, bad])

    assert report["written"] == 1 and report["failed"] == 1
    assert len(report["errors"]) == 1
    assert report["errors"][0]["pdb_id"] == "pdb_00003abc"
    assert report["errors"][0]["code"] == "MISSING_CHAIN_ID"
    assert "manifest" not in report
    assert not (store / "manifest.json").exists()


def test_all_unusable_files_keep_all_diagnostics(setup_run):
    root = setup_run
    files = {"a.cif": "data_empty\n", "b.pdb": "END\n", "c.pdb": b"\xff",
             "d.pdb": _missing_chain_input("pdb"), "e.cif": _missing_chain_input("cif", "4ABC")}
    with pytest.raises(InputValidationError):
        _run(root, "bad", files)
    summary = json.loads((root / "bad-web/summary.json").read_text())
    assert summary["status"] == "failed"
    assert summary["resources"]["input"]["files"] == len(files)
    assert len(summary["warnings"]) == len(files)
    assert all(name in entry["message"] for name, entry in zip(files, summary["warnings"]))
    assert sum("MISSING_CHAIN_ID" in entry["message"] for entry in summary["warnings"]) == 2
    assert summary["networks"] == summary["interactions"] == []


def test_late_bad_file_is_parsed_once_and_does_not_shift_identity(setup_run, monkeypatch):
    root = setup_run
    calls = []
    process = pipeline.process_single_file

    def counted(path, *args, **kwargs):
        calls.append(Path(path).name)
        return process(path, *args, **kwargs)

    monkeypatch.setattr(pipeline, "process_single_file", counted)
    monkeypatch.setattr(pipeline, "create_processing_batches", lambda paths, _inventory: ([p] for p in paths))
    summary, _, _ = _run(root, "late", {"a.pdb": _pdb(), "b.pdb": _missing_chain_input("pdb"), "c.pdb": _pdb("2XYZ")})
    assert calls == ["a.pdb", "b.pdb", "c.pdb"]
    assert summary["input_files"] == ["a.pdb", "c.pdb"]
    assert [identity["display_id"] for identity in summary["identities"]] == ["1ABC", "2XYZ"]
    assert summary["counts"]["structures"] == 2


@pytest.mark.parametrize("error_type", [RuntimeError, ValueError])
def test_unknown_parser_error_is_fatal_and_keeps_prior_diagnostics(setup_run, monkeypatch, error_type):
    root = setup_run
    process = pipeline.process_single_file

    def broken(path, *args, **kwargs):
        if Path(path).name == "b.pdb":
            raise error_type("unexpected parser implementation failure")
        return process(path, *args, **kwargs)

    monkeypatch.setattr(pipeline, "process_single_file", broken)
    monkeypatch.setattr(pipeline, "create_processing_batches", lambda paths, _inventory: ([p] for p in paths))
    with pytest.raises(error_type, match="unexpected parser"):
        _run(root, "fatal", {"a.pdb": "END\n", "b.pdb": _pdb()})
    summary = json.loads((root / "fatal-web/summary.json").read_text())
    assert summary["status"] == "failed"
    assert len(summary["warnings"]) == 1 and "a.pdb" in summary["warnings"][0]["message"]
    assert summary["networks"] == summary["interactions"] == []


@pytest.mark.parametrize("code", ["INPUT_CHANGED_DURING_PROCESSING", "INPUT_FILE_EXPANDED_BYTES_LIMIT_EXCEEDED"])
def test_integrity_and_resource_errors_remain_fatal(setup_run, monkeypatch, code):
    def fail(*_args, **_kwargs):
        raise InputValidationError(code, "fatal test failure")

    monkeypatch.setattr(pipeline, "process_single_file", fail)
    with pytest.raises(InputValidationError) as error:
        _run(setup_run, "fatal", {"good.pdb": _pdb()})
    assert error.value.code == code


def test_omitted_file_still_counts_towards_input_limit(setup_run, monkeypatch):
    monkeypatch.setitem(pipeline.config["resource_limits"], "max_input_files", 1)
    with pytest.raises(InputValidationError) as error:
        _run(setup_run, "limited", {"good.pdb": _pdb(), "empty.pdb": ""})
    assert error.value.code == "INPUT_FILE_COUNT_LIMIT_EXCEEDED"


def test_file_diagnostics_are_not_truncated_by_general_warning_limit():
    from pdb2net.outputs import _public_diagnostics

    ordinary = [{"code": "OTHER", "message": "other"}] * 105
    excluded = [{"code": "INPUT_FILE_SKIPPED", "message": f"file-{i}.pdb: NO_VALID_CHAINS"} for i in range(150)]
    assert _public_diagnostics(ordinary + excluded) == ordinary[:100] + excluded


@pytest.mark.parametrize("include_good", [True, False])
@pytest.mark.parametrize(("contents", "code"), [
    ("END\n", "NO_VALID_CHAINS"),
    (_missing_chain_input("pdb"), "MISSING_CHAIN_ID"),
])
def test_actual_cli_reports_exclusions_and_exit_status(setup_run, include_good, contents, code):
    root = setup_run
    inputs = root / "cli-inputs"
    inputs.mkdir()
    (inputs / "bad.pdb").write_text(contents)
    if include_good:
        (inputs / "good.pdb").write_text(_pdb())
    config_path = root / "config.json"
    config_path.write_text(json.dumps(config_loader.get_config()))
    defaults = root / "defaults"
    defaults.mkdir()
    (defaults / "config.base.json").write_text(json.dumps(config_loader.get_config()))
    env = {k: v for k, v in os.environ.items() if not k.startswith("PDB2NET_")}
    env["PDB2NET_CONFIG_DIR"] = str(defaults)
    web = root / "cli-web"
    result = subprocess.run(
        [sys.executable, "-m", "pdb2net", "run", "--input-dir", str(inputs), "--output-dir", str(root / "cli-work"),
         "--web-output-dir", str(web), "--config", str(config_path), "--headless"],
        env=env, capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == (0 if include_good else 1), result.stderr
    summary = json.loads((web / "summary.json").read_text())
    assert summary["status"] == ("success" if include_good else "failed")
    assert "bad.pdb" in summary["warnings"][0]["message"]
    assert code in summary["warnings"][0]["message"]
    if include_good:
        assert "Completed with exclusions: 1 of 2" in result.stderr
        assert "bad.pdb" in result.stderr


@pytest.mark.parametrize("failure", ["memory", "annotation"])
def test_parser_does_not_disguise_resource_or_annotation_bugs_as_bad_files(setup_run, monkeypatch, failure):
    from pdb2net import file_parser

    error_type = MemoryError if failure == "memory" else RuntimeError

    def fail(*_args):
        raise error_type("unexpected internal failure")

    if failure == "memory":
        monkeypatch.setattr(file_parser.gemmi.cif, "read_string", fail)
    else:
        monkeypatch.setattr(file_parser, "extract_embedded_annotations", fail)
    with pytest.raises(error_type, match="unexpected internal failure"):
        _run(setup_run, "fatal", {"problem.cif": "data_local\n"})


def test_mutation_followed_by_content_error_is_fatal_not_an_exclusion(setup_run, monkeypatch):
    from pdb2net import file_parser

    parse = file_parser.gemmi.read_pdb_string

    def mutate_then_fail(text):
        if "2XYZ" in text:
            (setup_run / "mutating/bad.pdb").write_text("END\n")
            raise RuntimeError("content error after the input changed")
        return parse(text)

    monkeypatch.setattr(file_parser.gemmi, "read_pdb_string", mutate_then_fail)
    with pytest.raises(InputValidationError) as error:
        _run(setup_run, "mutating", {"bad.pdb": _pdb("2XYZ"), "good.pdb": _pdb()})
    assert error.value.code == "INPUT_CHANGED_DURING_PROCESSING"
    summary = json.loads((setup_run / "mutating-web/summary.json").read_text())
    assert summary["status"] == "failed"
    assert summary["warnings"] == []
    assert summary["networks"] == summary["interactions"] == []


def test_user_interrupt_preserves_failed_summary(setup_run, monkeypatch):
    def interrupt(*args, **kwargs):
        raise KeyboardInterrupt
    monkeypatch.setattr(pipeline, "_run_blast_annotation", interrupt)
    with pytest.raises(KeyboardInterrupt):
        _run(setup_run, "interrupted", {"1abc.pdb": _pdb()})
    report = json.loads((setup_run / "interrupted-web/summary.json").read_text())
    assert report["status"] == "failed"
    assert report["errors"][0]["code"] == "RUN_INTERRUPTED"
    assert report["errors"][0]["message"] == "RUN_INTERRUPTED: Run interrupted by user."
    assert report["networks"] == report["interactions"] == []


@pytest.mark.parametrize("limit", [None, 30])
def test_search_limit_reaches_public_runtime_snapshot(setup_run, monkeypatch, limit):
    monkeypatch.setitem(pipeline.config, "sequence_search_timeout_seconds", limit)
    report, _, paths = _run(setup_run, "limit-snapshot", {"1abc.pdb": _pdb()})
    assert report["config"]["sequence_search_timeout_seconds"] == limit
    assert json.loads(Path(paths.summary_file).read_text())["config"]["sequence_search_timeout_seconds"] == limit


@pytest.mark.parametrize("reporting_step", ["write_failed_run_manifest", "collect_web_outputs"])
def test_unavailable_interrupt_report_cannot_replace_keyboard_interrupt(setup_run, monkeypatch, reporting_step):
    def interrupt(*args, **kwargs):
        raise KeyboardInterrupt
    def unavailable(*args, **kwargs):
        raise OSError("simulated report destination failure")
    monkeypatch.setattr(pipeline, "_run_blast_annotation", interrupt)
    monkeypatch.setattr(pipeline, reporting_step, unavailable)
    with pytest.raises(KeyboardInterrupt):
        _run(setup_run, "unreportable", {"1abc.pdb": _pdb()})


@pytest.mark.parametrize("command", ["run", "precompute"])
def test_cli_user_interrupt_returns_130(setup_run, monkeypatch, capsys, command):
    from pdb2net.__main__ import main
    from pdb2net import precomputed

    def interrupt(*args, **kwargs):
        raise KeyboardInterrupt
    if command == "run":
        monkeypatch.setattr(pipeline, "run_pipeline", interrupt)
        args = ["--output-dir", str(setup_run / "out")]
    else:
        monkeypatch.setattr(precomputed, "precompute_directory", interrupt)
        args = ["--store", str(setup_run / "store")]
    assert main([command, "--input-dir", str(setup_run), *args, "--headless"]) == 130
    assert "interrupted by user" in capsys.readouterr().err
