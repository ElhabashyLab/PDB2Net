"""Opt-in real-file integration tests; normal runs need no external dataset."""

import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess

import pytest


STRUCTURE_FILES = (
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
    "123456.cif",
    "1a1t.pdb",
    "1a1t.cif",
    "1a8o.pdb",
    "1a8o.cif",
    "8aly.cif",
    "pdb_00008aly_xyz-enrich.cif.gz",
    "AF-P69905-F1-model_v6.pdb",
    "AF-P69905-F1-model_v6.cif",
)
REFERENCE_FILES = ("pdb_seqres.txt", "sifts.tsv", "swissprot.fasta")


def pytest_addoption(parser):
    group = parser.getgroup("PDB2Net integration tests")
    group.addoption(
        "--run-integration",
        action="store_true",
        help="Enable real-file integration tests.",
    )
    group.addoption(
        "--integration-data", help="Directory containing the fixed structure fixtures."
    )
    group.addoption(
        "--integration-references",
        help="Directory containing the small reference FASTA/TSV fixtures.",
    )


def _readable_file(path, maximum_bytes):
    try:
        if not path.is_file() or not 0 < path.stat().st_size <= maximum_bytes:
            raise ValueError(
                f"expected a nonempty file of at most {maximum_bytes} bytes"
            )
        with path.open("rb") as stream:
            stream.read(1)
    except (OSError, ValueError) as exc:
        raise pytest.UsageError(
            f"Integration fixture unavailable: {path}: {exc}"
        ) from exc
    return path


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "integration: opt-in checks using external real structures and references",
    )
    enabled = config.getoption("--run-integration")
    data = config.getoption("--integration-data")
    references = config.getoption("--integration-references")
    suite = Path(__file__).resolve().parent / "integration"
    selected_directly = any(
        Path(str(arg).split("::", 1)[0]).resolve().is_relative_to(suite)
        for arg in config.args
    )
    if not enabled:
        if data or references or selected_directly:
            raise pytest.UsageError(
                "Real-file integration tests require --run-integration."
            )
        return
    if not data or not references:
        raise pytest.UsageError(
            "--run-integration requires --integration-data and --integration-references."
        )
    data, references = (
        Path(data).expanduser().resolve(),
        Path(references).expanduser().resolve(),
    )
    if not data.is_dir() or not references.is_dir():
        raise pytest.UsageError(
            "--integration-data and --integration-references must be existing directories."
        )
    files = {}
    for name in STRUCTURE_FILES:
        path = data / name
        if not path.exists():
            path = data / "review-ergaenzungen" / name
        files[name] = _readable_file(path, 12_000_000)
    for name in REFERENCE_FILES:
        _readable_file(references / name, 5_000_000)
    tools = {name: shutil.which(name) for name in ("blastp", "makeblastdb")}
    missing = [name for name, path in tools.items() if path is None]
    if missing:
        raise pytest.UsageError(
            "Integration tools missing from PATH: " + ", ".join(missing)
        )
    config._pdb2net_integration = {
        "files": files,
        "references": references,
        "tools": tools,
    }


@pytest.hookimpl(trylast=True)
def pytest_collection_modifyitems(config, items):
    if not config.getoption("--run-integration"):
        # Pytest has already evaluated -m, including parentheses and exclusions.
        if re.search(r"\bintegration\b", config.getoption("markexpr")) and any(
            item.get_closest_marker("integration") for item in items
        ):
            raise pytest.UsageError(
                "Real-file integration tests require --run-integration."
            )
        skip = pytest.mark.skip(
            reason="real-file integration tests require --run-integration and external fixtures"
        )
        for item in items:
            if item.get_closest_marker("integration"):
                item.add_marker(skip)


@pytest.fixture(scope="session")
def data_files(pytestconfig):
    return pytestconfig._pdb2net_integration["files"]


@pytest.fixture(scope="session")
def integration_config(pytestconfig, tmp_path_factory):
    settings = pytestconfig._pdb2net_integration
    references, tools = settings["references"], settings["tools"]
    work = tmp_path_factory.mktemp("integration-reference")
    database = work / "blast_db"
    database.mkdir()
    result = subprocess.run(
        [
            tools["makeblastdb"],
            "-in",
            str(references / "swissprot.fasta"),
            "-dbtype",
            "prot",
            "-out",
            str(database / "uniprot_db"),
        ],
        capture_output=True,
        text=True,
        timeout=60,
    )
    (work / "makeblastdb.log").write_text(
        result.stdout + result.stderr, encoding="utf-8"
    )
    if result.returncode:
        pytest.fail(
            "Cannot build the temporary integration BLAST database: " + result.stderr
        )
    hashes = {}
    for name in REFERENCE_FILES:
        with (references / name).open("rb") as stream:
            hashes[name] = hashlib.file_digest(stream, "sha256").hexdigest()
    versions = {}
    for name, executable in tools.items():
        probe = subprocess.run(
            [executable, "-version"],
            capture_output=True,
            text=True,
            timeout=10,
            check=True,
        )
        versions[name] = probe.stdout.strip()
    identity = hashlib.sha256(
        json.dumps([hashes, versions], sort_keys=True).encode()
    ).hexdigest()
    (work / "reference-identity.json").write_text(
        json.dumps(
            {"sha256": hashes, "tools": versions, "reference_manifest_id": identity},
            indent=2,
        ),
        encoding="utf-8",
    )
    repository = Path(__file__).resolve().parents[1]
    cfg = json.loads(
        (repository / "pdb2net/configs/config.base.json").read_text(encoding="utf-8")
    )
    cfg.update(
        {
            "input_folder_path": str(work / "unused-inputs"),
            "output_path": str(work / "unused-output"),
            "pdb_fasta_path": str(references / "pdb_seqres.txt"),
            "sifts_tsv_path": str(references / "sifts.tsv"),
            "uniprot_fasta_path": str(references / "swissprot.fasta"),
            "blast_db_path": str(database),
            "blastp_executable": tools["blastp"],
            "blast_cache_path": str(work / "cache.sqlite3"),
            "reference_manifest_id": identity,
            "open_in_cytoscape": False,
            "workers": {"parsing": 2, "blast_threads": 1},
        }
    )
    cfg["diamond"]["enabled"] = False
    cfg["resource_limits"].update(
        {
            "max_input_files": 8,
            "max_single_input_bytes": 12_000_000,
            "max_total_input_bytes": 30_000_000,
            "max_single_input_expanded_bytes": 12_000_000,
            "max_total_input_expanded_bytes": 30_000_000,
            "max_processing_batch_bytes": 12_000_000,
            "max_detailed_interaction_rows": 200_000,
            "max_detailed_interaction_bytes": 30_000_000,
        }
    )
    return cfg
