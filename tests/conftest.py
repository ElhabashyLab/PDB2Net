"""Frozen real-file tests use repository data and lazy, isolated tool setup."""

import hashlib
import json
from pathlib import Path
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
    "1bna.cif",
    "1ehz.cif",
    "2erk.cif",
    "1ao7.cif",
)
REFERENCE_FILES = ("pdb_seqres.txt", "sifts.tsv", "swissprot.fasta")
REFERENCE_DATASET = Path(__file__).parent / "fixtures/core_reference"


def pytest_addoption(parser):
    group = parser.getgroup("PDB2Net integration tests")
    group.addoption(
        "--run-integration",
        action="store_true",
        help="Explicitly select the default real-file setup; required for external path overrides.",
    )
    group.addoption(
        "--integration-data", help="Directory containing the fixed structure fixtures."
    )
    group.addoption(
        "--integration-references",
        help="Directory containing the small reference FASTA/TSV fixtures.",
    )


def _readable_file(path, maximum_bytes, expected):
    try:
        if not path.is_file() or not 0 < path.stat().st_size <= maximum_bytes:
            raise ValueError(
                f"expected a nonempty file of at most {maximum_bytes} bytes"
            )
        with path.open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        if path.stat().st_size != expected["bytes"] or digest != expected["sha256"]:
            raise ValueError("Frozen fixture checksum mismatch; use the pinned dataset bytes")
    except (OSError, ValueError) as exc:
        raise pytest.UsageError(
            f"Integration fixture unavailable: {path}: {exc}"
        ) from exc
    return path


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "integration: real structures with frozen repository references and local search tools",
    )
    enabled = config.getoption("--run-integration")
    data = config.getoption("--integration-data")
    references = config.getoption("--integration-references")
    if (data or references) and not enabled:
        raise pytest.UsageError("External fixture overrides require --run-integration.")
    if bool(data) != bool(references):
        raise pytest.UsageError(
            "External overrides require both --integration-data and --integration-references."
        )


@pytest.fixture(scope="session")
def integration_settings(pytestconfig):
    data = pytestconfig.getoption("--integration-data") or REFERENCE_DATASET / "structures"
    references = pytestconfig.getoption("--integration-references") or REFERENCE_DATASET / "references"
    data, references = (
        Path(data).expanduser().resolve(),
        Path(references).expanduser().resolve(),
    )
    if not data.is_dir() or not references.is_dir():
        raise pytest.UsageError(
            "--integration-data and --integration-references must be existing directories."
        )
    files = {}
    provenance = json.loads((REFERENCE_DATASET / "provenance.json").read_text(encoding="utf-8"))
    assert set(provenance["structures"]) == set(STRUCTURE_FILES)
    assert set(provenance["references"]) == set(REFERENCE_FILES)
    for name in STRUCTURE_FILES:
        path = data / name
        if not path.exists():
            path = data / "review-ergaenzungen" / name
        files[name] = _readable_file(path, 12_000_000, provenance["structures"][name])
    for name in REFERENCE_FILES:
        _readable_file(references / name, 5_000_000, provenance["references"][name])
    return {"files": files, "references": references}


@pytest.fixture(scope="session")
def data_files(integration_settings):
    return integration_settings["files"]


@pytest.fixture(scope="session")
def integration_config(integration_settings, tmp_path_factory):
    references = integration_settings["references"]
    tools = {name: shutil.which(name) for name in ("blastp", "makeblastdb")}
    missing = [name for name, path in tools.items() if path is None]
    if missing:
        pytest.fail(
            "Integration tools missing from PATH: " + ", ".join(missing)
        )
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
    cfg = json.loads(
        (REFERENCE_DATASET / "config.json").read_text(encoding="utf-8")
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
    return cfg
