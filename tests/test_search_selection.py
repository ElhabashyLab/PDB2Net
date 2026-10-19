"""Independent boundary cases for Swiss-Prot selection and its format retry."""

import copy
import subprocess
from pathlib import Path

import pytest

from pdb2net import config_loader, pipeline, uniprot_matcher as matcher
from pdb2net.input_contract import InputValidationError
from pdb2net.structure_facts import profile_id


def row(*, length=100, coverage=1.0, pident=100.0, evalue=1e-80, bitscore=250.0, qcovs=True):
    columns = ["query", "sp|P12345|TEST", str(pident), str(round(length * coverage)),
               "1", str(round(length * coverage)), str(length)]
    if qcovs:
        columns.append(str(coverage * 100))
    return "\t".join([*columns, str(evalue), str(bitscore), "Test protein"])


@pytest.mark.parametrize("coverage", [1.0, 0.8, 0.5])
@pytest.mark.parametrize("format", ["qcovs", "span", "invalid_qcovs"])
def test_alignment_coverage_has_the_same_meaning_in_each_supported_format(coverage, format):
    record = row(coverage=coverage, qcovs=format != "span")
    if format == "invalid_qcovs":
        columns = record.split("\t")
        columns[7] = "unavailable"
        record = "\t".join(columns)
    hit = matcher._select_blastp_swissprot_hit("A" * 100, [record], use_qcovs=format != "span")
    if coverage < 0.65:  # Strong hits get the established 0.05 coverage relaxation.
        assert hit is None
    else:
        assert hit is not None and hit.qcov == coverage


@pytest.mark.parametrize(("length", "limits"), [
    (40, (1e-6, .85, 35., 60.)),
    (100, (1e-10, .70, 30., 70.)),
    (250, (1e-20, .60, 25., 80.)),
])
@pytest.mark.parametrize("metric", ["evalue", "coverage", "pident", "bitscore"])
def test_each_quality_threshold_rejects_only_the_wrong_side(length, limits, metric):
    baseline = dict(zip(("evalue", "coverage", "pident", "bitscore"), limits))
    boundary = baseline[metric]
    step = boundary * 1e-6 if metric == "evalue" else .001
    for direction in (-1, 0, 1):
        values = {**baseline, metric: boundary + direction * step}
        hit = matcher._select_blastp_swissprot_hit(
            "A" * length, [row(length=length, **values)], use_qcovs=True,
        )
        expected = direction <= 0 if metric == "evalue" else direction >= 0
        assert (hit is not None) == expected, (length, metric, direction)


@pytest.mark.parametrize(("length", "accepted"), [(79, False), (80, True), (199, False), (200, True)])
def test_query_length_boundaries_select_the_declared_quality_limits(length, accepted):
    values = ({"evalue": 1e-10, "coverage": .7, "pident": 30., "bitscore": 70.}
              if length < 100 else {"evalue": 1e-20, "coverage": .6, "pident": 25., "bitscore": 80.})
    hit = matcher._select_blastp_swissprot_hit("A" * length, [row(length=length, **values)], use_qcovs=True)
    assert (hit is not None) == accepted


@pytest.mark.parametrize(("coverage", "pident", "bitscore", "accepted"), [
    (.65, 35., 110., True), (.649, 35., 110., False),
    (.65, 34.999, 110., False), (.65, 35., 109.999, False),
])
def test_coverage_relaxation_requires_both_stronger_identity_and_bitscore(coverage, pident, bitscore, accepted):
    hit = matcher._select_blastp_swissprot_hit(
        "A" * 100, [row(coverage=coverage, pident=pident, bitscore=bitscore)], use_qcovs=True,
    )
    assert (hit is not None) == accepted


def test_backend_format_retry_preserves_a_full_alignment(tmp_path, monkeypatch):
    formats = []

    def run(_self, command, *, tool):
        assert tool == "BLASTP"
        format = command[command.index("-outfmt") + 1]
        formats.append(format)
        if "qcovs" in format:
            return subprocess.CompletedProcess(command, 1, "", "qcovs unsupported")
        Path(command[command.index("-out") + 1]).write_text(row(length=40, qcovs=False) + "\n")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(matcher.SearchProcesses, "run", run)
    result = matcher._run_blastp_swissprot_batch([("candidate", "A" * 40, "candidate")])["candidate"]
    assert len(formats) == 2 and "qcovs" in formats[0] and "qcovs" not in formats[1]
    assert result.status == "hit" and result.hit.accession == "P12345" and result.hit.qcov == 1.0


@pytest.mark.parametrize(("known", "maximum", "eligible"), [(3, .2, False), (4, .2, True), (3, .5, True)])
def test_preflight_and_search_agree_at_the_configured_x_limit(tmp_path, monkeypatch, known, maximum, eligible):
    monkeypatch.setitem(matcher.config, "blast", {"max_x_fraction": maximum})
    monkeypatch.setitem(matcher.config, "diamond", {"enabled": False})
    monkeypatch.setitem(matcher.config, "blastp_executable", str(tmp_path / "missing"))
    monkeypatch.setattr(matcher, "_CACHE_ENABLED", False)
    queries = []

    def search(batch, **_kwargs):
        queries.extend(batch)
        return {key: matcher._no_hit_outcome() for key, _sequence, _label in batch}

    monkeypatch.setattr(matcher, "_run_blastp_swissprot_batch", search)
    data = [{"pdb_id": "LOCAL", "file_path": "local.cif", "atom_data": [{
        "chain_id": "A", "molecule_type": "Protein", "uniprot_id": None,
        "residues": [{"residue_name": name} for name in ["ALA"] * known + ["UNK"] * (5 - known)],
    }]}]
    assert pipeline._blast_fallback_needed(data) is eligible
    if eligible:
        with pytest.raises(InputValidationError, match="BLASTP_NOT_FOUND"):
            pipeline._validate_blast_ready_if_needed(data)
    else:
        pipeline._validate_blast_ready_if_needed(data)
    matcher.parallel_blast_search(data, max_workers=1)
    assert [sequence for _key, sequence, _label in queries] == (["A" * known + "X" * (5 - known)] if eligible else [])


def test_search_policy_invalidates_negative_cache_without_invalidating_geometry(tmp_path, monkeypatch):
    original_config = copy.deepcopy(config_loader.get_config())
    cfg = copy.deepcopy(original_config)
    cfg.update(blast_cache_path=str(tmp_path / "cache.sqlite3"), blast_db_path=str(tmp_path),
               uniprot_fasta_path=str(tmp_path / "uniprot.fasta"), reference_manifest_id="fixed-generation",
               structure_model_policy="first", diamond={"enabled": False},
               network_annotations={"use_embedded_sifts": True})
    monkeypatch.setattr(matcher, "_CACHE_ENABLED", True)
    policy = matcher._search_policy_payload

    def old_policy():
        value = policy()
        value["swissprot"]["selection_policy"]["version"] = 2
        return value

    try:
        monkeypatch.setattr(matcher, "_search_policy_payload", old_policy)
        config_loader.activate_config(cfg)
        old_profile = profile_id()
        key = matcher._seq_cache_key("A" * 40)
        matcher._cache_put(key, None)
        assert matcher._cache_get(key) == (True, None)

        monkeypatch.setattr(matcher, "_search_policy_payload", policy)
        config_loader.activate_config(cfg)
        assert matcher._cache_get(key) == (False, None)
        assert profile_id() == old_profile
    finally:
        if matcher._CACHE_CONN is not None:
            matcher._CACHE_CONN.close()
        matcher._CACHE_CONN = None
        matcher._CACHE_DB_SIG = None
        config_loader.activate_config(original_config)
