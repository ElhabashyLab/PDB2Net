"""One export must preserve raw scientific results when only some inputs are reused."""
import copy
import gzip
import json
from pathlib import Path

import pytest

from pdb2net import pipeline
from pdb2net.config_loader import activate_config, config
from pdb2net.input_contract import InputValidationError
from pdb2net.structure_analysis import analyze_mixed, prepare_structure, reusable_upload
from scripts.compare_cx2_outputs import compare_directories

DATA = Path(__file__).parent / "fixtures/core_reference/structures"


@pytest.fixture
def settings(tmp_path, integration_config):
    previous = copy.deepcopy(dict(config))
    value = copy.deepcopy(integration_config)
    value.update(output_path=str(tmp_path / "results"), export_detailed_interactions=True)
    value["workers"]["parsing"] = 1
    activate_config(value)
    yield value
    activate_config(previous)


@pytest.mark.integration
@pytest.mark.parametrize("replacement", [False, True])
def test_mixed_outputs_equal_raw_and_reuse_has_no_coordinate_parsing(tmp_path, settings, monkeypatch, replacement):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    for name in ["1bna.cif", "4hhb.pdb"]:
        (inputs / name).write_bytes((DATA / name).read_bytes())
    raw = pipeline.run_pipeline(str(inputs))
    rows = []
    facts = prepare_structure(inputs / "1bna.cif", contact_sink=rows.append)
    unchanged = copy.deepcopy(facts)
    assert len(rows) == 652
    mapping = {"pdb_00001bna": str(inputs / "1bna.cif")} if replacement else None
    if not replacement:
        (inputs / "1bna.cif").unlink()
    actual_parse = pipeline._parse_input_files
    def parse(paths, **kwargs):
        assert all(not name.endswith('1bna.cif') for name in paths)
        return actual_parse(paths, **kwargs)
    monkeypatch.setattr(pipeline, '_parse_input_files', parse)
    mixed = analyze_mixed(inputs, [facts], settings, {"pdb_00001bna": lambda: iter(rows)},
                          upload_sources=mapping, web_output_dir=str(tmp_path / "web"))
    comparison = compare_directories(Path(mixed.run_output_path), Path(raw.run_output_path), 1.)
    assert comparison['status'] == 'PASS', comparison
    assert comparison['files'] and not comparison['missing_expected_files'] and not comparison['extra_actual_files']
    assert len(list(Path(mixed.distances_dir).glob('*.csv'))) == 2
    assert facts == unchanged
    assert {p.name: p.read_bytes() for p in Path(raw.distances_dir).glob('*.csv')} == {
        p.name: p.read_bytes() for p in Path(mixed.distances_dir).glob('*.csv')}
    summary = json.loads((tmp_path / 'web/summary.json').read_text())
    assert summary['counts']['structures'] == summary['resources']['input']['files'] == 2
    assert [entry['pdb_id'] for entry in summary['prepared_inputs']['entries']] == ['pdb_00001bna']
    assert summary['output_contract_version'] == '3.1'


@pytest.mark.integration
def test_upload_match_requires_exact_content_filename_format_and_profile(tmp_path, settings):
    facts = prepare_structure(DATA / '1bna.cif')
    source = tmp_path / '1bna.cif.gz'
    source.write_bytes(gzip.compress((DATA / '1bna.cif').read_bytes()))
    assert reusable_upload(source, facts, settings)
    renamed = tmp_path / 'custom.cif.gz'
    renamed.write_bytes(source.read_bytes())
    assert not reusable_upload(renamed, facts, settings)
    wrong_format = tmp_path / '1bna.pdb'
    wrong_format.write_bytes((DATA / '1bna.cif').read_bytes())
    assert not reusable_upload(wrong_format, facts, settings)
    source.write_bytes(gzip.compress((DATA / '1bna.cif').read_bytes() + b'\n# changed\n'))
    assert not reusable_upload(source, facts, settings)
    source.write_bytes(gzip.compress((DATA / '1bna.cif').read_bytes()))
    custom = copy.deepcopy(settings)
    custom['distance_thresholds']['all_atoms_radius'] = 4.
    assert not reusable_upload(source, facts, custom)


@pytest.mark.integration
def test_mixed_rejects_duplicate_identity_and_combined_count_limit(tmp_path, settings):
    inputs = tmp_path / 'inputs'
    inputs.mkdir()
    (inputs / '1bna.cif').write_bytes((DATA / '1bna.cif').read_bytes())
    facts = prepare_structure(DATA / '1bna.cif')
    with pytest.raises(InputValidationError, match='DUPLICATE_STRUCTURE_IDENTITY'):
        analyze_mixed(inputs, [facts], settings)
    (inputs / '1bna.cif').unlink()
    (inputs / '4hhb.pdb').write_bytes((DATA / '4hhb.pdb').read_bytes())
    settings['resource_limits']['max_input_files'] = 1
    with pytest.raises(InputValidationError, match='INPUT_FILE_COUNT_LIMIT_EXCEEDED'):
        analyze_mixed(inputs, [facts], settings)


@pytest.mark.integration
def test_upload_changed_after_preflight_cannot_be_replaced(tmp_path, settings, monkeypatch):
    inputs = tmp_path / 'inputs'
    inputs.mkdir()
    source = inputs / '1bna.cif'
    changed = (DATA / '1bna.cif').read_bytes() + b'\n# changed\n'
    source.write_bytes(changed)
    facts = prepare_structure(source)
    source.write_bytes((DATA / '1bna.cif').read_bytes())
    preflight = pipeline._preflight_structure_identities
    def mutate(*args, **kwargs):
        identities = preflight(*args, **kwargs)
        source.write_bytes(changed)
        return identities
    monkeypatch.setattr(pipeline, '_preflight_structure_identities', mutate)
    with pytest.raises(InputValidationError, match='INPUT_CHANGED_DURING_PROCESSING'):
        analyze_mixed(inputs, [facts], settings, upload_sources={'pdb_00001bna': str(source)})


@pytest.mark.integration
@pytest.mark.parametrize('reason', ['radius', 'renamed'])
def test_incompatible_upload_candidate_retains_raw_outputs_and_diagnostics(tmp_path, settings, reason):
    inputs = tmp_path / 'inputs'
    inputs.mkdir()
    source = inputs / ('2abc.cif' if reason == 'renamed' else '1bna.cif')
    source.write_bytes((DATA / '1bna.cif').read_bytes())
    facts = prepare_structure(DATA / '1bna.cif')
    if reason == 'radius':
        settings['distance_thresholds']['all_atoms_radius'] = 4.
    activate_config(settings)
    raw = pipeline.run_pipeline(str(inputs))
    mixed = analyze_mixed(inputs, [facts], settings, upload_sources={'pdb_00001bna': str(source)},
                          web_output_dir=str(tmp_path / 'web'))
    comparison = compare_directories(Path(mixed.run_output_path), Path(raw.run_output_path), 1.)
    assert comparison['status'] == 'PASS' and comparison['files']
    report = json.loads((tmp_path / 'web/summary.json').read_text())
    assert report['prepared_inputs'] is None
    assert report['input_files'] == [source.name]
    if reason == 'renamed':
        assert any(w['code'] == 'STRUCTURE_FILENAME_ID_MISMATCH' for w in report['warnings'])


@pytest.mark.integration
@pytest.mark.parametrize('detailed', [False, True])
def test_all_reused_network_only_and_csv_only(tmp_path, settings, monkeypatch, detailed):
    inputs = tmp_path / 'inputs'
    inputs.mkdir()
    source = inputs / '1bna.cif'
    source.write_bytes((DATA / '1bna.cif').read_bytes())
    rows = []
    facts = prepare_structure(source, contact_sink=rows.append)
    settings['export_detailed_interactions'] = detailed
    if detailed:
        settings['networks'] = dict.fromkeys(settings['networks'], False)
    monkeypatch.setattr(pipeline, '_parse_input_files', lambda *_a, **_kw: pytest.fail('Reused upload was parsed'))
    streams = {'pdb_00001bna': lambda: iter(rows)} if detailed else {'pdb_00001bna': lambda: pytest.fail('Network-only read contacts')}
    result = analyze_mixed(inputs, [facts], settings, streams,
                          upload_sources={'pdb_00001bna': str(source)})
    report = json.loads(Path(result.summary_file).read_text())
    assert report['status'] == 'success' and report['counts']['structures'] == 1
    assert bool(list(Path(result.distances_dir).glob('*.csv'))) == detailed
    assert bool(list(Path(result.run_output_path).rglob('*.cx2'))) != detailed


@pytest.mark.integration
@pytest.mark.parametrize('failure', ['budget', 'corrupt'])
def test_prepared_failure_after_raw_export_publishes_only_diagnostics(tmp_path, settings, failure):
    inputs = tmp_path / 'inputs'
    inputs.mkdir()
    (inputs / '4hhb.pdb').write_bytes((DATA / '4hhb.pdb').read_bytes())
    rows = []
    facts = prepare_structure(DATA / '1bna.cif', contact_sink=rows.append)
    if failure == 'budget':
        settings['resource_limits']['max_detailed_interaction_rows'] = 2000
    else:
        rows = rows[:-1]
    with pytest.raises(InputValidationError):
        analyze_mixed(inputs, [facts], settings, {'pdb_00001bna': lambda: iter(rows)},
                      web_output_dir=str(tmp_path / 'web'))
    report = json.loads((tmp_path / 'web/summary.json').read_text())
    assert report['status'] == 'failed' and report['errors']
    assert report['errors'][0]['code'] == ('DETAILED_INTERACTION_ROW_LIMIT_EXCEEDED' if failure == 'budget' else 'INVALID_STRUCTURE_FACTS')
    assert any('4HHB' in path.name for path in (tmp_path / 'results').rglob('*.csv'))
    assert report['networks'] == report['interactions'] == []
    assert not list((tmp_path / 'web').rglob('*.csv'))
    assert not list((tmp_path / 'web').rglob('*.cx2'))


@pytest.mark.integration
def test_raw_reused_and_explicit_inputs_share_one_combined_network_run(tmp_path, settings):
    inputs = tmp_path / 'inputs'
    inputs.mkdir()
    names = ['1bna.cif', '4hhb.pdb', 'AF-P69905-F1-model_v6.cif']
    for name in names:
        (inputs / name).write_bytes((DATA / name).read_bytes())
    raw = pipeline.run_pipeline(str(inputs))
    facts, streams = [], {}
    for name in names[:2]:
        rows = []
        payload = prepare_structure(inputs / name, contact_sink=rows.append)
        facts.append(payload)
        streams[payload['structure_identity']['canonical_id']] = lambda rows=rows: iter(rows)
    (inputs / '4hhb.pdb').unlink()
    mixed = analyze_mixed(inputs, facts, settings, streams,
                          upload_sources={'pdb_00001bna': str(inputs / '1bna.cif')},
                          web_output_dir=str(tmp_path / 'web'))
    comparison = compare_directories(Path(mixed.run_output_path), Path(raw.run_output_path), 1.)
    assert comparison['status'] == 'PASS' and comparison['files']
    report = json.loads((tmp_path / 'web/summary.json').read_text())
    assert report['counts']['structures'] == report['resources']['input']['files'] == 3
    assert len(report['prepared_inputs']['entries']) == 2
