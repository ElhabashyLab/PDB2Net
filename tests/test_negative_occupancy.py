"""Independent boundary and geometric expectations for the 2% acceptance rule."""
import copy
import csv
import json
from pathlib import Path

import pytest

from pdb2net import pipeline
from pdb2net.config_loader import activate_config, config
from pdb2net.contact_atoms import ContactView, OccupancyAudit, contact_positions
from pdb2net.input_contract import InputValidationError
from pdb2net.structure_analysis import analyze_prepared, prepare_structure
from pdb2net.structure_facts import canonical_json, profile_id, validate_facts
from scripts.compare_cx2_outputs import compare_directories


@pytest.fixture(autouse=True)
def configuration():
    before = copy.deepcopy(dict(config))
    activate_config({'structure_model_policy': 'first', 'open_in_cytoscape': False})
    yield
    activate_config(before)


def atom(serial, number, occupancy, *, chain='A', x=0., element='P', component='U', record='ATOM', name='P'):
    token = '      ' if occupancy is None else f'{occupancy:6.2f}' if isinstance(occupancy, (int, float)) else f'{occupancy:>6}'
    return (f'{record:<6}{serial:5d} {name:^4} {component:>3} {chain}{number:4d}    '
            f'{x:8.3f}{0.:8.3f}{0.:8.3f}{token} 20.00          {element:>2}  \n')


def source(tmp_path, models, *, extra='', name='1abc.pdb'):
    text = 'HEADER'.ljust(62) + name[:4].upper() + '\n'
    for model, occupancies in enumerate(models, 1):
        if len(models) > 1:
            text += f'MODEL     {model:4d}\n'
        text += ''.join(atom(i, i, value) for i, value in enumerate(occupancies, 1))
        text += extra + 'TER\n'
        if len(models) > 1:
            text += 'ENDMDL\n'
    path = tmp_path / name
    path.write_text(text + 'END\n')
    return path


@pytest.mark.parametrize('total,accepted', [(49, False), (50, True), (51, True)])
def test_exact_boundary_keeps_source_and_raw_atom_metadata(tmp_path, total, accepted):
    path = source(tmp_path, [[-.5] + [1.] * (total - 1)])
    before = path.read_bytes()
    if not accepted:
        with pytest.raises(InputValidationError, match='exceeding the 2% limit'):
            pipeline.process_single_file(str(path))
    else:
        result = pipeline.process_single_file(str(path))
        chain = result['atom_data'][0]
        assert chain['residues'][0]['atoms'][0]['occupancy'] == -.5
        assert len(ContactView(chain).positions) == total - 1
        warnings = [w for w in result['input_warnings'] if w['code'] == 'NEGATIVE_OCCUPANCY_ATOMS_SKIPPED']
        assert len(warnings) == 1 and 'excluded 1 atom positions' in warnings[0]['message']
    assert path.read_bytes() == before


def test_unknown_occupancy_counts_but_zero_does_not_dilute(tmp_path):
    path = source(tmp_path, [[-.5] + [None] * 49])
    assert len(ContactView(pipeline.process_single_file(str(path))['atom_data'][0]).positions) == 49
    path = source(tmp_path, [[-.5] + [1.] * 48 + [0.] * 200])
    with pytest.raises(InputValidationError, match='1/49'):
        pipeline.process_single_file(str(path))


@pytest.mark.parametrize('element,component,record,name', [('H', 'U', 'ATOM', 'H'), ('O', 'HOH', 'HETATM', 'O')])
def test_hydrogen_and_water_do_not_dilute_denominator(tmp_path, element, component, record, name):
    extra = ''.join(atom(100+i, 100+i, 1., element=element, component=component, record=record, name=name) for i in range(100))
    path = source(tmp_path, [[-.5] + [1.] * 48], extra=extra)
    with pytest.raises(InputValidationError, match='1/49'):
        pipeline.process_single_file(str(path))


def test_unsupported_polymer_and_nonpolymer_positions_do_not_count():
    audit = OccupancyAudit()
    for polymer, kind in [(True, 'Pna'), (True, 'SaccharideD'), (False, 'PeptideL')]:
        residue = {'is_polymer': polymer, 'polymer_type': kind, 'atoms': [{'occupancy': -.5}] * 100}
        assert list(contact_positions({'residues': [residue]}, occupancy_audit=audit)) == []
    assert audit.finish('irrelevant.pdb') == []
    assert audit.models[1] == [0, 0]


def test_model_limit_is_not_averaged_across_models(tmp_path):
    path = source(tmp_path, [[1.] * 1000, [-.5] + [1.] * 48])
    assert not pipeline.process_single_file(str(path))['input_warnings']
    config['structure_model_policy'] = 'all'
    with pytest.raises(InputValidationError, match='model 2 has 1/49'):
        pipeline.process_single_file(str(path))


def test_counts_are_shared_across_chains_of_one_model(tmp_path):
    # Each chain alone would have a different fraction: apply the agreed model boundary.
    path = source(tmp_path, [[-.5] + [1.] * 24], extra='TER\n' + ''.join(atom(100+i, i, 1., chain='B') for i in range(1, 26)))
    result = pipeline.process_single_file(str(path))
    assert [len(ContactView(c).positions) for c in result['atom_data']] == [24, 25]
    assert '2%' in result['input_warnings'][0]['message']


@pytest.mark.parametrize('invalid', ['xyz', float('nan'), float('inf'), float('-inf')])
def test_malformed_or_nonfinite_after_negative_remains_error(tmp_path, invalid):
    path = source(tmp_path, [[-.5] + [1.] * 100 + [invalid]])
    with pytest.raises(InputValidationError, match='INVALID_ATOM_OCCUPANCY'):
        pipeline.process_single_file(str(path))


def test_warning_has_at_most_ten_examples_and_one_message_per_file(tmp_path):
    path = source(tmp_path, [[-.1] * 11 + [1.] * 989, [-.1] + [1.] * 49])
    config['structure_model_policy'] = 'all'
    result = pipeline.process_single_file(str(path))
    warning, = [w for w in result['input_warnings'] if w['code'] == 'NEGATIVE_OCCUPANCY_ATOMS_SKIPPED']
    assert 'excluded 12 atom positions' in warning['message']
    assert 'largest affected-model fraction 2%' in warning['message']
    assert warning['message'].count('atom P') == 10
    assert '2 further excluded positions not listed' in warning['message']
    assert len(warning['message']) < 4000


def test_all_zero_remains_excluded_without_false_negative_warning(tmp_path):
    path = source(tmp_path, [[0.] * 50])
    with pytest.raises(InputValidationError, match='NO_VALID_CHAINS'):
        prepare_structure(path)


def test_profile_stays_equal_to_existing_published_profile():
    assert profile_id() == '4470c9971dfecf118c832a80b62f0dbab3e4d37c1a33a02a9cca9cd3c3c4fc7c'


def test_raw_and_prepared_contacts_csv_and_warning_match_independent_geometry(tmp_path, integration_config):
    cfg = copy.deepcopy(integration_config)
    cfg.update(output_path=str(tmp_path / 'runs'), export_detailed_interactions=True)
    cfg['networks'].update(protein_per_pdb=False, combined_protein_network=False)
    cfg['workers']['parsing'] = 1
    activate_config(cfg)
    inputs = tmp_path / 'inputs'
    inputs.mkdir()
    path = source(inputs, [[-.5] + [1.] * 49], extra='TER\n' + atom(51, 1, 1., chain='B', x=3.))
    rows = []
    facts = prepare_structure(path, contact_sink=rows.append)
    immutable = canonical_json(facts)
    validate_facts(facts)
    # 49 eligible A positions at x=0; one B position at x=3. Every distance is exactly 3.
    assert facts['pairs'] == [{'chain_a': '1ABC:A', 'chain_b': '1ABC:B', 'all_atoms_count': 49, 'ca_neighbors': 0}]
    assert len(rows) == 49 and all(r[-1] == 3. for r in rows)
    assert {r[2] for r in rows} == {str(i) for i in range(2, 51)}
    assert [c['nt_len'] for c in facts['chains']] == [50, 1]
    raw = pipeline.run_pipeline(str(inputs), web_output_dir=str(tmp_path / 'raw-web'))
    path.unlink()  # Prepared analysis must use the saved facts, never the source.
    prepared = analyze_prepared([facts], cfg, contact_sources={'pdb_00001abc': lambda: iter(rows)}, web_output_dir=str(tmp_path / 'prepared-web'))
    comparison = compare_directories(Path(prepared.run_output_path), Path(raw.run_output_path), 1.)
    assert comparison['status'] == 'PASS', comparison
    assert canonical_json(facts) == immutable
    raw_csv, = Path(raw.distances_dir).glob('*.csv')
    prepared_csv, = Path(prepared.distances_dir).glob('*.csv')
    assert raw_csv.read_bytes() == prepared_csv.read_bytes()
    with raw_csv.open() as stream:
        assert len(list(csv.reader(stream))) - 1 == 49
    for directory in ('raw-web', 'prepared-web'):
        summary = json.loads((tmp_path / directory / 'summary.json').read_text())
        warnings = [w for w in summary['warnings'] if w['code'] == 'NEGATIVE_OCCUPANCY_ATOMS_SKIPPED']
        assert len(warnings) == 1 and 'excluded 1 atom positions' in warnings[0]['message']
        assert not any(w['code'] == 'INPUT_FILE_SKIPPED' for w in summary['warnings'])
