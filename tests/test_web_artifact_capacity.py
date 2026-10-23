"""CSV capacity is independent of the bounded network JSON capacity."""
import pytest
from pathlib import Path

from pdb2net import outputs
from pdb2net.input_contract import InputValidationError
from pdb2net.structure_identity import identity_from_official_id


def test_csv_can_exceed_network_cap_but_not_its_own_cap(tmp_path, monkeypatch):
    monkeypatch.setattr(outputs, 'MAX_WEB_ARTIFACT_BYTES', 8)
    monkeypatch.setattr(outputs, 'MAX_WEB_CSV_BYTES', 32)
    source = tmp_path / 'contacts.csv'
    source.write_text('chain,distance\nA,1.0\n')
    record = outputs._artifact_metadata(source, kind='interaction')
    assert record['rows'] == 1 and record['columns'] == ['chain', 'distance']
    with pytest.raises(InputValidationError, match='WEB_ARTIFACT_TOO_LARGE'):
        outputs._artifact_metadata(source, kind='network')
    source.write_text('chain,distance\n' + 'A,1.0\n' * 10)
    with pytest.raises(InputValidationError, match='WEB_ARTIFACT_TOO_LARGE'):
        outputs._artifact_metadata(source, kind='interaction')


@pytest.mark.parametrize('aggregate_failure', [False, True])
def test_public_csv_publication_uses_csv_and_aggregate_budgets(tmp_path, monkeypatch, aggregate_failure):
    monkeypatch.setattr(outputs, 'MAX_WEB_ARTIFACT_BYTES', 8)
    monkeypatch.setattr(outputs, 'MAX_WEB_CSV_BYTES', 32)
    monkeypatch.setattr(outputs, 'MAX_WEB_ARTIFACT_TOTAL_BYTES', 10 if aggregate_failure else 100)
    paths = outputs.create_run_output_paths(str(tmp_path / 'internal'))
    source = Path(paths.distances_dir) / 'contacts.csv'
    source.write_text('chain,distance\nA,1.0\n')
    identity = identity_from_official_id('1abc').as_dict()
    outputs.write_run_manifest(paths.manifest_file, input_files=['/inputs/1abc.pdb'], output_paths=paths,
        config_snapshot={'networks': {'chain_per_pdb': False}, 'structure_model_policy': 'first', 'export_detailed_interactions': True},
        status='success', started_at='2026-01-02T03:04:05', finished_at='2026-01-02T03:04:06', total_time=1,
        references={'manifest_id': 'fixture'}, extra_counts={'structures': 1}, identities=[identity],
        structure_inputs=[{'file': '/inputs/1abc.pdb', 'identity': identity, 'format': 'pdb', 'kind': 'pdb',
                           'embedded_annotation_counts': dict.fromkeys(['uniprot', 'pfam', 'cath', 'scop2'], 0)}])
    target = tmp_path / 'web'
    if aggregate_failure:
        with pytest.raises(InputValidationError, match='WEB_ARTIFACT_TOTAL_TOO_LARGE'):
            outputs.collect_web_outputs(paths, str(target))
        assert not list(target.rglob('*.csv')) and not (target / 'summary.json').exists()
    else:
        report = outputs.collect_web_outputs(paths, str(target))
        assert report['status'] == 'success'
        assert (target / 'interactions/contacts.csv').read_bytes() == source.read_bytes()
        assert report['artifacts']['interactions'][0]['rows'] == 1
