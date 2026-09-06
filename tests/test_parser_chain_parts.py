"""Minimal independent reproduction of atom loss from disconnected chains."""
from collections import Counter
import gzip

import gemmi
import pytest

from pdb2net import file_parser


@pytest.mark.parametrize('extension', ['pdb', 'cif', 'pdb.gz', 'cif.gz'])
@pytest.mark.parametrize('model_count', [1, 2])
def test_disconnected_author_chains_retain_every_atom(tmp_path, extension, model_count):
    lines = ['HEADER    TEST                                    01-JAN-00   1ABC\n']
    for model in range(1, model_count + 1):
        lines.append(f'MODEL     {model:4d}\n')
        for serial, (group, atom, residue, chain, seq, element) in enumerate([
            ('ATOM', 'CA', 'ALA', 'A', 1, 'C'),
            ('ATOM', 'CA', 'GLY', 'B', 1, 'C'),
            ('HETATM', 'CA', 'LEU', 'A', 2, 'C'),
            ('HETATM', 'O', 'HOH', 'B', 2, 'O'),
        ], start=1):
            lines.append(f'{group:6s}{serial:5d} {atom:^4s} {residue:3s} {chain}{seq:4d}    {float(serial):8.3f}{float(model):8.3f}{0.:8.3f}  1.00 20.00          {element:>2s}  \n')
        lines.append('ENDMDL\n')
    lines.append('END\n')
    content = ''.join(lines)
    if extension.startswith('cif'):
        source = gemmi.read_pdb_string(content)
        source.name = '1ABC'
        source.setup_entities()
        content = source.make_mmcif_document().as_string()
    path = tmp_path / f'1abc.{extension}'
    payload = content.encode()
    path.write_bytes(gzip.compress(payload) if extension.endswith('.gz') else payload)
    reference = gemmi.read_structure(str(path))
    actual = file_parser.parse_structure_input(str(path))['structure']
    assert len(actual) == model_count
    for actual_model, expected_model in zip(actual, reference):
        assert [chain.name for chain in actual_model] == ['A', 'B']
        def atoms(model):
            return Counter((c.name, r.name, str(r.seqid), a.name, tuple(a.pos)) for c in model for r in c for a in r)
        assert atoms(actual_model) == atoms(expected_model)
        assert sum(atoms(actual_model).values()) == 4
