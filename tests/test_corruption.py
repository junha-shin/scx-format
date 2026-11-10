"""Damaged transport files must fail with a file name and no CLI traceback."""
import gzip

import pytest

from scx.anndata import from_anndata
from scx.cli import main
from scx.validate import validate
from test_stability import sample


@pytest.mark.parametrize('name,damage', [
    ('cells.tsv.gz', 'truncate'),
    ('features.tsv.gz', 'crc'),
    ('cells.tsv.gz', 'deflate'),
    ('matrices/counts.mtx.gz', 'truncate'),
    ('matrices/normalized.mtx.gz', 'crc'),
    ('matrices/counts.mtx.gz', 'missing'),
    ('matrices/counts.mtx.gz', 'mtx_short'),
    ('matrices/counts.mtx.gz', 'mtx_index'),
    ('cells.tsv.gz', 'utf8'),
    ('cells.tsv.gz', 'table_width'),
    ('manifest.tsv', 'manifest_width'),
    ('manifest.tsv', 'manifest_empty'),
    ('metadata.schema.json', 'json'),
    ('metadata.schema.json', 'schema_type'),
])
def test_corrupt_files_are_identified(tmp_path, capsys, name, damage):
    root = tmp_path/'out.scx'
    from_anndata(sample(), root)
    path = root/name
    if damage == 'truncate':
        path.write_bytes(path.read_bytes()[:-5])
    elif damage == 'deflate':
        path.write_bytes(bytes.fromhex('1f8b08000000000002ff') + b'\xff' * 20)
    elif damage == 'crc':
        data = bytearray(path.read_bytes())
        data[-8] ^= 1
        path.write_bytes(data)
    elif damage == 'missing':
        path.unlink()
    elif damage.startswith('mtx_'):
        text = '%%MatrixMarket matrix coordinate real general\n2 3 2\n1 1 1\n'
        if damage == 'mtx_index':
            text += '3 1 1\n'
        with gzip.open(path, 'wt') as handle:
            handle.write(text)
    elif damage in {'utf8', 'table_width'}:
        with gzip.open(path, 'wb') as handle:
            handle.write(b'cell_id\n\xff\n' if damage == 'utf8' else b'cell_id\tx\n0\n')
    elif damage == 'manifest_width':
        path.write_text(path.read_text() + 'extra\tembedding\tx\tunexpected\n')
    elif damage == 'manifest_empty':
        path.write_text('')
    elif damage == 'json':
        path.write_text('{')
    else:
        path.write_text('[1, 2]')
    with pytest.raises(ValueError) as error:
        validate(root)
    assert name in str(error.value).replace('\\', '/')
    with pytest.raises(SystemExit) as exc:
        main(['validate', str(root)])
    output = capsys.readouterr().err
    assert exc.value.code == 1
    assert name in output.replace('\\', '/')
    assert 'Restore or regenerate' in output
    assert 'Traceback' not in output


def test_corrupt_source_does_not_replace_h5ad(tmp_path, capsys):
    root, target = tmp_path/'out.scx', tmp_path/'precious.h5ad'
    from_anndata(sample(), root)
    (root/'matrices/counts.mtx.gz').write_bytes(b'bad gzip')
    target.write_bytes(b'existing output')
    with pytest.raises(SystemExit):
        main(['to-h5ad', str(root), str(target), '--overwrite'])
    assert target.read_bytes() == b'existing output'
    assert 'counts.mtx.gz' in capsys.readouterr().err
