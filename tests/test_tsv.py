import csv
import gzip

import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from anndata import AnnData

import scx
from scx.core import write_matrix
from test_scx import sample


def read(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8', newline='') as handle:
        return list(csv.reader(handle, delimiter='\t'))


@pytest.mark.parametrize('matrix', ['counts', 'normalized'])
@pytest.mark.parametrize('suffix', ['.tsv', '.tsv.gz'])
def test_values_selection_order_and_protection(tmp_path, matrix, suffix):
    a = sample()
    a.obs_names = ['cell"α', 'c2']; a.var_names = ['gene"β', 'g2']
    scx.write_scx(a, tmp_path/'in.scx', analysis=True)
    out = tmp_path / ('out'+suffix)
    scx.to_tsv(tmp_path/'in.scx', out, matrix=matrix,
               cells=['c2', 'cell"α'], features=['g2', 'gene"β'])
    table = read(out)
    assert table[0] == ['cell_id', 'g2', 'gene"β']
    assert [row[0] for row in table[1:]] == ['c2', 'cell"α']
    expected = a.layers['counts'] if matrix == 'counts' else a.X
    np.testing.assert_array_equal(np.array([r[1:] for r in table[1:]], dtype=float), expected.toarray()[::-1, ::-1])
    before = out.read_bytes()
    with pytest.raises(FileExistsError):
        scx.to_tsv(tmp_path/'in.scx', out)
    assert out.read_bytes() == before


@pytest.mark.parametrize('kwargs', [dict(matrix='X'), dict(cells=[]), dict(features=[]),
    dict(cells=['missing']), dict(features=['missing']), dict(cells=['c1','c1']),
    dict(features=['g1','g1']), dict(cells='c1'), dict(cells=[0])])
def test_bad_selection_has_no_output(tmp_path, kwargs):
    scx.write_scx(sample(), tmp_path/'in.scx', analysis=True)
    with pytest.raises(ValueError):
        scx.to_tsv(tmp_path/'in.scx', tmp_path/'out.tsv', **kwargs)
    assert not (tmp_path/'out.tsv').exists()
    assert not list(tmp_path.glob('.scx-lock-*'))


def test_missing_normalized_corruption_and_input_protection(tmp_path):
    scx.write_scx(sample(), tmp_path/'in.scx', normalized=None, analysis=True)
    with pytest.raises(ValueError, match='no normalized'):
        scx.to_tsv(tmp_path/'in.scx', tmp_path/'out.tsv', matrix='normalized')
    with pytest.raises(ValueError, match='outside'):
        scx.to_tsv(tmp_path/'in.scx', tmp_path/'in.scx/new/out.tsv')
    assert not (tmp_path/'in.scx/new').exists()
    # Even an unselected invalid graph must fail full input validation.
    write_matrix(tmp_path/'in.scx/graphs/connectivities.mtx.gz', sparse.eye(3))
    with pytest.raises(ValueError, match='graph'):
        scx.to_tsv(tmp_path/'in.scx', tmp_path/'out.tsv')
    assert not (tmp_path/'out.tsv').exists()


def test_one_dense_row_and_exact_integer_output(tmp_path, monkeypatch):
    a = AnnData(sparse.csr_matrix((3,20000)), obs=pd.DataFrame(index=['c1','c2','c3']),
                var=pd.DataFrame(index=[f'g{i}' for i in range(20000)]))
    a.layers['counts'] = sparse.csr_matrix(([2**53+1,7], ([0,2],[0,19999])),shape=a.shape,dtype=np.int64)
    scx.write_scx(a, tmp_path/'in.scx', normalized=None)
    original = sparse.csr_matrix.toarray
    shapes = []
    def guarded(x, *args, **kwargs):
        assert x.shape[0] == 1, 'Full matrix densification is forbidden'
        shapes.append(x.shape)
        return original(x, *args, **kwargs)
    monkeypatch.setattr(sparse.csr_matrix,'toarray',guarded)
    scx.to_tsv(tmp_path/'in.scx', tmp_path/'out.tsv')
    table = read(tmp_path/'out.tsv')
    assert table[1][1] == '9007199254740993'
    assert table[3][-1] == '7'
    assert shapes == [(1,20000)]*3


def test_failed_writer_leaves_no_output(tmp_path, monkeypatch):
    import scx.tsv
    scx.write_scx(sample(), tmp_path/'in.scx', analysis=True)
    def fail(*args, **kwargs):
        raise OSError('injected disk failure')
    monkeypatch.setattr(scx.tsv.csv, 'writer', fail)
    with pytest.raises(OSError, match='injected'):
        scx.to_tsv(tmp_path/'in.scx', tmp_path/'out.tsv.gz')
    assert not (tmp_path/'out.tsv.gz').exists()
    assert not list(tmp_path.glob('.scx-stage-*'))
    assert not list(tmp_path.glob('.scx-lock-*'))
