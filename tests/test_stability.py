import gzip
import numpy as np
import pandas as pd
import pytest
from anndata import AnnData, read_h5ad
from scipy import sparse
from scx.anndata import from_anndata
from scx.cli import main
from scx.core import read_matrix, write_table
from scx.validate import validate


def sample():
    a = AnnData(sparse.csr_matrix([[1., 0., 2.], [0., 3., 0.]]))
    a.layers['counts'] = a.X.copy()
    return a


def test_h5ad_cli_round_trip(tmp_path):
    source, folder, target = tmp_path/'in.h5ad', tmp_path/'data.scx', tmp_path/'out.h5ad'
    sample().write_h5ad(source)
    assert main(['from-h5ad', str(source), str(folder)]) == 0
    assert main(['validate', str(folder)]) == 0
    assert main(['to-h5ad', str(folder), str(target)]) == 0
    np.testing.assert_array_equal(read_h5ad(target).X.toarray(), sample().X.toarray())
    before = target.read_bytes()
    with pytest.raises(SystemExit):
        main(['to-h5ad', str(folder), str(target)])
    assert target.read_bytes() == before


def test_failed_export_leaves_no_output(tmp_path):
    a = sample()
    a.layers['counts'][0, 0] = -1
    root = tmp_path/'bad.scx'
    with pytest.raises(ValueError):
        from_anndata(a, root)
    assert not root.exists()
    assert not list(tmp_path.glob('.scx-*'))


def test_existing_output_preserved(tmp_path):
    with pytest.raises(FileExistsError):
        from_anndata(sample(), tmp_path)


@pytest.mark.parametrize('mutation', ['version', 'missing', 'duplicate', 'path', 'ids', 'embedding'])
def test_invalid_dataset(tmp_path, mutation):
    root = tmp_path/'data.scx'
    from_anndata(sample(), root)
    manifest = root/'manifest.tsv'
    text = manifest.read_text()
    if mutation == 'version':
        manifest.write_text(text.replace('0.1-alpha', '99'))
    elif mutation == 'missing':
        (root/'matrices/counts.mtx.gz').unlink()
    elif mutation == 'duplicate':
        manifest.write_text(text + text.splitlines()[-1] + '\n')
    elif mutation == 'path':
        manifest.write_text(text + '../escape\tembedding\tx\n')
    elif mutation == 'ids':
        write_table(root/'cells.tsv.gz', pd.DataFrame({'cell_id': ['same', 'same']}))
    else:
        write_table(root/'embeddings/x.tsv.gz', pd.DataFrame({'cell_id': ['1', '0'], 'dim_1': [0, 1]}))
        manifest.write_text(text + 'embeddings/x.tsv.gz\tembedding\tx\n')
    with pytest.raises((ValueError, OSError)):
        validate(root)


def test_unsupported_matrix_and_table(tmp_path):
    path = tmp_path/'matrix.mtx.gz'
    with gzip.open(path, 'wt') as f:
        f.write('%%MatrixMarket matrix array real general\n1 1\n1\n')
    with pytest.raises(ValueError):
        read_matrix(path)
    with pytest.raises(ValueError):
        write_table(tmp_path/'x.tsv.gz', pd.DataFrame({'x': ['bad\tvalue']}))


def test_loss_warning(tmp_path):
    a = sample()
    a.uns['unsupported'] = 'analysis settings'
    with pytest.warns(UserWarning, match='SCX omits'):
        from_anndata(a, tmp_path/'data.scx')


def test_nullable_metadata_round_trip(tmp_path):
    from scx.anndata import to_anndata
    a = sample()
    a.obs['group'] = pd.Categorical(['a', None])
    a.obs['number'] = pd.array([1, None], dtype='Int64')
    a.obs['flag'] = pd.array([True, None], dtype='boolean')
    a.var['score'] = np.array([0.12345678, 1.2345678, 2.3456789], dtype=np.float32)
    from_anndata(a, tmp_path/'nullable.scx')
    restored = to_anndata(tmp_path/'nullable.scx')
    pd.testing.assert_frame_equal(restored.obs, a.obs)
    pd.testing.assert_frame_equal(restored.var, a.var)


def test_float32_text_precision(tmp_path):
    from scx.anndata import to_anndata
    a = sample()
    a.X = sparse.csr_matrix(np.array([[1.2345678, 0, 2.3456789], [0, 3.456789, 0]], dtype=np.float32))
    a.obsm['X_umap'] = np.array([[1.2345678, 3.456789], [2.3456789, 0.12345678]], dtype=np.float32)
    from_anndata(a, tmp_path/'precision.scx')
    restored = to_anndata(tmp_path/'precision.scx')
    np.testing.assert_array_equal(restored.X.toarray(), a.X.toarray().astype(np.float64))
    np.testing.assert_allclose(restored.obsm['X_umap'], a.obsm['X_umap'].astype(np.float64), rtol=0, atol=1e-15)


@pytest.mark.parametrize('storage', [np.asarray, sparse.csr_matrix, sparse.csc_matrix])
def test_counts_only_dense_and_sparse(tmp_path, storage):
    from scx.anndata import to_anndata
    a = sample()
    a.layers['counts'] = storage(np.array([[1., 0., 2.], [0., 3., 0.]]))
    a.X = None
    root = tmp_path/'counts.scx'
    from_anndata(a, root)
    restored = to_anndata(root)
    assert restored.shape == (2, 3)
    assert not (root/'matrices/normalized.mtx.gz').exists()
    np.testing.assert_array_equal(restored.X.toarray(), [[1, 0, 2], [0, 3, 0]])
    np.testing.assert_array_equal(restored.layers['counts'].toarray(), restored.X.toarray())
