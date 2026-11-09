import json
import numpy as np
import pandas as pd
import pytest
from scx.anndata import from_anndata, to_anndata
from scx.validate import validate
from scx.core import read_table, write_table
from test_stability import sample


@pytest.mark.parametrize('dtype,values', [
    ('int8', [-128, 127]), ('uint64', [0, 2**64-1]),
    ('Int64', [2**63-1, None]), ('UInt64', [2**64-1, None]),
    ('float32', [0.12345678, np.nan]), ('float64', [1.2345678901234567, np.nan]),
    ('Float32', [0.12345678, None]), ('Float64', [0.12345678, None]),
    ('bool', [False, True]), ('boolean', [False, None]),
    ('object', ['', None]), ('string', ['null', None]),
    ('object', ['tab\tline\nquote"slash\\', 'null']),
])
def test_dtype_round_trip(tmp_path, dtype, values):
    a = sample()
    a.obs['value'] = pd.array(values, dtype=dtype)
    from_anndata(a, tmp_path/'data.scx')
    pd.testing.assert_frame_equal(to_anndata(tmp_path/'data.scx').obs, a.obs)


@pytest.mark.parametrize('categories', [['z', '', 'unused'], [5, 2, 99]])
def test_category_order_and_unused(tmp_path, categories):
    a = sample()
    a.obs['value'] = pd.Categorical([categories[1], None], categories=categories, ordered=True)
    from_anndata(a, tmp_path/'data.scx')
    pd.testing.assert_frame_equal(to_anndata(tmp_path/'data.scx').obs, a.obs)


@pytest.mark.parametrize('mutation', ['dtype', 'version', 'columns', 'code', 'boolean', 'overflow'])
def test_corrupt_schema_rejected(tmp_path, mutation):
    a = sample()
    a.obs['category'] = pd.Categorical(['a', 'b'])
    a.obs['flag'] = [True, False]
    a.obs['tiny'] = np.array([1, 2], dtype='int8')
    root = tmp_path/'data.scx'
    from_anndata(a, root)
    path = root/'metadata.schema.json'
    schema = json.loads(path.read_text())
    table = read_table(root/'cells.tsv.gz')
    if mutation == 'dtype':
        schema['cells']['columns'][0]['dtype'] = 'pickle'
    elif mutation == 'version':
        schema['schema_version'] = '99'
    elif mutation == 'columns':
        schema['cells']['columns'].pop()
    elif mutation == 'code':
        table.loc[0, 'category'] = '99'
    elif mutation == 'boolean':
        table.loc[0, 'flag'] = '1'
    else:
        table.loc[0, 'tiny'] = '999'
    path.write_text(json.dumps(schema))
    write_table(root/'cells.tsv.gz', table)
    with pytest.raises(ValueError):
        validate(root)


def test_legacy_read(tmp_path):
    a = sample()
    root = tmp_path/'data.scx'
    from_anndata(a, root)
    path = root/'manifest.tsv'
    path.write_text('\n'.join(x for x in path.read_text().splitlines() if 'metadata_schema' not in x)+'\n')
    write_table(root/'cells.tsv.gz', pd.DataFrame({'cell_id': ['0', '1'], 'legacy': ['null', '']}))
    with pytest.warns(UserWarning, match='Legacy'):
        restored = to_anndata(root)
    assert restored.shape == a.shape
    assert restored.obs['legacy'].tolist() == ['null', '']


def test_unsupported_objects_fail_without_output(tmp_path):
    a = sample()
    a.obs['objects'] = [[1], {'x': 2}]
    with pytest.raises(ValueError, match='only strings'):
        from_anndata(a, tmp_path/'data.scx')
    assert not (tmp_path/'data.scx').exists()


def test_typed_h5ad_cli(tmp_path):
    from anndata import read_h5ad
    from scx.cli import main
    a = sample()
    a.obs['integer'] = pd.array([2**63-1, None], dtype='Int64')
    a.obs['flag'] = pd.array([True, None], dtype='boolean')
    a.obs['category'] = pd.Categorical(['b', None], categories=['b', 'a', 'unused'], ordered=True)
    a.obs['float'] = np.array([0.12345678, np.nan], dtype='float32')
    source, target, folder = tmp_path/'source.h5ad', tmp_path/'target.h5ad', tmp_path/'data.scx'
    a.write_h5ad(source)
    main(['from-h5ad', str(source), str(folder)])
    main(['to-h5ad', str(folder), str(target)])
    pd.testing.assert_frame_equal(read_h5ad(source).obs, read_h5ad(target).obs, check_exact=True)
