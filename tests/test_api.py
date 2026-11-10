"""Public object API: defaults, explicit selection and safety contracts."""
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

import scx
from scx.anndata import from_anndata
from scx.cli import main
from test_scx import sample


def test_import_does_not_load_analysis_stack():
    subprocess.run([sys.executable, '-c',
        "import scx, sys; assert not {'numpy', 'pandas', 'scipy', 'anndata'} & sys.modules.keys(); assert callable(scx.write_scx)"], check=True)


def test_basic_default_and_analysis_opt_in(tmp_path):
    a = sample()
    with pytest.warns(UserWarning, match='analysis=True'):
        scx.write_scx(a, tmp_path/'basic.scx')
    basic = scx.read_scx(tmp_path/'basic.scx')
    assert not basic.obsm and not basic.obsp
    pd.testing.assert_frame_equal(basic.obs, a.obs)
    pd.testing.assert_frame_equal(basic.var, a.var)
    np.testing.assert_array_equal(basic.layers['counts'].toarray(), a.layers['counts'].toarray())
    np.testing.assert_array_equal(basic.X.toarray(), a.X.toarray())
    scx.write_scx(a, tmp_path/'analysis.scx', analysis=True)
    restored = scx.read_scx(tmp_path/'analysis.scx')
    np.testing.assert_array_equal(restored.obsm['X_pca'], a.obsm['X_pca'])
    np.testing.assert_array_equal(restored.obsp['connectivities'].toarray(), a.obsp['connectivities'].toarray())
    assert scx.validate_scx(tmp_path/'analysis.scx')['cells'] == a.n_obs


def test_explicit_counts_source_and_no_raw_count_guessing(tmp_path):
    a = sample()
    del a.layers['counts']
    with pytest.raises(ValueError, match="counts='X'"):
        scx.write_scx(a, tmp_path/'missing.scx')
    assert not (tmp_path/'missing.scx').exists()
    with pytest.warns(UserWarning):
        scx.write_scx(a, tmp_path/'counts.scx', counts='X')
    assert not (tmp_path/'counts.scx/matrices/normalized.mtx.gz').exists()
    np.testing.assert_array_equal(scx.read_scx(tmp_path/'counts.scx').layers['counts'].toarray(), a.X.toarray())
    a.layers['raw_counts'] = a.X.copy()
    a.layers['log'] = a.X.log1p()
    with pytest.warns(UserWarning):
        scx.write_scx(a, tmp_path/'selected.scx', counts='raw_counts', normalized='log')
    np.testing.assert_array_equal(scx.read_scx(tmp_path/'selected.scx').X.toarray(), a.layers['log'].toarray())


def test_public_api_retains_safety_and_old_adapter_defaults(tmp_path):
    a = sample()
    scx.write_scx(a, tmp_path/'data.scx', analysis=True)
    before = (tmp_path/'data.scx/matrices/counts.mtx.gz').read_bytes()
    with pytest.raises(FileExistsError):
        scx.write_scx(a, tmp_path/'data.scx')
    bad = a.copy()
    bad.layers['counts'].data[0] = -1
    with pytest.raises(ValueError):
        scx.write_scx(bad, tmp_path/'data.scx', overwrite=True)
    assert (tmp_path/'data.scx/matrices/counts.mtx.gz').read_bytes() == before
    from_anndata(a, tmp_path/'legacy.scx')
    assert scx.read_scx(tmp_path/'legacy.scx').obsm.keys() == a.obsm.keys()
    with pytest.raises(ValueError, match='analysis must'):
        scx.write_scx(a, tmp_path/'invalid.scx', analysis='false')


def test_cli_uses_same_basic_and_analysis_profiles(tmp_path):
    a = sample()
    a.write_h5ad(tmp_path/'in.h5ad')
    with pytest.warns(UserWarning, match='analysis=True'):
        main(['from-h5ad', str(tmp_path/'in.h5ad'), str(tmp_path/'basic.scx')])
    assert not scx.read_scx(tmp_path/'basic.scx').obsm
    main(['from-h5ad', str(tmp_path/'in.h5ad'), str(tmp_path/'full.scx'), '--analysis'])
    assert 'X_pca' in scx.read_scx(tmp_path/'full.scx').obsm
