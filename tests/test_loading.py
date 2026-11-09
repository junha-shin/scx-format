"""Regression checks for validated single-pass imports and bounded validation."""
import importlib
import weakref

import numpy as np
import pytest

from scx.anndata import from_anndata, to_anndata
from scx.core import write_matrix
from test_scx import sample


def test_import_reads_each_resource_once(tmp_path, monkeypatch):
    root = tmp_path/'data.scx'
    a = sample()
    from_anndata(a, root)
    validation = importlib.import_module('scx.validate')
    metadata = importlib.import_module('scx.metadata')
    seen = []
    for module, name in [(validation, 'read_matrix'), (validation, 'read_table'), (metadata, 'read_table')]:
        original = getattr(module, name)
        def counted(path, original=original):
            seen.append(str(path))
            return original(path)
        monkeypatch.setattr(module, name, counted)
    restored = to_anndata(root)
    assert len(seen) == len(set(seen)) == 6  # counts, normalized, graph, cells, features, embedding
    np.testing.assert_array_equal(restored.X.toarray(), a.X.toarray())
    np.testing.assert_array_equal(restored.obsm['X_pca'], a.obsm['X_pca'])


def test_standalone_validation_releases_prior_matrix(tmp_path, monkeypatch):
    root = tmp_path/'data.scx'
    from_anndata(sample(), root)
    validation = importlib.import_module('scx.validate')
    original = validation.read_matrix
    previous = []
    def bounded(path):
        assert all(ref() is None for ref in previous)
        matrix = original(path)
        previous.append(weakref.ref(matrix))
        return matrix
    monkeypatch.setattr(validation, 'read_matrix', bounded)
    assert validation.validate(root)['cells'] == 2
    assert len(previous) == 3


@pytest.mark.parametrize('file, matrix', [
    ('matrices/counts.mtx.gz', [[-1, 0], [0, 1]]),
    ('matrices/normalized.mtx.gz', [[np.nan, 0], [0, 1]]),
    ('graphs/connectivities.mtx.gz', np.eye(3)),
])
def test_single_pass_import_still_rejects_invalid_matrices(tmp_path, file, matrix):
    root = tmp_path/'data.scx'
    from_anndata(sample(), root)
    write_matrix(root/file, matrix)
    with pytest.raises(ValueError):
        to_anndata(root)
