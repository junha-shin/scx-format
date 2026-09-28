import numpy as np
import pandas as pd
import pytest
from scipy import sparse
from anndata import AnnData

from scx.anndata import from_anndata, to_anndata
from scx.core import write_matrix
from scx.validate import validate


def sample():
    adata = AnnData(
        X=sparse.csr_matrix([[0.0, 1.2], [2.3, 0.0]]),
        obs=pd.DataFrame({"batch": ["a", "b"]}, index=["c1", "c2"]),
        var=pd.DataFrame({"symbol": ["A", "B"]}, index=["g1", "g2"]),
    )
    adata.layers["counts"] = sparse.csr_matrix([[0, 2], [3, 0]])
    adata.obsm["X_pca"] = np.array([[1.0, 2.0], [3.0, 4.0]])
    adata.obsp["connectivities"] = sparse.eye(2, format="csr")
    return adata


def test_round_trip(tmp_path):
    root = tmp_path / "sample.scx"
    from_anndata(sample(), root)
    assert validate(root)["cells"] == 2
    restored = to_anndata(root)
    assert restored.obs_names.tolist() == ["c1", "c2"]
    assert restored.var_names.tolist() == ["g1", "g2"]
    np.testing.assert_array_equal(restored.layers["counts"].toarray(), [[0, 2], [3, 0]])
    np.testing.assert_allclose(restored.X.toarray(), [[0, 1.2], [2.3, 0]])
    np.testing.assert_array_equal(restored.obsm["X_pca"], [[1, 2], [3, 4]])
    np.testing.assert_array_equal(restored.obsp["connectivities"].toarray(), np.eye(2))


def test_invalid_counts_shape(tmp_path):
    root = tmp_path / "sample.scx"
    from_anndata(sample(), root)
    write_matrix(root / "matrices/counts.mtx.gz", sparse.eye(3))
    with pytest.raises(ValueError, match="counts shape"):
        validate(root)


def test_negative_counts(tmp_path):
    root = tmp_path / "sample.scx"
    from_anndata(sample(), root)
    write_matrix(root / "matrices/counts.mtx.gz", [[-1, 0], [0, 1]])
    with pytest.raises(ValueError, match="nonnegative"):
        validate(root)
