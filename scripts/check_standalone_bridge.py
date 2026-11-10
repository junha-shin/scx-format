"""Read the standalone R fixtures with only the downloaded Python script."""
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "standalone"))
import scx
import numpy as np
import pandas as pd

assert Path(scx.__file__).resolve() == ROOT / "standalone/scx.py"
folder = Path(sys.argv[1])
for kind in ("Assay", "Assay5"):
    a = scx.read_scx(folder / f"{kind}-full.scx")
    assert a.shape == (2, 3)
    counts = np.array([[1, 0, 3], [0, 2, 4]])
    np.testing.assert_array_equal(a.layers["counts"].toarray(), counts)
    np.testing.assert_allclose(a.X.toarray(), np.log1p(counts), rtol=0, atol=1e-12)
    assert str(a.obs["integer"].dtype) == "Int32"
    assert str(a.obs["logical"].dtype) == "boolean"
    assert a.obs["category"].cat.ordered
    assert list(a.obs["category"].cat.categories) == ["b", "a", "unused"]
    assert pd.isna(a.obs["integer"].iloc[1])
    np.testing.assert_allclose(a.obsm["pca"], [[1.23456789, 3], [2, 4]], rtol=0, atol=1e-12)
    np.testing.assert_array_equal(a.obsp["RNA_nn"].toarray(), [[0, 1], [2, 0]])
report = {"passed": True, "standalone_R_to_Python": ["Assay", "Assay5"],
          "expression_orientation": "CELL x FEATURE", "typed_metadata": True, "analysis": True}
(folder / "bridge-report.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report))
