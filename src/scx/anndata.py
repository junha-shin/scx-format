"""AnnData adapters for the SCX alpha profile."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

from .core import REQUIRED, read_manifest, read_matrix, read_table, write_manifest, write_matrix, write_table, NAME
from .validate import validate


def _table(frame: pd.DataFrame, key: str) -> pd.DataFrame:
    table = frame.copy()
    table.insert(0, key, frame.index.astype(str))
    if table.columns.duplicated().any():
        raise ValueError(f"duplicate metadata column: {key}")
    return table.fillna("").astype(str)


def from_anndata(adata, root: str | Path) -> None:
    root = Path(root)
    if "counts" not in adata.layers:
        raise ValueError("AnnData layers['counts'] is required")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("cell and feature names must be unique")
    entries = [(path, role, name) for role, (path, name) in REQUIRED.items()]
    write_table(root / "cells.tsv.gz", _table(adata.obs, "cell_id"))
    write_table(root / "features.tsv.gz", _table(adata.var, "feature_id"))
    write_matrix(root / "matrices/counts.mtx.gz", adata.layers["counts"])
    if adata.X is not None:
        write_matrix(root / "matrices/normalized.mtx.gz", adata.X)
        entries.append(("matrices/normalized.mtx.gz", "normalized", "normalized"))
    for name, values in adata.obsm.items():
        if not NAME.fullmatch(name):
            continue
        array = np.asarray(values)
        if array.ndim != 2 or not np.issubdtype(array.dtype, np.number):
            continue
        table = pd.DataFrame(array, columns=[f"dim_{i+1}" for i in range(array.shape[1])])
        table.insert(0, "cell_id", adata.obs_names.astype(str))
        path = f"embeddings/{name}.tsv.gz"
        write_table(root / path, table)
        entries.append((path, "embedding", name))
    for name, values in adata.obsp.items():
        if not NAME.fullmatch(name):
            continue
        path = f"graphs/{name}.mtx.gz"
        write_matrix(root / path, values)
        entries.append((path, "cell_graph", name))
    write_manifest(root, entries)
    validate(root)


def to_anndata(root: str | Path):
    from anndata import AnnData

    root = Path(root)
    validate(root)
    manifest = read_manifest(root)
    cells = read_table(root / "cells.tsv.gz").set_index("cell_id")
    features = read_table(root / "features.tsv.gz").set_index("feature_id")
    counts = read_matrix(root / "matrices/counts.mtx.gz")
    normalized = manifest[manifest.role == "normalized"]
    x = read_matrix(root / normalized.iloc[0]["path"]) if len(normalized) else counts.copy()
    adata = AnnData(X=x, obs=cells, var=features)
    adata.layers["counts"] = counts
    for row in manifest.itertuples(index=False):
        if row.role == "embedding":
            adata.obsm[row.name] = read_table(root / row.path).iloc[:, 1:].to_numpy(dtype=float)
        elif row.role == "cell_graph":
            adata.obsp[row.name] = read_matrix(root / row.path)
    return adata
