"""Structural validation of SCX datasets."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

from .core import read_manifest, read_matrix, read_table, safe_path


def _ids(table: pd.DataFrame, column: str) -> list[str]:
    if table.empty or len(table.columns) == 0 or table.columns[0] != column:
        raise ValueError(f"first column must be {column} and table must be nonempty")
    values = table[column].tolist()
    if any(not value for value in values) or len(set(values)) != len(values):
        raise ValueError(f"{column} values must be nonempty and unique")
    return values


def validate(root: str | Path) -> dict[str, int]:
    root = Path(root)
    manifest = read_manifest(root)
    for row in manifest.itertuples(index=False):
        if not safe_path(root, row.path).is_file():
            raise ValueError(f"missing file: {row.path}")
    cells = _ids(read_table(root / "cells.tsv.gz"), "cell_id")
    features = _ids(read_table(root / "features.tsv.gz"), "feature_id")
    counts = read_matrix(root / "matrices/counts.mtx.gz")
    if counts.shape != (len(cells), len(features)):
        raise ValueError("counts shape does not match cell and feature tables")
    if not np.isfinite(counts.data).all() or (counts.data < 0).any():
        raise ValueError("counts must be finite and nonnegative")
    for row in manifest.itertuples(index=False):
        path = safe_path(root, row.path)
        if row.role == "normalized":
            matrix = read_matrix(path)
            if matrix.shape != counts.shape or not np.isfinite(matrix.data).all():
                raise ValueError("normalized matrix has invalid shape or values")
        elif row.role == "cell_graph":
            matrix = read_matrix(path)
            if matrix.shape != (len(cells), len(cells)) or not np.isfinite(matrix.data).all():
                raise ValueError(f"cell graph {row.name} has invalid shape or values")
        elif row.role == "embedding":
            table = read_table(path)
            if len(table.columns) < 2 or table.columns[0] != "cell_id" or table.cell_id.tolist() != cells:
                raise ValueError(f"embedding {row.name} has invalid cell axis")
            try:
                values = table.iloc[:, 1:].apply(pd.to_numeric, errors="raise").to_numpy(dtype=float)
            except (ValueError, TypeError) as exc:
                raise ValueError(f"embedding {row.name} must be numeric") from exc
            if not np.isfinite(values).all():
                raise ValueError(f"embedding {row.name} has nonfinite values")
    return {"cells": len(cells), "features": len(features), "files": len(manifest)}
