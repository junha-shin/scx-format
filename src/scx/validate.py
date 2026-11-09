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


def _inspect(root: str | Path, *, retain=False):
    """Check every resource once; retain it only for an importing caller."""
    root = Path(root)
    manifest = read_manifest(root)
    for row in manifest.itertuples(index=False):
        if not safe_path(root, row.path).is_file():
            raise ValueError(f"missing file: {row.path}")
    cell_table = read_table(root / "cells.tsv.gz")
    feature_table = read_table(root / "features.tsv.gz")
    cells = _ids(cell_table, "cell_id")
    features = _ids(feature_table, "feature_id")
    if (manifest.role == 'metadata_schema').any():
        from .metadata import read_schema, read_metadata
        schema = read_schema(root / 'metadata.schema.json')
        try:
            cell_metadata = read_metadata(root / 'cells.tsv.gz', schema['cells'], 'cell_id', table=cell_table)
            feature_metadata = read_metadata(root / 'features.tsv.gz', schema['features'], 'feature_id', table=feature_table)
        except (KeyError, TypeError) as exc:
            raise ValueError('invalid metadata schema') from exc
    else:
        cell_metadata = cell_table.set_index('cell_id')
        feature_metadata = feature_table.set_index('feature_id')
    shape = (len(cells), len(features))
    payload = {'cells': cell_metadata, 'features': feature_metadata, 'normalized': None,
               'embeddings': {}, 'graphs': {}, 'typed_metadata': (manifest.role == 'metadata_schema').any()} if retain else None
    del cell_table, feature_table, cell_metadata, feature_metadata
    matrix = read_matrix(root / "matrices/counts.mtx.gz")
    if matrix.shape != shape:
        raise ValueError("counts shape does not match cell and feature tables")
    if not np.isfinite(matrix.data).all() or (matrix.data < 0).any():
        raise ValueError("counts must be finite and nonnegative")
    if retain:
        payload['counts'] = matrix
    del matrix
    for row in manifest.itertuples(index=False):
        path = safe_path(root, row.path)
        if row.role == "normalized":
            matrix = read_matrix(path)
            if matrix.shape != shape or not np.isfinite(matrix.data).all():
                raise ValueError("normalized matrix has invalid shape or values")
            if retain:
                payload['normalized'] = matrix
            del matrix
        elif row.role == "cell_graph":
            matrix = read_matrix(path)
            if matrix.shape != (len(cells), len(cells)) or not np.isfinite(matrix.data).all():
                raise ValueError(f"cell graph {row.name} has invalid shape or values")
            if retain:
                payload['graphs'][row.name] = matrix
            del matrix
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
            if retain:
                # Keep the original importer's round-trip float parser after
                # the strict pandas numeric check, without rereading the TSV.
                payload['embeddings'][row.name] = table.iloc[:, 1:].to_numpy(dtype=float)
            del table, values
    return {"cells": len(cells), "features": len(features), "files": len(manifest)}, payload


def validate(root: str | Path) -> dict[str, int]:
    return _inspect(root)[0]


def load_validated(root: str | Path):
    """Return checked import resources without rereading them or skipping checks."""
    return _inspect(root, retain=True)[1]
