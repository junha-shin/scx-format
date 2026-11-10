"""Export a validated sparse expression matrix one dense output row at a time."""
import csv
import gzip
from pathlib import Path

from .publish import staged_output
from .validate import _inspect


def _select(ids, requested, label):
    if requested is None:
        return list(range(len(ids))), list(ids)
    if isinstance(requested, str):
        raise ValueError(f"{label} must be a nonempty sequence of IDs, not one string")
    requested = list(requested)
    if not requested or any(not isinstance(x, str) for x in requested):
        raise ValueError(f"{label} must be a nonempty sequence of string IDs")
    if len(set(requested)) != len(requested):
        raise ValueError(f"{label} contains duplicate IDs")
    positions = {name: i for i, name in enumerate(ids)}
    missing = [x for x in requested if x not in positions]
    if missing:
        raise ValueError(f"unknown {label} IDs: {missing[:5]}")
    return [positions[x] for x in requested], requested


def to_tsv(path, output, *, matrix="counts", cells=None, features=None):
    """Export CELL x FEATURE TSV[.gz]; IDs select/order rows and columns.

    Existing output is always refused. Full SCX validation retains only the
    selected expression matrix, not analysis objects. Dense working space is
    one selected output row; sparse input and validation still require memory.
    """
    if matrix not in ("counts", "normalized"):
        raise ValueError("matrix must be counts or normalized")
    root, output = Path(path).resolve(), Path(output).absolute()
    if output.resolve().is_relative_to(root):
        raise ValueError("TSV output must be outside the input SCX directory")
    with staged_output(output, overwrite=False) as stage:
        _, data = _inspect(root, retain_expression=matrix)
        values = data.get(matrix)
        if values is None:
            raise ValueError("SCX has no normalized matrix")
        rows, cell_ids = _select(data['cells'].index, cells, "cells")
        columns, feature_ids = _select(data['features'].index, features, "features")
        # Compress according to the requested name, not the private staging name.
        opener = gzip.open if str(output).endswith('.gz') else open
        with opener(stage, 'wt', encoding='utf-8', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(['cell_id', *feature_ids])
            for row, cell_id in zip(rows, cell_ids):
                dense = values[row:row+1, columns].toarray().ravel()
                writer.writerow([cell_id, *(str(int(x)) if dense.dtype.kind in 'iu' else format(x, '.17g') for x in dense)])
