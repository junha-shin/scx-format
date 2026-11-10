"""Core SCX file and manifest handling."""
from __future__ import annotations

import gzip
import csv
from pathlib import Path, PurePosixPath
import re

import pandas as pd
import numpy as np
from scipy import sparse
from scipy.io import mmread, mmwrite
from .errors import resource

HEADER = ("path", "role", "name")
VERSION = "0.1-alpha"
REQUIRED = {
    "cells": ("cells.tsv.gz", "cells"),
    "features": ("features.tsv.gz", "features"),
    "counts": ("matrices/counts.mtx.gz", "counts"),
}
ROLES = set(REQUIRED) | {"normalized", "embedding", "cell_graph", "metadata_schema"}
NAME = re.compile(r"^[A-Za-z0-9_.-]+$")


def safe_path(root: Path, name: str) -> Path:
    path = PurePosixPath(name)
    if (not name or path.is_absolute() or "\\" in name or any(c in name for c in '\t\r\n\x00')
            or any(p in ("", ".", "..") for p in name.split('/')) or ":" in name):
        raise ValueError(f"unsafe manifest path: {name}")
    target = root.joinpath(*path.parts)
    if not target.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"path escapes dataset: {name}")
    return target


def read_manifest(root: Path) -> pd.DataFrame:
    root = Path(root)
    with resource(root / 'manifest.tsv'):
        return _read_manifest(root)


def _read_manifest(root):
    with open(safe_path(root, 'manifest.tsv'), encoding="utf-8", newline='') as handle:
        if handle.readline().rstrip("\n\r") != f"# scx_version={VERSION}":
            raise ValueError("missing or unsupported SCX manifest version")
        rows = list(csv.reader(handle, delimiter='\t', strict=True))
    if not rows or tuple(rows[0]) != HEADER:
        raise ValueError("manifest header must be path, role, name")
    if any(len(row) != len(HEADER) for row in rows[1:]):
        raise ValueError('manifest rows must contain exactly path, role, name')
    manifest = pd.DataFrame(rows[1:], columns=HEADER)
    if manifest.duplicated("path").any() or manifest.duplicated(["role", "name"]).any():
        raise ValueError("duplicate manifest path or role/name")
    for row in manifest.itertuples(index=False):
        safe_path(root, row.path)
        if row.role not in ROLES or not NAME.fullmatch(row.name):
            raise ValueError(f"invalid role or name: {row.role}/{row.name}")
    for role, (path, name) in REQUIRED.items():
        rows = manifest[manifest.role == role]
        if len(rows) != 1 or rows.iloc[0]["path"] != path or rows.iloc[0]["name"] != name:
            raise ValueError(f"required manifest entry invalid: {role}")
    if (manifest.role == "normalized").sum() > 1:
        raise ValueError("only one normalized matrix is supported")
    schemas = manifest[manifest.role == 'metadata_schema']
    if len(schemas) > 1 or (len(schemas) and (schemas.iloc[0]['path'] != 'metadata.schema.json' or schemas.iloc[0]['name'] != 'metadata')):
        raise ValueError('invalid metadata schema manifest entry')
    return manifest


def write_manifest(root: Path, entries: list[tuple[str, str, str]]) -> None:
    root.mkdir(parents=True, exist_ok=True)
    with open(root / "manifest.tsv", "w", encoding="utf-8", newline="") as handle:
        handle.write(f"# scx_version={VERSION}\n")
        pd.DataFrame(entries, columns=HEADER).to_csv(handle, sep="\t", index=False)
    read_manifest(root)


def read_table(path: Path) -> pd.DataFrame:
    with resource(path):
        return _read_table(path)


def _read_table(path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t", strict=True))
    if not rows or not rows[0] or any(not h for h in rows[0]) or len(set(rows[0])) != len(rows[0]):
        raise ValueError("table headers must be nonempty and unique")
    if any(len(row) != len(rows[0]) for row in rows) or any(any(c in value for c in "\t\r\n") for row in rows for value in row):
        raise ValueError("invalid table width or tab/newline in value")
    return pd.DataFrame(rows[1:], columns=rows[0])


def write_table(path: Path, table: pd.DataFrame) -> None:
    def invalid(values):
        return any(any(c in str(value) for c in "\t\r\n") for value in values)
    if invalid(table.columns) or any(invalid(column.array) for _, column in table.items()
                                     if not pd.api.types.is_numeric_dtype(column.dtype)):
        raise ValueError("tabs and newlines are unsupported in tables")
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, sep="\t", index=False, compression={"method": "gzip", "compresslevel": 6})


def read_matrix(path: Path) -> sparse.csr_matrix:
    with resource(path):
        return _read_matrix(path)


def _read_matrix(path):
    with gzip.open(path, "rb") as handle:
        header = handle.readline().decode("ascii").strip().lower().split()
        if len(header) != 5 or header[:3] != ["%%matrixmarket", "matrix", "coordinate"] or header[3] not in {"real", "integer"} or header[4] != "general":
            raise ValueError("matrix must use Matrix Market coordinate real/integer general")
        handle.seek(0)
        matrix = sparse.csr_matrix(mmread(handle, spmatrix=True))
        # Readers can stop after declared entries. Consume the gzip trailer too
        # so truncated streams and CRC failures cannot validate successfully.
        while handle.read(1024 * 1024):
            pass
        return matrix


def write_matrix(path: Path, matrix) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    matrix = sparse.coo_matrix(matrix)
    # Promote float32 before requesting 17 digits. Some SciPy writers emit
    # embedded NUL bytes when precision exceeds the native float32 buffer.
    if np.issubdtype(matrix.dtype, np.floating):
        matrix = matrix.astype(np.float64)
    with gzip.open(path, "wb", compresslevel=6) as handle:
        mmwrite(handle, matrix, symmetry="general", precision=17)
