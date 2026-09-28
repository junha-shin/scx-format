"""Core SCX file and manifest handling."""
from __future__ import annotations

import gzip
from pathlib import Path, PurePosixPath
import re

import pandas as pd
from scipy import sparse
from scipy.io import mmread, mmwrite

HEADER = ("path", "role", "name")
REQUIRED = {
    "cells": ("cells.tsv.gz", "cells"),
    "features": ("features.tsv.gz", "features"),
    "counts": ("matrices/counts.mtx.gz", "counts"),
}
ROLES = set(REQUIRED) | {"normalized", "embedding", "cell_graph"}
NAME = re.compile(r"^[A-Za-z0-9_.-]+$")


def safe_path(root: Path, name: str) -> Path:
    path = PurePosixPath(name)
    if not name or path.is_absolute() or "\\" in name or any(p in ("", ".", "..") for p in path.parts) or ":" in name:
        raise ValueError(f"unsafe manifest path: {name}")
    target = root.joinpath(*path.parts)
    if not target.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"path escapes dataset: {name}")
    return target


def read_manifest(root: Path) -> pd.DataFrame:
    manifest = pd.read_csv(root / "manifest.tsv", sep="\t", dtype=str, keep_default_na=False)
    if tuple(manifest.columns) != HEADER:
        raise ValueError("manifest header must be path, role, name")
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
    return manifest


def write_manifest(root: Path, entries: list[tuple[str, str, str]]) -> None:
    root.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(entries, columns=HEADER).to_csv(root / "manifest.tsv", sep="\t", index=False)
    read_manifest(root)


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def write_table(path: Path, table: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, sep="\t", index=False, compression="gzip")


def read_matrix(path: Path) -> sparse.csr_matrix:
    with gzip.open(path, "rb") as handle:
        return sparse.csr_matrix(mmread(handle))


def write_matrix(path: Path, matrix) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb") as handle:
        mmwrite(handle, sparse.coo_matrix(matrix))
