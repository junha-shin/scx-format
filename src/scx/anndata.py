"""AnnData adapters for the SCX alpha profile."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import warnings
import tempfile
import shutil
import json
from .metadata import write_metadata

from .core import REQUIRED, write_manifest, write_matrix, write_table, NAME
from .validate import validate, load_validated


def from_anndata(adata, root: str | Path) -> None:
    """Validate a staged dataset and publish only to a new directory."""
    root = Path(root)
    if root.exists():
        raise FileExistsError(f"output already exists: {root}")
    root.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=".scx-", dir=root.parent))
    try:
        _from_anndata(adata, stage)
        stage.rename(root)
    finally:
        if stage.exists():
            shutil.rmtree(stage)


def _from_anndata(adata, root: str | Path) -> None:
    root = Path(root)
    if "counts" not in adata.layers:
        raise ValueError("AnnData layers['counts'] is required")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("cell and feature names must be unique")
    if set(adata.layers.keys()) - {"counts", None} or len(adata.uns) or list(adata.varm.keys()) or list(adata.varp.keys()) or adata.raw is not None:
        warnings.warn("SCX omits extra layers, uns, varm, varp and raw", UserWarning)
    entries = [(path, role, name) for role, (path, name) in REQUIRED.items()]
    schema = {'schema_version': '1.0',
              'cells': write_metadata(root / 'cells.tsv.gz', adata.obs, 'cell_id'),
              'features': write_metadata(root / 'features.tsv.gz', adata.var, 'feature_id')}
    (root / 'metadata.schema.json').write_text(json.dumps(schema, ensure_ascii=False, allow_nan=False, indent=2), encoding='utf-8')
    entries.append(('metadata.schema.json', 'metadata_schema', 'metadata'))
    write_matrix(root / "matrices/counts.mtx.gz", adata.layers["counts"])
    if adata.X is not None:
        write_matrix(root / "matrices/normalized.mtx.gz", adata.X)
        entries.append(("matrices/normalized.mtx.gz", "normalized", "normalized"))
    for name, values in adata.obsm.items():
        if not NAME.fullmatch(name):
            warnings.warn(f"Skipping unsupported embedding name: {name}", UserWarning)
            continue
        array = np.asarray(values)
        if array.ndim != 2 or not np.issubdtype(array.dtype, np.number):
            warnings.warn(f"Skipping unsupported embedding: {name}", UserWarning)
            continue
        if np.issubdtype(array.dtype, np.floating):
            array = array.astype(np.float64)
        table = pd.DataFrame(array, columns=[f"dim_{i+1}" for i in range(array.shape[1])])
        table.insert(0, "cell_id", adata.obs_names.astype(str))
        path = f"embeddings/{name}.tsv.gz"
        write_table(root / path, table)
        entries.append((path, "embedding", name))
    for name, values in adata.obsp.items():
        if not NAME.fullmatch(name):
            warnings.warn(f"Skipping unsupported graph name: {name}", UserWarning)
            continue
        path = f"graphs/{name}.mtx.gz"
        write_matrix(root / path, values)
        entries.append((path, "cell_graph", name))
    write_manifest(root, entries)
    validate(root)


def to_anndata(root: str | Path):
    from anndata import AnnData

    data = load_validated(root)
    if not data['typed_metadata']:
        warnings.warn('Legacy SCX metadata have no dtype schema; loading as text', UserWarning)
    counts = data['counts']
    x = data['normalized'] if data['normalized'] is not None else counts.copy()
    adata = AnnData(X=x, obs=data['cells'], var=data['features'])
    adata.layers["counts"] = counts
    for name, values in data['embeddings'].items():
        adata.obsm[name] = values
    for name, graph in data['graphs'].items():
        adata.obsp[name] = graph
    return adata
