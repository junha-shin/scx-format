# SCX — Single Cell eXchange

SCX v0.1-alpha is a directory format for one scRNA-seq feature space. It uses gzipped TSV metadata and Matrix Market expression data. Every expression matrix is **CELL × FEATURE**: `cells.tsv.gz` row `i` is matrix row `i`, and `features.tsv.gz` row `j` is matrix column `j`.

Required files: `manifest.tsv`, `cells.tsv.gz`, `features.tsv.gz`, `matrices/counts.mtx.gz`. Optional files include `matrices/normalized.mtx.gz`, `embeddings/*.tsv.gz`, and `graphs/*.mtx.gz`. See [the specification](specification/SCX-v0.1.md).

## Install and use

```sh
pip install -e '.[anndata]'
scx from-h5ad input.h5ad output.scx
scx validate output.scx
scx to-h5ad output.scx restored.h5ad
```

`from-h5ad` requires `adata.layers['counts']`; it stores `adata.X` as `normalized` by default. `to-h5ad` restores counts to `layers['counts']`, and sets `X` to normalized when present, otherwise counts. Supported embeddings and cell graphs map to `obsm` and `obsp`. AnnData metadata are serialized as text; specialized pandas dtypes and arbitrary `uns` values are outside this alpha round trip.

Seurat import/export is planned, not implemented. A Seurat adapter will transpose its FEATURE × CELL assay layers at the boundary.

## Development

```sh
pip install -e '.[test]'
pytest
```
