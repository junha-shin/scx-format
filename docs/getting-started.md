# Read, write and convert

Start with the [PBMC3k examples](../tutorials/README.md).
This page shows how to use your own objects.

Use the [no-install files](install.md) in your existing environment.
Save scx.py beside your notebook, or load scx.R with `source("scx.R")`.
Download v0.1.0 and [check its SHA-256](release-verification.md) before loading it.
Both files expose `SCX_VERSION`.

## Seurat → SCX → AnnData

```r
source("scx.R")
write_scx(seurat_object, "sample.scx")
```

```python
import scx
adata = scx.read_scx("sample.scx")
```

Save H5AD safely with `python scx.py to-h5ad sample.scx sample.h5ad`.

## AnnData → SCX → Seurat

```python
import scx
scx.write_scx(adata, "sample.scx")
```

```r
source("scx.R")
seurat_object <- read_scx("sample.scx")
```

Both languages use `write_scx(object, path)` and `read_scx(path)`.
These calls check the data. Use `validate_scx(path)` for a separate check.
No SCX package install is needed.

## Default scope

| Resource | Default write | Option |
|---|---|---|
| Counts | Required; CELL × FEATURE | Explicit source selection |
| Normalized | If present | Select source or omit |
| Cell/feature metadata | Typed numeric, boolean, categorical, missing values | [dtype limits](../compatibility/README.md) apply |
| Embeddings / cell graphs | Omitted with a warning if present | `analysis=True` / `analysis=TRUE` |
| Extra assays/layers, raw, loadings, analysis internals | Outside profile | Keep original object |

```python
scx.write_scx(adata, "analysis.scx", analysis=True)
# Only when X actually contains raw counts:
scx.write_scx(adata, "counts-only.scx", counts="X", normalized=None)
scx.write_scx(adata, "selected.scx", normalized="logcounts")
```

```r
write_scx(seurat_object, "analysis.scx", analysis = TRUE)
write_scx(seurat_object, "counts-only.scx", normalized = NULL)
write_scx(seurat_object, "selected.scx", assay = "RNA", normalized = "data")
```

Python uses `layers['counts']` / `X`. R uses counts/data in the selected assay.
Counts are required. Normalized values are optional.
On read, counts go to `layers['counts']` / `counts`.
Normalized values go to AnnData `X` / Seurat `data`.
Without normalized values, AnnData `X` uses counts.
Reads also restore supported embeddings and graphs found in the input.

## Safety, errors and migration

The R file refuses existing outputs and `overwrite=TRUE`. Choose a new folder.
Python supports explicit `overwrite=True` (CLI: `--overwrite`).
The optional R package also supports replacement and rollback.
An interrupted R file copy may leave a partial folder and lock.
Only recover the lock after the writer stops. See [recovery](safe-io.md).

For an existing output, choose a new path or use supported overwrite options.
For missing counts, select the actual raw counts.
For type errors, use supported scalar metadata types.
For damaged files or wrong axes, export again from the original object.
SCX does not guess counts or silently change unsupported metadata.

```sh
python scx.py from-h5ad input.h5ad sample.scx
python scx.py from-h5ad input.h5ad analysis.scx --analysis
python scx.py validate sample.scx
python scx.py to-h5ad sample.scx restored.h5ad
```

The API and CLI write the basic profile by default.
Add `--analysis` to include embeddings and graphs.
Legacy `from_anndata` and `write_scx_seurat` keep their full defaults.
The file format and metadata schema are unchanged. Existing SCX files remain readable.
Current R tools do not need Python. Historical alpha.5 files keep their original behavior.

Python: `scx.to_tsv("sample.scx", "counts.tsv.gz")`.
R: `to_tsv("sample.scx", "counts.tsv.gz")`.
Use `matrix="normalized"` for normalized values. Use IDs to select cells or genes.
Existing output is refused. See [TSV options and memory needs](tsv-export.md).
