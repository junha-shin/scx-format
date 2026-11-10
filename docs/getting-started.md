# Read, write and convert

Start with [installation](install.md) in your existing environment.

## Seurat → SCX → AnnData

```r
library(scxSeurat)
write_scx(seurat_object, "sample.scx")
```

```python
import scx
adata = scx.read_scx("sample.scx")
```

To save H5AD with SCX's protected publication, use
`scx to-h5ad sample.scx sample.h5ad`.

## AnnData → SCX → Seurat

```python
import scx
scx.write_scx(adata, "sample.scx")
```

```r
library(scxSeurat)
seurat_object <- read_scx("sample.scx")
```

Both languages use `write_scx(object, path)` and `read_scx(path)`.
`validate_scx(path)` is optional diagnosis; read/write already validate.

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

Python defaults to `layers['counts']` / `X`; R to the selected assay's
`counts` / `data`. Missing counts fail explicitly. Normalized expression is
optional. Reads restore normalized to AnnData `X` / Seurat `data` and counts
to `layers['counts']` / `counts`; AnnData `X` falls back to counts when needed.
Reads restore supported embeddings/graphs present in existing datasets.

## Safety, errors and migration

Existing outputs require `overwrite=True` (R: `overwrite=TRUE`; CLI:
`--overwrite`). A failed write preserves the previous dataset or retains a
recoverable backup and lock. Only recover a lock after the writer stops.
See [recovery](safe-io.md).

For “already exists”, choose a new output or deliberate overwrite. For missing
counts, supply actual raw counts. For dtype errors, choose supported scalar
metadata types. For corruption/axis mismatch, re-export from the original
object. SCX never guesses counts or silently coerces unsupported metadata.

```sh
scx from-h5ad input.h5ad sample.scx
scx from-h5ad input.h5ad analysis.scx --analysis
scx validate sample.scx
scx to-h5ad sample.scx restored.h5ad
```

**New API and CLI exports use the basic profile by default.** Add `--analysis`
to existing CLI workflows requiring embeddings/graphs. Legacy `from_anndata`
and `write_scx_seurat` retain their full defaults. The format/metadata schema
are unchanged, and existing SCX datasets remain readable. R no longer needs
Python; alpha.5 release assets remain unchanged.