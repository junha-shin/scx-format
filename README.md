# SCX — Single Cell eXchange

Move one scRNA-seq assay between Seurat and AnnData using an inspectable directory.
Expression matrices always use **CELL × FEATURE** orientation.

## Install in your existing environment

R / Seurat:

```r
source("https://raw.githubusercontent.com/junha-shin/scx-format/v0.1.0-alpha.6/install.R")
library(scxSeurat)
```

Python / scverse:

```sh
pip install "https://github.com/junha-shin/scx-format/releases/download/v0.1.0-alpha.6/scx_format-0.1.0a6-py3-none-any.whl"
```

Compatible dependencies are reused. R does not need Python. Prebuilt R packages
target R 4.5 on Windows x86_64, macOS arm64 and Ubuntu x86_64. Other environments
can build from source. See [installation](docs/install.md) for exact coverage.
Packages are distributed through GitHub, not PyPI or CRAN.

## Read and write

```r
write_scx(seurat_object, "sample.scx")
restored <- read_scx("sample.scx")
```

```python
import scx
scx.write_scx(adata, "sample.scx")
restored = scx.read_scx("sample.scx")
```

The default preserves **counts, optional normalized expression, and typed
cell/feature metadata**. Python uses `layers['counts']` and `X`; R uses the
selected assay's `counts` and `data`. Counts are required and never inferred.
Add `analysis=True` (R: `analysis=TRUE`) to write embeddings and cell graphs.
Reads restore supported resources present in the input. Omissions produce a warning.

Read and write validate automatically. Existing destinations require explicit
`overwrite=True` / `overwrite=TRUE`; replacement is staged and validated.
See [safe I/O and recovery](docs/safe-io.md).

- [Usage and Seurat ↔ AnnData examples](docs/getting-started.md)
- [Support, dtype limits and errors](compatibility/README.md)
- [Format](specification/SCX-v0.1.md) and [metadata schema](specification/metadata-schema-v1.md)
- [Reproducible package verification](docs/release-verification.md)
- [PBMC3k](examples/pbmc3k/README.md), [Seurat validation](examples/seurat/README.md), [benchmarks](benchmarks/README.md)
- [Changes](CHANGELOG.md) and [alpha.6 release notes](docs/releases/v0.1.0-alpha.6.md)

One feature space / selected assay is supported. Extra assays, arbitrary analysis
internals, loadings, spatial and multiome data remain outside the profile.
Keep the original analysis object.

## Development

```sh
pip install -e '.[test]'
pytest
```

See the [R development guide](docs/r-only.md). Legacy adapters
`from_anndata` / `to_anndata` and `write_scx_seurat` / `read_scx_seurat`
remain available with their original defaults.