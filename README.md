# SCX — Single Cell eXchange

Move one scRNA-seq assay between Seurat and AnnData using an inspectable directory.
Expression matrices always use **CELL × FEATURE** orientation.

## Use in your existing Scanpy / Seurat environment

**No SCX package installation is needed.** Save the single
[Python file](https://raw.githubusercontent.com/junha-shin/scx-format/main/standalone/scx.py)
next to your notebook/script:

```python
import scx
scx.write_scx(adata, "sample.scx")
restored = scx.read_scx("sample.scx")
```

In an existing Seurat environment:

```r
source("https://raw.githubusercontent.com/junha-shin/scx-format/main/standalone/scx.R")
write_scx(seurat_object, "sample.scx")
restored <- read_scx("sample.scx")
```

Or save the [R file](https://raw.githubusercontent.com/junha-shin/scx-format/main/standalone/scx.R)
and use `source("scx.R")` offline. The scripts reuse analysis dependencies already
present; they never install packages. R does not need Python or an SCX compiler.
See [no-install setup, dependencies and limits](docs/install.md).

The default preserves **counts, optional normalized expression, and typed
cell/feature metadata**. Python uses `layers['counts']` and `X`; R uses the
selected assay's `counts` and `data`. Counts are required and never inferred.
Add `analysis=True` (R: `analysis=TRUE`) to write embeddings and cell graphs.
Reads restore supported resources present in the input. Omissions produce a warning.

Read and write validate automatically. Existing destinations are protected.
**Standalone R only writes new directories**; it verifies staged files and
publishes the manifest last. Python retains staged replacement via
`overwrite=True`. The [optional R package](docs/package-install.md) provides
native validation and replacement with rollback. See [safety differences and
recovery](docs/install.md#shared-default-and-safety).

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

Run `python scripts/build_standalone.py` after changing adapter sources;
`--check` verifies the generated one-file scripts. CI tests standalone Python/R
without installed SCX on Linux, Windows and macOS, including both conversion
directions. See the [R package development guide](docs/r-only.md). Legacy adapters
`from_anndata` / `to_anndata` and `write_scx_seurat` / `read_scx_seurat`
remain available with their original defaults.