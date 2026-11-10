# SCX — Single Cell eXchange

SCX v0.1-alpha is a directory format for one scRNA-seq feature space. It uses gzipped TSV metadata and Matrix Market expression data. Every expression matrix is **CELL × FEATURE**: `cells.tsv.gz` row `i` is matrix row `i`, and `features.tsv.gz` row `j` is matrix column `j`.

Required files: `manifest.tsv`, `cells.tsv.gz`, `features.tsv.gz`, `matrices/counts.mtx.gz`. Optional files include `matrices/normalized.mtx.gz`, `embeddings/*.tsv.gz`, and `graphs/*.mtx.gz`. See [the specification](specification/SCX-v0.1.md).

## First public alpha

[v0.1.0-alpha.5](https://github.com/junha-shin/scx-format/releases/tag/v0.1.0-alpha.5)
packages Python **0.1.0a5** and R **0.1.3**. Start with the
[combined Python/R installation guide](docs/getting-started.md), including
Seurat → AnnData and AnnData → Seurat examples, or the
[external-user verification procedure](docs/release-verification.md).
See [changes](CHANGELOG.md) and [known limits / release notes](docs/releases/v0.1.0-alpha.5.md).
This is a GitHub prerelease; PyPI and CRAN publication are not included.

**Current main:** R adapter **0.1.4** supports Seurat ↔ SCX entirely in R,
including validation and safe output publication. See [R-only installation and
validation](docs/r-only.md). The existing alpha.5 release contains R 0.1.3 and
keeps its original Python requirement; its published assets are unchanged.

## Use

After installing the Python package in the guide:

```sh
scx from-h5ad input.h5ad output.scx
scx validate output.scx
scx to-h5ad output.scx restored.h5ad
# Explicit replacement after a successful staged write:
scx from-h5ad input.h5ad output.scx --overwrite
scx to-h5ad output.scx restored.h5ad --overwrite
```

`from-h5ad` requires `adata.layers['counts']`; it stores `adata.X` as `normalized` by default. `to-h5ad` restores counts to `layers['counts']`, and sets `X` to normalized when present, otherwise counts. Supported embeddings and cell graphs map to `obsm` and `obsp`. New exports include `metadata.schema.json` to preserve numeric widths, boolean, nullable values, strings and categorical order/unused levels in `obs`/`var`. See [metadata schema 1.0](specification/metadata-schema-v1.md). Arbitrary `uns` values remain outside this alpha round trip.

The [Seurat adapter MVP](adapters/seurat/README.md) supports single-assay R import/export, typed metadata, embeddings and cell graphs. It transposes FEATURE × CELL expression layers at the boundary. See the [Seurat PBMC3k validation](examples/seurat/README.md) for real-data reproduction and results.

See [compatibility and support tables](compatibility/README.md) for verified
AnnData environments, legacy Assay/Assay5, dtype mappings, omissions and errors.

See [scaling benchmarks](benchmarks/README.md) for reproducible Python/Seurat
time, peak-memory and file-size measurements at 1k, 10k and 50k cells, including
before/after results and the gzip level 6 storage tradeoff.

## Development

From a checkout, `pip install -e '.[test]'` installs the development environment.
Published-package checks instead install the release wheel and R tarball; run
`python examples/quickstart/verify.py --output quickstart-output-bridge` after
following the installation guide. Add `--python-only` for a no-R check.

See [PBMC3k validation](examples/pbmc3k/README.md) for a reproducible real-data
H5AD → SCX → H5AD example, pinned source checksum and recorded results.

Exports refuse existing destinations by default. `--overwrite` (Python APIs:
`overwrite=True`; R: `overwrite = TRUE`) enables staged replacement. SCX is
validated before publication; H5AD is reopened with AnnData and its axes checked.
Failed directory replacements restore the previous dataset, or retain its backup
and lock with recovery paths if restoration fails. Replacing an SCX directory
briefly makes its name unavailable to readers. See [safe I/O and recovery](docs/safe-io.md)
for interruption behavior, disk-space needs, errors and filesystem limits.

The manifest begins with `# scx_version=0.1-alpha`; unversioned draft datasets need regeneration. Versioned legacy datasets without a metadata schema load as text with a warning. Unsupported mixed objects/dtypes fail explicitly; omitted analysis objects produce a warning. GitHub Actions runs Python 3.10/3.12 compatibility tests, Seurat bridge tests, benchmarks and Windows safety tests.

```sh
pip install -e '.[test]'
pytest
```
