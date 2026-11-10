# Changelog

## Unreleased — R adapter 0.1.4

- Seurat ↔ SCX runs without Python, including the public `validate_scx()` API.
- R-native structural, typed metadata, axis and sparse matrix validation, with
  registered C helpers for complete gzip CRC/trailer validation and atomic
  no-replace directory rename. Source installation needs R build tools/zlib.
- Staging, explicit overwrite protection, rollback/backup retention and lock
  identity compatible with Python 0.1.0a5. Legacy `python` arguments are ignored.
- Linux/Windows R-only CI hides Python and forbids subprocess calls during the
  installed-package round trip; corruption/fault/termination tests and separate
  Python/R acceptance and cross-language lock checks remain in CI.
- Wide metadata integer and numeric-category restrictions remain explicit; the
  R validator checks R-readable data. Integer MTX values exceeding native R's
  exact integer range are refused. The alpha.5 tag/assets are unchanged.

## v0.1.0-alpha.5 — first public alpha

Python package: **0.1.0a5**. R adapter: **0.1.3**. Format: **0.1-alpha**.
Metadata schema: **1.0**. Alpha.1–alpha.4 were internal development iterations;
this is the first tagged GitHub prerelease, with no earlier public release promised.

### Included

- CELL × FEATURE counts and optional normalized expression, typed cell/feature
  metadata, embeddings and cell graphs in gzipped TSV/Matrix Market files.
- Manifest validation, AnnData import/export and `scx validate`, `from-h5ad`,
  `to-h5ad`; single-assay Seurat Assay/Assay5 import/export.
- Numeric widths, boolean, categories, category order/unused levels and missing
  values preserved by metadata schema 1.0, with explicit R/H5AD transport limits.
- Staged output publication, default overwrite protection, explicit replacement,
  directory rollback/retained backups and shared Python validation for R export.
- Contextual corrupt-input diagnostics, safe-I/O recovery documentation and
  fault-injection/hard-interruption tests on Linux and Windows.
- Single-pass validated import, sparse graph finite-value checks, faster metadata
  coding and gzip level 6; reproducible scaling measurements through 50k cells.
- Pinned AnnData compatibility profiles, real PBMC3k checks, unified installation
  guide and installed-package quickstart for both conversion directions.
- Seurat import restores the exact source metadata columns after assay setters;
  automatically generated QC columns are not added to an AnnData-origin table.

### Limits and migration

- One RNA feature space / selected assay only. Additional layers, loadings,
  arbitrary analysis settings, identities/history, spatial and multiome objects
  are outside this alpha profile. Keep the original analysis object.
- R export **and** import require the Python package. Native R cannot safely
  represent every pandas dtype. See the [support table](compatibility/README.md).
- CLI conversion needs actual counts in AnnData `layers['counts']`; normalized
  or scaled matrices are not substituted. Legacy schema-less metadata load as
  text with a warning; unversioned draft manifests must be regenerated.
- Directory replacement briefly removes the final name; hard termination may
  require manual recovery. Local-filesystem tests do not certify network or
  cloud-sync filesystems. See [safe I/O](docs/safe-io.md).
- This is an in-memory alpha, not an out-of-core or million-cell qualification.
  No PyPI/CRAN publication or production stability guarantee is implied.

See [release verification](docs/release-verification.md) for reproducible checks
and the [release notes](docs/releases/v0.1.0-alpha.5.md) for downloadable assets.
