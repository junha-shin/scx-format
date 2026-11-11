# Changelog

## Unreleased — faster R matrix TSV rendering

- Profiling identified repeated formatting of implicit zeros and sparse row
  extraction as the main costs in the baseline R exporter.
- Convert to CSR once and format only stored entries into a zero-filled text row;
  keep explicit stored/signed zeros and the existing 17-digit numeric format.
- No new packages or numeric dense matrices. Full validation, gzip CRC, output
  protection, staging and verified exclusive publication remain unchanged.
- Add sparse edge-case, no-densification and compatible-run comparison tests.
- PBMC3k full R exports measured 3.64–4.09× faster with 0 numeric error; plain
  TSV hashes and all file sizes unchanged. Peak memory increased by 0.5–21.8 MiB
  in these observations; complete before/after results and profiles are recorded.


## Unreleased — PBMC3k TSV verification and measurements

- Checksum-pinned 10x PBMC3k preparation and no-install R/Python TSV benchmark.
- Counts/normalized, plain/gzip and ordered subsets checked against both the
  source sparse matrix and the other language, one row at a time.
- 24 isolated measured exports / 12 paired comparisons; all values and axes
  matched, with maximum absolute error 0 in the recorded Windows run.
- Recorded wall/API time, worker peak memory, output sizes, implementation/input
  hashes, versions and reproducible methods; three-platform synthetic CI smoke.


## Unreleased — matrix TSV export

- Shared no-install `to_tsv(path, output, matrix, cells, features)` API for Python
  and R. Counts/normalized selection, ordered ID subsets, UTF-8 TSV and TSV.gz.
- Full SCX validation retains only the selected expression matrix; output is
  densified one row at a time. Sparse input/validation memory remains required.
- Existing outputs and writes inside the SCX input are refused. Python uses
  staged exclusive publication; pure R stages, exclusively creates, copies and
  verifies bytes, retaining partial files/locks on publication failure.
- Tests cover output values, ID ordering/quoting, gzip, invalid/corrupt inputs,
  single-row memory behavior, failed copies and R/Python TSV agreement.

## Unreleased — no-install scripts

- Self-contained `standalone/scx.py` and `standalone/scx.R`, generated from the
  reference adapters. Reuse existing Scanpy/Seurat dependencies; no SCX package
  installation, SCX compiler, Python-from-R bridge or automatic dependency installs.
- Same read/write API and basic/analysis profiles. Typed metadata and canonical
  CELL × FEATURE axes remain validated. Python retains reference safe publication.
- Pure R single-member gzip CRC/length checks using scratch disk, exclusive shared
  locks, new-directory reservation, SHA-256 copy verification and manifest-last
  publication. Existing outputs and overwrite=TRUE are refused in standalone R.
- Installation-free usage is the documentation default. Alpha.6 installed packages
  remain optional and unchanged, including native R replacement/rollback.
- Three-platform no-SCX CI verifies Assay/Assay5, types, optional analysis, both
  conversion directions, corrupt gzip and interrupted-copy rejection.


## v0.1.0-alpha.6 — simpler installation and object I/O

Python **0.1.0a6**, R **0.1.5**. Format **0.1-alpha** and schema **1.0** unchanged.

- Same `read_scx` / `write_scx` names in R and Python; validate automatically.
- Basic default: counts, optional normalized matrix and typed metadata. Add
  `analysis=True` / `analysis=TRUE` for embeddings and cell graphs. CLI exports
  now use the basic default: existing workflows need `--analysis` to retain them.
- Explicit source selection via `counts` / `normalized`. Counts are never guessed.
  Legacy adapters preserve their original broader defaults.
- One install in the existing environment: release wheel and checksum-verified
  R binary installer; only missing R dependencies are installed. No Python for R.
- SciPy minimum lowered to 1.11.4 with a dedicated compatibility CI profile.
- R binary build/install round trips on Windows, macOS and Ubuntu, plus existing
  corruption/interruption/overwrite and cross-language tests. Published alpha.5
  packages are unchanged. Single-assay and native R dtype limits still apply.

### Included R-only groundwork (previously main-only 0.1.4)

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
