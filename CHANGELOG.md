# Changelog

## v0.1.0 - first stable utilities

- Existing-object quickstarts replace synthetic data setup with four API calls.
- First-use dependency, counts-source and recovery guidance; actionable errors.
- README files use plain ASCII, enforced by repository instructions and CI.
- Analysis-only CI environments execute both tutorial levels, round-trip and
  failure scenarios without installing SCX. Optional R package is 0.1.7.
- On-disk format 0.1-alpha and metadata schema 1.0 remain compatible.

### Guided tutorials

- Executed Python notebook and sectioned RStudio script use the same tiny data
  and pinned alpha.7 files, with verified downloads and no SCX installation.
- Walk through typed metadata, SCX round-trip, counts/normalized TSV, ordered
  subsets, overwrite refusal and an optional own-data example.
- Test both lessons on three operating systems, including offline full reruns,
  preservation of earlier output, damaged-cache refusal and R/Python TSV equality.
  Published alpha.7 assets and tags are unchanged.

## v0.1.0-alpha.7 — versioned no-install files

Standalone **v0.1.0-alpha.7**; optional Python **0.1.0a7** / R **0.1.6**.
Format **0.1-alpha** and metadata schema **1.0** unchanged.

- Primary release assets are self-contained `scx.py` and `scx.R`, canonical
  UTF-8/LF, with SHA-256 checksums, offline `USAGE.md` and verification evidence.
  Both expose `SCX_VERSION`; docs/downloads are pinned to this tag.
- Existing Scanpy/Seurat dependencies are reused without automatic installation.
  R scripts need no SCX package, Python, compiler or compiled SCX helper.
- Unified read/write/validate APIs plus matrix `to_tsv`: counts/normalized,
  ordered ID subsets, UTF-8 quoting, gzip and refusal of existing TSV files.
- Pure R validates single-member gzip CRC/length using scratch disk. New SCX
  outputs use shared exclusive locks, verified copies and manifest-last
  publication. Existing SCX and `overwrite=TRUE` are refused in standalone R.
  Optional native R packages retain stronger publication/replacement features.
- R TSV rendering converts to CSR once and formats only stored entries into a
  zero-filled text row. No numeric dense matrix or added package is needed.
- Recorded PBMC3k optimization: full R outputs 3.64–4.09× faster in two local
  repeats, 0 numeric error, identical plain TSV hashes and unchanged file sizes.
  Observed peak memory increased by 0.5–21.8 MiB; complete profiles/tradeoffs
  remain recorded rather than claiming memory savings.
- Three-platform CI tests the staged release files with SCX absent, including
  Assay/Assay5, both conversion directions, dtype/axes, corruption, interrupted
  copies, sparse edge cases and TSV agreement. Checksum tampering fails.
- Previous alpha.5/alpha.6 tags and assets remain unchanged. Numeric R dtype
  limits, one feature space, optional analysis scope and filesystem limits apply.

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
