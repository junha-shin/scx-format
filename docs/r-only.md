# Seurat ↔ SCX without Python

R adapter **scxSeurat 0.1.6** performs export, validation, safe publication
and import entirely in R plus its registered native helpers. No Python process,
reticulate, NumPy, pandas, SciPy or AnnData installation is required. Python is
needed separately for AnnData/H5AD conversions.

The published `v0.1.0-alpha.5` assets contain R 0.1.3 and retain their original
Python requirement. Alpha.6 includes the R-only improvement;
existing release assets and tags have not been replaced.

## Install

In your existing R environment:

```r
source("https://raw.githubusercontent.com/junha-shin/scx-format/v0.1.0-alpha.7/install.R")
library(scxSeurat)
packageVersion("scxSeurat")  # 0.1.6
```

See [installation and binary coverage](install.md). To build from a checkout, use
`install.packages("adapters/seurat", repos=NULL, type="source")`.

Matrix is a dependency of SeuratObject. Existing compatible R dependencies are
reused. Source installation requires a C toolchain and zlib development files:
Windows Rtools matching R (Rtools45 for R 4.5); Linux build tools and zlib headers;
macOS command-line developer tools and zlib. Those tools build the R package;
they are not Python dependencies. Windows R-only CI also creates an R 4.5 binary
ZIP in its artifact; the release uses a platform-labelled ZIP, installable with `install.packages("scxSeurat_0.1.6-windows-x86_64.zip",
repos = NULL, type = "win.binary")` in a compatible Windows R environment.

Use `library(scxSeurat)` after installation. Development `source("adapters/seurat/R/scx.R")`
loads the R support files but still needs the installed native helpers. A lone
copied R file cannot provide the operating-system rename and gzip guarantees.

## Export, validate and import

```r
write_scx(object, "result.scx", assay = "RNA")
validate_scx("result.scx")  # cells, features, manifest file count
restored <- read_scx("result.scx")
legacy <- read_scx("result.scx", assay_class = "Assay")
# Deliberate whole-directory replacement:
write_scx(object, "result.scx", overwrite = TRUE)
```

Both Assay5 and legacy Assay under SeuratObject >=5 are supported. Exact counts
are required, normalized is optional, and typed metadata are retained by default. Add `analysis=TRUE` for embeddings
and cell graphs. Expression is CELL × FEATURE in SCX and
transposed at the Seurat boundary. Cell graphs remain CELL × CELL.
Legacy `write_scx_seurat` / `read_scx_seurat` calls with `python = ...`
continue to work; that argument is ignored.

## Validation and dtype scope

`validate_scx()` checks the **R-readable** SCX profile before export publication
and before import. It checks versioned manifest structure, required entries,
unique roles/names/paths, path traversal and escaping symlinks, required files,
strict UTF-8 TSV headers/widths/IDs, schema structure and JSON scalar types,
nullability, native integer bounds, float32 exactness, categories/order/codes,
matrix header/count/index/shape/finite values and embedding cell order.

Gzip is checked as a stream using zlib before parsing, including complete
trailers, CRC, concatenated members and trailing garbage. Matrix checks operate
on sparse stored values, not a dense cell graph. Import retains validated
resources, so it does not reopen unchecked matrices for assembly.

Native R limitations are unchanged: wide integer metadata (int64/uint32/uint64),
INT_MIN and numeric/boolean category labels are rejected explicitly. A dataset
may therefore pass the general Python validator while being outside R's readable
subset. Integer Matrix Market values outside ±(2^53−1), including duplicate sums,
are refused rather than rounded. Sparse dimensions/entry counts must fit R/Matrix
integer indexing. See [dtype mappings](../compatibility/README.md#dtype-support).
Duplicate integer groups whose absolute term sum exceeds that range are also
refused, even if cancellation would yield a small final value; intermediate
double arithmetic could otherwise silently lose integer precision.

## Safe publication and recovery

Staging and validation complete before the previous output is moved. An exclusive
lock file uses the **same destination hash as Python 0.1.0a5**, so cooperating R
and Python writers recognize each other's locks. The native rename refuses any
existing destination, even a racing empty directory. There is no unsafe fallback
to a replacing rename. Directory replacement retains `stage/previous` until
publication succeeds; errors and interrupts roll back, and failed restoration
keeps both backup and lock with recovery paths.
Windows output reparse points, including junctions, are refused using native
file attributes. Lock keys expand short (8.3) parent paths before hashing, so
recovery and cross-language locking use the same canonical destination.

Directory replacement still has a brief name gap. Hard process termination and
power loss bypass normal cleanup; neither implementation promises durable
power-loss transactions or network/cloud-sync filesystem semantics. Follow
[manual recovery](safe-io.md#manual-recovery-after-a-hard-interruption), inspecting
the lock's `pid`, `output`, `stage` and `backup`. Do not remove another active
writer's lock or replace an independently created output.

## Reproduce checks

After installing 0.1.4, from the repository root:

```sh
Rscript adapters/seurat/scripts/run-tests.R
Rscript adapters/seurat/scripts/r-only.R .r-tests/my-r-only-run
```

The second command uses the installed package, removes Python from PATH and
blocks both `system()` and `system2()` while testing both assay forms, typed
metadata, counts, normalized, embedding, asymmetric graph and replacement.
It writes a small `report.json`; use a fresh output directory each time.
Windows and Linux CI run both the fault tests and this no-subprocess proof.

Optional Python interoperability checks remain separate:

```sh
python scripts/check_r_validation.py --output .compat-results/my-native-validation
python examples/quickstart/verify.py --output quickstart-output-r-native
```

These test corruption acceptance, intentional R dtype limits, locks in both
language directions and Seurat/AnnData round trips. They do not make Python a
runtime dependency of the R adapter. Historical performance measurements and
the first alpha's installed-package reports describe their recorded versions;
they are not measurements of this new R-native validation path.

For a precomputed single-assay Seurat object saved as RDS:

```sh
Rscript adapters/seurat/scripts/validate-rds.R /path/to/source.rds .r-tests/my-rds-run
```

This also hides Python and blocks subprocesses, verifies exact supported values
and metadata, and records the source SHA-256. The existing public Seurat PBMC3k
snapshot (2,700 cells × 13,714 features) passed this path with matrix errors 0
and exact metadata. It reuses the saved analysis; analysis itself is not rerun.
The [recorded result](../examples/seurat/r-native-report.json) includes its
source SHA-256 and dependency versions.
