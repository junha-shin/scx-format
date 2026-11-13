# Seurat ↔ SCX without Python

scxSeurat 0.1.7 reads, writes, and checks SCX using R and native R helpers.
It does not need Python, reticulate, NumPy, pandas, SciPy, or AnnData.
Use Python separately for AnnData and H5AD conversion.

Historical alpha.5 includes R 0.1.3, which needs Python.
R-only use started in alpha.6. Old release files and tags are unchanged.

## Install the optional package

In your existing R environment:

```r
source("https://raw.githubusercontent.com/junha-shin/scx-format/v0.1.0/install.R")
library(scxSeurat)
packageVersion("scxSeurat")  # 0.1.7
```

See [package setup and binary support](package-install.md).
To build from a checkout, use `install.packages("adapters/seurat", repos=NULL, type="source")`.

Matrix comes with SeuratObject. Compatible installed packages are reused.
Source builds need a C compiler and zlib development files.
Use Rtools matching R on Windows, such as Rtools45 for R 4.5.
Linux needs build tools and zlib headers. macOS needs command-line developer tools and zlib.
Windows R-only CI also builds an R 4.5 ZIP.
In a compatible environment, install the release ZIP with
`install.packages("scxSeurat_0.1.7-windows-x86_64.zip", repos=NULL, type="win.binary")`.

Load the installed package with `library(scxSeurat)`.
Development `source("adapters/seurat/R/scx.R")` still needs installed native helpers.
It differs from the generated one-file scx.R, which needs no native SCX package.
The development file alone cannot provide native rename and gzip checks.

## Export, validate and import

```r
write_scx(object, "result.scx", assay = "RNA")
validate_scx("result.scx")  # cells, features, manifest file count
restored <- read_scx("result.scx")
legacy <- read_scx("result.scx", assay_class = "Assay")
# Deliberate whole-directory replacement:
write_scx(object, "result.scx", overwrite = TRUE)
```

Assay5 and legacy Assay are supported under SeuratObject >=5.
Counts are required. Normalized values are optional. Supported metadata types are kept.
Add `analysis=TRUE` for embeddings and graphs.
SCX stores cells x features; the adapter transposes Seurat matrices.
Cell graphs stay cells x cells. Legacy calls still accept `python=...` but ignore it.

## Validation and dtype scope

`validate_scx()` checks the R-readable subset before a write is published or data are read.
It checks the manifest version, required rows, unique roles/names/paths, and safe paths and symlinks.
It also checks UTF-8 tables, headers, row widths, IDs, schema types, missing values,
integer bounds, exact float32 values, categories and their codes and order.
Matrix checks cover headers, entry counts, indices, shapes, finite values, and embedding cell order.

Before parsing, zlib checks gzip trailers, CRC, multiple members, and trailing garbage.
Matrix checks use sparse values. They do not build a dense graph.
The reader keeps checked resources and uses them to build the object.

R rejects int64/uint32/uint64 metadata, INT_MIN, and numeric or boolean category labels.
Python may accept a file that R cannot read.
Integer MTX values and duplicate sums must fit within +/-(2^53-1).
Duplicate groups are also rejected if the sum of absolute terms exceeds that limit,
even if cancellation would give a small result. This avoids lost integer precision.
Sparse dimensions and entry counts must fit R/Matrix integer indexing.
See [type mappings](../compatibility/README.md#dtype-support).

## Safe publication and recovery

SCX writes and checks the new data before moving the old output.
R and Python use the same destination hash for their lock files.
Native rename refuses any existing target, even a folder created during the write.
There is no fallback that replaces that target.
Replacement keeps the old output in stage/previous until the new output is published.
Errors roll back the change. If restore fails, the backup and lock remain, with paths in the error.
Windows reparse points, including junctions, are refused.
Lock keys expand short 8.3 parent paths before hashing, so both languages use the same path.

Folder replacement has a brief gap when the final name is absent.
A hard process stop or power loss can bypass cleanup.
There is no power-loss transaction or network/cloud-sync guarantee.
See [manual recovery](safe-io.md#manual-recovery-after-a-hard-interruption).
Check pid, output, stage, and backup in the lock.
Do not remove an active writer's lock or replace an output made by another program.

## Reproduce checks

After installing the current package, run from the repo root:

```sh
Rscript adapters/seurat/scripts/run-tests.R
Rscript adapters/seurat/scripts/r-only.R .r-tests/my-r-only-run
```

The second command hides Python and blocks `system()` and `system2()`.
It tests both assay forms, metadata, counts, normalized values, embeddings, graphs, and replacement.
It writes report.json. Choose a new output folder each time.
Windows and Linux CI also run failure tests.

To test Python/R conversion separately:

```sh
python scripts/check_r_validation.py --output .compat-results/my-native-validation
python examples/quickstart/verify.py --output quickstart-output-r-native
```

These checks test damaged files, R type limits, locks in both directions, and round trips.
Python is needed for these extra checks, not for the R adapter.
Historical timings and first-alpha reports apply to their recorded versions.
They do not measure the newer native R validation path.

To test a saved single-assay Seurat object:

```sh
Rscript adapters/seurat/scripts/validate-rds.R /path/to/source.rds .r-tests/my-rds-run
```

This command also hides Python and blocks subprocesses.
It checks values and metadata and records the source SHA-256.
The saved Seurat PBMC3k object passed with 2700 cells x 13714 features,
zero matrix error, and exact metadata. It reuses the saved analysis rather than rerunning it.
The [report](../examples/seurat/r-native-report.json) lists the source hash and dependency versions.
