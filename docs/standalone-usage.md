# SCX v0.1.0 — no-install files

Release: https://github.com/junha-shin/scx-format/releases/tag/v0.1.0

Download scx.py or scx.R and SHA256SUMS.
Use your existing Scanpy or Seurat environment. No SCX package or compiler is needed.
The files do not install dependencies. They use UTF-8 and LF line endings.

## Verify before running

Compute SHA-256. Compare it with the entry for that exact file in SHA256SUMS.

```sh
# Linux; all listed files must be present for this combined command
sha256sum -c SHA256SUMS
# macOS: hash one downloaded file, then compare with its entry
shasum -a 256 scx.py
```

```powershell
Get-FileHash .\scx.py -Algorithm SHA256
Get-FileHash .\scx.R -Algorithm SHA256
```

For one language file, check just that file's entry.
The full list also covers this guide and optional packages.
Use the fixed release URLs below in saved workflows, not latest or main.
Keep the file and checksum with your analysis.
If a hash differs, do not load the file. Download the release file again.

## Python / existing Scanpy environment

Save scx.py beside your script or notebook:

```python
import scx
assert scx.SCX_VERSION == "v0.1.0"
scx.write_scx(adata, "sample.scx")
restored = scx.read_scx("sample.scx")
scx.to_tsv("sample.scx", "counts.tsv.gz")
scx.to_tsv("sample.scx", "selected.tsv.gz", matrix="normalized",
           cells=["cell3", "cell1"], features=["geneB", "geneA"])
```

Replace the example IDs with IDs from your data.
You need Python >=3.10, NumPy >=1.24, pandas >=2, and SciPy >=1.11.4.
Object and H5AD conversion need AnnData >=0.11.4.
Missing packages are not installed for you.
If another scx module is loaded, restart the kernel and check `scx.__file__`.

```sh
python scx.py from-h5ad input.h5ad sample.scx
python scx.py validate sample.scx
python scx.py to-h5ad sample.scx restored.h5ad
```

## R / existing Seurat environment

Save scx.R, then run:

```r
source("scx.R")
stopifnot(SCX_VERSION == "v0.1.0")
write_scx(seurat_object, "sample.scx")
restored <- read_scx("sample.scx")
to_tsv("sample.scx", "counts.tsv.gz")
to_tsv("sample.scx", "selected.tsv.gz", matrix="normalized",
       cells=c("cell3", "cell1"), features=c("geneB", "geneA"))
```

You need SeuratObject >=5, Matrix, jsonlite, and digest.
Assay and Assay5 are tested under SeuratObject >=5, not an older Seurat runtime.
You do not need scxSeurat, Python, reticulate, or compiled SCX helpers.
Missing packages cause an error with their names. Nothing is installed or upgraded.
Use a UTF-8 R locale for non-ASCII IDs.

## Scope and safety

- The canonical orientation is CELL × FEATURE; counts are required. Defaults:
  AnnData `layers['counts']` / `X`, Seurat selected-assay `counts` / `data`.
- Basic writes preserve counts, optional normalized expression and typed cell /
  feature metadata. Add `analysis=True` (R `TRUE`) for embeddings/cell graphs.
  Other assays/layers, raw, loadings and arbitrary analysis internals are omitted
  with warnings. Keep the original object.
- Reads/writes validate. TSV uses ordered unique IDs and includes zeros; existing
  TSV output is always refused. It is an expression view, not a metadata round trip.
- Standalone R writes new SCX directories only; `overwrite=TRUE` is refused.
  It reserves a directory, verifies copied bytes and publishes the manifest last.
  TSV is staged, exclusively created and copied with SHA-256 verification.
  Wait for a successful return before opening an R TSV. Interrupted copies can
  retain a partial output plus lock; inspect only after the writer stops.
- Python retains staged publication and explicit SCX replacement with rollback.
  Replacement is not an atomic reader snapshot. Neither path certifies network /
  cloud-sync filesystems or hostile concurrent filesystem modification.
- R validates single-member gzip CRC/length using scratch disk. Native R integer
  bounds and character-category constraints apply; unsupported dtypes fail.
- Matrix input remains sparse. Python TSV densifies one selected row; R uses CSR
  and one text row, formatting only stored entries. Full input validation and
  transient sparse conversion still need RAM; this is not out-of-core parsing or
  a fixed RAM cap. Budget scratch/staging space and space for dense TSV output.

Optional packages are also available, but are not needed for these examples.
SCX_VERSION is the script release ID.
The optional Python package is 0.1.0; the R package is 0.1.7.
The format remains 0.1-alpha and metadata schema 1.0.

Fixed release files:

- https://github.com/junha-shin/scx-format/releases/download/v0.1.0/scx.py
- https://github.com/junha-shin/scx-format/releases/download/v0.1.0/scx.R
- https://github.com/junha-shin/scx-format/releases/download/v0.1.0/SHA256SUMS

Documentation for the tag:
https://github.com/junha-shin/scx-format/blob/v0.1.0/docs/release-verification.md
