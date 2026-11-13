# SCX v0.1.0 — no-install files

Release: https://github.com/junha-shin/scx-format/releases/tag/v0.1.0

Download only the language file you need (`scx.py` or `scx.R`) and `SHA256SUMS`.
Use the Scanpy/Seurat environment you already have. No SCX installation,
compiler or automatic dependency installation is used. Files are canonical UTF-8/LF.

## Verify before running

Compare the computed SHA-256 with the exact file-name entry in `SHA256SUMS`.

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

If downloading only one language file, compare that file's entry individually.
The complete checksum list also covers this guide and optional package assets.
Use the release-specific URLs below, not `latest` or `main`, in saved workflows.
Keep the script and checksum alongside the analysis. A mismatch means the file
must not be imported/sourced; download the release asset again.

## Python / existing Scanpy environment

Save `scx.py` beside your notebook/script:

```python
import scx
assert scx.SCX_VERSION == "v0.1.0"
scx.write_scx(adata, "sample.scx")
restored = scx.read_scx("sample.scx")
scx.to_tsv("sample.scx", "counts.tsv.gz")
scx.to_tsv("sample.scx", "selected.tsv.gz", matrix="normalized",
           cells=["cell3", "cell1"], features=["geneB", "geneA"])
```

The selection example uses placeholder IDs: substitute IDs actually in your data.
Existing NumPy >=1.24, pandas >=2, SciPy >=1.11.4 and Python >=3.10 are needed.
Object/H5AD conversion also needs AnnData >=0.11.4. No missing package is installed
automatically. If a different `scx` is already imported, restart the kernel and
check `scx.__file__` points to this saved file.

```sh
python scx.py from-h5ad input.h5ad sample.scx
python scx.py validate sample.scx
python scx.py to-h5ad sample.scx restored.h5ad
```

## R / existing Seurat environment

Save `scx.R`, then:

```r
source("scx.R")
stopifnot(SCX_VERSION == "v0.1.0")
write_scx(seurat_object, "sample.scx")
restored <- read_scx("sample.scx")
to_tsv("sample.scx", "counts.tsv.gz")
to_tsv("sample.scx", "selected.tsv.gz", matrix="normalized",
       cells=c("cell3", "cell1"), features=c("geneB", "geneA"))
```

SeuratObject >=5, Matrix, jsonlite and digest must already be present. Assay and
Assay5 under SeuratObject >=5 are supported; an old Seurat runtime is not implied.
No scxSeurat package, Python, reticulate or compiled SCX helper is needed.
Missing dependencies stop with their names; nothing is installed/upgraded.
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

Optional installed packages in this release are for users who prefer them;
they are not required for these examples. Script release version is
`SCX_VERSION`; optional package versions are Python 0.1.0 / R 0.1.7.
Format remains 0.1-alpha and metadata schema 1.0.

Pinned files:

- https://github.com/junha-shin/scx-format/releases/download/v0.1.0/scx.py
- https://github.com/junha-shin/scx-format/releases/download/v0.1.0/scx.R
- https://github.com/junha-shin/scx-format/releases/download/v0.1.0/SHA256SUMS

Tagged documentation:
https://github.com/junha-shin/scx-format/blob/v0.1.0/docs/release-verification.md
