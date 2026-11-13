# Use SCX 0.1.0 without installation

Download [scx.py](https://github.com/junha-shin/scx-format/releases/download/v0.1.0/scx.py)
or [scx.R](https://github.com/junha-shin/scx-format/releases/download/v0.1.0/scx.R)
and [SHA256SUMS](https://github.com/junha-shin/scx-format/releases/download/v0.1.0/SHA256SUMS).
[Check the hash](release-verification.md) before loading the file.
An [offline guide](https://github.com/junha-shin/scx-format/releases/download/v0.1.0/USAGE.md) is also available.


Use your existing Scanpy or Seurat environment.
You need one SCX file, not an installed SCX package or a compiler.
The files use your existing dependencies. They do not install packages.

## Python / Scanpy

Save [scx.py](https://github.com/junha-shin/scx-format/releases/download/v0.1.0/scx.py)
beside your script or notebook. Then run:

```python
import scx
scx.write_scx(adata, "sample.scx")
restored = scx.read_scx("sample.scx")
```

You can also run the file from a shell:

```sh
python scx.py from-h5ad input.h5ad sample.scx
python scx.py validate sample.scx
python scx.py to-h5ad sample.scx restored.h5ad
```

The file includes the SCX code. You do not need to clone the repository.
It needs Python >=3.10, NumPy >=1.24, pandas >=2, SciPy >=1.11.4, and AnnData >=0.11.4.
See the [tested versions](../compatibility/README.md).
Older environments may fail. Missing dependencies are not installed for you.
If another scx module is loaded, restart the kernel and put the saved file beside the notebook.

## R / Seurat

Save scx.R and check its hash before loading it:

```r
source("scx.R")
write_scx(seurat_object, "sample.scx")
restored <- read_scx("sample.scx")
```

Save [scx.R](https://github.com/junha-shin/scx-format/releases/download/v0.1.0/scx.R)
once, then use `source("scx.R")` offline.
It needs SeuratObject >=5, Matrix, jsonlite, and digest.
If a package is missing, the file names it and stops. It does not install or upgrade packages.
You do not need scxSeurat, Python, reticulate, or compiled SCX helpers.
Assay and Assay5 are supported. Native R type limits still apply.

## Shared default and safety

`write_scx(object, path)` writes counts, optional normalized values, and typed cell/feature metadata.
Add `analysis=True` (R: `analysis=TRUE`; CLI: `--analysis`) for embeddings and cell graphs.
Counts are required. Reads and writes check types, IDs, order, and cells x features shape.
See [examples](getting-started.md).

| Guarantee / behavior | Python file | R source file |
|---|---|---|
| Schema, axes, shape, finite/count values, typed metadata validation | Yes | Yes, within native R dtype limits |
| Gzip integrity | Python reference checks | CRC + decompressed length; single-member gzip profile |
| Existing output protection | Yes | Yes; `overwrite=TRUE` always refused |
| New output publication | Staged directory rename | Atomic directory reservation; verified copies; manifest published last |
| Explicit replacement with rollback | `overwrite=True` | Optional installed package only |
| Concurrent cooperating writers | Shared exclusive lock | Same lock identity as Python/package |

R uses temporary disk space to check gzip CRC and length without a compiler.
Allow room for staging, a copied dataset, and the largest uncompressed resource.
The R file accepts single-member gzip files, such as SCX writes.
It rejects concatenated gzip files. The native R package supports more gzip forms and is faster.
Use ordinary local storage. Network drives, cloud sync, and changes by other programs are not certified.

If R stops before writing manifest.tsv, the partial folder cannot be read.
Its lock remains. Confirm the writer has stopped, then inspect or rename the partial folder.
Remove only that lock and retry to a new path.
If R stops just after writing the manifest, the folder may be valid. Check it first.
No existing dataset is moved or deleted. See [safe writes](safe-io.md).

## Reproducibility and optional packages

These links point to v0.1.0. Keep the file and checksum with your analysis.
Both files expose that tag as `SCX_VERSION`.
Do not silently replace a file in a saved analysis.
CI checks the generated files, versions, and exact release bytes.

[Packages are optional](package-install.md).
Use them for an installed CLI, reusable imports, or faster native R checks and replacement.
The stable release has Python package 0.1.0 and R package 0.1.7.
The R package's native checks and replacement differ from the no-install R file.
Historical release files are unchanged.

## Ordinary TSV output

The same files provide `to_tsv()` for counts, normalized values, ID subsets, and .tsv.gz.
No extra dependency is needed. See [TSV options and limits](tsv-export.md).
