# Use SCX without installation

Use the Scanpy or Seurat environment you already have. No SCX package, compiler,
Python-from-R bridge or dependency installer is used. Downloaded scripts are
SCX code, but they are ordinary files rather than installed packages.

## Python / Scanpy

Save [scx.py](https://raw.githubusercontent.com/junha-shin/scx-format/main/standalone/scx.py)
next to your script or notebook. In that same environment:

```python
import scx
scx.write_scx(adata, "sample.scx")
restored = scx.read_scx("sample.scx")
```

Or run the file directly:

```sh
python scx.py from-h5ad input.h5ad sample.scx
python scx.py validate sample.scx
python scx.py to-h5ad sample.scx restored.h5ad
```

Only this one file is needed; the repository and an installed `scx-format` are
unnecessary. It bundles the reference Python code and imports NumPy, pandas,
SciPy and AnnData already present in your analysis environment. Tested dependency
profiles are in [compatibility](../compatibility/README.md): Python >=3.10,
pandas >=2, SciPy >=1.11.4, NumPy >=1.24 and AnnData >=0.11.4. Older environments
are not promised to work. Missing dependencies fail; nothing installs itself.
If another `scx` module is already imported, restart the kernel and put the
saved file in the notebook directory before importing it.

## R / Seurat

```r
source("https://raw.githubusercontent.com/junha-shin/scx-format/main/standalone/scx.R")
write_scx(seurat_object, "sample.scx")
restored <- read_scx("sample.scx")
```

For offline/repeated use, save [scx.R](https://raw.githubusercontent.com/junha-shin/scx-format/main/standalone/scx.R)
once and use `source("scx.R")`. No other repository files are needed. This uses
SeuratObject >=5, Matrix, jsonlite and digest from the existing Seurat environment.
If any are absent the script names them and stops; it never installs or upgrades
packages. `scxSeurat`, Python, reticulate and compiled SCX helpers are unnecessary.
The same Assay and Assay5 forms and native R metadata limits apply.

## Shared default and safety

`write_scx(object, path)` writes counts, optional normalized expression and typed
cell/feature metadata. `analysis=True` (R: `analysis=TRUE`; CLI: `--analysis`)
adds embeddings/cell graphs. Counts are required; CELL × FEATURE orientation,
canonical axis order and supported metadata types are validated on read/write.
See [usage examples](getting-started.md).

| Guarantee / behavior | Python file | R source file |
|---|---|---|
| Schema, axes, shape, finite/count values, typed metadata validation | Yes | Yes, within native R dtype limits |
| Gzip integrity | Python reference checks | CRC + decompressed length; single-member gzip profile |
| Existing output protection | Yes | Yes; `overwrite=TRUE` always refused |
| New output publication | Staged directory rename | Atomic directory reservation; verified copies; manifest published last |
| Explicit replacement with rollback | `overwrite=True` | Optional installed package only |
| Concurrent cooperating writers | Shared exclusive lock | Same lock identity as Python/package |

R streams each gzip member to temporary scratch disk to verify CRC without a
compiler. Budget space for the staged dataset, its copy and the largest
uncompressed gzip resource. It rejects ordinary concatenated gzip files; use
single-member files written by SCX. Native package validation remains faster and
supports the wider gzip profile. Scripts target ordinary local filesystems;
network/cloud-sync filesystems and hostile concurrent filesystem modification
are not certified.

If R stops during publication, its incomplete folder has no `manifest.tsv` and
cannot be imported. Its lock is retained. After confirming the writer stopped,
inspect/rename the partial folder and remove only the identified lock, then write
to a fresh destination. A stop just after manifest publication can leave a valid
folder plus a lock: validate before recovery. No existing dataset is moved or
deleted. See also [safe I/O](safe-io.md).

## Reproducibility and optional packages

The links above follow main. Save the file with your workflow; for an immutable
copy, replace `main` in the URL with the full Git commit SHA you verified. Never
silently refresh a script in a reproducible analysis. The two bundled scripts
are generated from checked-in sources and CI rejects stale generated files.

Package installation is [optional](package-install.md): useful for an installed
CLI, reusable imports, or R's faster native validator and staged replacement.
The alpha.6 packages are unchanged and do not contain the standalone files.
