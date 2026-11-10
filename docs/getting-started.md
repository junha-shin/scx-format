# Install and convert Seurat / AnnData

For **main-branch R adapter 0.1.4 without Python**, use the newer
[R-only installation guide](r-only.md). The version-pinned instructions below
continue to describe the original alpha.5 release and its R 0.1.3 assets.

This guide targets the first public GitHub prerelease, **v0.1.0-alpha.5**:
Python `scx-format 0.1.0a5`, R `scxSeurat 0.1.3`, SCX format `0.1-alpha` and
metadata schema `1.0`. Earlier alpha numbers were internal development builds.
This is a GitHub release; PyPI and CRAN publication are not part of this release.

## 1. Get the release and create a Python environment

Use Python **3.12** for this guide. Verified dependency profiles cover AnnData
0.11.4, 0.12.3 and 0.13.4; the default below uses 0.13.4. Download and extract
GitHub's complete tagged **Source code** archive from the [release page](https://github.com/junha-shin/scx-format/releases/tag/v0.1.0-alpha.5)
(the Python sdist asset is a Python package, not the whole repository),
or clone it:

```sh
git clone --branch v0.1.0-alpha.5 --depth 1 https://github.com/junha-shin/scx-format.git
cd scx-format
python -m venv .venv
```

Activate it in your shell:

| Shell | Command |
|---|---|
| Linux/macOS bash or zsh | `source .venv/bin/activate` |
| Windows PowerShell | `.venv\Scripts\Activate.ps1` |
| Windows Command Prompt | `.venv\Scripts\activate.bat` |

If PowerShell blocks activation, use `.venv\Scripts\python.exe` explicitly for
each Python command; there is no need to change a system execution policy.

Install the **released wheel**, with the verified dependency profile from the
tagged checkout (keep the quotation marks around the URL requirement):

```sh
python -m pip install -r compatibility/anndata-0.13.4.txt "scx-format[anndata] @ https://github.com/junha-shin/scx-format/releases/download/v0.1.0-alpha.5/scx_format-0.1.0a5-py3-none-any.whl"
python -c "import sys, scx; print(sys.executable); print(scx.__version__)"
python -m scx.cli --help
```

Record the printed **absolute Python executable path** for the R examples.
Alternatively download the wheel from the release page and install that local
file with `python -m pip install "scx_format-0.1.0a5-py3-none-any.whl[anndata]"`
plus the same constraints. To build from the tagged source instead, run
`python -m pip install -r compatibility/anndata-0.13.4.txt ".[anndata]"`.
Avoid mixing this environment with an unrelated existing analysis environment.
The attached `requirements-verified.txt` records all 26 runtime dependency
versions from the clean Python 3.12 release check. Use it in place of the five
direct profile pins when reproducing that exact dependency set; other platforms
and R dependency versions still need their own checks.

Run the no-download Python check now:

```sh
python examples/quickstart/verify.py --python-only --output quickstart-output-python
```

Success prints `"passed": true` and writes `report.json` and `commands.json`.
Use a fresh output directory on each run. This fixture is 3 cells × 4 features;
it is a transport check, not a biological analysis.

## 2. Install the R adapter

Use a current R installation compatible with SeuratObject >= 5. This release
was locally checked with R 4.5.2 / SeuratObject 5.4.0; see the full
[support matrix](../compatibility/README.md). Legacy Assay and Assay5 objects
are supported under SeuratObject 5; an old Seurat 3/4 runtime is not qualified.

In R or RStudio, install dependencies and then the **release source package**:

```r
install.packages(c("SeuratObject", "jsonlite"), repos = "https://cloud.r-project.org")
install.packages(
  "https://github.com/junha-shin/scx-format/releases/download/v0.1.0-alpha.5/scxSeurat_0.1.3.tar.gz",
  repos = NULL, type = "source")
library(scxSeurat)
packageVersion("scxSeurat")  # 0.1.3
python <- "/absolute/path/printed/by/Python"
# Windows example: "C:/path/to/scx-format/.venv/Scripts/python.exe"
```

Use a writable R user library when R asks where to install. The adapter itself
is pure R; compiling dependencies from source can require the usual R build
tools. Both import and export require Python with scx-format installed. RStudio
does not have to inherit the shell's activated environment when `python` is
passed explicitly.

## 3. Seurat → SCX → AnnData

Given your Seurat object `object`, export one assay in R:

```r
write_scx_seurat(object, "from-seurat.scx", assay = "RNA", python = python)
```

In the activated Python shell:

```sh
python -m scx.cli validate from-seurat.scx
python -m scx.cli to-h5ad from-seurat.scx from-seurat.h5ad
```

```python
import anndata as ad
adata = ad.read_h5ad("from-seurat.h5ad")
print(adata.shape)  # cells × features
```

The selected Seurat `counts` layer becomes `layers['counts']`, `data` becomes
AnnData `X`, metadata become `obs`/`var`, reductions become `obsm`, and cell
graphs become `obsp`. Expression is transposed from FEATURE × CELL at the R
boundary; graphs are CELL × CELL throughout. Missing `data` warns and exports
counts-only. Split layers must be joined explicitly or selected by exact name.
See [adapter options](../adapters/seurat/README.md) for counts-only selection.

## 4. AnnData → SCX → Seurat

The source H5AD must have `adata.layers['counts']`. Supply actual raw counts;
SCX cannot reconstruct them from log-normalized or scaled expression. `adata.X`
is stored as `normalized` when present. Export in the Python shell:

```sh
python -m scx.cli from-h5ad input.h5ad from-anndata.scx
python -m scx.cli validate from-anndata.scx
```

Then in R:

```r
object <- read_scx_seurat("from-anndata.scx", python = python, assay_class = "Assay5")
# To request a legacy assay class under SeuratObject 5:
legacy <- read_scx_seurat("from-anndata.scx", python = python, assay_class = "Assay")
```

Seurat has narrower native metadata types than pandas. For example, R import
rejects int64/uint32/uint64 and numeric categorical labels. Choose a supported
type before export only when the values fit; do not cast to int32 blindly.
[Dtype mappings and omissions](../compatibility/README.md#dtype-support) define
the supported round trip. Metadata schema dtypes are preserved by SCX, with
documented native-R and H5AD transport conversions at the boundaries.

## 5. Verify both installed packages

Ensure the R package is in the user library visible to `Rscript` and the Python
environment is activated. From the tagged checkout:

```sh
python examples/quickstart/verify.py --output quickstart-output-bridge
```

If Rscript is not on PATH, pass its executable with
`--rscript "/absolute/path/to/Rscript"`. This runs Python → SCX → H5AD,
AnnData → Seurat → AnnData and Seurat → AnnData → Seurat. It imports the
**installed** Python/R packages, exercises public CLI commands, and compares
axes, counts, normalized values, embedding, an asymmetric graph, nullable
integer/boolean/text metadata and ordered categories with an unused level.
Matrix tolerance is 1e-12; metadata comparisons are exact within documented
H5AD category-label storage inference. It needs no external dataset.

For further validation, see [release verification](release-verification.md),
the [Scanpy PBMC3k example](../examples/pbmc3k/README.md) and
[Seurat PBMC3k example](../examples/seurat/README.md). Real-data reproduction
downloads the separately attributed dataset and installs additional analysis
dependencies; the quickstart does neither.

## 6. Existing outputs and failures

Existing destinations are refused. Prefer a fresh output name. For deliberate
whole-output replacement, use `--overwrite`, Python `overwrite=True` or R
`overwrite = TRUE`. Unlisted files inside a replaced SCX folder are removed too.

| Message / symptom | Next step |
|---|---|
| `layers['counts'] is required` | Supply real raw counts in that layer; normalized values are not a substitute. |
| `No module named scx` / R Python validation fails | Use the exact Python path printed in step 1 and verify `python -m scx.cli --help`. |
| R cannot find `scxSeurat` | Install the R asset into the user library visible to this R/Rscript process. |
| output already exists | Choose a new output or explicitly opt in to replacement. |
| output is locked | Confirm the prior writer has stopped, then follow manual recovery below. |
| missing file, gzip CRC/EOF, invalid MTX or schema | Restore/regenerate the named source file with its matching tables/schema. |
| unsupported metadata dtype | Check the mapping; convert explicitly only if safe or omit the column intentionally. |
| permission denied / no disk space | Close software holding files, check write access and make space for staging. |

See [safe I/O and manual recovery](safe-io.md) for retained backups, stale locks,
directory replacement's brief name gap, Ctrl-C and hard termination. Do not
delete lock/staging paths while a writer is active. Ordinary conversion does not
require network access after packages and inputs are installed.
