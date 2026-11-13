# First use in your existing environment

Download scx.py or scx.R from the v0.1.0 release into your working directory.
The quickstarts load PBMC3k (reusing local data when present); substitute your
own adata or Seurat object after trying the lesson. No SCX installation is needed.
The [Python notebook](../tutorials/quickstart.ipynb) and [R walkthrough](../tutorials/quickstart.md)
show input/restored shapes, save/read calls, and a small TSV preview. The detailed lessons retain their pinned alpha.7
examples for reproducibility.

| Environment | Existing dependencies | Tested scope |
|---|---|---|
| Scanpy / AnnData | Python >=3.10, NumPy >=1.24, pandas >=2, SciPy >=1.11.4, AnnData >=0.11.4 | AnnData 0.11.4, 0.12.3, 0.13.4; Linux, Windows, macOS standalone checks |
| Seurat | R, SeuratObject >=5, Matrix, jsonlite, digest | R 4.5.2, existing Assay and Assay5; Linux, Windows, macOS standalone checks |

These dependencies must already be present. An arbitrary older Scanpy or Seurat
installation is not guaranteed to contain compatible versions. SCX never installs
or upgrades dependencies. Native R SCX helpers and Python are unnecessary for R.

## Counts choices

Python defaults to raw counts in layers['counts'] and normalized values in X.
If X contains raw counts, use `write_scx(adata, "sample.scx", counts="X", normalized=None)`.
If your counts layer has another name, select it explicitly with `counts="raw_counts"`.
For counts-only objects use `normalized=None`. raw.X is never guessed to be counts.

R defaults to counts/data in the active assay. Select a different assay or exact
layer name explicitly; counts-only objects use `normalized=NULL`. Split Assay5
layers must be joined explicitly when appropriate for the analysis.

## Common errors

| Problem | Recovery |
|---|---|
| Cannot import scx or source scx.R | Save the release file in the current working folder; restart the Python kernel if another scx was loaded |
| Missing dependency | Choose your existing analysis environment with the dependencies above |
| Missing counts | Inspect adata.layers or SeuratObject::Layers(object[[assay]], search=NA), then select the actual raw-count source |
| Missing normalized layer | Select its actual name, or set normalized=None (R: NULL) for counts-only export |
| Input does not exist | Check the working folder and full input path; read_scx expects an SCX directory |
| Output already exists | Choose a new output name; existing results are protected |
| Corrupt or incomplete SCX | Restore/regenerate the named file from the original analysis object; see safe-io.md |

## Verification

CI installs only the analysis dependencies, never SCX, for standalone jobs.
It checks Assay/Assay5 and Python/R conversions on three operating systems;
installed SCX imports and R subprocess calls are forbidden in standalone checks.
The quickstart runner executes the exact PBMC3k lessons, checks complete sparse
counts and axes, verifies dataset hashes, and compares the R/Python TSV preview.
Standalone and core tests cover missing inputs and output protection. No dependency is installed by
any lesson or converter. Detailed lessons additionally test metadata, offline
reruns, corruption checks and ordered subsets.

Run `python scripts/check_tutorials.py --output .compat-results/first-use` in the
prepared development environment. Reports include both quickstart and detailed
results. See compatibility/README.md for dtype limits and examples/ for PBMC3k.
