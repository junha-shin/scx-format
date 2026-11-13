# Quickstart with PBMC3k

- [Python notebook](quickstart.ipynb): code, shapes, and actual table outputs.
- [R walkthrough](quickstart.md): rendered Markdown with code and actual outputs.
  [R Markdown source](quickstart.Rmd) and [plain R script](quickstart.R) are also available.

Use your existing Scanpy/AnnData or Seurat environment. Run from the repository
root, tutorials/, or any working folder. The lesson downloads scx.py or scx.R
only when absent; no SCX package installation is needed. Existing local files
are reused, so use the published stable file when starting in a new folder.

Python reuses examples/pbmc3k/data/pbmc3k_raw.h5ad; R reuses
examples/seurat/data/hg19 or the previously downloaded pbmc3k.tar.gz archive.
Raw data are ignored by Git and are not shipped in the repository. In a fresh
checkout, each lesson downloads its input once (about 6 MB Python, 8 MB R).
The download/cache steps are part of the runnable example.

Both lessons write/read all 2700 cells x 32738 genes and display the dimensions.
R's in-memory Seurat matrix is 32738 features x 2700 cells; SCX uses cells x
features. A 5-cell x 5-gene TSV preview uses the genes with the highest totals,
so only a tiny table is converted to dense form. Counts are raw; this transport
lesson does not perform normalization or a biological analysis.

Every full run creates a fresh scx-pbmc3k-* output folder. Rerunning only a
write block refuses an existing output; rerun the folder-creation block for a
new destination. Keep original analysis objects. To use your own data, replace
the loading section with your existing adata or Seurat object.

R Markdown rendering is optional: the published .md already contains results.
You can run the script or copy the R blocks without installing knitr, rmarkdown,
Jupyter or Python. CI renders the Markdown directly from the executed Rmd code
blocks using base R and the existing SCX dependencies.

[Counts choices and common errors](../docs/first-use.md).
Dataset provenance: [Python](../examples/pbmc3k/README.md),
[R](../examples/seurat/README.md). Upstream dataset terms apply.

## Detailed walkthrough

The previous full lessons are preserved as [detailed.ipynb](detailed.ipynb) and
[detailed.R](detailed.R). They cover verified downloads, typed metadata, normalized
TSV, subsets, offline reruns and recovery. The instructions below describe these
**detailed** lessons.

# SCX in five minutes

Choose **one** lesson for the environment you already use. Neither lesson
installs SCX or any other package. They use the published **v0.1.0-alpha.7** files;
you do not need to clone this repository.

| You use | Download | Run |
|---|---|---|
| Scanpy / AnnData | [detailed.ipynb](https://raw.githubusercontent.com/junha-shin/scx-format/main/tutorials/detailed.ipynb) (save as `.ipynb`) | Open in your existing Jupyter environment; select your Scanpy kernel, then **Restart Kernel and Run All Cells** |
| Seurat | [detailed.R](https://raw.githubusercontent.com/junha-shin/scx-format/main/tutorials/detailed.R) (save as `.R`) | Open in RStudio; select a numbered section and Ctrl/Cmd+Enter, or run `source("detailed.R")` |

The [notebook preview](detailed.ipynb) includes executed results; you can read
it on GitHub before downloading. R uses an ordinary script so that you do not
need Jupyter, IRkernel or Python. The notebook assumes that you already have
Jupyter for your Python environment; a plain-script reference remains available
in [getting started](../docs/getting-started.md).

## What happens

1. Check your existing dependencies: Python 3.10+, AnnData, NumPy, pandas, SciPy;
   or R with SeuratObject >=5, Matrix, jsonlite and digest. A missing dependency
   stops with a message; nothing is installed automatically.
2. Download the fixed alpha.7 script and checksum file. Both are verified against
   embedded SHA-256 hashes **before** loading code.
3. Generate 3 cells and 4 features, with counts, log1p values, nullable numbers,
   boolean values and an ordered categorical column with an unused level.
4. Write SCX and read it back. Assertions check counts, normalized values,
   metadata types and identities. Expect `PASS`.
5. Export `counts.tsv`, `normalized.tsv.gz`, and `subset.tsv.gz`. See actual tables.
   Ordered subsets select cell3/cell1 and geneB/geneA, yielding rows `[1,5]`, `[0,1]`.
6. Attempt to overwrite an existing TSV. Expect refusal and unchanged bytes,
   then `SUCCESS: all demo checks passed`.

The counts table in both languages is:

| cell_id | geneA | geneB | geneC | geneD |
|---|---:|---:|---:|---:|
| cell1 | 1 | 0 | 3 | 0 |
| cell2 | 0 | 2 | 0 | 4 |
| cell3 | 5 | 1 | 0 | 0 |

SCX and TSV are **CELL x FEATURE**. Seurat's in-memory matrix is FEATURE x CELL;
the R adapter handles that transpose. Seurat also adds QC metadata, so the two
objects' complete metadata tables are not expected to be identical across languages.
The lessons check exact metadata round-trip within each environment.

## Files, reruns and offline use

Outputs are relative to your **current working folder**, not necessarily the
downloaded lesson's folder. The lessons do not change that folder:

```text
scx-tutorial/
  files/                  # verified scx.py or scx.R, and SHA256SUMS
  output/run-*/           # a fresh directory for each full run
    demo.scx/
    counts.tsv
    normalized.tsv.gz
    subset.tsv.gz
```

A full rerun keeps prior output directories. Rerunning only a write cell/section
refuses the existing destination: rerun the whole lesson, or the setup section
that creates a new run directory. No old output is deleted automatically.

First use needs internet. Later runs reuse verified cache files without network.
For first-run offline use, download the matching script and `SHA256SUMS` from the
[alpha.7 release](https://github.com/junha-shin/scx-format/releases/tag/v0.1.0-alpha.7)
on another machine and copy them to `scx-tutorial/files/` in your working folder.
Keep the lesson's pinned version and hashes unchanged.

| Message / symptom | Recovery |
|---|---|
| Missing dependency | Select your existing Scanpy kernel or Seurat environment; see [dependency list](../docs/install.md) |
| Download failed | Check internet access, or use the offline cache procedure above |
| Checksum mismatch | Stop; inspect the cached file, move it aside, and obtain the exact alpha.7 asset again; the lesson never runs or silently replaces mismatched code |
| Another standalone SCX is loaded | Restart the Python kernel and Run All |
| Destination already exists | Run setup again for a new output directory; keep the existing result |
| Missing counts / split Seurat layers in own-data example | Select the actual raw-count source explicitly; join split layers explicitly if appropriate; see [conversion guide](../docs/getting-started.md) |

## Your own data and next steps

The last example is **disabled** initially. After the demo passes, provide your
actual H5AD/RDS path and enable its switch. Inspect available layers and choose
raw counts explicitly. For raw counts in AnnData `X`, use `counts="X"` and
`normalized=None`; normalized-only `X` cannot supply counts. `raw.X` is not
selected automatically. For Seurat counts-only data, set `normalized_layer <- NULL`.
Neither lesson overwrites the input analysis object.

Default transport preserves counts, optional normalized values and typed
cell/feature metadata. Supported embeddings/graphs require `analysis=True` in
Python or `analysis=TRUE` in R; the tiny demo does not contain these. Additional
assays and arbitrary analysis internals are omitted. Keep original objects.
Dense TSV can be very large: use subsets and gzip; do not densify real input
matrices for display. TSV is an expression table, not a metadata-preserving
replacement for SCX.

- [Seurat <-> AnnData conversion recipes](../docs/getting-started.md)
- [Real PBMC3k lesson and validation](../examples/pbmc3k/README.md)
- [PBMC3k R/Python TSV correctness and performance](../examples/tsv-pbmc3k/README.md)
- [TSV memory and format limits](../docs/tsv-export.md)

## Reproduce the tutorial checks (developers)

The [recorded local result](verification.json) passed on Windows with Python
3.12.14, AnnData 0.13.4, pandas 3.0.6, NumPy 2.5.3, SciPy 1.18.1, R 4.5.2 and
SeuratObject 5.4.0. Counts matched exactly; normalized maximum absolute error
was zero in both lessons. This tiny transport check complements the real-data
validation linked above.

In an environment with the analysis dependencies and R dependencies above:

```sh
python -m pip install -r tutorials/requirements-test.txt
python scripts/check_tutorials.py --output .compat-results/tutorials
```

These are development-only execution tools, not SCX dependencies. Use a **new**
output path each time. The runner also executes both PBMC3k quickstarts and checks full sparse counts, axis order, dataset hashes, and cross-language TSV preview values. It uses a private kernel configuration and executes
both lessons twice, forbids network on the second run, checks prior outputs and
damaged-cache rejection, and compares actual R/Python TSV values and ID order.
Reports and the executed notebook are saved under the output path. CI repeats
this on Linux, Windows and macOS; installed SCX imports and R subprocess calls
are forbidden during the lessons.
