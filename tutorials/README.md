# Try SCX with PBMC3k

Use your existing Scanpy or Seurat environment. No SCX package install is needed.

- [Python notebook](quickstart.ipynb): code, matrix sizes, and real results.
- [R guide](quickstart.md): code and real results on one page.
- [R script](quickstart.R) | [R Markdown source](quickstart.Rmd)

## Run the example

Run from the repo root, tutorials/, or your own folder.
The example uses local PBMC3k data if it finds them. Otherwise, it downloads them once.
The Python file is about 6 MB. The R data file is about 8 MB.
Raw data are not included in Git.

Both examples save and read **2700 cells x 32738 genes**.
They show a **5 x 5** TSV table with the same values in Python and R.
These are raw counts. The example does not normalize them.

Each full run makes a new `scx-pbmc3k-*` folder. Earlier results are kept.
To repeat a write step, first make a new output folder.
For your own data, replace the loading step with your object.
Keep your original analysis object.

You can read the R guide or run its code without extra notebook tools.
Use R Markdown only if you already have it.

[Requirements and common errors](../docs/first-use.md).
Dataset sources: [Python](../examples/pbmc3k/README.md) | [R](../examples/seurat/README.md).
The dataset has its own terms of use.

## Detailed walkthrough

Use [detailed.ipynb](detailed.ipynb) or [detailed.R](detailed.R) for more options.
They cover checksums, metadata types, normalized values, subsets, and common errors.
They use small sample data and fixed alpha.7 files so the checks can be repeated.

Run the Python notebook from top to bottom in your analysis kernel.
In RStudio, run each section or use `source("detailed.R")`.
No SCX package is installed.

Outputs go under `scx-tutorial/` in your working folder.
Each full run makes a new output folder. Earlier results are kept.
After the first run, the saved files work offline.
For a first run offline, copy the alpha.7 script and SHA256SUMS into
`scx-tutorial/files/`. See [the detailed notebook](detailed.ipynb) for the steps.

The own-data section is off at first. Try the sample data before turning it on.
Check your counts source and input path first.
[See counts choices and errors](../docs/first-use.md).

## Check the examples

[Saved results](verification.json) cover PBMC3k and the detailed examples.
The checks compare counts, cell and gene order, metadata, and TSV values.
CI runs them on Linux, Windows and macOS.

For development, use an environment with the Python and R requirements:

```sh
python -m pip install -r tutorials/requirements-test.txt
python scripts/check_tutorials.py --output .compat-results/tutorials
```

These extra tools run the tests. They are not needed to use SCX.
Choose a new output folder for each test run.
