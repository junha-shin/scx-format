# Export expression as a matrix TSV

Use the [no-install files](install.md). No extra dependency is needed.
TSV shows expression values. It does not keep SCX metadata or replace an SCX round trip.
The SCX format is unchanged.

## Same API in Python and R

```python
import scx  # saved standalone/scx.py
scx.to_tsv("sample.scx", "counts.tsv")
scx.to_tsv("sample.scx", "selected.tsv.gz", matrix="normalized",
           cells=["cell3", "cell1"], features=["geneB", "geneA"])
```

```r
source("scx.R")  # saved standalone/scx.R
to_tsv("sample.scx", "counts.tsv")
to_tsv("sample.scx", "selected.tsv.gz", matrix="normalized",
        cells=c("cell3", "cell1"), features=c("geneB", "geneA"))
```

| Argument | Meaning |
|---|---|
| `path` | SCX directory; all resources are validated |
| `output` | File outside the input directory; `.gz` suffix enables gzip |
| `matrix` | `"counts"` (default) or `"normalized"`; missing normalized is an error |
| `cells` | `None` / `NULL` for all, or a nonempty list/vector of unique cell IDs |
| `features` | `None` / `NULL` for all, or a nonempty list/vector of unique feature IDs |

Select by ID, not position. This avoids Python/R indexing differences.
The output follows your requested order. Unselected axes keep their SCX order.
Unknown, duplicate, or empty ID lists fail. Use a list or vector even for one ID.

The first column is cell_id. Each later column is a selected feature.
Each selected cell has a row, including all zero values.
IDs use UTF-8 and normal TSV double-quote escaping.
R needs a UTF-8 locale for non-ASCII IDs. On Windows, use
`Sys.setlocale("LC_CTYPE", ".UTF-8")` if needed. An invalid locale may cause an input error.
Numbers use up to 17 significant digits. Python keeps integer values beyond 2^53.
R's native limits still apply; unsupported MTX integers are rejected.
A feature named cell_id keeps that name. If headers repeat, read the first column by position.

```text
cell_id geneB geneA
cell3   0     7
cell1   2     0
```

This example uses spaces for display. The file uses tabs.

## Memory and file size

The exporter does not build a full dense matrix or an AnnData/Seurat object.
It checks all SCX files and keeps only the selected expression matrix.
Python makes one selected row dense at a time.
R converts to CSR and fills a text row from stored entries; it builds no dense numeric row.
Other matrices and graphs are checked one at a time, then released.

Input parsing still needs memory for sparse matrices, metadata, and checks.
There is no fixed RAM limit or out-of-core MTX reader.
Peak memory may include the selected matrix plus the largest other resource being checked.
The row buffer grows with selected feature count, not total output size.
R may briefly keep both sparse layouts during CSR conversion (O(nnz)).
Subsets reduce output and row-buffer size, but all input is still checked.
The R file also needs scratch disk for gzip checks. See [setup](install.md).

Dense TSV includes zeros and can be much larger than MTX.
Use a small subset for viewing and .tsv.gz for larger output.
R needs room for the staged TSV and its final copy. Gzip is applied while writing.

## Existing files and failures

Both APIs refuse existing outputs. There is no overwrite option.
Choose a new path outside the input SCX folder.
Bad selections, missing normalized values, and corrupt inputs produce no output.

Python stages the TSV, then creates the output only if the path is free.
R stages it, creates a new output, copies the data, and checks SHA-256.
Both use the SCX lock for the destination.
Wait until the R call succeeds before reading the TSV; its final copy is not atomic.
A failed or killed copy can leave a partial file and lock.
After the writer stops, inspect and rename the partial file, remove only its lock,
and retry to a new path. Existing files are not truncated or deleted.

CI tests counts, normalized values, subset order, gzip, quoted/Unicode IDs,
bad input, wide Python integers, row memory use, failed writes, and R/Python agreement.
No-install checks run on Linux, Windows, and macOS.
to_tsv is available in alpha.7 and stable releases. Historical alpha.6 packages do not have it.

See the [PBMC3k TSV example](../examples/tsv-pbmc3k/README.md)
for value checks and time, memory, and file-size measurements.

The [R speed study](../examples/tsv-pbmc3k/optimization.md) records profiles,
exact value/byte checks, speed changes, and peak-memory tradeoffs.
