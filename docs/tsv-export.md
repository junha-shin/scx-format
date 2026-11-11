# Export expression as a matrix TSV

Use the current [no-install scripts](install.md). No new dependencies are needed.
The TSV is an expression view, not a replacement for SCX metadata/schema or an
object round-trip format. The SCX format itself is unchanged.

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

Selections use IDs, not numeric positions, so R/Python indexing conventions do
not differ. Output follows the requested order; unspecified axes retain SCX order.
Unknown, duplicate or empty selections fail. Use a list/vector even for one ID.

The first column is `cell_id`, followed by one column per selected feature ID.
Every selected cell has one row; **all zero values are included**. UTF-8 IDs are
preserved with ordinary TSV double-quote escaping. R needs a UTF-8 locale for
non-ASCII IDs, as in the SCX reader; an invalid inherited locale can cause an
explicit input error. On Windows, `Sys.setlocale("LC_CTYPE", ".UTF-8")` selects
the native UTF-8 locale. Numeric output uses up to 17
significant digits; Python integer matrices retain integer values beyond 2^53.
Native R numeric/integer limitations still apply and unsupported integer MTX
values are rejected. A feature literally named `cell_id` retains that label;
interpret the first column by position if this creates duplicate header names.

```text
cell_id geneB geneA
cell3   0     7
cell1   2     0
```

(The example above is visually spaced; actual fields are separated by tabs.)

## Memory and file size

The exporter never builds a dense CELL × FEATURE array or an AnnData/Seurat
object. It validates all SCX resources, retains only the selected expression
matrix, and builds **one selected output row at a time**. Python densifies only
that row. R converts to CSR and formats stored entries into a zero-filled text
row, without a numeric dense matrix or dense row. Other expression matrices and
graphs are validated individually and discarded.

This is not out-of-core MTX parsing or a fixed RAM cap: sparse matrix input,
metadata and validation still need memory. Peak memory can include the selected
sparse matrix plus the largest other resource being validated. Row conversion
and formatting use space proportional to selected feature count, not the total
output cell × feature count. R CSR conversion can transiently hold both sparse
layouts (O(nnz)); it does not change the bounded row-buffer contract. Subsets
reduce output/row-buffer size; they do not
avoid full SCX input validation. Existing R standalone gzip validation also
uses scratch disk, as documented in [setup](install.md).

Dense TSV can be much larger than MTX because it includes zeros. Prefer a small
selection for inspection and `.tsv.gz` for larger outputs. Allow space for a
complete staged TSV plus the published file in R; gzip is applied during writing.

## Existing files and failures

**Both APIs always refuse an existing output.** There is no overwrite option;
choose a new path. Writing inside the input SCX directory is also refused.
Invalid selections, absent normalized data and corrupt inputs produce no output.

Python stages the TSV and publishes via exclusive file creation after successful
writing, retaining its reference locking and cleanup behavior. R stages the TSV,
then exclusively creates the output and verifies the copied bytes with SHA-256.
Cooperating writers share the SCX output lock identity. R's final copy is not an
atomic reader snapshot: wait for successful return before opening the TSV.
If that copy fails or is killed, a partial file and lock can remain; inspect only
after the writer stops, rename the partial file and remove the identified lock,
then retry to a new path. Existing files are never truncated or deleted.

CI tests counts/normalized, subsets/order, gzip, Unicode/quoted IDs, missing and
invalid inputs, exact wide Python integers, Python single-row densification, R no-dense sparse rendering, failed
publication and cross-language TSV agreement. No-install checks run on Linux,
Windows and macOS. Optional installed-package APIs receive the same addition on
main; immutable alpha.6 packages do not yet expose `to_tsv`.

For real-data equivalence checks and isolated-process time, memory and file-size
measurements, see the [reproducible PBMC3k TSV example](../examples/tsv-pbmc3k/README.md).

The [R optimization study](../examples/tsv-pbmc3k/optimization.md) records
profiling, exact-value/byte checks, speedups and peak-memory tradeoffs on PBMC3k.
