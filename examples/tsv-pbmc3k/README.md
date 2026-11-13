# PBMC3k TSV checks and speed

Check that Python and R write the same expression values. Then measure time,
memory, and file size. These steps need no SCX package install.

This example uses the same public 10x PBMC3k hg19 archive as the existing
[Seurat real-data example](../seurat/README.md). It measures the current
**no-install Python and R scripts**, not an installed SCX package.

The archive SHA-256 is
`847d6ebd9a1ec9a768f2be7e40ca42cbfe75ebeb6d76a4c24167041699dc28b5`.
Preparation reads named archive members without extracting filesystem paths,
keeps genes detected in at least three cells, and uses Ensembl feature IDs.
The result is **2,700 cells x 13,714 features**, **2,282,976 nonzero counts**.
Normalized expression is float64 `log1p(counts / filtered-cell total * 10000)`.
Both exporters receive this exact same expression-only SCX directory. This is
an export/transport test, not a comparison of Scanpy and Seurat analysis methods.

## Current R optimization

[Before/after performance, memory and safety](optimization.md). Full R outputs
are about 3.6-4.1x faster on this PBMC3k input.
output values and sizes are
unchanged. Small subset validation remains the time floor. The report includes
The small measured peak-memory increases instead of claiming memory savings.

## Recorded baseline real-data results

[Measured results and interpretation](results.md) | [all measurements and hashes](results.json).
Counts and normalized matched exactly in both repeats, including all zeros and
requested subset order. The results are local Windows measurements; CI smoke
checks are separately described below.

## Reproduce

From a checkout, use your existing AnnData/scverse and Seurat environments.
The measurement controller needs `psutil` (`pip install psutil` only if absent).
There is no SCX installation. R reuses SeuratObject, Matrix, jsonlite and digest.
set `R_LIBS_USER` if those packages live in a custom library. The R worker uses
`--vanilla` and does not read an R profile.

Download the [10x archive](https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz)
(or reuse the archive from the Seurat example), then run:

```sh
python examples/tsv-pbmc3k/benchmark.py prepare --archive examples/seurat/data/pbmc3k.tar.gz --output benchmarks/work/tsv-pbmc3k/input.scx
python examples/tsv-pbmc3k/benchmark.py run --input benchmarks/work/tsv-pbmc3k/input.scx --output benchmarks/work/tsv-pbmc3k-run --repeats 2 --rscript Rscript
```

Use a fresh output directory. On Windows, `--rscript` can be the full path to
`Rscript.exe`, including spaces. Preparation refuses a wrong archive checksum.
Work files, input SCX, TSVs and worker logs remain under ignored
`benchmarks/work/`.
No downloaded expression data are committed. The public
dataset keeps its upstream terms; MIT covers this repository's code.

The runner writes `report.json` and `summary.md`. A worker or numeric/axis check
failure stops the run without a success report. Raw logs can contain local
paths.
They are not part of the checked-in results.

## What is compared

For **counts and normalized**, each engine writes full `.tsv`, full `.tsv.gz`,
And an ID-selected `.tsv.gz`. The subset uses up to 100 evenly spaced cells and
200 evenly spaced features, then reverses both axes to exercise requested order.
Exact IDs are kept in each run's case configuration and report.

Every output header, cell ID, row/column count and numeric field is checked
against the original sparse SCX matrix, including every zero. R and Python are
also compared directly. Counts must match exactly.
Normalized values have an
absolute tolerance of `1e-12`. TSVs are read one row at a time.
comparisons do not
build a full dense matrix. Gzip and platform newline/float-format differences
can produce different file hashes despite equivalent tables.

## Measurement contract

Each export runs in its own fresh process, sequentially. Two repeats reverse
R/Python order on the second repeat. Summary time is the median.
memory is the
maximum recorded peak across repeats.
every observation remains in JSON.

- **Wall time** includes process startup, code imports / R source, complete SCX
  validation, expression export, compression and protected publication.
- **API time** starts after code dependencies are loaded, before the `to_tsv`
  call. Input preparation and streamed result comparison are excluded.
- **Windows memory** is the largest worker/descendant process's OS peak working
  set, queried every 10 ms. This follows venv/R launchers to the actual worker.
  Concurrent process-tree sampled RSS is also recorded. It is not Python/R heap
  size or total system RAM. Very short processes/startup/exit transients can be
  missed by polling.
- **Linux/macOS memory** is process-tree RSS sampled every 10 ms, an measured
  lower bound instead of a guaranteed OS high-water value.
- Requested BLAS/OpenMP threads are one; SciPy Matrix Market threading remains
  automatic. Default API gzip settings are kept and recorded (Python 9,
  R 6).
The compression settings are not controlled to be identical.
- The source is read by the controller before timing and the OS cache is not
  reset. These are warm-cache local measurements. Disk, background load and
  runtime/library versions affect results.
two repeats are not a broad hardware
  or statistical performance qualification.

The exporters keep input sparse and build one output row. Python densifies one
row.
optimized R formats stored CSR entries into a zero-filled text row. Full input
validation still loads sparse resources.
subset export therefore does not imply
subset-sized input memory. R also uses scratch disk for CRC validation and
stages/copies output. Full dense TSV writes zeros and can be much larger than
MTX.
allow space for both engines and both repeats.

## CI and tests

The same runner and comparison code run on a small nonsquare synthetic fixture
on Linux, Windows and macOS, with SCX not installed. CI covers all six cases and
uploads JSON/Markdown results. It does not download PBMC3k or turn shared-runner
smoke timings into real-data performance claims:

```sh
python examples/tsv-pbmc3k/benchmark.py smoke --output benchmarks/work/tsv-smoke --rscript Rscript
python -m pytest tests/test_tsv_benchmark.py
```

Comparator unit tests reject matching-but-wrong values, reordered axes,
truncation, extra rows and nonfinite values. Existing TSV, corruption, overwrite,
limited-densification and adapter tests remain in CI.
## Profiling and optimized-run comparison

`profile.R` is a separate diagnostic run with Rprof and elapsed-time wrappers.
It measures input validation, row extraction/densification, number formatting
And writing, and saves top inclusive/self samples. Sampling/wrapper times are
not additive and are not used as performance benchmark measurements.

```sh
Rscript examples/tsv-pbmc3k/profile.R standalone/scx.R benchmarks/work/tsv-pbmc3k/input.scx benchmarks/work/profile.tsv benchmarks/work/profile
python examples/tsv-pbmc3k/compare_optimization.py examples/tsv-pbmc3k/results.json benchmarks/work/optimized/report.json benchmarks/work/comparison.json
```

The run command accepts `--r-script /path/to/saved-baseline.R` to measure an
immutable historical script without installing it. Save `standalone/scx.R` from
commit `3fc185059fdef35f03da2acace95a250796c227f` for the pre-optimization version.
Use new output paths for every diagnostic and benchmark run. The comparator
refuses changed input hashes, axes/selections, methods, R versions/dependencies
or gzip defaults, and checks uncompressed R TSV hashes against the baseline.
Historical baseline and optimized measurements are separate sessions on the
same machine.
cache and background load are not reset or controlled.
