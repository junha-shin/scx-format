# R TSV rendering optimization — PBMC3k

## Result

On the same checksum-pinned **2,700 × 13,714** PBMC3k input (2,282,976 stored
nonzero counts), full R TSV exports became **3.64–4.09× faster**. Both runs used
two repeats, the same R/dependency versions, gzip defaults, selections, thread
requests, input file hashes and measurement method. The before run is the
[historical baseline](results.json); the after run is a separate session on the
same Windows machine, not a tightly controlled same-session experiment. OS cache
and background activity were not reset or controlled.

All **12 optimized paired comparisons passed** against the source sparse matrix
and Python, with **maximum absolute error 0**, including zeros and every axis ID.
Uncompressed R TSV SHA-256 hashes are identical to the baseline. File sizes of
all R exports are unchanged. Gzip byte identity is not required by the comparator.

Time is median wall time including process start, code loading, complete SCX
validation, export, gzip and protected publication. Peak is the maximum observed
worker OS peak working set over two repeats; launchers/children are followed.
Preparation and streamed correctness checks are outside the timing. API times,
all Python/R observations and sampled process-tree RSS are also in JSON.

| R export | Before s | After s | Speedup | Before peak MiB | After peak MiB | Output MiB (unchanged) |
|---|---:|---:|---:|---:|---:|---:|
| counts-full-tsv | 149.52 | 37.82 | 3.95× | 449.2 | 466.8 | 71.0067 |
| counts-full-gzip | 154.41 | 37.77 | 4.09× | 415.8 | 423.6 | 3.7978 |
| counts-subset-gzip | 19.39 | 20.12 | 0.96× | 388.5 | 410.3 | 0.0041 |
| normalized-full-tsv | 152.86 | 39.40 | 3.88× | 423.5 | 434.2 | 107.6818 |
| normalized-full-gzip | 144.21 | 39.65 | 3.64× | 376.8 | 389.2 | 5.1946 |
| normalized-subset-gzip | 20.58 | 20.58 | 1.00× | 377.1 | 377.5 | 0.0086 |

Small subset exports did not materially improve: complete input validation
remains their main cost. No speedup is claimed for those cases.

## Why the baseline was slow

A separate full-counts diagnostic used base R `Rprof` at 10 ms plus elapsed
wrappers. It is **not** one of the benchmark observations. Inclusive samples and
wrapper measurements are not additive. The saved [profile report](optimization-profile.json)
and top [before self](profile-before-self.csv) / [after self](profile-after-self.csv)
and [before inclusive](profile-before-total.csv) / [after inclusive](profile-after-total.csv)
samples retain the evidence.

| Diagnostic component | Before seconds | After seconds |
|---|---:|---:|
| Full API call (instrumented) | 140.90 | 33.73 |
| Complete input validation | 15.76 | 16.00 |
| Number formatting | 91.57 | 6.15 |
| Old row extraction + dense conversion wrapper | 21.71 | 0 calls |
| `writeLines` after forcing formatted text | 5.46 | 5.04 |

The baseline formatted all 37,027,800 output values and repeatedly extracted
rows from column-oriented storage. About 93.8% of the counts matrix is implicit
zero. Sampling attributed about 65% of the baseline diagnostic to number
formatting; the elapsed wrappers confirm that this was the largest cost.
The old dense-conversion wrapper is no longer called; CSR traversal/conversion
still has a cost in the remaining total, not literally zero CPU time.

## Change and safety

- Convert the validated selected expression matrix to existing Matrix CSR
  storage once, then access each row's stored entries directly.
- Start one selected text row with shared `"0"` strings. Format only stored
  entries using the same `%.17g` conversion and map them into requested feature
  order. Explicit stored/signed zeros are still formatted, preserving behavior.
- No new packages, compiler, Python-from-R bridge, numeric dense row or full
  dense matrix. Existing output naming, ID selection/quoting and API stay the same.
- Full manifest/schema/axis/value/gzip CRC validation, refusal of existing outputs,
  staging, exclusive creation and SHA-256-verified copy are unchanged.
- Tests cover zero rows/columns, cancellation/duplicates, signed zero, extreme
  finite values, reordered selections, every output field, no matrix densification,
  formatting only stored entries, corrupt input and interrupted copies.

## Memory limits and tradeoff

The bound remains **O(nnz) sparse input/validation + O(selected feature count)
text-row workspace**, not O(cells × features) numeric dense storage. CSR
conversion can transiently hold both sparse layouts. It does not provide a fixed
RAM cap or out-of-core MTX parsing.

This is a speed optimization, **not a peak-memory reduction**. Observed R peaks
increased by about **0.5–21.8 MiB** across these cases (at most about 5.6%). Full
counts peak was **449.2 → 466.8 MiB**. Sparse conversion and allocator/runtime
variation affect the peak; two sessions cannot isolate their individual causes.
The complete before/after peak values are retained above, including regressions.

## Reproduce and artifacts

Follow [the PBMC3k benchmark guide](README.md). Preparation verifies the public
10x archive SHA-256 and the same input file hashes are required for comparison.
Save the baseline `standalone/scx.R` from commit
`3fc185059fdef35f03da2acace95a250796c227f`; its bytes were verified against the
historical result's R source hash. The optimized renderer is in code commit
`15e226fe35f617bdc9d778aa42f4de8d80dd1b9b` and the source hashes in the reports.

```sh
# Measure an immutable baseline again if a same-session A/B study is desired:
python examples/tsv-pbmc3k/benchmark.py run --input benchmarks/work/tsv-pbmc3k/input.scx --output benchmarks/work/before --r-script baseline.R --rscript Rscript --repeats 2
python examples/tsv-pbmc3k/benchmark.py run --input benchmarks/work/tsv-pbmc3k/input.scx --output benchmarks/work/after --rscript Rscript --repeats 2
python examples/tsv-pbmc3k/compare_optimization.py benchmarks/work/before/report.json benchmarks/work/after/report.json benchmarks/work/comparison.json
```

Use fresh paths. `profile.R` reproduces the separate diagnostic; see its guide.
Reports: [all optimized observations](optimization-results.json),
[validated before/after comparison](optimization-comparison.json),
[diagnostic and source/input hashes](optimization-profile.json).
The comparator refuses missing/failed comparisons, changed input/methods/R
versions/compression, or changed uncompressed R TSV hashes. Tests inject each
failure condition. CI runs no-install export/edge tests and benchmark smoke on
Linux, Windows and macOS; the real-data performance measurements above are local.
