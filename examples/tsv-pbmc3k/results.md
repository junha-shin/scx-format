# Recorded PBMC3k TSV results

Measured on Windows x86_64 using the standalone scripts identified by SHA-256
in [results.json](results.json). Input: **2,700 × 13,714**, 2,282,976 nonzero
counts. Six cases × two languages × two repeats = **24 isolated export runs**.
Time below is the median of two runs; memory is the largest observed worker
OS peak working set across the repeats. See [methods and reproduction](README.md).

## Equivalence

All **12 paired comparisons passed**. Every header, ID, row/column count and
numeric field (including zeros) was compared against the source sparse matrix
and between languages. Counts and normalized both had **maximum absolute error
0** in this run. Normalized acceptance tolerance was 1e-12, not a promise of
bitwise identity on every platform. Each full-table comparison checked
37,027,800 numeric fields; each reversed-order subset checked 100 × 200.

## Time, peak memory and file size

| Case | Engine | Wall s (median) | API s (median) | Peak MiB (max) | File MiB |
|---|---|---:|---:|---:|---:|
| counts-full-tsv | python | 31.950 | 30.224 | 170.8 | 71.004 |
| counts-full-tsv | R | 149.516 | 144.805 | 449.2 | 71.007 |
| counts-full-gzip | python | 54.693 | 52.752 | 170.8 | 3.176 |
| counts-full-gzip | R | 154.411 | 151.105 | 415.8 | 3.798 |
| counts-subset-gzip | python | 2.674 | 0.952 | 170.9 | 0.004 |
| counts-subset-gzip | R | 19.387 | 16.250 | 388.5 | 0.004 |
| normalized-full-tsv | python | 34.359 | 32.583 | 144.3 | 107.679 |
| normalized-full-tsv | R | 152.865 | 149.615 | 423.5 | 107.682 |
| normalized-full-gzip | python | 43.615 | 41.963 | 144.5 | 4.881 |
| normalized-full-gzip | R | 144.211 | 140.945 | 376.8 | 5.195 |
| normalized-subset-gzip | python | 2.643 | 0.928 | 144.4 | 0.009 |
| normalized-subset-gzip | R | 20.581 | 17.390 | 377.1 | 0.009 |

## Interpretation

- R and Python transported the same values. Plain TSV byte hashes differ on
  Windows because R writes CRLF and Python writes LF. Gzip headers/default
  compression and floating-point text formatting can also differ; byte equality
  is not the correctness criterion.
- Python full plain TSV took about 32–34 seconds; R took about 150–153 seconds
  on this machine. These measurements identify a reason to profile R's large
  output path; they do not isolate sparse row extraction, formatting or I/O as
  the cause. No exporter optimization was mixed into this measurement.
- Subsets reduced output to a few KiB and wall time to about 2.6–2.7 seconds in
  Python / 19–21 seconds in R. Peak memory remained substantial because the
  complete sparse input is validated before selecting output rows/columns.
- Full counts TSV was about 71 MiB and normalized about 108 MiB. Gzip reduced
  these to roughly 3–5 MiB. Python default gzip level was 9; R default was 6,
  so compressed time/size is not a controlled equal-compression comparison.
- Preparation and streamed comparison are excluded from timing. Imports/source
  are included in wall time but excluded from API time. The unmeasured controller
  retained reference matrices; its memory is not included in worker peaks.
- The source was read before timing and the OS cache was not reset. Two local
  repeats do not characterize cold-cache, all machines or all software versions.
  Runtime versions, hardware, source/input hashes, every observation, and exact
  subset IDs are retained in the JSON. No raw TSVs or expression files are committed.

## Validation evidence

The reproducible runner and comparator passed the three-platform CI smoke
checks in [run 42893730427](https://github.com/junha-shin/scx-format/actions/runs/42893730427).
Those CI checks use a small synthetic fixture, independently of this measured
real-data run. Comparator tests reject matching-but-wrong values, misaligned
axes, truncation, extra rows and nonfinite values.
