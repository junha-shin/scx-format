# Measured scaling results

Wall time includes process startup, imports and I/O. Python: median of two runs; R: one exploratory run. Peaks: maximum observed per case. R memory is sampled process-tree RSS; Python memory is OS process peak. See README for methods and limits.

## Python CLI

| Cells | Operation | Before s | After s | Speedup | Before peak MiB | After peak MiB | Before output MiB | After output MiB |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | from-h5ad | 5.92 | 3.90 | 1.52 | 121.54 | 119.92 | 1.22 | 1.28 |
| 1,000 | to-h5ad | 3.80 | 3.18 | 1.19 | 116.27 | 114.75 | 2.87 | 2.87 |
| 1,000 | validate | 1.97 | 2.09 | 0.94 | 85.99 | 84.91 | 1.22 | 1.28 |
| 10,000 | from-h5ad | 20.86 | 8.91 | 2.34 | 189.39 | 189.90 | 12.33 | 12.88 |
| 10,000 | to-h5ad | 6.47 | 4.34 | 1.49 | 159.15 | 158.28 | 27.37 | 27.37 |
| 10,000 | validate | 3.55 | 3.15 | 1.13 | 131.64 | 111.98 | 12.33 | 12.88 |
| 50,000 | from-h5ad | 80.08 | 32.68 | 2.45 | 498.86 | 498.49 | 62.00 | 64.42 |
| 50,000 | to-h5ad | 20.65 | 10.67 | 1.94 | 366.19 | 368.25 | 136.19 | 136.19 |
| 50,000 | validate | 14.99 | 8.00 | 1.87 | 338.61 | 227.77 | 62.00 | 64.42 |

## Seurat R

| Cells | Operation | Before s | After s | Speedup | Before peak MiB | After peak MiB | Before output MiB | After output MiB |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 1,000 | seurat-export | 6.71 | 4.41 | 1.52 | 415.15 | 304.90 | 1.23 | 1.23 |
| 1,000 | seurat-import | 5.62 | 4.91 | 1.14 | 301.04 | 293.83 | 2.86 | 2.86 |
| 10,000 | seurat-export | 35.24 | 19.52 | 1.80 | 1196.10 | 476.93 | 12.43 | 12.43 |
| 10,000 | seurat-import | 10.23 | 9.28 | 1.10 | 382.49 | 405.88 | 28.21 | 28.21 |
| 50,000 | seurat-export | not measured | 87.61 | not measured | not measured | 971.04 | not measured | 63.00 |
| 50,000 | seurat-import | not measured | 33.00 | not measured | not measured | 800.44 | not measured | 141.08 |

R baseline 50,000-cell cases were not completed: the old finite check can allocate a dense 50,000 × 50,000 logical graph (~10 GB). No measured time, peak or speedup is inferred for those cases.

SCX gzip level changed from 9 to 6 in Python. File size differences are measured, not assumed; H5AD and RDS outputs use their documented uncompressed defaults.
