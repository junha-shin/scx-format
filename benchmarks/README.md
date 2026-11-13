# Scaling benchmark and bottleneck fixes

Measured on 2026-11-09. [All result tables](results/results.md), [raw Python
before](results/python-before.json)/[after](results/python-after.json), [raw R
before](results/r-before.json)/[after](results/r-after.json), and
[CSV](results/summary.csv) are checked in. Large inputs/outputs stay under ignored
`benchmarks/work/`. Baseline source is commit `600dcba`; optimized source hashes
are recorded in the reports (Python 0.1.0a4, R adapter 0.1.2).

## Observed results

| Case | Before | After | Tradeoff / scope |
|---|---:|---:|---|
| Python H5AD -> SCX, 50k cells | 80.08 s | 32.68 s | 2.45x observed speedup; SCX output grows 3.9% |
| Python SCX -> H5AD, 50k | 20.65 s | 10.67 s | 1.94x observed speedup; peak memory essentially unchanged |
| Python standalone validation peak, 50k | 338.61 MiB | 227.77 MiB | 32.7% lower OS process peak |
| R Seurat -> SCX, 10k | 35.24 s / 1196.10 MiB | 19.52 s / 476.93 MiB | 1.80x observed speedup / 60.1% lower sampled peak |
| R Seurat -> SCX, 50k | not completed | 87.61 s / 971.04 MiB | Old quadratic graph check was stopped; no baseline speedup is claimed |
| R SCX -> Seurat, 50k | not completed | 33.00 s / 800.44 MiB | Improved sparse path completed and round-trip checks passed |

Memory figures are defined below. Gains are not uniform: the small Python
validation case did not speed up, Python export/import peaks barely changed,
and R import at 10k had a slightly higher sampled peak. All figures, including
those cases, are retained in the raw reports. These are descriptive observations
from one machine, not universal speed or capacity guarantees.

## Workload and measurement

- Sizes: 1,000, 10,000 and 50,000 cells; 2,000 features; 100 random count draws
  per cell. Duplicate positions are summed, giving about 98 stored features per
  cell (~5% density). Seed 42. No dense CELL x FEATURE allocation is used by the
  generators. This synthetic scaling workload is not a biological analysis.
- Counts + log1p-normalized expression, 30D PCA-shaped and 2D UMAP-shaped random
  embeddings, a directed 15-neighbor sparse graph, numeric/nullable boolean and
  categorical metadata. R also includes native generated Seurat metadata and a
  repeated nullable text column. Embeddings are synthetic; PCA/UMAP analysis
  time is not part of conversion time.
- Python and R use separate seeded generators and native numeric widths. Compare
  before/after **within** each engine; their absolute timings are not a direct
  Python-vs-R language comparison.
- Inputs are generated once, reused byte-for-byte, and their SHA-256 hashes are
  checked by the comparison script. Preparation and correctness checks occur
  outside timed conversions. Source H5AD uses float32 sparse expression; Matrix
  Market import produces float64, so uncompressed restored H5AD can be larger
  while numeric values remain identical. Native R expression is double.
- Each command runs in a fresh process, sequentially. Wall time includes startup,
  imports, input/output and the production validation gates. Python CLI workers
  call the actual `scx.cli.main`; instrumentation adds small common overhead.
  BLAS/OpenMP limits are requested as one through threadpoolctl/environment
  settings; R's Python validator receives the same setup through a benchmark-only
  site hook. Verification showed OpenBLAS at one, while SciPy Matrix Market's
  lazily loaded pool retained its default/automatic parallelism (PARALLELISM=0).
  Matrix Market threads were not explicitly limited; this is not a universal
  single-core benchmark. The same setup is used before and after.
- Python: two repeats per size/command, reported median time and largest peak.
  With two repeats the median equals their mean; no outlier rejection or formal
  confidence interval is applied. Raw timing spread remains available.
- R: one exploratory repeat per case. Memory is sampled every 10 ms as the sum
  of process-tree RSS, including R launcher/worker and the Python validator.
  Sampling may miss short peaks and shared pages may be counted more than once.
  This differs from Python's OS-reported peak process RSS (Windows peak working
  set; Linux `ru_maxrss`). Do not compare the two memory columns as identical metrics.
- Environment: Windows 11, Python 3.12.14, 4 reported cores and about 16 GB RAM.
  R 4.5.2, SeuratObject 5.4.0, Matrix 1.7-4, jsonlite 2.0.0. Full hardware and
  Python versions are in [environment.json](results/environment.json); Python
  installation pins are in [python-requirements.txt](results/python-requirements.txt).
  Disk type and OS cache state were not characterized or flushed. Inputs and
  earlier outputs may be cached; background activity caused visible timing spread.
- Output bytes count regular files in SCX directories or the final H5AD/RDS.
  H5AD and RDS use uncompressed storage in this benchmark. Profiling is a separate
  run and is excluded from all timing tables. This is an in-memory reference
  implementation, not an out-of-core or million-cell capacity qualification.

## Changes supported by measurements

1. A validated import now checks and retains resources in one pass. At 10k cells,
   the separate [profile](results/profile.json) shows matrix reads **6 -> 3**, TSV
   reads **12 -> 4**, and typed metadata decodes **4 -> 2**. There is no global cache,
   stale-file shortcut or skipped validation.
2. Standalone validation releases each matrix before opening the next, instead
   of retaining counts/previous layers through the whole validation loop.
3. Table output checks text without making a full numeric table/string-list copy.
   Python gzip defaults change from level 9 to 6. At 50k, SCX grows from 62.00 to
   64.42 MiB; restored H5AD remains 136.19 MiB. This is a measured CPU/storage tradeoff.
4. R sparse finite/nonnegative checks inspect stored values. The old graph check
   could allocate a dense-sized CELL x CELL logical intermediate. At 50k its
   logical payload alone is about 10 GB, before other temporaries. The preliminary
   baseline run was stopped before completion of that case; the retained baseline
   R report is a clean rerun at 1k/10k only. No OOM result or 50k baseline metric
   is fabricated. The improved 50k case completed with sparse checks.
5. R metadata uses native vector numeric/logical formatting and encodes repeated
   text once per unique value, while jsonlite still handles string escaping.
   Import parses one JSON array per column instead of one parser call per cell.
   Numeric formatting remains 17-digit and category/missing semantics are tested.

Both Python before/after outputs were checked outside timing against input
counts, normalized values, embeddings, graph, axes and exact typed metadata;
all numeric errors were zero in this run. Reports are in
[before verification](results/python-before-verification.json) and
[after verification](results/python-after-verification.json). Each R run also
checks counts/data, metadata, embeddings and graph against its original RDS.

## Reproduce

From the repository root, use a separate Python 3.12 environment and install the
recorded package pins plus the checkout:

```sh
python -m venv .venv-bench
# Activate the environment for your shell.
pip install -r benchmarks/results/python-requirements.txt
pip install -e . --no-deps
```

Use the recorded R environment above, or record your own versions when comparing
on another system. R only needs SeuratObject, Matrix and jsonlite for this
synthetic benchmark. Put Rscript on PATH or pass `--rscript /path/to/Rscript`.
Create an immutable baseline source snapshot without switching the working tree:

```sh
python -c "from pathlib import Path; Path('benchmarks/work').mkdir(parents=True, exist_ok=True)"
git archive --format=zip --output=benchmarks/work/baseline.zip 600dcba src adapters/seurat/R
python -m zipfile -e benchmarks/work/baseline.zip benchmarks/work/baseline
python benchmarks/scale.py run --implementation benchmarks/work/baseline/src --output benchmarks/work/python-before --label 600dcba
python benchmarks/scale.py run --implementation src --output benchmarks/work/python-after --label current
python benchmarks/check.py --output benchmarks/work/python-before
python benchmarks/check.py --output benchmarks/work/python-after
```

R before/after must also run sequentially. The baseline is deliberately limited
to 10k cells to avoid the known quadratic graph check:

```sh
python benchmarks/scale.py r-run --implementation benchmarks/work/baseline/src --adapter benchmarks/work/baseline/adapters/seurat/R/scx.R --output benchmarks/work/r-before --label 600dcba --cells 1000 10000
python benchmarks/scale.py r-run --implementation src --adapter adapters/seurat/R/scx.R --output benchmarks/work/r-after --label current
python benchmarks/summarize.py --baseline benchmarks/work/python-before/report.json --improved benchmarks/work/python-after/report.json --r-baseline benchmarks/work/r-before/report.json --r-improved benchmarks/work/r-after/report.json --output benchmarks/work/summary
```

Choose fresh output directories. Do not run competing benchmarks concurrently.
`--cells`, `--genes`, `--per-cell`, `--seed`, and `--repeats` control workload size;
changing them requires matching before/after inputs. Do not extend the old R
baseline above 10k without removing its dense graph check or constraining memory.
For diagnostic Python profiling, use a separate `run ... --profile` output:
profiled reports are explicitly rejected by the timing summarizer.

CI runs tiny Python and R benchmark smoke cases with correctness/positive-metric
checks. Performance regressions are guarded structurally by single-read import,
matrix-release and sparse-value tests; CI does not enforce hardware-dependent
latency or RSS thresholds. The full scaling measurements are a separate manual run.
