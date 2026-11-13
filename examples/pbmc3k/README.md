# PBMC3k real-data validation

This example validates the supported SCX alpha profile using real single-cell
RNA-seq counts, not synthetic expression data. It exercises H5AD -> SCX -> H5AD
through all three CLI commands and fails on any comparison error.

## Source and source details

- Dataset: 10x Genomics 3k PBMCs from a healthy donor, as distributed by Scanpy.
- [Scanpy dataset documentation](https://scanpy.readthedocs.io/en/stable/generated/scanpy.datasets.pbmc3k.html)
- Download: https://falexwolf.de/data/pbmc3k_raw.h5ad
- SHA-256: `89a96f1beaa2dd83a687666d3f19a4513ac27a2a2d12581fcd77afed7ea653a1`
- Original shape: 2,700 cells x 32,738 genes.

Raw and generated datasets are downloaded locally and excluded from Git. The
repository MIT license covers our code, not the upstream dataset.
consult the
linked dataset documentation and 10x source for upstream terms and attribution.

## Reproduce

From the repository root, using Python 3.12:

```sh
python -m venv .venv
# Activate the virtual environment for your shell.
python -m pip install -r examples/pbmc3k/requirements.txt
python -m pip install -e . --no-deps
python examples/pbmc3k/validate_roundtrip.py
```

The input is downloaded automatically, checked against its pinned SHA-256, and
cached under `examples/pbmc3k/data/`. For offline use, supply the same file with
`--data /path/to/pbmc3k_raw.h5ad`. Each run needs a fresh `--output` directory.
The script writes input/restored H5AD files, an SCX directory, environment pins,
And `report.json`. Dependencies are pinned from the recorded validation
environment.
different operating systems and numerical libraries can produce
different PCA/UMAP coordinates. This test compares each run's analysis against
its own restored copy, not against reference coordinates.

## Analysis and assertions

1. Keep genes detected in at least three cells.
keep all 2,700 cells.
2. Preserve raw counts in `layers['counts']`.
normalize to 10,000 per cell and
   apply log1p. No scaled matrix is presented as raw counts.
3. Select 1,000 highly variable genes, compute 30 PCs, a 15-neighbor graph and
   UMAP, using random seed 0.
4. Add derived QC metadata: numeric totals, boolean, nullable integer/boolean,
   and ordered library-size categories with an unused category and missing value.
   These are QC annotations, not inferred biological cell types.
5. Save H5AD, invoke `from-h5ad`, `validate`, and `to-h5ad`, then read both H5ADs.
6. Assert exact cell/feature order and metadata values, dtypes, categorical order,
   unused levels and missing masks with `assert_frame_equal(check_exact=True)`.
compare counts,
   normalized expression, PCA, UMAP, connectivities and distances with maximum
   absolute error <= 1e-6. Sparse comparisons do not densify expression matrices.

## Scope and issues found

The checked-in `report.json` records the actual results and dependency versions.
The final Windows/Python 3.12.14 run on 2026-11-09 passed all comparisons for
2,700 cells x 13,714 kept genes (2,282,976 nonzero counts). Counts,
Normalized expression, PCA, UMAP, connectivities and distances all had maximum
absolute error **0.0**. Cell/feature order and typed metadata matched exactly.

This is one real dataset and one analysis pipeline, not a general compatibility
or performance benchmark. `uns` analysis settings and `varm` PCA loadings remain
outside the supported alpha profile. Metadata dtypes are now preserved through
[schema 1.0](../../specification/metadata-schema-v1.md). AnnData may itself convert
repeated text columns to categorical on H5AD save, so the baseline is the source
H5AD after saving.

Validation work found and fixed:

- Missing values in categorical and nullable numeric/boolean metadata previously
  raised an exception. Schema 1.0 now preserves their dtype and missing mask,
  distinguishing empty strings from missing values. Regression tests cover this.
- Default float32 text formatting introduced errors up to about 4.9e-7
  in the initial run. Matrix and embedding exports now promote floating data to
  float64 before writing sufficient digits. This also avoids malformed output
  measured when requesting 17-digit Matrix Market formatting directly on
  float32 in the installed SciPy writer. A float32 regression test covers it.

Small regression tests run in ordinary CI. The separate quickstarts test raw
PBMC3k transport in CI. Run this example to repeat the full biological analysis.
