# Seurat -> SCX -> AnnData PBMC3k validation

Main's R adapter 0.1.4 also passed a **Python-free RDS -> SCX -> R**
round trip using this example's saved source object. See the
[recorded native report](r-native-report.json): 2,700 cells x 13,714 features,
exact metadata and zero errors for counts, normalized data, PCA/UMAP and NN/SNN
graphs. Python was hidden and external processes were forbidden. Reproduce
after installing 0.1.4 with `Rscript adapters/seurat/scripts/validate-rds.R
examples/seurat/output-verified/source.rds .r-tests/native-pbmc-check`, substituting
The location of your generated `source.rds`. See [R-only setup](../../docs/r-only.md).

This example begins with the official 10x PBMC3k raw count archive, creates a
real Seurat object, and runs Seurat analysis before exporting to SCX. Independent
native Seurat reference files are compared against an actual restored H5AD.
AnnData is then exported back to SCX and read into Seurat for a reverse check.

## Input and repeat the check

- [Dataset context and 10x source links](https://scanpy.readthedocs.io/en/stable/generated/scanpy.datasets.pbmc3k.html)
- Archive: https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
- SHA-256: `847d6ebd9a1ec9a768f2be7e40ca42cbfe75ebeb6d76a4c24167041699dc28b5`
- Original: 2,700 cells x 32,738 genes. Ensembl gene IDs are used to avoid gene-symbol duplication or Seurat name rewriting.

From the repository root, with R and Python available:

```r
install.packages(c("Seurat", "jsonlite", "data.table", "testthat"))
```

```sh
pip install -e '.[test]'
python examples/seurat/validate_roundtrip.py --rscript /path/to/Rscript
```

The script downloads and verifies the archive. For offline use, pass the same
archive with `--data /path/to/pbmc3k.tar.gz`. Use a fresh `--output` directory for
Each run. The default is `examples/seurat/output`. Raw data, source RDS, matrices
And generated SCX/H5AD files stay local and are excluded from Git. The dataset
keeps its upstream terms.
The repository MIT license covers adapter code.

## Analysis and comparisons

Seurat keeps genes expressed in at least three cells, runs LogNormalize
(scale factor 10,000), selects 1,000 variable genes, scales these genes, computes
30 PCs, 15-neighbor NN/SNN graphs and UMAP with seed 42. Derived QC annotations
exercise integer, double, logical, ordered factor with an unused level, character
empty/literal-null strings, and missing values. These are not biological labels.

R writes independent reference expression matrices in native FEATURE x CELL
orientation. Python compares their transpose to AnnData counts and X, compares
PCA/UMAP and both graphs, and checks cell/feature IDs and order. Counts need
zero error.
floating matrices need maximum absolute error <= 1e-12. Metadata
values, native dtype mappings, factor levels/order and missing masks are checked
exactly, both in memory and after H5AD saving. Pandas may infer string
storage for categorical labels on H5AD read.
category values/order are kept.

The reverse check compares the restored Seurat object with the saved original
RDS for the same supported data. Generated Seurat reduction component keys and
analysis internals outside the adapter profile are not compared.

For the versioned Assay/Assay5 matrix, see [compatibility results](../../compatibility/README.md).
Across AnnData versions, missing character metadata are stored in H5AD using
nullable strings with an explicit mask.
The native R character values and NA
positions are kept on reverse import.

## Recorded results and environment

The 2026-11-09 validation used Seurat 5.6.0 / SeuratObject 5.4.0 with R 4.5.2
And Python 3.12.14. It passed for 2,700 cells x 13,714 kept genes. Counts,
Normalized expression, PCA, UMAP and NN/SNN graphs all had maximum absolute
error **0.0**. Cell/feature order and typed metadata matched, and the reverse
AnnData -> SCX -> Seurat comparison passed.

See `report.json` for measured errors and R/Python versions. `R-packages.tsv`
And `R-session.txt` record the actual R environment; `python-requirements.txt`
records Python dependency versions. UMAP/PCA can vary across numeric libraries
And package versions, so each run compares its own original analysis with its
restored copy instead of expecting identical coordinates across machines.

Issues detected and fixed during validation were Matrix's automatic
pattern/symmetric graph formats, jsonlite's default double rounding, and
AnnData's automatic string-to-categorical conversion on H5AD save. The adapter
writes numeric general graphs and 17-digit floating metadata, and the CLI
preserves string columns when writing H5AD. Regression tests cover these cases.

CI runs the R unit tests plus a small independent Seurat/AnnData bridge and the
reverse import on every push. This full PBMC3k analysis is a separate explicit
example.
this analysis is not rerun by ordinary CI. The separate raw-count quickstarts
download PBMC3k and test transport only.
