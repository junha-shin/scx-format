# Compatibility and support matrix

The [one-file tools](../docs/install.md) need no SCX install.
They use the same adapters and metadata mappings as the packages.
R writes new outputs only and checks single-member gzip with base R and digest.
Python uses the reference safety code.
CI runs without SCX installed on Linux, Windows, and macOS.
It checks damaged files, partial copies, and both conversion directions.
The older reports below keep their recorded versions.


The optional R package supports [R-only checks and safe writes](../docs/r-only.md).
Historical alpha.7 used R package 0.1.6. Alpha.5 used R 0.1.3 and still needs Python.

This page covers format 0.1-alpha and metadata schema 1.0.
Minimum requirements include Python >=3.10 and AnnData >=0.11.4.
Meeting a minimum does not mean every dependency combination was tested.
The tables below list the tested environments.

## Versions and evidence

| AnnData | Python in CI | pandas | NumPy | SciPy | Zarr | Assay forms tested | Local PBMC3k result |
|---|---|---|---|---|---|---|---|
| 0.11.4 | 3.10 (unit) | 2.2.3 | 1.26.4 | 1.15.3 | 2.18.3 / numcodecs 0.13.1 | Unit suite | PBMC3k uses 3.12 profile below |
| 0.11.4 | 3.10 (minimum SciPy) | 2.2.3 | 1.26.4 | 1.11.4 | 2.18.3 / numcodecs 0.13.1 | Unit suite | Synthetic CI fixture |
| 0.11.4 | 3.12 (bridge) | 2.2.3 | 2.2.6 | 1.15.3 | 2.18.7 | legacy Assay + Assay5 | Both pass; matrix error 0 |
| 0.12.3 | 3.12 | 2.2.3 | 2.2.6 | 1.15.3 | 2.18.7 | legacy Assay + Assay5 | Both pass; matrix error 0 |
| 0.13.4 | 3.12 | 3.0.6 | 2.5.3 | 1.18.1 | 3.4.0 | legacy Assay + Assay5 | Both pass; matrix error 0 |

| R / Seurat combination | Verification / scope |
|---|---|
| R 4.5.2, SeuratObject 5.4.0, Matrix 1.7.4, jsonlite 2.0.0 | Local Windows Assay/Assay5 unit tests and PBMC3k bridge across all three AnnData environments |
| Seurat 5.6.0 | Generated the recorded PBMC3k normalized matrix, PCA/UMAP and NN/SNN graphs; both assay forms retain those supported values |
| Ubuntu CI R/SeuratObject | CI resolves installed R dependencies and records their exact versions in each bridge artifact |
| Seurat 3/4 runtime or SeuratObject < 5 | Not supported/tested; legacy **Assay class** support is tested under SeuratObject 5, not an old Seurat runtime |
| Other AnnData/dependency versions | Not qualified by this matrix; no guarantee is inferred from an installation satisfying the minimum |

The local pbmc-anndata-*.json reports cover 2700 cells x 13714 genes.
They record versions, source RDS hashes, errors, axes, metadata, and reverse import checks.
Both assay forms come from the same PBMC3k analysis.
Before export, conversion to legacy Assay checks that counts, normalized values,
and feature metadata are unchanged. See [the dataset source](../examples/seurat/README.md).

CI tests four Python environments, three Seurat bridges, Windows safety,
benchmarks, two R-only environments, three binary builds, and three standalone environments.
Default writes keep expression and metadata. Add analysis=True / analysis=TRUE
(CLI: --analysis) for embeddings and graphs.
Each bridge checks both assay classes, H5AD, non-square axes, directed graphs,
types, missing values, categories, and reverse import.
Small fixtures run in CI. The listed PBMC3k analysis reports are separate local runs.
CI also runs the raw PBMC3k quickstarts without rerunning biological analysis.
Reports are evidence for those runs, not general speed or biological accuracy claims.
The Python suite has 125 cases; two symlink cases may skip on Windows.
R tests cover both assay classes, validation, corruption, rollback, races, and hard stops.
Counts-only and dense/CSR/CSC Python input are tested.
See [safe writes](../docs/safe-io.md).

## Dtype support

| Metadata dtype | Python / SCX | Native R export | R import |
|---|---|---|---|
| Signed/unsigned integer widths 8/16/32/64 | Preserved; nullable Int/UInt variants retain missing mask | integer -> int32 or Int32 | signed <=32, unsigned <=16 -> integer; INT_MIN rejected; wider/uint32/uint64 rejected |
| float32/64; nullable Float32/64 | Width, values and missing mask preserved in SCX | double -> float64; NA/NaN -> null | native double; original width/nullable-wrapper provenance is not retained |
| bool / boolean | Values and missing mask preserved | logical -> bool/boolean | logical |
| object/string/str text | Empty, literal `null` and missing distinct; Python/Arrow storage backend not portable | character -> object text | character |
| categorical | Values, category dtype, order and unused categories preserved | factor/ordered factor -> character-category schema | character categories -> factor; numeric/boolean category labels rejected |
| datetime, complex, bytes, arbitrary/mixed objects | Unsupported metadata types; export rejected | Date/POSIXt, list, custom numeric/bit64 classes rejected | Unknown or nonrepresentable schema types rejected |

H5AD has its own limits. AnnData 0.11/0.12 cannot directly save object columns
with missing strings. The CLI writes them with nullable-string encoding and warns.
The AnnData setting is changed only for that write.
Values and missing masks stay intact, but object/str storage may become pandas string.
Numbers are not changed to strings; strings are not changed to categories.
In-memory SCX reads keep the declared logical type.
H5AD readers may infer another string storage type for category labels.
Category values and order remain intact. NaN payload bits and Python-specific
missing sentinels are outside the schema.

## Exchanged and omitted data

| Item | Supported mapping | Limit / omission |
|---|---|---|
| counts | Seurat FEATURE x CELL transposed to SCX CELL x FEATURE; AnnData `layers['counts']` | One exact layer, nonnegative finite values |
| normalized | Seurat `data` transposed; SCX normalized -> AnnData X | Optional; otherwise X uses counts |
| obs / var | Typed cell/feature metadata, canonical ID order | Single feature space; metadata index IDs are text |
| embedding | Named Seurat reduction -> CELL x COMPONENT -> `obsm[name]` | Component labels/Seurat keys regenerated; loadings omitted |
| cell graph | CELL x CELL -> `obsp[name]`; directed entries retained | No transpose, no graph recomputation |
| additional assays/layers, scale.data, raw | Not in profile | Explicit assay/layer selection required; unsupported objects omitted with warning |
| identities, command history, misc / uns, varm, varp, images/spatial | Not in profile | Keep source Seurat/H5AD if these internals are needed |

## Error conditions

| Trigger | Result |
|---|---|
| Existing export directory or H5AD file | Refuse by default; explicit --overwrite / overwrite=True / overwrite=TRUE enables staged replacement with rollback/backup retention |
| R adapter 0.1.4 without Python | Supported for Seurat <-> SCX; native R validation/publication is always enabled |
| R adapter 0.1.3 (alpha.5) without Python | Its original validated import/export requires Python; release assets are unchanged |
| Output lock, output symlink, unrelated directory replacement | Refuse; inspect the named lock/recovery paths or choose a new destination |
| Missing counts, split layers without an exact selected layer | Fail; explicitly JoinLayers/select a layer first |
| Missing default normalized layer | R exporter warns and writes counts-only; normalized_layer=NULL selects counts-only intentionally |
| Duplicate/empty IDs, metadata headers, reserved ID column | Fail |
| Partial/extra cell or feature coverage, reduction/graph from another assay | Fail; matching complete axes are reordered explicitly |
| Negative/nonfinite counts, nonfinite normalized/embedding/graph values | Fail |
| Unknown manifest/schema version, unsafe/missing paths, wrong dimensions or malformed typed values | Validator rejects before import |
| Numeric overflow, invalid bool/category codes, null in nonnullable integer/bool | Fail instead of coercion |
| Unsupported metadata class, missing factor level, unsafe object name | Export fails; no output directory published |
| R cannot represent dtype or Seurat rewrites feature IDs | Import fails explicitly |
| Seurat assay class requested for import | `assay_class='Assay5'` (default) or `'Assay'`; independent of global Seurat options |

## Reproduce and extend the matrix

Run from the repo root with the R dependencies listed in the adapter guide:

```sh
pip install -e '.[test]' -r compatibility/anndata-0.11.4.txt
pytest -q
Rscript adapters/seurat/scripts/compat-fixtures.R .r-tests/fixtures
python scripts/check_compatibility.py --fixtures .r-tests/fixtures --output .compat-results/0.11.4
```

Use a separate environment and new output folder for each constraints file.
For real data, run the Seurat PBMC3k example first, then:

```sh
Rscript adapters/seurat/scripts/compat-fixtures.R .r-tests/pbmc-fixtures examples/seurat/output/source.rds
python scripts/check_compatibility.py --fixtures .r-tests/pbmc-fixtures --output .compat-results/pbmc-0.11.4 --dataset PBMC3k
```

Repeat the Python command in each fixed environment.
Each run copies the R fixtures to a new output folder.
Before adding supported versions, record new results and pass CI.

to_tsv checks the same input and type rules. It exports expression only.
See [TSV options](../docs/tsv-export.md) for row buffers, sparse input memory,
subsets, and safe writes.
