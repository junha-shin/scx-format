# Compatibility and support matrix

The default [standalone scripts](../docs/install.md) require no SCX installation.
They reuse the same Python/R adapters and metadata mappings. R standalone writes
new outputs only and validates single-member gzip via base R + existing digest;
Python retains reference safety. Dedicated Linux/Windows/macOS CI runs with SCX
absent, tests corruption/partial copies, and checks both cross-language directions.
Existing package/real-data reports below retain their original versions.


Alpha.6 R adapter 0.1.5 provides [R-only validation and safe
publication](../docs/r-only.md), independently of Python. The alpha.5 release
contains R 0.1.3 and retains its original requirements.

This document defines the tested SCX v0.1-alpha / metadata-schema 1.0 profile.
The Python package declares Python >= 3.10 and AnnData >= 0.11.4; that minimum
does not imply every newer release and dependency combination has been tested.
Exact environment constraints below are the reproducible supported test profiles.

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

`pbmc-anndata-*.json` are measured local reports for 2,700 cells × 13,714 genes.
They include dependency versions, source RDS hashes, matrix errors, axis order,
metadata checks and the reverse Seurat assay-class check. Both forms are exported
from the same genuine PBMC3k Seurat analysis; conversion to legacy Assay first
asserts that counts, normalized values and feature metadata have not changed.
The original 10x archive provenance is in [the real-data example](../examples/seurat/README.md).

CI has four Python environment jobs, three bridge jobs, a Windows safety job
and a benchmark smoke job, two R-only jobs and three binary-distribution jobs.
Default public writes preserve expression and typed metadata; embeddings/graphs
require `analysis=True` / `analysis=TRUE` (CLI `--analysis`). Read scope and
dtype safety are unchanged. Recorded PBMC3k reports predate the new wrapper API.
Every bridge job
tests both assay classes, H5AD serialization, CELL × FEATURE orientation on a
nonsquare matrix, an asymmetric graph, dtype/missing/category preservation and
AnnData → SCX → Seurat import into the requested assay class. Small fixtures run
in CI; PBMC3k reports record the separate local real-data validation. CI uploads
JUnit and per-environment JSON artifacts. Recorded success is evidence for the
listed runs, not a performance or general biological correctness claim.
The Python suite includes 109 cases per AnnData environment (two symlink cases
may skip on Windows without link-creation privileges); R tests cover both assay
classes plus native validation, corruption, rollback, races and hard interruption.
Counts-only selection and dense/CSR/CSC
Python input are included. Safety cases cover failed writes, interrupted processes,
writer races, rollback and corrupt transport; see [safe I/O](../docs/safe-io.md).

## Dtype support

| Metadata dtype | Python / SCX | Native R export | R import |
|---|---|---|---|
| Signed/unsigned integer widths 8/16/32/64 | Preserved; nullable Int/UInt variants retain missing mask | integer → int32 or Int32 | signed <=32, unsigned <=16 → integer; INT_MIN rejected; wider/uint32/uint64 rejected |
| float32/64; nullable Float32/64 | Width, values and missing mask preserved in SCX | double → float64; NA/NaN → null | native double; original width/nullable-wrapper provenance is not retained |
| bool / boolean | Values and missing mask preserved | logical → bool/boolean | logical |
| object/string/str text | Empty, literal `null` and missing distinct; Python/Arrow storage backend not portable | character → object text | character |
| categorical | Values, category dtype, order and unused categories preserved | factor/ordered factor → character-category schema | character categories → factor; numeric/boolean category labels rejected |
| datetime, complex, bytes, arbitrary/mixed objects | Unsupported metadata types; export rejected | Date/POSIXt, list, custom numeric/bit64 classes rejected | Unknown or nonrepresentable schema types rejected |

H5AD is a separate transport constraint: AnnData 0.11/0.12 cannot serialize an
object column with missing strings directly. The CLI stores missing text using
the native H5AD nullable-string encoding, enables the AnnData setting only within
the write, and warns about the storage dtype conversion. Values and missing masks
are preserved; pandas object/str storage may become pandas `string` after H5AD.
No numeric-to-string conversion or string-to-categorical conversion is applied.
In-memory SCX import retains the declared logical schema dtype. H5AD readers may
also infer a different string storage dtype for category labels, while retaining
category values/order. NaN payload bits and distinct Python missing sentinels are
not part of the schema.

## Exchanged and omitted data

| Item | Supported mapping | Limit / omission |
|---|---|---|
| counts | Seurat FEATURE × CELL transposed to SCX CELL × FEATURE; AnnData `layers['counts']` | One exact layer, nonnegative finite values |
| normalized | Seurat `data` transposed; SCX normalized → AnnData X | Optional; otherwise X uses counts |
| obs / var | Typed cell/feature metadata, canonical ID order | Single feature space; metadata index IDs are text |
| embedding | Named Seurat reduction → CELL × COMPONENT → `obsm[name]` | Component labels/Seurat keys regenerated; loadings omitted |
| cell graph | CELL × CELL → `obsp[name]`; directed entries retained | No transpose, no graph recomputation |
| additional assays/layers, scale.data, raw | Not in profile | Explicit assay/layer selection required; unsupported objects omitted with warning |
| identities, command history, misc / uns, varm, varp, images/spatial | Not in profile | Keep source Seurat/H5AD if these internals are needed |

## Error conditions

| Trigger | Result |
|---|---|
| Existing export directory or H5AD file | Refuse by default; explicit --overwrite / overwrite=True / overwrite=TRUE enables staged replacement with rollback/backup retention |
| R adapter 0.1.4 without Python | Supported for Seurat ↔ SCX; native R validation/publication is always enabled |
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

Run from the repository root with R dependencies from the adapter README:

```sh
pip install -e '.[test]' -r compatibility/anndata-0.11.4.txt
pytest -q
Rscript adapters/seurat/scripts/compat-fixtures.R .r-tests/fixtures
python scripts/check_compatibility.py --fixtures .r-tests/fixtures --output .compat-results/0.11.4
```

Use separate environments and fresh output directories for each constraints
file. For the real-data matrix, first run the Seurat PBMC3k example, then:

```sh
Rscript adapters/seurat/scripts/compat-fixtures.R .r-tests/pbmc-fixtures examples/seurat/output/source.rds
python scripts/check_compatibility.py --fixtures .r-tests/pbmc-fixtures --output .compat-results/pbmc-0.11.4 --dataset PBMC3k
```

Repeat the Python command under each pinned environment. The R fixtures are
immutable inputs; each run copies them to a fresh output. Changing dependencies
requires a new report and CI check before expanding this table's claims.

`to_tsv` uses the same input validation/dtype rules and exports expression only.
See [the TSV guide](../docs/tsv-export.md) for single-row densification, sparse
input memory, subset semantics and protected file publication.
