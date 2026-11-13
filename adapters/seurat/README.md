# Seurat <-> SCX

Use one R file in your existing SeuratObject >=5 environment.
No SCX package, compiler, or Python is needed:

```r
source("scx.R")
write_scx(pbmc, "pbmc.scx")
restored <- read_scx("pbmc.scx")
write_scx(pbmc, "analysis.scx", analysis = TRUE)
```

Save scx.R once and source it offline.
It uses existing Matrix, jsonlite, and digest packages.
Missing packages cause an error, not an install.
The R file refuses replacement and accepts single-member gzip.
See [setup and safety](../../docs/install.md).
The sections below also cover the optional installed package.

Default writes keep counts, optional normalized values, and typed cell/feature metadata.
Python and R use the same read_scx / write_scx names.
R uses counts/data by default. See [setup](../../docs/install.md) and
[conversion examples](../../docs/getting-started.md).
Legacy write_scx_seurat / read_scx_seurat keep their full defaults,
including reductions/graphs selection.
Development source() still needs installed native helpers. See [R-only use](../../docs/r-only.md).

The installed package stages and checks all files before publishing them with native helpers.
It needs explicit overwrite=TRUE to replace existing output.
If replacement fails, it restores or keeps the old data.
See [recovery and storage limits](../../docs/safe-io.md).
Reads reuse checked resources. They fail if Seurat changes feature IDs rather than risking wrong alignment.

validate_scx() checks the subset R can read, including native type limits.
Python also handles wide integers and numeric categories that R cannot read.
The legacy python argument is accepted but ignored.
Python and AnnData are needed only for H5AD conversion.

## Mapping

| Seurat | SCX | AnnData |
|---|---|---|
| selected assay `counts` (FEATURE x CELL) | transposed counts (CELL x FEATURE) | `layers['counts']` |
| selected assay `data` | transposed normalized (CELL x FEATURE) | `X` |
| `object[[]]` | typed cells table | `obs` |
| `object[[assay]][[]]` | typed features table | `var` |
| `Embeddings(object[[name]])` | CELL x COMPONENT | `obsm[name]` |
| cell graph `object[[name]]` | CELL x CELL; no transpose | `obsp[name]` |

Reduction and graph names are kept.
Embedding component labels are not portable in v0.1; R creates valid Seurat keys.
Graphs use numeric general MTX, with both directions for symmetric graphs
and explicit values for binary graphs.

Use exact layer names. Join split layers such as counts.batch1 with JoinLayers when needed.
Use counts and normalized (legacy counts_layer and normalized_layer) to select layers.
Use normalized=NULL for counts only. A missing default data layer warns and exports counts only.
Legacy reductions and graphs select named objects; character() omits them.
Selected objects must cover the full axis. Matching IDs are reordered.
Extra or missing cells/features and reductions/graphs from another assay are rejected.

## Metadata types

Uses [metadata schema 1.0](../../specification/metadata-schema-v1.md).

| Native R | SCX / pandas |
|---|---|
| integer without NA | `int32` |
| integer with NA | `Int32` |
| double, NA/NaN allowed | `float64` (missing -> NaN) |
| logical without NA | `bool` |
| logical with NA | `boolean` |
| character | `object` strings; NA -> null |
| factor / ordered factor | categorical with all levels and order retained |

Empty text, literal "null", and NA stay distinct.
Dates, list columns, custom numeric classes (including bit64 integer64),
complex values, and infinity are rejected. Unsupported classes are not silently converted.

R reads smaller integers as native integer and float32/nullable floats as double.
The original width cannot always be kept in native R types.
It rejects int64/uint32/uint64 rather than rounding them.
It also rejects INT_MIN, whose bit pattern R uses for NA.
Numeric and boolean category labels are rejected because R factor levels are text.
These limits do not change native R -> SCX -> AnnData type preservation.

The [compatibility tables](../../compatibility/README.md) separate legacy Assay
under SeuratObject 5 from the untested Seurat 4 runtime.

Other assays/layers, identities, loadings, scaled values, history, images,
and arbitrary analysis settings are outside the profile.
The exporter warns about omitted data. Keep the original Seurat object.

## Tests and real data

```sh
Rscript adapters/seurat/scripts/run-tests.R
Rscript adapters/seurat/scripts/bridge.R .r-tests/bridge
python adapters/seurat/scripts/check_bridge.py .r-tests/bridge
Rscript examples/seurat/check_reverse.R .r-tests/bridge /path/to/python
```

The bridge uses a non-square expression matrix and a directed graph.
A wrong transpose cannot pass just because a matrix is square or symmetric.
[PBMC3k checks](../../examples/seurat/README.md) use real Seurat normalization,
PCA/UMAP, and cell graphs.

API references: [LayerData](https://satijalab.github.io/seurat-object/reference/Layers.html),
[CreateSeuratObject](https://satijalab.github.io/seurat-object/reference/CreateSeuratObject.html),
[CreateDimReducObject](https://satijalab.github.io/seurat-object/reference/CreateDimReducObject.html).
