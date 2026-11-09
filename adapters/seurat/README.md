# Seurat adapter MVP

Requires R, SeuratObject >= 5, Matrix and jsonlite. The adapter uses public
`LayerData` APIs and accepts a single Seurat assay (Assay5 or Assay). Source it
from the repository root, or install the R package in this directory.

```r
install.packages(c("SeuratObject", "jsonlite"))
source("adapters/seurat/R/scx.R")
write_scx_seurat(pbmc, "pbmc.scx", assay = "RNA")
# Import also requires Python with scx-format[anndata] installed.
restored <- read_scx_seurat("pbmc.scx", python = "/path/to/python")
# Explicit legacy output; independent of Seurat.object.assay.version option.
legacy <- read_scx_seurat("pbmc.scx", python = "/path/to/python", assay_class = "Assay")
```

Exports require a new destination and stage all files before publishing. Run
`scx validate pbmc.scx` for independent Python validation. R import always calls
that validator before reading paths and data. It refuses feature IDs changed by
Seurat rather than silently changing their alignment.

## Mapping

| Seurat | SCX | AnnData |
|---|---|---|
| selected assay `counts` (FEATURE × CELL) | transposed counts (CELL × FEATURE) | `layers['counts']` |
| selected assay `data` | transposed normalized (CELL × FEATURE) | `X` |
| `object[[]]` | typed cells table | `obs` |
| `object[[assay]][[]]` | typed features table | `var` |
| `Embeddings(object[[name]])` | CELL × COMPONENT | `obsm[name]` |
| cell graph `object[[name]]` | CELL × CELL; no transpose | `obsp[name]` |

SCX stores reduction and graph names verbatim. Embedding component labels are
not portable in v0.1; R import generates valid Seurat component keys. Graphs are
written as explicit numeric Matrix Market `general` matrices, including both
directions for symmetric graphs and values for binary graphs.

Exact layer names are required: split layers such as `counts.batch1` must first
be explicitly joined using Seurat's `JoinLayers`. Use `counts_layer` and
`normalized_layer` to select different exact layer names; set
`normalized_layer = NULL` for counts-only. A missing default `data` layer warns
and exports counts-only. Use `reductions` and `graphs` to choose named objects
or `character()` to omit them. Every selected object must cover the canonical
axis exactly; matching IDs are reordered, while partial/extra cells or features
are rejected. Reductions/graphs from another assay are rejected.

## Metadata types

Uses [SCX metadata schema 1.0](../../specification/metadata-schema-v1.md).

| Native R | SCX / pandas |
|---|---|
| integer without NA | `int32` |
| integer with NA | `Int32` |
| double, NA/NaN allowed | `float64` (missing → NaN) |
| logical without NA | `bool` |
| logical with NA | `boolean` |
| character | `object` strings; NA → null |
| factor / ordered factor | categorical with all levels and order retained |

Empty strings, the literal `"null"`, and NA are distinct. Dates, list columns,
custom numeric classes (including bit64 integer64), complex values and infinity
are rejected. No implicit conversion is applied to unsupported classes.

R import maps smaller integer widths to native integer and float32/nullable
float types to native double; width provenance cannot be retained by R's native
types. It rejects int64/uint32/uint64 types rather than silently round them, and
rejects INT_MIN (R reserves its bit pattern for NA). Numeric/boolean category
labels are rejected on import because R factor levels must be character.
These import limits do not affect native R → SCX → AnnData dtype preservation.

The [compatibility tables](../../compatibility/README.md) distinguish legacy
Assay class support under SeuratObject 5 from the untested Seurat 4 runtime.

Active identities, additional assays/layers, feature loadings, scaled values,
command history, images and arbitrary analysis settings are not in this MVP;
the exporter emits a warning identifying its scope. Existing datasets should
retain their original Seurat objects if those analysis internals are needed.

## Tests and real data

```sh
Rscript adapters/seurat/scripts/run-tests.R
Rscript adapters/seurat/scripts/bridge.R .r-tests/bridge
python adapters/seurat/scripts/check_bridge.py .r-tests/bridge
Rscript examples/seurat/check_reverse.R .r-tests/bridge /path/to/python
```

The bridge tests include a nonsquare expression matrix and an asymmetric cell
graph, so incorrect transposition cannot pass by shape or symmetry alone.
[PBMC3k validation](../../examples/seurat/README.md) exercises genuine Seurat
normalization, PCA/UMAP and cell graphs on a public dataset.

API references: [LayerData](https://satijalab.github.io/seurat-object/reference/Layers.html),
[CreateSeuratObject](https://satijalab.github.io/seurat-object/reference/CreateSeuratObject.html),
[CreateDimReducObject](https://satijalab.github.io/seurat-object/reference/CreateDimReducObject.html).
