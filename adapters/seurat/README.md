# Seurat <-> SCX

Use one R file in your existing SeuratObject >=5 environment; no SCX package,
compiler or Python is required:

```r
source("scx.R")
write_scx(pbmc, "pbmc.scx")
restored <- read_scx("pbmc.scx")
write_scx(pbmc, "analysis.scx", analysis = TRUE)
```

Save the file and `source("scx.R")` for offline use. Existing Matrix, jsonlite and
digest are reused; missing dependencies cause an error, never an installation.
Standalone R refuses replacement and supports single-member gzip. See
[no-install setup and safety](../../docs/install.md). The rest of this document
also describes optional installed-package features, which remain unchanged.

The basic profile preserves counts, optional normalized expression and typed
cell/feature metadata. Both languages use the same `read_scx` / `write_scx`
names; default R sources are `counts` / `data`. See [installation](../../docs/install.md)
and [cross-language examples](../../docs/getting-started.md). Legacy
`write_scx_seurat` / `read_scx_seurat` retain their full defaults, including
selective `reductions` / `graphs` controls. Development `source()` still needs
installed native helpers; see [R-only details](../../docs/r-only.md).

Installed-package exports stage all files, run the R-native validator, and publish using
R's transaction code and registered native helpers. Existing destinations require explicit
`overwrite = TRUE`; failed replacement restores or retains the previous output.
See [safe I/O and manual recovery](../../docs/safe-io.md) for filesystem limits,
the brief directory-replacement gap and hard-interruption behavior. R import
retains resources from its validation pass, avoiding unvalidated rereads.
It refuses feature IDs changed by
Seurat rather than silently changing their alignment.

`validate_scx()` checks the R-readable subset, including native R dtype limits
below; the Python validator additionally handles wide integer metadata and
numeric categories that R cannot import. The optional `python` argument remains
accepted and ignored for source compatibility. Python/AnnData are needed only
when actually converting to or from H5AD, not for these R operations.

## Mapping

| Seurat | SCX | AnnData |
|---|---|---|
| selected assay `counts` (FEATURE x CELL) | transposed counts (CELL x FEATURE) | `layers['counts']` |
| selected assay `data` | transposed normalized (CELL x FEATURE) | `X` |
| `object[[]]` | typed cells table | `obs` |
| `object[[assay]][[]]` | typed features table | `var` |
| `Embeddings(object[[name]])` | CELL x COMPONENT | `obsm[name]` |
| cell graph `object[[name]]` | CELL x CELL; no transpose | `obsp[name]` |

SCX stores reduction and graph names verbatim. Embedding component labels are
not portable in v0.1; R import generates valid Seurat component keys. Graphs are
written as explicit numeric Matrix Market `general` matrices, including both
directions for symmetric graphs and values for binary graphs.

Exact layer names are required: split layers such as `counts.batch1` must first
be explicitly joined using Seurat's `JoinLayers`. Use public `counts` and `normalized` (legacy `counts_layer` and
`normalized_layer`) to select different exact layer names; set
`normalized = NULL` for counts-only. A missing default `data` layer warns
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
| double, NA/NaN allowed | `float64` (missing -> NaN) |
| logical without NA | `bool` |
| logical with NA | `boolean` |
| character | `object` strings; NA -> null |
| factor / ordered factor | categorical with all levels and order retained |

Empty strings, the literal `"null"`, and NA are distinct. Dates, list columns,
custom numeric classes (including bit64 integer64), complex values and infinity
are rejected. No implicit conversion is applied to unsupported classes.

R import maps smaller integer widths to native integer and float32/nullable
float types to native double; width provenance cannot be retained by R's native
types. It rejects int64/uint32/uint64 types rather than silently round them, and
rejects INT_MIN (R reserves its bit pattern for NA). Numeric/boolean category
labels are rejected on import because R factor levels must be character.
These import limits do not affect native R -> SCX -> AnnData dtype preservation.

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
