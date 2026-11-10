\name{write_scx_seurat}
\alias{write_scx_seurat}
\alias{read_scx_seurat}
\alias{validate_scx}
\title{Exchange Single-Assay Seurat Objects with SCX}
\description{Export counts, normalized data, typed metadata, cell embeddings and graphs in SCX canonical cell-by-feature orientation. Export, validation and import run in R without Python.}
\usage{
write_scx_seurat(object, path, assay = SeuratObject::DefaultAssay(object),
  counts_layer = "counts", normalized_layer = "data",
  reductions = SeuratObject::Reductions(object),
  graphs = SeuratObject::Graphs(object), overwrite = FALSE,
  python = NULL)
read_scx_seurat(path, assay = "RNA", python = NULL,
  assay_class = c("Assay5", "Assay"))
validate_scx(path)
}
\arguments{
\item{object}{A Seurat object using SeuratObject version 5 or later.}
\item{path}{An output directory for export; an existing SCX directory for import. Export refuses existing paths by default.}
\item{assay}{One assay name.}
\item{counts_layer}{Exact counts layer name. Split layers must be joined explicitly.}
\item{normalized_layer}{Exact normalized layer name, or NULL for counts-only export.}
\item{reductions}{Names of cell reductions to export, or character() to omit.}
\item{graphs}{Names of cell graphs to export, or character() to omit.}
\item{overwrite}{Logical, default FALSE. Explicitly replace an existing SCX directory after staging and validation. Replaces the whole directory, including unlisted files.}
\item{python}{Unused compatibility argument. Python is not inspected or invoked.}
\item{assay_class}{Explicit output assay class: Assay5 (default) or legacy Assay.}
}
\value{Export invisibly returns the output path; import returns a Seurat object. Validation returns a list with cells, features and files counts.}
\details{Metadata use schema 1.0 with native R type mappings. The validator checks the R-readable subset: wide integer metadata, INT_MIN and numeric/boolean categorical levels are rejected explicitly. Matrix integer values beyond native R's exact range are rejected. Additional assays, identities, scaled layers, feature loadings and analysis history are not exchanged.

Native helpers check complete gzip streams including CRC/trailers and provide atomic no-replace directory rename. Installing from source requires an R compiler toolchain and zlib; Python is not required. Import retains checked resources from the validation pass.

Export stages and validates files before publication. Existing outputs require explicit overwrite consent. Failed directory replacement restores the previous output or retains its backup and lock with recovery paths. Replacement briefly makes the final directory name unavailable to readers. Hard termination and power loss require manual recovery; locks are not automatically removed. Output symlinks and unrelated directory replacement are refused. See the repository adapter README and docs/safe-io.md for recovery instructions and filesystem limits.}
