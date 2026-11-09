\name{write_scx_seurat}
\alias{write_scx_seurat}
\alias{read_scx_seurat}
\title{Exchange Single-Assay Seurat Objects with SCX}
\description{Export counts, normalized data, typed metadata, cell embeddings and graphs in SCX canonical cell-by-feature orientation. Import validates SCX using the Python reference validator.}
\usage{
write_scx_seurat(object, path, assay = SeuratObject::DefaultAssay(object),
  counts_layer = "counts", normalized_layer = "data",
  reductions = SeuratObject::Reductions(object),
  graphs = SeuratObject::Graphs(object))
read_scx_seurat(path, assay = "RNA", python = Sys.which("python"))
}
\arguments{
\item{object}{A Seurat object using SeuratObject version 5 or later.}
\item{path}{A new output directory for export; an existing SCX directory for import.}
\item{assay}{One assay name.}
\item{counts_layer}{Exact counts layer name. Split layers must be joined explicitly.}
\item{normalized_layer}{Exact normalized layer name, or NULL for counts-only export.}
\item{reductions}{Names of cell reductions to export, or character() to omit.}
\item{graphs}{Names of cell graphs to export, or character() to omit.}
\item{python}{Python executable with scx-format installed, used for validation.}
}
\value{Export invisibly returns the output path; import returns a Seurat object.}
\details{Metadata use schema 1.0 with native R type mappings. Unsupported classes and unsafe integer imports fail explicitly. Additional assays, identities, scaled layers, feature loadings and analysis history are not exchanged. See the repository adapter README for the mapping and limits.}
