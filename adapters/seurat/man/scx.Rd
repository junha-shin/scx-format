\name{write_scx}
\alias{write_scx}
\alias{read_scx}
\title{Read and Write the Basic SCX Profile}
\description{Safe object I/O with the same read_scx and write_scx names as Python. Python is not required for Seurat I/O.}
\usage{
write_scx(object, path, counts = "counts", normalized = "data",
  analysis = FALSE, overwrite = FALSE,
  assay = SeuratObject::DefaultAssay(object))
read_scx(path, assay = "RNA", assay_class = c("Assay5", "Assay"))
}
\arguments{
\item{object}{A Seurat object.}
\item{path}{An SCX directory.}
\item{counts}{Exact source counts layer name. Raw counts are required.}
\item{normalized}{Exact normalized layer name, or NULL. Defaults to data; missing data exports counts only with a warning.}
\item{analysis}{Include supported embeddings and cell graphs. Default FALSE.}
\item{overwrite}{Explicitly replace an existing whole SCX output after validation. Default FALSE.}
\item{assay}{One selected assay, defaulting to the object's active assay on write and RNA on read.}
\item{assay_class}{Assay5 (default) or legacy Assay under SeuratObject >=5.}
}
\value{Writing invisibly returns the output path. Reading returns a Seurat object.}
\details{The basic profile contains counts, optional normalized expression and typed cell/feature metadata. Set analysis=TRUE for coordinates and cell graphs. Reading restores every supported field present in the directory. Omitted embeddings/graphs and extra assays are reported. Additional analysis internals are outside the profile. Validation, staging, overwrite protection, locking and rollback are automatic. validate_scx is an optional diagnostic. Legacy read_scx_seurat/write_scx_seurat remain compatible and retain their broader defaults.}
