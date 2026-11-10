\name{to_tsv}
\alias{to_tsv}
\title{Export SCX Expression as a Matrix TSV}
\description{Validate SCX and export CELL by FEATURE expression one dense output row at a time. No Seurat object is constructed. Existing output is always refused.}
\usage{to_tsv(path, output, matrix = "counts", cells = NULL, features = NULL)}
\arguments{
\item{path}{Input SCX directory.}
\item{output}{TSV file outside the input directory. A .gz suffix enables gzip.}
\item{matrix}{counts or normalized. Missing normalized is an error.}
\item{cells}{NULL for all cells, or nonempty unique character IDs in output order.}
\item{features}{NULL for all features, or nonempty unique character IDs in output order.}
}
\value{Invisibly returns the output path.}
\details{Sparse input and validation require memory; dense working space is limited to one selected row. The first column is cell_id. Feature IDs form the remaining headers. Duplicate or unknown selections fail explicitly. No replacement is supported. R stages the complete TSV, exclusively creates the destination, then copies and verifies bytes. A copy failure or hard interruption can leave a partial TSV and lock; wait until the writer stops and inspect before recovery. This is an expression view, not a metadata-preserving interchange format.}
