# CI fixture and R unit tests, run from repository root.
source("adapters/seurat/R/scx.R")
source("adapters/seurat/tests/testthat/test-scx.R")
source("examples/seurat/reference.R")
args <- commandArgs(trailingOnly = TRUE)
output <- args[1]
stopifnot(!dir.exists(output))
dir.create(output, recursive = TRUE)
object <- tiny_seurat(if (length(args) >= 2L) args[2] else "Assay5")
save_reference(object, output)
write_scx_seurat(object, file.path(output, "seurat.scx"))
