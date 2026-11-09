# CI fixture and R unit tests, run from repository root.
source("adapters/seurat/R/scx.R")
source("adapters/seurat/tests/testthat/test-scx.R")
source("examples/seurat/reference.R")
output <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(!dir.exists(output))
dir.create(output, recursive = TRUE)
object <- tiny_seurat()
save_reference(object, output)
write_scx_seurat(object, file.path(output, "seurat.scx"))
