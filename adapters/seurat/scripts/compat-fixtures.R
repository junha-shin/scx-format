# Rscript .../compat-fixtures.R <new fixture root> [existing PBMC3k source.rds]
source("adapters/seurat/R/scx.R")
source("adapters/seurat/tests/testthat/helper-fixtures.R")
source("adapters/seurat/tests/testthat/test-scx.R")
source("examples/seurat/reference.R")
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(1L, 2L), !dir.exists(args[1]))
dir.create(args[1], recursive = TRUE)
for (assay_class in c("Assay", "Assay5")) {
  object <- if (length(args) == 2L) readRDS(args[2]) else tiny_seurat(assay_class)
  if (!inherits(object[["RNA"]], assay_class)) {
    original <- object[["RNA"]]
    object[["RNA"]] <- methods::as(original, assay_class)
    for (layer in c("counts", "data")) stopifnot(isTRUE(all.equal(
      SeuratObject::LayerData(original, layer = layer),
      SeuratObject::LayerData(object[["RNA"]], layer = layer), tolerance = 0)))
    stopifnot(identical(original[[]], object[["RNA"]][[]]))
  }
  output <- file.path(args[1], assay_class); dir.create(output)
  save_reference(object, output)
  write_scx_seurat(object, file.path(output, "seurat.scx"))
}
versions <- list(R = R.version.string, SeuratObject = as.character(packageVersion("SeuratObject")),
                 Matrix = as.character(packageVersion("Matrix")), jsonlite = as.character(packageVersion("jsonlite")))
writeLines(.scx_json(versions), file.path(args[1], "R-versions.json"))
