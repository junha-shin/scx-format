args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
source("adapters/seurat/R/scx.R")
source_object <- readRDS(file.path(args[1], "source.rds"))
restored <- read_scx_seurat(file.path(args[1], "ann.scx"), python = args[2])
stopifnot(identical(colnames(source_object), colnames(restored)), identical(rownames(source_object), rownames(restored)))
for (pair in list(list(source_object[[]], restored[[]]), list(source_object[["RNA"]][[]], restored[["RNA"]][[]]))) {
  comparison <- all.equal(pair[[1]], pair[[2]], tolerance = 0)
  if (!isTRUE(comparison)) stop("metadata comparison failed: ", paste(comparison, collapse = "; "))
}
for (layer in c("counts", "data")) {
  stopifnot(max(abs(SeuratObject::LayerData(source_object, layer = layer) - SeuratObject::LayerData(restored, layer = layer))) <= 1e-12)
}
for (name in SeuratObject::Reductions(source_object)) {
  stopifnot(max(abs(SeuratObject::Embeddings(source_object[[name]]) - SeuratObject::Embeddings(restored[[name]]))) <= 1e-12)
}
for (name in SeuratObject::Graphs(source_object)) stopifnot(max(abs(source_object[[name]] - restored[[name]])) <= 1e-12)
cat("Reverse AnnData -> SCX -> Seurat comparison passed\n")
