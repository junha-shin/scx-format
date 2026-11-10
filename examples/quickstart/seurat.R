# Uses the installed R package, never source() from the repository.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
library(scxSeurat)
mode <- args[1]; input <- args[2]; output <- args[3]; python <- args[4]
if (mode == "create") {
  counts <- Matrix::Matrix(matrix(c(1, 0, 2, 0, 0, 3, 0, 1, 4, 0, 0, 2), 4, 3), sparse = TRUE)
  dimnames(counts) <- list(paste0("gene", 1:4), paste0("cell", 1:3))
  object <- SeuratObject::CreateSeuratObject(SeuratObject::CreateAssay5Object(counts = counts))
  SeuratObject::LayerData(object, layer = "data") <- log1p(counts)
  object@meta.data <- data.frame(
    batch = ordered(c("b", "a", NA), levels = c("b", "a", "unused")),
    passing = c(TRUE, NA, FALSE), number = c(1L, NA_integer_, 3L),
    note = c("", NA_character_, "null"), row.names = colnames(counts))
  object[["RNA"]] <- SeuratObject::AddMetaData(object[["RNA"]],
    data.frame(score = c(0.5, NA, 1.5, 2.5), row.names = rownames(counts)))
  embedding <- matrix(c(1, 3, 5, 2, 4, 6), 3,
    dimnames = list(colnames(counts), c("PC_1", "PC_2")))
  object[["pca"]] <- SeuratObject::CreateDimReducObject(embedding, assay = "RNA", key = "PC_")
  graph <- Matrix::Matrix(matrix(c(0, 0, 3, 1, 0, 0, 0, 2, 0), 3), sparse = TRUE)
  dimnames(graph) <- list(colnames(counts), colnames(counts))
  graph <- SeuratObject::as.Graph(graph); SeuratObject::DefaultAssay(graph) <- "RNA"
  object[["RNA_nn"]] <- graph
  writeLines(jsonlite::toJSON(list(R = R.version.string,
    scxSeurat = as.character(packageVersion("scxSeurat")),
    SeuratObject = as.character(packageVersion("SeuratObject"))), auto_unbox = TRUE),
    file.path(dirname(input), "R-versions.json"))
  saveRDS(object, input)
  write_scx_seurat(object, output, python = python)
} else if (mode == "convert") {
  object <- read_scx_seurat(input, python = python)
  write_scx_seurat(object, output, python = python)
} else if (mode == "verify") {
  before <- readRDS(input)
  after <- read_scx_seurat(output, python = python)
  stopifnot(identical(colnames(before), colnames(after)), identical(rownames(before), rownames(after)))
  stopifnot(identical(before[[]], after[[]]), identical(before[["RNA"]][[]], after[["RNA"]][[]]))
  for (layer in c("counts", "data")) stopifnot(
    max(abs(SeuratObject::LayerData(before, layer = layer) - SeuratObject::LayerData(after, layer = layer))) <= 1e-12)
  for (name in SeuratObject::Reductions(before)) stopifnot(
    max(abs(SeuratObject::Embeddings(before[[name]]) - SeuratObject::Embeddings(after[[name]]))) <= 1e-12)
  for (name in SeuratObject::Graphs(before)) stopifnot(max(abs(before[[name]] - after[[name]])) <= 1e-12)
  cat("Installed Seurat package round trip passed\n")
} else stop("unknown quickstart mode")
