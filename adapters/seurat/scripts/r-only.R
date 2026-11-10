# Installed-package proof: no Python on PATH and no external process allowed.
library(scxSeurat)
library(SeuratObject)
library(Matrix)
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, !dir.exists(args[1]))
dir.create(args[1], recursive = TRUE)
.run_r_only <- function() {
original_path <- Sys.getenv("PATH")
on.exit({Sys.setenv(PATH = original_path); untrace(base::system); untrace(base::system2)}, add = TRUE)
trace(base::system, tracer = quote(stop("external process forbidden")), print = FALSE)
trace(base::system2, tracer = quote(stop("external process forbidden")), print = FALSE)
Sys.setenv(PATH = "", RETICULATE_PYTHON = "/nonexistent/python", PYTHON = "/nonexistent/python")
stopifnot(identical(Sys.getenv("PATH"), ""))
for (assay_class in c("Assay", "Assay5")) {
  counts <- Matrix(matrix(c(1, 0, 3, 0, 2, 4), 3, 2), sparse = TRUE)
  dimnames(counts) <- list(c("g1", "g2", "g3"), c("c1", "c2"))
  assay <- if (assay_class == "Assay") CreateAssayObject(counts = counts) else CreateAssay5Object(counts = counts)
  before <- CreateSeuratObject(assay)
  LayerData(before, layer = "data") <- log1p(counts)
  before$number <- c(1L, NA_integer_); before$flag <- c(TRUE, NA)
  before$group <- ordered(c("b", NA), levels = c("b", "a", "unused"))
  before$text <- c("", NA_character_)
  before[["RNA"]] <- AddMetaData(before[["RNA"]], data.frame(score = c(0.5, NA, 2.5), row.names = rownames(counts)))
  before[["pca"]] <- CreateDimReducObject(matrix(c(1, 2, 3, 4), 2,
    dimnames = list(colnames(counts), c("PC_1", "PC_2"))), assay = "RNA", key = "PC_")
  graph <- Matrix(matrix(c(0, 2, 1, 0), 2), sparse = TRUE)
  dimnames(graph) <- list(colnames(counts), colnames(counts))
  graph <- as.Graph(graph); DefaultAssay(graph) <- "RNA"
  before[["RNA_nn"]] <- graph
  root <- file.path(args[1], assay_class)
  suppressWarnings(write_scx_seurat(before, root))
  suppressWarnings(write_scx_seurat(before, root, overwrite = TRUE, python = "/nonexistent/python"))
  stopifnot(validate_scx(root)$cells == 2L)
  after <- read_scx_seurat(root, assay_class = assay_class, python = "/nonexistent/python")
  stopifnot(identical(before[[]], after[[]]), identical(before[["RNA"]][[]], after[["RNA"]][[]]))
  for (layer in c("counts", "data")) stopifnot(max(abs(LayerData(before, layer = layer) - LayerData(after, layer = layer))) == 0)
  stopifnot(max(abs(Embeddings(before[["pca"]]) - Embeddings(after[["pca"]]))) == 0)
  stopifnot(max(abs(before[["RNA_nn"]] - after[["RNA_nn"]])) == 0)
}
writeLines(jsonlite::toJSON(list(passed = TRUE, python_on_path = FALSE, external_process_calls = 0,
  assays = c("Assay", "Assay5"), R = R.version.string, scxSeurat = as.character(packageVersion("scxSeurat"))),
  auto_unbox = TRUE), file.path(args[1], "report.json"))
cat("R-only installed-package validation passed for Assay and Assay5\n")
}
.run_r_only()
