# Reproduce a supported single-assay RDS round trip without Python/subprocesses.
library(scxSeurat)
library(SeuratObject)
args <- commandArgs(TRUE)
stopifnot(length(args) == 2L, !dir.exists(args[2]))
before <- readRDS(args[1]); assay <- DefaultAssay(before)
trace(base::system, tracer = quote(stop("external process forbidden")), print = FALSE)
trace(base::system2, tracer = quote(stop("external process forbidden")), print = FALSE)
Sys.setenv(PATH = "", PYTHON = "/nonexistent/python", RETICULATE_PYTHON = "/nonexistent/python")
dir.create(args[2], recursive = TRUE)
root <- file.path(args[2], "dataset.scx")
suppressWarnings(write_scx_seurat(before, root, assay = assay))
suppressWarnings(write_scx_seurat(before, root, assay = assay, overwrite = TRUE))
after <- read_scx_seurat(root, assay = assay,
  assay_class = if (inherits(before[[assay]], "Assay5")) "Assay5" else "Assay")
stopifnot(identical(before[[]], after[[]]), identical(before[[assay]][[]], after[[assay]][[]]))
stopifnot(identical(colnames(before), colnames(after)), identical(rownames(before), rownames(after)))
errors <- list()
for (layer in intersect(c("counts", "data"), Layers(before[[assay]], search = NA)))
  errors[[layer]] <- max(abs(LayerData(before, assay = assay, layer = layer) - LayerData(after, assay = assay, layer = layer)))
for (name in Reductions(before)) errors[[name]] <- max(abs(Embeddings(before[[name]]) - Embeddings(after[[name]])))
for (name in Graphs(before)) errors[[name]] <- max(abs(before[[name]] - after[[name]]))
stopifnot(all(unlist(errors) == 0))
report <- list(passed = TRUE, python_on_path = FALSE, external_process_calls = 0,
  source_sha256 = digest::digest(args[1], algo = "sha256", file = TRUE),
  cells = ncol(after), features = nrow(after), errors = errors, metadata = "exact",
  R = R.version.string, scxSeurat = as.character(packageVersion("scxSeurat")),
  SeuratObject = as.character(packageVersion("SeuratObject")))
writeLines(jsonlite::toJSON(report, auto_unbox = TRUE, pretty = TRUE), file.path(args[2], "report.json"))
cat(jsonlite::toJSON(report, auto_unbox = TRUE), "\n")
