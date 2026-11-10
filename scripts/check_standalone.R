# Test the one-file source() route; no scxSeurat, native helpers, or subprocesses.
args <- commandArgs(TRUE)
stopifnot(length(args) %in% c(1L, 2L), !dir.exists(args[1]))
dir.create(args[1], recursive = TRUE)
library(SeuratObject)
library(Matrix)
loadNamespace("jsonlite"); loadNamespace("digest")
if (Sys.getenv("SCX_REQUIRE_UNINSTALLED") == "1") stopifnot(!nzchar(system.file(package="scxSeurat")))
trace(base::loadNamespace, tracer = quote(if (identical(package, "scxSeurat")) stop("SCX package loading forbidden")), print = FALSE)
trace(base::system, tracer = quote(stop("subprocess forbidden")), print = FALSE)
trace(base::system2, tracer = quote(stop("subprocess forbidden")), print = FALSE)
Sys.setenv(PATH = "", PYTHON = "/nonexistent/python", RETICULATE_PYTHON = "/nonexistent/python")
source("standalone/scx.R")
stopifnot(!"scxSeurat" %in% loadedNamespaces(), !is.loaded("scx_check_gzip", PACKAGE="scxSeurat"))
source("adapters/seurat/tests/testthat/helper-fixtures.R")
fails <- function(code, pattern = NULL) {
  error <- tryCatch({force(code); NULL}, error = conditionMessage)
  stopifnot(!is.null(error))
  if (!is.null(pattern)) stopifnot(grepl(pattern, error, fixed=TRUE))
}
for (type in c("Assay", "Assay5")) {
  before <- tiny_seurat(type)
  basic <- file.path(args[1], paste0(type, "-basic.scx"))
  full <- file.path(args[1], paste0(type, "-full.scx"))
  suppressWarnings(write_scx(before, basic))
  restored <- read_scx(basic, assay_class=type)
  stopifnot(identical(restored[[]], before[[]]), identical(restored[["RNA"]][[]], before[["RNA"]][[]]),
    !length(Reductions(restored)), !length(Graphs(restored)))
  write_scx(before, full, analysis=TRUE)
  restored <- read_scx(full, assay_class=type)
  stopifnot(identical(before[[]], restored[[]]), identical(before[["RNA"]][[]], restored[["RNA"]][[]]))
  stopifnot(max(abs(LayerData(before, layer="counts")-LayerData(restored, layer="counts"))) == 0,
    max(abs(LayerData(before, layer="data")-LayerData(restored, layer="data"))) < 1e-12,
    max(abs(Embeddings(before[["pca"]])-Embeddings(restored[["pca"]]))) < 1e-12,
    max(abs(before[["RNA_nn"]]-restored[["RNA_nn"]])) == 0)
  hash <- digest::digest(file.path(full, "matrices/counts.mtx.gz"), file=TRUE, algo="sha256")
  fails(write_scx(before, full), "already exists")
  fails(write_scx(before, full, overwrite=TRUE), "overwrite=FALSE")
  stopifnot(hash == digest::digest(file.path(full, "matrices/counts.mtx.gz"), file=TRUE, algo="sha256"))
  # CRC mutation and missing trailer must fail despite otherwise readable values.
  path <- file.path(full, "matrices/counts.mtx.gz")
  original <- readBin(path, "raw", file.info(path)$size)
  changed <- original; changed[length(changed)-7L] <- as.raw(bitwXor(as.integer(changed[length(changed)-7L]), 1L))
  writeBin(changed, path); fails(read_scx(full))
  writeBin(head(original, -8L), path); fails(read_scx(full))
  writeBin(original, path)
  # A bad canonical axis is rejected by structural validation.
  cells <- file.path(full, "cells.tsv.gz")
  original_cells <- readBin(cells, "raw", file.info(cells)$size)
  con <- gzfile(cells, "wt"); writeLines(c("cell_id", "duplicate", "duplicate"), con); close(con)
  fails(read_scx(full))
  writeBin(original_cells, cells)
}
if (length(args) == 2L) for (type in c("Assay", "Assay5")) {
  python <- read_scx(file.path(args[2], "cli.scx"), assay_class=type)
  metadata <- python[[]]  # Seurat's $ accessor itself drops ordered-factor class.
  stopifnot(all(dim(python) == c(2L,3L)), identical(colnames(python),c("c1","c2","c3")),
    identical(rownames(python),c("g1","g2")), identical(unname(python$number),c(1L,NA_integer_,3L)),
    identical(unname(python$flag),c(TRUE,FALSE,NA)), is.ordered(metadata$group),
    identical(levels(metadata$group),c("b","a","unused")),
    max(abs(as.matrix(LayerData(python,layer="counts"))-matrix(c(0,2,3,0,4,1),2))) == 0,
    max(abs(as.matrix(python[["nn"]])-matrix(c(0,1,0,2,0,0,0,3,0),3))) == 0)
}
# Inject a failed copy: the partial output has no manifest and cannot be read.
env <- environment(write_scx)
original_copy <- env$.scx_copy
env$.scx_copy <- function(...) stop("injected disk/copy failure")
partial <- file.path(args[1], "interrupted.scx")
fails(write_scx(tiny_seurat(), partial, analysis=TRUE), "Incomplete output retained")
stopifnot(dir.exists(partial), !file.exists(file.path(partial,"manifest.tsv")))
fails(read_scx(partial))
fails(write_scx(tiny_seurat(), partial), "already exists")
env$.scx_copy <- original_copy
stopifnot(!"scxSeurat" %in% loadedNamespaces())
report <- list(passed=TRUE, installed_scx_loads=0, package_installations=0,
  external_process_calls=0, assays=c("Assay","Assay5"), shape=c(2,3), typed_metadata=TRUE,
  overwrite_protection=TRUE, interrupted_copy_rejected=TRUE, gzip_crc_and_truncation_rejected=TRUE)
writeLines(jsonlite::toJSON(report,auto_unbox=TRUE,pretty=TRUE), file.path(args[1],"report.json"))
cat("Standalone R validation passed\n")
