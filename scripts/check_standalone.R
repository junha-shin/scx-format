# Test the one-file source() route; no scxSeurat, native helpers, or subprocesses.
args <- commandArgs(TRUE)
# Rscript can inherit an invalid Unix locale name on Windows. Use its native
# UTF-8 locale for Unicode fixture creation (CI already uses this locale).
if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]])) Sys.setlocale("LC_CTYPE", ".UTF-8")
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
  tsv <- file.path(args[1],paste0(type,"-counts.tsv"))
  to_tsv(full,tsv)
  table <- read.delim(tsv,check.names=FALSE,stringsAsFactors=FALSE)
  stopifnot(identical(names(table),c("cell_id","g1","g2","g3")),
    identical(table$cell_id,c("c1","c2")), max(abs(as.matrix(table[-1])-matrix(c(1,0,0,2,3,4),2))) == 0)
  # The optimized R exporter never densifies a sparse matrix, even one row.
  tsv_env <- environment(to_tsv)
  original_as_matrix <- as.matrix
  tsv_env$as.matrix <- function(x,...) {
    if (inherits(x,"Matrix")) stop("sparse matrix densification is forbidden")
    original_as_matrix(x,...)
  }
  gz <- file.path(args[1],paste0(type,"-normalized.tsv.gz"))
  to_tsv(full,gz,matrix="normalized",cells=c("c2","c1"),features=c("g3","g1"))
  rm("as.matrix",envir=tsv_env)
  con <- gzfile(gz,"rt"); table <- read.delim(con,check.names=FALSE); close(con)
  stopifnot(identical(names(table),c("cell_id","g3","g1")), identical(table$cell_id,c("c2","c1")),
    max(abs(as.matrix(table[-1])-log1p(matrix(c(4,3,0,1),2)))) < 1e-12)
  hash <- digest::digest(tsv,file=TRUE,algo="sha256")
  fails(to_tsv(full,tsv),"already exists")
  stopifnot(hash == digest::digest(tsv,file=TRUE,algo="sha256"))
  for (selection in list(character(),c("c1","c1"),"missing",NA_character_))
    fails(to_tsv(full,file.path(args[1],"bad.tsv"),cells=selection))
  fails(to_tsv(full,file.path(full,"new/out.tsv")),"outside")
  stopifnot(!dir.exists(file.path(full,"new")), !file.exists(file.path(args[1],"bad.tsv")))
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
# IDs containing quotes/Unicode must survive ordinary TSV quoting.
counts <- Matrix(matrix(c(1,0,3,0,2,4),3),sparse=TRUE)
dimnames(counts) <- list(c('gene"\u03b2','g2','g3'),c('cell"\u03b1','c2'))
quoted <- CreateSeuratObject(CreateAssay5Object(counts=counts))
only <- file.path(args[1],"quoted.scx")
write_scx(quoted,only,normalized=NULL)
quoted_tsv <- file.path(args[1],"quoted.tsv")
to_tsv(only,quoted_tsv)
table <- read.delim(quoted_tsv,check.names=FALSE,stringsAsFactors=FALSE,fileEncoding="UTF-8")
stopifnot(names(table)[2] == 'gene"\u03b2', table$cell_id[1] == 'cell"\u03b1')
fails(to_tsv(only,file.path(args[1],"no-normalized.tsv"),matrix="normalized"),"no normalized")
# A failed raw copy must retain the partial file and lock, never overwrite.
failed_tsv <- env$.scx_tsv_target(file.path(args[1],"failed-copy.tsv.gz"))
env$writeBin <- function(object,con,...) {
  if (identical(summary(con)$description,failed_tsv)) stop("injected TSV copy failure")
  base::writeBin(object,con,...)
}
fails(to_tsv(only,failed_tsv),"Incomplete TSV retained")
rm("writeBin",envir=env)
stopifnot(file.exists(failed_tsv), file.exists(env$.scx_lock_path(failed_tsv)))
# Sparse edge cases: empty rows/columns, cancellation, explicit signed zero,
# very small/large normalized values, and reordered selections.
edge_counts <- Matrix::sparseMatrix(i=c(1,4,2,2),j=c(2,2,4,4),x=c(1,5,2,-2),
  dims=c(4,5),dimnames=list(paste0("g",1:4),paste0("c",1:5)))
edge_norm <- Matrix::sparseMatrix(i=c(1,4,2,3,4),j=c(2,2,4,4,5),
  x=c(1e-300,-1e300,pi,-0,.Machine$double.xmin),dims=dim(edge_counts),dimnames=dimnames(edge_counts))
edge <- CreateSeuratObject(CreateAssay5Object(counts=edge_counts))
LayerData(edge,layer="data") <- edge_norm
edge_root <- file.path(args[1],"sparse-edge.scx")
write_scx(edge,edge_root)
row_ids <- c("c5","c2","c1","c4","c3"); col_ids <- c("g4","g1","g3","g2")
for (role in c("counts","normalized")) {
  checked <- env$.scx_inspect(edge_root,retain_expression=role)
  resource <- checked$manifest$path[checked$manifest$role==role]
  sparse_source <- checked$matrices[[resource]]
  oracle <- as.matrix(sparse_source)[match(row_ids,checked$cells),match(col_ids,checked$features),drop=FALSE]
  oracle <- matrix(sprintf("%.17g",oracle),nrow=length(row_ids))
  formatted <- 0L
  env$as.matrix <- function(x,...) {
    if (inherits(x,"Matrix")) stop("sparse matrix densification is forbidden")
    original_as_matrix(x,...)
  }
  env$sprintf <- function(fmt,...) {
    values <- list(...)
    if (identical(fmt,"%.17g")) formatted <<- formatted + length(values[[1]])
    do.call(base::sprintf,c(list(fmt),values))
  }
  edge_out <- file.path(args[1],paste0("sparse-edge-",role,".tsv.gz"))
  to_tsv(edge_root,edge_out,matrix=role,cells=row_ids,features=col_ids)
  rm(list=c("as.matrix","sprintf"),envir=env)
  con <- gzfile(edge_out,"rt"); table <- read.delim(con,check.names=FALSE,colClasses="character"); close(con)
  stopifnot(identical(table$cell_id,row_ids),identical(names(table)[-1],col_ids),
    identical(unname(as.matrix(table[-1])),unname(oracle)),
    formatted == length(sparse_source@x),formatted < length(row_ids)*length(col_ids))
}
zero_out <- file.path(args[1],"sparse-zero-column.tsv")
to_tsv(edge_root,zero_out,cells=c("c2","c1"),features="g3")
zero_table <- read.delim(zero_out,check.names=FALSE,colClasses="character")
stopifnot(identical(zero_table$cell_id,c("c2","c1")),identical(zero_table$g3,c("0","0")))
stopifnot(!"scxSeurat" %in% loadedNamespaces())
report <- list(passed=TRUE, installed_scx_loads=0, package_installations=0,
  external_process_calls=0, assays=c("Assay","Assay5"), shape=c(2,3), typed_metadata=TRUE,
  overwrite_protection=TRUE, interrupted_copy_rejected=TRUE, gzip_crc_and_truncation_rejected=TRUE,
  TSV=TRUE, TSV_no_dense_matrix=TRUE, TSV_interrupted_copy_rejected=TRUE)
writeLines(jsonlite::toJSON(report,auto_unbox=TRUE,pretty=TRUE), file.path(args[1],"report.json"))
cat("Standalone R validation passed\n")
