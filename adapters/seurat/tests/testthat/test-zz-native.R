native_export <- function(root) suppressWarnings(write_scx_seurat(tiny_seurat(), root))
native_bytes <- function(path) readBin(path, "raw", n = file.info(path)$size)
native_snapshot <- function(root) {
  paths <- list.files(root, recursive = TRUE, full.names = TRUE)
  setNames(lapply(paths, native_bytes), substring(paths, nchar(root) + 2L))
}
native_json <- function(root, edit) {
  path <- file.path(root, "metadata.schema.json")
  schema <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  writeLines(.scx_json(edit(schema)), path)
}

for (class in c("Assay", "Assay5")) testthat::test_that(paste("native R round trip", class), {
  before <- tiny_seurat(class)
  root <- tempfile(); on.exit(unlink(root, recursive = TRUE))
  suppressWarnings(write_scx_seurat(before, root, python = "/nonexistent/python"))
  testthat::expect_equal(validate_scx(root), list(cells = 2L, features = 3L, files = 7L))
  after <- read_scx_seurat(root, assay_class = class, python = "/nonexistent/python")
  testthat::expect_identical(before[[]], after[[]])
  testthat::expect_identical(before[["RNA"]][[]], after[["RNA"]][[]])
  for (layer in c("counts", "data")) testthat::expect_equal(
    SeuratObject::LayerData(before, layer = layer), SeuratObject::LayerData(after, layer = layer), tolerance = 0)
  testthat::expect_equal(unname(SeuratObject::Embeddings(before[["pca"]])), unname(SeuratObject::Embeddings(after[["pca"]])), tolerance = 0)
  testthat::expect_equal(as.matrix(before[["RNA_nn"]]), as.matrix(after[["RNA_nn"]]), tolerance = 0)
})

for (damage in c("truncate", "crc", "deflate", "plain", "extra", "short", "index", "shape", "nonfinite", "header"))
testthat::test_that(paste("native validation rejects damaged MTX", damage), {
  root <- tempfile(); on.exit(unlink(root, recursive = TRUE)); native_export(root)
  path <- file.path(root, "matrices/counts.mtx.gz")
  bytes <- native_bytes(path)
  if (damage == "truncate") bytes <- head(bytes, -5L)
  if (damage == "crc") bytes[length(bytes) - 7L] <- as.raw(bitwXor(as.integer(bytes[length(bytes) - 7L]), 1L))
  if (damage == "deflate") bytes <- c(as.raw(c(31, 139, 8, 0, 0, 0, 0, 0, 2, 255)), rep(as.raw(255), 20))
  if (damage == "plain") bytes <- charToRaw("not gzip")
  if (damage %in% c("truncate", "crc", "deflate", "plain")) writeBin(bytes, path) else {
    text <- switch(damage,
      extra = "2 3 1\n1 1 1\n2 2 2",
      short = "2 3 2\n1 1 1",
      index = "2 3 1\n3 1 1",
      shape = "3 2 1\n1 1 1",
      nonfinite = "2 3 1\n1 1 NaN",
      header = "2 3\n1\n")
    con <- gzfile(path, "wt")
    writeLines(c(if (damage == "header") "%%MatrixMarket matrix array real general" else "%%MatrixMarket matrix coordinate real general", text), con)
    close(con)
  }
  testthat::expect_error(validate_scx(root), "counts.mtx.gz")
  testthat::expect_error(read_scx_seurat(root), "counts.mtx.gz")
})

for (damage in c("version", "width", "duplicate", "unsafe", "missing", "ids", "embedding", "schema", "bool", "overflow", "category"))
testthat::test_that(paste("native structural/typed validation", damage), {
  root <- tempfile(); on.exit(unlink(root, recursive = TRUE)); native_export(root)
  manifest <- file.path(root, "manifest.tsv")
  text <- readLines(manifest)
  if (damage == "version") writeLines(sub("0.1-alpha", "99", text), manifest)
  if (damage == "width") writeLines(c(text, "bad\tembedding\tx\textra"), manifest)
  if (damage == "duplicate") writeLines(c(text, text[3]), manifest)
  if (damage == "unsafe") writeLines(c(text, "../escape\tembedding\tx"), manifest)
  if (damage == "missing") unlink(file.path(root, "features.tsv.gz"))
  if (damage == "ids") .scx_tsv(data.frame(cell_id = c("same", "same")), file.path(root, "cells.tsv.gz"))
  if (damage == "embedding") .scx_tsv(data.frame(cell_id = c("c2", "c1"), dim_1 = c(1, 2)), file.path(root, "embeddings/pca.tsv.gz"))
  if (damage == "schema") writeLines("{", file.path(root, "metadata.schema.json"))
  if (damage %in% c("bool", "overflow", "category")) {
    path <- file.path(root, "cells.tsv.gz"); table <- .scx_read_tsv(path)
    if (damage == "bool") table$logical[1] <- "1"
    if (damage == "overflow") table$integer[1] <- "2147483648"
    if (damage == "category") table$category[1] <- "99"
    .scx_tsv(table, path)
  }
  testthat::expect_error(validate_scx(root))
})

testthat::test_that("gzip members, trailing garbage and malformed table quoting are checked", {
  path <- tempfile(fileext = ".gz"); on.exit(unlink(path))
  writeBin(c(memCompress(charToRaw("cell_id\n"), "gzip"), memCompress(charToRaw("c1\n"), "gzip")), path)
  # memCompress(type=gzip) uses zlib, not a gzip member: it must be rejected.
  testthat::expect_error(.scx_read_tsv(path), "gzip")
  con <- gzfile(path, "wt"); writeLines(c("cell_id", "c1"), con); close(con)
  member <- native_bytes(path)
  writeBin(c(member, member), path)
  testthat::expect_equal(.scx_read_tsv(path)$cell_id, c("c1", "cell_id", "c1"))
  writeBin(c(member, charToRaw("garbage")), path)
  testthat::expect_error(.scx_read_tsv(path), "gzip")
  con <- gzfile(path, "wt"); writeLines(c("cell_id\tx", "c1\t\"unterminated"), con); close(con)
  testthat::expect_error(.scx_read_tsv(path))
})

testthat::test_that("failed native publication restores old output", {
  root <- tempfile(); on.exit(unlink(root, recursive = TRUE)); native_export(root)
  before <- native_snapshot(root)
  original <- .scx_rename_new
  assign(".scx_rename_new", function(source, target) {
    if (basename(source) == "payload") stop("publication denied")
    original(source, target)
  }, envir = environment(write_scx_seurat))
  on.exit(assign(".scx_rename_new", original, envir = environment(write_scx_seurat)), add = TRUE)
  testthat::expect_error(suppressWarnings(write_scx_seurat(tiny_seurat(), root, overwrite = TRUE)), "publication denied")
  testthat::expect_identical(native_snapshot(root), before)
  testthat::expect_false(file.exists(.scx_lock_path(root)))
})

testthat::test_that("failed restoration retains backup and lock", {
  parent <- tempfile(); dir.create(parent); on.exit(unlink(parent, recursive = TRUE))
  root <- file.path(parent, "out"); native_export(root); before <- native_snapshot(root)
  original <- .scx_rename_new
  assign(".scx_rename_new", function(source, target) {
    if (basename(source) %in% c("payload", "previous")) stop("destination unavailable")
    original(source, target)
  }, envir = environment(write_scx_seurat))
  on.exit(assign(".scx_rename_new", original, envir = environment(write_scx_seurat)), add = TRUE)
  testthat::expect_error(suppressWarnings(write_scx_seurat(tiny_seurat(), root, overwrite = TRUE)), "backup retained")
  record <- jsonlite::fromJSON(.scx_lock_path(root))
  testthat::expect_identical(native_snapshot(record$backup), before)
  testthat::expect_error(suppressWarnings(write_scx_seurat(tiny_seurat(), root)), "locked")
})

testthat::test_that("a racing empty directory is not replaced", {
  parent <- tempfile(); dir.create(parent); on.exit(unlink(parent, recursive = TRUE))
  root <- file.path(parent, "out"); original <- .scx_rename_new
  assign(".scx_rename_new", function(source, target) {
    if (basename(source) == "payload") dir.create(target)
    original(source, target)
  }, envir = environment(write_scx_seurat))
  on.exit(assign(".scx_rename_new", original, envir = environment(write_scx_seurat)), add = TRUE)
  testthat::expect_error(suppressWarnings(write_scx_seurat(tiny_seurat(), root)), "rename failed")
  testthat::expect_true(dir.exists(root)); testthat::expect_length(list.files(root), 0L)
})

testthat::test_that("interrupt during publication restores previous output", {
  root <- tempfile(); on.exit(unlink(root, recursive = TRUE)); native_export(root)
  before <- native_snapshot(root); original <- .scx_rename_new
  assign(".scx_rename_new", function(source, target) {
    if (basename(source) == "payload") stop(structure(list(message = "interrupted"), class = c("interrupt", "condition")))
    original(source, target)
  }, envir = environment(write_scx_seurat))
  on.exit(assign(".scx_rename_new", original, envir = environment(write_scx_seurat)), add = TRUE)
  result <- tryCatch(suppressWarnings(write_scx_seurat(tiny_seurat(), root, overwrite = TRUE)), interrupt = function(e) "caught")
  testthat::expect_identical(result, "caught")
  testthat::expect_identical(native_snapshot(root), before)
})

testthat::test_that("process termination in replacement gap keeps recoverable backup", {
  parent <- tempfile(); dir.create(parent); on.exit(unlink(parent, recursive = TRUE))
  root <- file.path(parent, "out"); native_export(root); before <- native_snapshot(root)
  script <- file.path(parent, "crash.R")
  writeLines(c('library(scxSeurat)', 'root <- commandArgs(TRUE)[1]',
    'object <- read_scx_seurat(root)',
    'trace(".scx_rename_new", where=asNamespace("scxSeurat"), print=FALSE, tracer=quote(if(basename(source)=="payload") quit(save="no",status=19L,runLast=FALSE)))',
    'write_scx_seurat(object,root,overwrite=TRUE)'), script)
  status <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    c(shQuote(script), shQuote(root)), stdout = file.path(parent, "log"), stderr = file.path(parent, "errors")))
  testthat::expect_equal(status, 19L)
  record <- jsonlite::fromJSON(.scx_lock_path(root))
  testthat::expect_false(dir.exists(root))
  testthat::expect_identical(native_snapshot(record$backup), before)
  testthat::expect_equal(validate_scx(record$backup)$cells, 2L)
})

testthat::test_that("invalid UTF-8 tables are rejected with file context", {
  root <- tempfile(); on.exit(unlink(root, recursive = TRUE)); native_export(root)
  path <- file.path(root, "cells.tsv.gz")
  con <- gzfile(path, "wb"); writeBin(c(charToRaw("cell_id\n"), as.raw(255), charToRaw("\n")), con); close(con)
  testthat::expect_error(validate_scx(root), "cells.tsv.gz")
})

testthat::test_that("escaping resource and output symlinks are refused", {
  parent <- tempfile(); dir.create(parent); on.exit(unlink(parent, recursive = TRUE))
  root <- file.path(parent, "out"); native_export(root)
  matrix <- file.path(root, "matrices/counts.mtx.gz")
  outside <- file.path(parent, "outside.gz"); file.copy(matrix, outside); unlink(matrix)
  linked <- suppressWarnings(file.symlink(outside, matrix))
  if (!linked) testthat::skip("account cannot create symlinks")
  testthat::expect_error(validate_scx(root), "escapes")
  output <- file.path(parent, "link")
  if (!suppressWarnings(file.symlink(root, output))) testthat::skip("account cannot create directory symlinks")
  testthat::expect_error(suppressWarnings(write_scx_seurat(tiny_seurat(), output, overwrite = TRUE)), "symlinks are refused")
})
