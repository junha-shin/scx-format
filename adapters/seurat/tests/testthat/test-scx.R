tiny_seurat <- function() {
  counts <- Matrix::Matrix(matrix(c(1, 0, 3, 0, 2, 4), 3, 2), sparse = TRUE)
  dimnames(counts) <- list(c("g1", "g2", "g3"), c("c1", "c2"))
  object <- SeuratObject::CreateSeuratObject(counts)
  SeuratObject::LayerData(object, layer = "data") <- log1p(counts)
  object$integer <- c(1L, NA_integer_)
  object$double <- c(0.12345678901234567, 1.234567890123456e-30)
  object$logical <- c(TRUE, NA)
  object$category <- ordered(c("b", NA), levels = c("b", "a", "unused"))
  object$text <- c("", "null")
  object[["RNA"]] <- SeuratObject::AddMetaData(object[["RNA"]], data.frame(score = c(1.2, NA, 3.4), row.names = rownames(counts)))
  embedding <- matrix(c(1.23456789, 2, 3, 4), 2, dimnames = list(colnames(counts), c("PC_1", "PC_2")))
  object[["pca"]] <- SeuratObject::CreateDimReducObject(embedding, assay = "RNA", key = "PC_")
  graph <- Matrix::Matrix(matrix(c(0, 2, 1, 0), 2), sparse = TRUE)
  dimnames(graph) <- list(colnames(counts), colnames(counts))
  graph <- SeuratObject::as.Graph(graph)
  SeuratObject::DefaultAssay(graph) <- "RNA"
  object[["RNA_nn"]] <- graph
  object
}

testthat::test_that("export is cells by features with typed schema", {
  object <- tiny_seurat()
  root <- tempfile(); on.exit(unlink(root, recursive = TRUE))
  testthat::expect_warning(write_scx_seurat(object, root), "SCX exports")
  con <- gzfile(file.path(root, "matrices/counts.mtx.gz"), "rt")
  x <- Matrix::readMM(con); close(con)
  testthat::expect_equal(as.matrix(x), unname(as.matrix(Matrix::t(SeuratObject::LayerData(object, layer = "counts")))))
  schema <- jsonlite::fromJSON(file.path(root, "metadata.schema.json"), simplifyVector = FALSE)
  types <- setNames(vapply(schema$cells$columns, `[[`, character(1), "dtype"), names(object[[]]))
  testthat::expect_equal(unname(types[c("integer", "logical", "category", "text")]), c("Int32", "boolean", "category", "object"))
  decoded <- .scx_decode_frame(.scx_read_tsv(file.path(root, "cells.tsv.gz")), schema$cells)
  testthat::expect_identical(decoded, object[[]])
})

testthat::test_that("existing output and partial axes fail without publishing", {
  object <- tiny_seurat()
  testthat::expect_error(write_scx_seurat(object, tempdir()), "already exists")
  root <- tempfile()
  testthat::expect_error(write_scx_seurat(object, root, counts_layer = "counts.missing"), "exact counts")
  testthat::expect_false(dir.exists(root))
  testthat::expect_error(.scx_aligned(matrix(1, 1, 1, dimnames = list("c1", "g1")), c("c1", "c2")), "row axis")
  object$bad <- as.Date(c("2020-01-01", "2020-01-02"))
  testthat::expect_error(write_scx_seurat(object, root), "unsupported R metadata")
  testthat::expect_false(dir.exists(root))
})

testthat::test_that("native R dtype limits fail explicitly", {
  frame <- data.frame(cell_id = c("c1", "c2"), x = c("1", "null"))
  spec <- list(columns = list(list(name = "x", dtype = "Int64")))
  testthat::expect_error(.scx_decode_frame(frame, spec), "not safely representable")
  object <- tiny_seurat()
  object$bad_factor <- factor(c("a", NA), exclude = NULL)
  root <- tempfile()
  testthat::expect_error(write_scx_seurat(object, root), "factor levels must be nonmissing")
  testthat::expect_false(dir.exists(root))
})

testthat::test_that("symmetric and binary graphs use explicit numeric general format", {
  for (x in list(Matrix::Matrix(matrix(1, 2, 2), sparse = TRUE), Matrix::Diagonal(2))) {
    path <- tempfile(fileext = ".mtx.gz"); on.exit(unlink(path), add = TRUE)
    .scx_mtx(x, path)
    con <- gzfile(path, "rt")
    testthat::expect_identical(readLines(con, n = 1), "%%MatrixMarket matrix coordinate real general")
    close(con)
    con <- gzfile(path, "rt"); restored <- Matrix::readMM(con); close(con)
    testthat::expect_equal(as.matrix(restored), as.matrix(x), ignore_attr = TRUE)
  }
})

testthat::test_that("JSON metadata preserves binary64 precision", {
  x <- c(0.12345678901234567, 1.234567890123456e-30, 1.234567890123456e30)
  testthat::expect_identical(jsonlite::fromJSON(.scx_json(x)), x)
})
