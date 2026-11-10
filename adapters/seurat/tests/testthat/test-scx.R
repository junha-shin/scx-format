tiny_seurat <- function(assay_class = "Assay5") {
  counts <- Matrix::Matrix(matrix(c(1, 0, 3, 0, 2, 4), 3, 2), sparse = TRUE)
  dimnames(counts) <- list(c("g1", "g2", "g3"), c("c1", "c2"))
  assay <- if (assay_class == "Assay") SeuratObject::CreateAssayObject(counts = counts) else SeuratObject::CreateAssay5Object(counts = counts)
  object <- SeuratObject::CreateSeuratObject(assay)
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

for (assay_class in c("Assay", "Assay5")) testthat::test_that(paste(assay_class, "exports the same canonical data"), {
  object <- tiny_seurat(assay_class)
  testthat::expect_s4_class(object[["RNA"]], assay_class)
  root <- tempfile(); on.exit(unlink(root, recursive = TRUE))
  testthat::expect_warning(write_scx_seurat(object, root), "SCX exports")
  con <- gzfile(file.path(root, "matrices/counts.mtx.gz"), "rt")
  x <- Matrix::readMM(con); close(con)
  testthat::expect_equal(as.matrix(x), unname(as.matrix(Matrix::t(SeuratObject::LayerData(object, layer = "counts")))))
})

for (assay_class in c("Assay", "Assay5")) testthat::test_that(paste(assay_class, "supports explicit counts-only selection"), {
  object <- tiny_seurat(assay_class)
  root <- tempfile(); on.exit(unlink(root, recursive = TRUE))
  testthat::expect_warning(write_scx_seurat(object, root, normalized_layer = NULL,
                                          reductions = character(), graphs = character()), "SCX exports")
  manifest <- read.delim(file.path(root, "manifest.tsv"), skip = 1)
  testthat::expect_setequal(manifest$role, c("cells", "features", "counts", "metadata_schema"))
})

testthat::test_that("bulk metadata coding matches scalar JSON semantics", {
  vectors <- list(c(1L, NA_integer_, 2147483647L), c(0.12345678901234567, NA, 1.234567890123456e-30),
                  c(TRUE, NA, FALSE), c("", "null", NA_character_, "quote\" slash\\ tab\t line\n", "null"))
  for (values in vectors) {
    reference <- lapply(seq_along(values), function(i) jsonlite::fromJSON(.scx_json(values[i]), simplifyVector = FALSE))
    encoded <- .scx_column_json(values)
    restored <- jsonlite::fromJSON(paste0("[", paste(encoded, collapse = ","), "]"), simplifyVector = FALSE)
    testthat::expect_identical(restored, reference)
  }
})

testthat::test_that("sparse finite checks inspect only stored values", {
  n <- 100000L
  x <- Matrix::sparseMatrix(i = seq_len(n), j = seq_len(n), x = 1, dims = c(n, n))
  testthat::expect_length(.scx_values(x), n)
  testthat::expect_true(all(is.finite(.scx_values(x))))
  x@x[1] <- NaN
  testthat::expect_true(any(!is.finite(.scx_values(x))))
})

testthat::test_that("validated replacement and failed writes preserve existing data", {
  object <- tiny_seurat()
  root <- tempfile(); on.exit(unlink(root, recursive = TRUE))
  testthat::expect_warning(write_scx_seurat(object, root), "SCX exports")
  snapshot <- function() {
    files <- list.files(root, recursive = TRUE, full.names = TRUE)
    lapply(files, function(path) readBin(path, "raw", n = file.info(path)$size))
  }
  before <- snapshot()
  testthat::expect_error(write_scx_seurat(object, root), "already exists")
  testthat::expect_identical(snapshot(), before)
  original <- .scx_mtx
  assign(".scx_mtx", function(x, path) stop("simulated disk full"), envir = environment(write_scx_seurat))
  on.exit(assign(".scx_mtx", original, envir = environment(write_scx_seurat)), add = TRUE)
  testthat::expect_error(write_scx_seurat(object, root, overwrite = TRUE), "simulated disk full")
  testthat::expect_identical(snapshot(), before)
  assign(".scx_mtx", original, envir = environment(write_scx_seurat))
  object$integer <- c(9L, NA_integer_)
  testthat::expect_warning(write_scx_seurat(object, root, overwrite = TRUE), "SCX exports")
  schema <- jsonlite::fromJSON(file.path(root, "metadata.schema.json"), simplifyVector = FALSE)
  decoded <- .scx_decode_frame(.scx_read_tsv(file.path(root, "cells.tsv.gz")), schema$cells)
  testthat::expect_identical(decoded$integer, c(9L, NA_integer_))
})

testthat::test_that("Python validation rejects corrupt R output before replacement", {
  object <- tiny_seurat()
  root <- tempfile(); on.exit(unlink(root, recursive = TRUE))
  testthat::expect_warning(write_scx_seurat(object, root), "SCX exports")
  path <- file.path(root, "matrices/counts.mtx.gz")
  before <- readBin(path, "raw", n = file.info(path)$size)
  original <- .scx_mtx
  assign(".scx_mtx", function(x, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines("broken gzip", path)
  }, envir = environment(write_scx_seurat))
  on.exit(assign(".scx_mtx", original, envir = environment(write_scx_seurat)), add = TRUE)
  testthat::expect_error(suppressWarnings(write_scx_seurat(object, root, overwrite = TRUE)), "counts.mtx.gz")
  testthat::expect_identical(readBin(path, "raw", n = file.info(path)$size), before)
})

testthat::test_that("import does not add QC columns absent from source metadata", {
  object <- tiny_seurat()
  object@meta.data <- object[[]][, c("integer", "logical", "category", "text"), drop = FALSE]
  root <- tempfile(); on.exit(unlink(root, recursive = TRUE))
  testthat::expect_warning(write_scx_seurat(object, root), "SCX exports")
  restored <- read_scx_seurat(root)
  testthat::expect_identical(restored[[]], object[[]])
})
