# SCX v0.1-alpha, single-assay SeuratObject >= 5 adapter.
# Development source() retains access to the installed package's native helpers.
.scx_ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (!isNamespace(environment()) && !is.null(.scx_ofile) && !exists(".scx_inspect", envir = environment(), inherits = FALSE)) {
  .scx_source_dir <- dirname(.scx_ofile)
  for (.scx_support in c("publication.R", "validation.R", "api.R"))
    sys.source(file.path(.scx_source_dir, .scx_support), envir = environment())
  rm(.scx_source_dir, .scx_support)
}
rm(.scx_ofile)
.scx_json <- function(x) as.character(jsonlite::toJSON(x, auto_unbox = TRUE, na = "null", null = "null", digits = 17))

.scx_column_json <- function(values) {
  if (is.integer(values)) {
    encoded <- as.character(values)
  } else if (is.double(values)) {
    encoded <- sprintf("%.17g", values)
  } else if (is.logical(values)) {
    encoded <- ifelse(values, "true", "false")
  } else {
    # Repeated text is encoded once per unique value; jsonlite still handles
    # escaping rather than treating R string escapes as JSON escapes.
    unique_values <- unique(values)
    encoded_unique <- vapply(seq_along(unique_values), function(i) .scx_json(unique_values[i]), character(1))
    encoded <- encoded_unique[match(values, unique_values)]
  }
  encoded[is.na(values)] <- "null"
  unname(encoded)
}

.scx_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, "wt", encoding = "UTF-8")
  on.exit(close(con))
  write.table(x, con, sep = "\t", row.names = FALSE, quote = TRUE, qmethod = "double", na = "")
}

.scx_read_tsv <- function(path) {
  .scx_table(path)
}

.scx_mtx <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  x <- methods::as(Matrix::Matrix(x, sparse = TRUE), "generalMatrix")
  x <- methods::as(methods::as(x, "dMatrix"), "TsparseMatrix")
  output <- gzfile(path, "wt", encoding = "UTF-8"); on.exit(close(output))
  # Matrix::writeMM auto-selects pattern/symmetric even for general matrices.
  # SCX requires numeric coordinate general, including all directed entries.
  writeLines(c("%%MatrixMarket matrix coordinate real general", "% SCX CELL x FEATURE or CELL x CELL",
               paste(nrow(x), ncol(x), length(x@x))), output)
  if (length(x@x)) for (start in seq.int(1L, length(x@x), by = 100000L)) {
    take <- start:min(start + 99999L, length(x@x))
    writeLines(paste(x@i[take] + 1L, x@j[take] + 1L, sprintf("%.17g", x@x[take])), output)
  }
}

.scx_ids <- function(ids) {
  if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) ||
      any(grepl("[\t\r\n]", ids))) stop("IDs must be unique nonempty text without tabs/newlines")
}

.scx_name <- function(name) {
  if (length(name) != 1L || !grepl("^[A-Za-z0-9_.-]+$", name)) stop("unsafe SCX object name")
}

.scx_meta <- function(frame, ids, key, path) {
  if (!is.data.frame(frame) || !identical(rownames(frame), ids)) stop("metadata cell/feature axis mismatch")
  if (anyDuplicated(names(frame)) || any(!nzchar(names(frame))) || key %in% names(frame) ||
      any(grepl("[\t\r\n]", names(frame)))) stop("invalid or reserved metadata column name")
  table <- data.frame(ids, check.names = FALSE); names(table) <- key
  columns <- lapply(names(frame), function(name) {
    x <- frame[[name]]
    if (is.factor(x)) {
      if (anyNA(levels(x))) stop("factor levels must be nonmissing")
      spec <- list(name = name, dtype = "category", ordered = is.ordered(x),
                   categories_dtype = "object", categories = as.list(levels(x)))
      values <- as.integer(x) - 1L
    } else {
      dtype <- if (is.integer(x) && !is.object(x)) {
        if (anyNA(x)) "Int32" else "int32"
      } else if (is.double(x) && !is.object(x)) {
        if (any(is.infinite(x))) stop("metadata infinity unsupported")
        "float64"
      } else if (is.logical(x) && !is.object(x)) {
        if (anyNA(x)) "boolean" else "bool"
      } else if (is.character(x) && !is.object(x)) "object" else stop("unsupported R metadata class: ", paste(class(x), collapse = "/"))
      spec <- list(name = name, dtype = dtype)
      values <- x
    }
    table[[name]] <<- .scx_column_json(values)
    spec
  })
  .scx_tsv(table, path)
  list(index_name = NULL, columns = columns)
}

.scx_aligned <- function(x, rows, cols = NULL) {
  if (!setequal(rownames(x), rows) || nrow(x) != length(rows)) stop("row axis does not cover canonical IDs")
  if (is.null(cols)) return(x[rows, , drop = FALSE])
  if (!setequal(colnames(x), cols) || ncol(x) != length(cols)) stop("column axis does not cover canonical IDs")
  x[rows, cols, drop = FALSE]
}

.scx_values <- function(x) {
  # Implicit sparse zeros are finite and nonnegative. Checking the whole sparse
  # matrix with is.finite() can materialize a dense CELL x CELL logical matrix.
  if (inherits(x, "sparseMatrix")) {
    if ("x" %in% methods::slotNames(x)) return(methods::slot(x, "x"))
    return(1)  # Pattern matrix: stored entries are implicit ones.
  }
  x
}

#' Export a single Seurat assay with validated, rollback-safe publication.
#' @export
write_scx_seurat <- function(object, path, assay = SeuratObject::DefaultAssay(object),
                             counts_layer = "counts", normalized_layer = "data",
                             reductions = SeuratObject::Reductions(object),
                             graphs = SeuratObject::Graphs(object), overwrite = FALSE,
                             python = NULL, warn_scope = TRUE) {
  if (!inherits(object, "Seurat")) stop("object must be Seurat")
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) stop("overwrite must be TRUE or FALSE")
  if ((file.exists(path) || dir.exists(path)) && !overwrite) stop("output already exists; choose a new path or use overwrite = TRUE")
  # Kept for source compatibility; no Python executable is inspected or invoked.
  if (!assay %in% SeuratObject::Assays(object)) stop("unknown assay")
  layers <- SeuratObject::Layers(object[[assay]], search = NA)
  if (!counts_layer %in% layers) stop("exact counts layer required; join split layers explicitly with JoinLayers first")
  counts <- SeuratObject::LayerData(object[[assay]], layer = counts_layer)
  features <- rownames(counts); cells <- colnames(counts)
  .scx_ids(cells); .scx_ids(features)
  if (!length(cells) || !length(features)) stop("empty expression matrix")
  values <- .scx_values(counts)
  if (is.complex(values) || any(!is.finite(values)) || any(values < 0)) stop("counts must be finite real and nonnegative")
  rm(values)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path <- file.path(normalizePath(dirname(path), winslash = "/", mustWork = TRUE), basename(path))
  work <- tempfile(".scx-stage-", tmpdir = dirname(path)); dir.create(work)
  stage <- file.path(work, "payload"); dir.create(stage)
  on.exit(if (!dir.exists(file.path(work, "previous"))) unlink(work, recursive = TRUE), add = TRUE)
  entries <- data.frame(path = c("cells.tsv.gz", "features.tsv.gz", "matrices/counts.mtx.gz", "metadata.schema.json"),
                        role = c("cells", "features", "counts", "metadata_schema"),
                        name = c("cells", "features", "counts", "metadata"))
  schema <- list(schema_version = "1.0",
                 cells = .scx_meta(object[[]][cells, , drop = FALSE], cells, "cell_id", file.path(stage, "cells.tsv.gz")),
                 features = .scx_meta(object[[assay]][[]][features, , drop = FALSE], features, "feature_id", file.path(stage, "features.tsv.gz")))
  writeLines(.scx_json(schema), file.path(stage, "metadata.schema.json"), useBytes = TRUE)
  .scx_mtx(Matrix::t(counts), file.path(stage, "matrices/counts.mtx.gz"))
  if (!is.null(normalized_layer)) {
    if (normalized_layer %in% layers) {
      normalized <- .scx_aligned(SeuratObject::LayerData(object[[assay]], layer = normalized_layer), features, cells)
      values <- .scx_values(normalized)
      if (is.complex(values) || any(!is.finite(values))) stop("normalized matrix must be finite real")
      rm(values)
      .scx_mtx(Matrix::t(normalized), file.path(stage, "matrices/normalized.mtx.gz"))
      entries <- rbind(entries, data.frame(path = "matrices/normalized.mtx.gz", role = "normalized", name = "normalized"))
    } else warning("normalized layer absent; exporting counts only")
  }
  for (name in reductions) {
    .scx_name(name)
    reduction <- object[[name]]
    if (SeuratObject::DefaultAssay(reduction) != assay) stop("embedding belongs to a different assay: ", name)
    x <- .scx_aligned(SeuratObject::Embeddings(reduction), cells)
    if (!ncol(x) || !is.numeric(x) || any(!is.finite(x))) stop("invalid embedding")
    table <- data.frame(cell_id = cells, x, check.names = FALSE)
    names(table) <- c("cell_id", paste0("dim_", seq_len(ncol(x))))
    relative <- paste0("embeddings/", name, ".tsv.gz")
    # Explicit 17-digit formatting before write.table (which formats doubles differently).
    table[-1] <- lapply(table[-1], function(v) sprintf("%.17g", v))
    .scx_tsv(table, file.path(stage, relative))
    entries <- rbind(entries, data.frame(path = relative, role = "embedding", name = name))
  }
  for (name in graphs) {
    .scx_name(name)
    graph <- object[[name]]
    if (SeuratObject::DefaultAssay(graph) != assay) stop("graph belongs to a different assay: ", name)
    x <- .scx_aligned(graph, cells, cells)
    values <- .scx_values(x)
    if (is.complex(values) || any(!is.finite(values))) stop("graph must be finite real")
    rm(values)
    relative <- paste0("graphs/", name, ".mtx.gz")
    .scx_mtx(x, file.path(stage, relative))
    entries <- rbind(entries, data.frame(path = relative, role = "cell_graph", name = name))
  }
  con <- file(file.path(stage, "manifest.tsv"), "wt", encoding = "UTF-8")
  writeLines("# scx_version=0.1-alpha", con)
  write.table(entries, con, sep = "\t", row.names = FALSE, quote = FALSE)
  close(con)
  if (warn_scope) warning("SCX exports one assay and selected layers, metadata, embeddings and graphs; identities, other layers, loadings and analysis history are not exported")
  .scx_publish(work, path, overwrite)
  invisible(path)
}

#' Import an SCX dataset via the native R validator and typed metadata decoder.
#' @export
read_scx_seurat <- function(path, assay = "RNA", python = NULL, assay_class = c("Assay5", "Assay")) {
  assay_class <- match.arg(assay_class)
  checked <- .scx_inspect(path, retain = TRUE)
  manifest <- checked$manifest
  ids <- checked$cells; genes <- checked$features
  meta <- checked$meta; feature_meta <- checked$feature_meta
  load_matrix <- function(relative) checked$matrices[[relative]]
  counts <- Matrix::t(load_matrix("matrices/counts.mtx.gz")); dimnames(counts) <- list(genes, ids)
  assay_object <- if (assay_class == "Assay5") SeuratObject::CreateAssay5Object(counts = counts) else SeuratObject::CreateAssayObject(counts = counts)
  object <- SeuratObject::CreateSeuratObject(assay_object, assay = assay, meta.data = meta)
  if (!identical(rownames(object[[assay]]), genes)) stop("Seurat altered feature IDs; import refused")
  if (ncol(feature_meta)) object[[assay]] <- SeuratObject::AddMetaData(object[[assay]], feature_meta)
  for (i in seq_len(nrow(manifest))) {
    entry <- manifest[i, ]
    if (entry$role == "normalized") {
      x <- Matrix::t(load_matrix(entry$path)); dimnames(x) <- list(genes, ids)
      SeuratObject::LayerData(object, assay = assay, layer = "data") <- x
    } else if (entry$role == "embedding") {
      x <- checked$embeddings[[entry$path]]; rownames(x) <- ids
      key <- paste0("SCX", i, "_"); colnames(x) <- paste0(key, seq_len(ncol(x)))
      object[[entry$name]] <- SeuratObject::CreateDimReducObject(x, assay = assay, key = key)
    } else if (entry$role == "cell_graph") {
      x <- load_matrix(entry$path); dimnames(x) <- list(ids, ids)
      graph <- SeuratObject::as.Graph(x); SeuratObject::DefaultAssay(graph) <- assay
      object[[entry$name]] <- graph
    }
  }
  # Assay/feature-metadata setters may regenerate QC columns. Restore the
  # validated SCX table after every setter to retain exactly the source columns.
  object@meta.data <- meta
  object
}

.scx_decode_frame <- function(table, spec) {
  frame <- data.frame(row.names = table[[1]])
  if (is.null(spec)) {
    for (name in names(table)[-1]) frame[[name]] <- table[[name]]
    return(frame)
  }
  convert <- function(values, dtype) {
    missing <- vapply(values, is.null, logical(1))
    if (dtype %in% c("object", "str", "string")) {
      x <- rep(NA_character_, length(values))
    } else if (dtype %in% c("bool", "boolean")) {
      x <- rep(NA, length(values))
    } else if (dtype %in% c("int8", "int16", "int32", "Int8", "Int16", "Int32", "uint8", "uint16", "UInt8", "UInt16")) {
      x <- rep(NA_integer_, length(values))
      if (any(vapply(values[!missing], function(v) v <= -.Machine$integer.max - 1, logical(1)))) stop("R cannot represent INT_MIN without losing missingness")
    } else if (dtype %in% c("float32", "float64", "Float32", "Float64")) {
      x <- rep(NA_real_, length(values))
    } else stop("SCX dtype not safely representable in native R: ", dtype)
    x[!missing] <- unlist(values[!missing], use.names = FALSE)
    if (is.integer(x)) x <- as.integer(x)
    x
  }
  for (column in spec$columns) {
    values <- jsonlite::fromJSON(paste0("[", paste(table[[column$name]], collapse = ","), "]"), simplifyVector = FALSE)
    if (column$dtype == "category") {
      levels <- convert(column$categories, column$categories_dtype)
      codes <- vapply(values, function(x) if (is.null(x)) NA_integer_ else as.integer(x) + 1L, integer(1))
      # R factors have character levels; numeric categories cannot retain their type.
      if (!is.character(levels)) stop("numeric/boolean categorical levels are unsupported by R factors")
      x <- structure(codes, levels = levels, class = if (column$ordered) c("ordered", "factor") else "factor")
    } else x <- convert(values, column$dtype)
    frame[[column$name]] <- x
  }
  frame
}
