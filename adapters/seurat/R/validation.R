.scx_resource <- function(path, code) {
  tryCatch(withCallingHandlers(force(code), warning = function(w) stop(conditionMessage(w), call. = FALSE)),
           error = function(e) stop(path, ": ", conditionMessage(e), call. = FALSE))
}

.scx_safe_path <- function(root, relative) {
  if (length(relative) != 1L || is.na(relative) || !nzchar(relative) ||
      grepl("[\\\\:\t\r\n]", relative) || startsWith(relative, "/") ||
      any(strsplit(relative, "/", fixed = TRUE)[[1]] %in% c("", ".", "..")) || endsWith(relative, "/"))
    stop("unsafe manifest path: ", relative)
  target <- file.path(root, relative)
  resolved <- normalizePath(target, winslash = "/", mustWork = FALSE)
  base <- normalizePath(root, winslash = "/", mustWork = TRUE)
  if (.Platform$OS.type == "windows") { resolved <- tolower(resolved); base <- tolower(base) }
  if (!startsWith(resolved, paste0(base, "/"))) stop("manifest path escapes dataset: ", relative)
  target
}

.scx_table_lines <- function(lines) {
  if (!length(lines) || anyNA(iconv(lines, from = "UTF-8", to = "UTF-8", sub = NA))) stop("empty or invalid UTF-8 table")
  con <- textConnection(lines); on.exit(close(con))
  width <- count.fields(con, sep = "\t", quote = '"', comment.char = "", blank.lines.skip = FALSE)
  if (anyNA(width) || any(width != width[1]) || width[1] < 1L) stop("invalid table width or quoted newline")
  table <- read.delim(text = lines, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
                      na.strings = NULL, comment.char = "", fill = FALSE, blank.lines.skip = FALSE)
  if (anyDuplicated(names(table)) || any(!nzchar(names(table))) ||
      any(grepl("[\t\r\n]", c(names(table), unlist(table, use.names = FALSE)))))
    stop("table headers must be nonempty/unique; literal tabs/newlines are unsupported")
  table
}

.scx_table <- function(path, compressed = TRUE) .scx_resource(path, {
  if (compressed) .scx_native("scx_check_gzip", path)
  con <- if (compressed) gzfile(path, "rt", encoding = "UTF-8") else file(path, "rt", encoding = "UTF-8")
  on.exit(close(con))
  .scx_table_lines(readLines(con, warn = FALSE))
})

.scx_manifest <- function(root) .scx_resource(file.path(root, "manifest.tsv"), {
  path <- .scx_safe_path(root, "manifest.tsv")
  con <- file(path, "rt", encoding = "UTF-8"); on.exit(close(con))
  if (!identical(readLines(con, 1L, warn = FALSE), "# scx_version=0.1-alpha")) stop("missing or unsupported manifest version")
  manifest <- .scx_table_lines(readLines(con, warn = FALSE))
  if (!identical(names(manifest), c("path", "role", "name"))) stop("manifest header must be path, role, name")
  if (anyDuplicated(manifest$path) || anyDuplicated(manifest[c("role", "name")])) stop("duplicate manifest path or role/name")
  roles <- c("cells", "features", "counts", "normalized", "embedding", "cell_graph", "metadata_schema")
  if (any(!manifest$role %in% roles) || any(!grepl("^[A-Za-z0-9_.-]+$", manifest$name))) stop("invalid role or name")
  for (p in manifest$path) .scx_safe_path(root, p)
  for (role in c("cells", "features", "counts")) {
    row <- manifest[manifest$role == role, , drop = FALSE]
    expected <- switch(role, cells = "cells.tsv.gz", features = "features.tsv.gz", counts = "matrices/counts.mtx.gz")
    if (nrow(row) != 1L || row$path != expected || row$name != role) stop("required manifest entry invalid: ", role)
  }
  if (sum(manifest$role == "normalized") > 1L) stop("only one normalized matrix supported")
  row <- manifest[manifest$role == "metadata_schema", , drop = FALSE]
  if (nrow(row) > 1L || (nrow(row) && (row$path != "metadata.schema.json" || row$name != "metadata")))
    stop("invalid metadata schema manifest entry")
  manifest
})

.scx_matrix <- function(path, shape) .scx_resource(path, {
  .scx_native("scx_check_gzip", path)
  con <- gzfile(path, "rt", encoding = "UTF-8"); on.exit(close(con))
  header <- strsplit(tolower(trimws(readLines(con, 1L, warn = FALSE))), "[[:space:]]+")[[1]]
  if (length(header) != 5L || !identical(header[c(1, 2, 3, 5)], c("%%matrixmarket", "matrix", "coordinate", "general")) ||
      !header[4] %in% c("integer", "real")) stop("matrix must use Matrix Market coordinate real/integer general")
  repeat {
    line <- readLines(con, 1L, warn = FALSE)
    if (!length(line)) stop("missing matrix dimensions")
    if (nzchar(trimws(line)) && !startsWith(trimws(line), "%")) break
  }
  fields <- strsplit(trimws(line), "[[:space:]]+")[[1]]
  if (length(fields) != 3L || any(!grepl("^[0-9]+$", fields))) stop("invalid matrix dimensions")
  dims <- as.numeric(fields)
  if (any(!is.finite(dims)) || any(dims > .Machine$integer.max) || any(dims[1:2] != shape)) stop("matrix dimensions do not match canonical axes or exceed R limits")
  n <- dims[3]
  entries <- if (n) scan(con, what = list(i = integer(), j = integer(), x = numeric()), nmax = n,
                         quiet = TRUE, multi.line = FALSE, comment.char = "%") else list(i = integer(), j = integer(), x = numeric())
  if (length(entries$x) != n) stop("truncated MTX: declared entry count does not match data")
  repeat {
    remaining <- readLines(con, 10000L, warn = FALSE)
    if (!length(remaining)) break
    if (any(nzchar(trimws(remaining)) & !startsWith(trimws(remaining), "%"))) stop("extra entries after declared MTX count")
  }
  if (anyNA(entries$i) || anyNA(entries$j) || any(entries$i < 1L | entries$i > dims[1]) ||
      any(entries$j < 1L | entries$j > dims[2])) stop("MTX index out of bounds")
  if (any(!is.finite(entries$x))) stop("matrix values must be finite")
  if (header[4] == "integer" && any(entries$x != trunc(entries$x) | abs(entries$x) > 2^53 - 1))
    stop("integer MTX value cannot be represented exactly in native R")
  if (header[4] == "integer" && sum(abs(entries$x)) > 2^53 - 1) {
    bounds <- Matrix::sparseMatrix(i = entries$i, j = entries$j, x = abs(entries$x), dims = as.integer(dims[1:2]))
    if (any(.scx_values(bounds) > 2^53 - 1)) stop("integer duplicate group can lose native R precision")
    rm(bounds)
  }
  x <- Matrix::sparseMatrix(i = entries$i, j = entries$j, x = entries$x, dims = as.integer(dims[1:2]))
  if (any(!is.finite(.scx_values(x)))) stop("nonfinite summed matrix values")
  if (header[4] == "integer" && any(abs(.scx_values(x)) > 2^53 - 1)) stop("summed integer MTX value exceeds native R exact range")
  x
})

.scx_keys <- function(x, keys) {
  is.list(x) && !anyDuplicated(names(x)) && setequal(names(x), keys)
}
.scx_text <- function(x) is.character(x) && length(x) == 1L && !is.na(x)

.scx_check_values <- function(values, dtype, tokens = NULL) {
  supported <- c("object", "str", "string", "bool", "boolean", "float32", "Float32", "float64", "Float64",
                 "int8", "Int8", "int16", "Int16", "int32", "Int32", "uint8", "UInt8", "uint16", "UInt16")
  if (!.scx_text(dtype) || !dtype %in% supported) stop("SCX dtype not safely representable in native R: ", paste(dtype, collapse = "/"))
  missing <- vapply(values, is.null, logical(1))
  present <- values[!missing]
  scalar <- function(v) length(v) == 1L && !is.list(v) && !is.na(v)
  if (any(!vapply(present, scalar, logical(1)))) stop("metadata must contain JSON scalars or null")
  if (dtype %in% c("object", "str", "string")) {
    if (any(!vapply(present, is.character, logical(1)))) stop("string metadata contains a non-string value")
  } else if (dtype %in% c("bool", "boolean")) {
    if (any(!vapply(present, is.logical, logical(1))) || (dtype == "bool" && any(missing))) stop("boolean metadata requires true/false and valid nullability")
  } else {
    if (any(!vapply(present, is.numeric, logical(1)))) stop("numeric metadata requires JSON numbers")
    numbers <- unlist(present, use.names = FALSE)
    if (any(!is.finite(numbers))) stop("nonfinite numeric metadata")
    if (grepl("[Ii]nt", dtype)) {
      if (!is.null(tokens) && any(!grepl("^-?(0|[1-9][0-9]*)$", tokens[!missing]))) stop("integer metadata requires JSON integer tokens")
      bits <- as.numeric(sub(".*?([0-9]+)$", "\\1", dtype))
      unsigned <- startsWith(tolower(dtype), "u")
      low <- if (unsigned) 0 else -2^(bits - 1)
      high <- if (unsigned) 2^bits - 1 else 2^(bits - 1) - 1
      if (any(numbers != trunc(numbers) | numbers < low | numbers > high)) stop("integer metadata overflow or fraction")
      if (any(numbers == -.Machine$integer.max - 1)) stop("R cannot represent INT_MIN without losing missingness")
      if (grepl("^(int|uint)", dtype) && any(missing)) stop("nonnullable integer contains null")
    } else {
      if (grepl("32$", dtype) && length(numbers)) {
        rounded <- readBin(writeBin(as.double(numbers), raw(), size = 4L), "double", n = length(numbers), size = 4L)
        if (any(!is.finite(rounded) | rounded != numbers)) stop("lossy float32 metadata cast")
      }
      if (!is.null(tokens)) {
        integer_token <- grepl("^-?(0|[1-9][0-9]*)$", tokens[!missing])
        if (any(integer_token & sprintf("%.0f", numbers) != tokens[!missing] & tokens[!missing] != "-0"))
          stop("integer JSON value loses precision in floating dtype")
      }
    }
  }
  invisible(NULL)
}

.scx_metadata <- function(table, spec, key) {
  if (is.null(spec)) return(.scx_decode_frame(table, NULL))
  if (!.scx_keys(spec, c("index_name", "columns")) ||
      (!is.null(spec$index_name) && !.scx_text(spec$index_name)) || !is.list(spec$columns)) stop("invalid metadata table schema")
  names <- vapply(spec$columns, function(x) {
    if (!is.list(x) || !.scx_text(x$name)) stop("invalid metadata column descriptor")
    x$name
  }, character(1))
  if (!identical(colnames(table), c(key, names))) stop("metadata schema columns do not match TSV")
  for (column in spec$columns) {
    tokens <- table[[column$name]]
    values <- jsonlite::fromJSON(paste0("[", paste(tokens, collapse = ","), "]"), simplifyVector = FALSE)
    if (identical(column$dtype, "category")) {
      if (!.scx_keys(column, c("name", "dtype", "ordered", "categories_dtype", "categories")) ||
          !is.logical(column$ordered) || length(column$ordered) != 1L || is.na(column$ordered) || !is.list(column$categories)) stop("invalid categorical schema")
      .scx_check_values(column$categories, column$categories_dtype)
      if (!column$categories_dtype %in% c("object", "str", "string")) stop("numeric/boolean categorical levels are unsupported by R factors")
      if (any(vapply(column$categories, is.null, logical(1))) || anyDuplicated(unlist(column$categories))) stop("categories must be nonmissing and unique")
      missing <- vapply(values, is.null, logical(1))
      codes <- unlist(values[!missing], use.names = FALSE)
      if (any(!grepl("^(0|[1-9][0-9]*)$", tokens[!missing])) ||
          any(!vapply(values[!missing], is.numeric, logical(1))) ||
          any(codes < 0 | codes >= length(column$categories))) stop("categorical code out of range or wrong type")
    } else {
      if (!.scx_keys(column, c("name", "dtype"))) stop("invalid column schema")
      .scx_check_values(values, column$dtype, tokens)
    }
  }
  .scx_decode_frame(table, spec)
}

.scx_inspect <- function(path, retain = FALSE, retain_expression = NULL) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  manifest <- .scx_manifest(path)
  for (relative in manifest$path) if (!file.exists(.scx_safe_path(path, relative))) stop("missing file: ", file.path(path, relative))
  cells <- .scx_table(file.path(path, "cells.tsv.gz"))
  features <- .scx_table(file.path(path, "features.tsv.gz"))
  for (pair in list(list(cells, "cell_id"), list(features, "feature_id"))) {
    .scx_resource(file.path(path, if (pair[[2]] == "cell_id") "cells.tsv.gz" else "features.tsv.gz"), {
      if (!nrow(pair[[1]]) || names(pair[[1]])[1] != pair[[2]]) stop("invalid nonempty canonical ID column: ", pair[[2]])
      .scx_ids(pair[[1]][[1]])
    })
  }
  schema <- NULL
  if (!"metadata_schema" %in% manifest$role) warning("Legacy SCX metadata have no dtype schema; loading as text")
  if ("metadata_schema" %in% manifest$role) schema <- .scx_resource(file.path(path, "metadata.schema.json"), {
    spec <- jsonlite::fromJSON(file.path(path, "metadata.schema.json"), simplifyVector = FALSE)
    if (!.scx_keys(spec, c("schema_version", "cells", "features")) || !identical(spec$schema_version, "1.0") ||
        !.scx_keys(spec$cells, c("index_name", "columns")) || !.scx_keys(spec$features, c("index_name", "columns")))
      stop("invalid or unsupported metadata schema")
    spec
  })
  meta <- .scx_resource(file.path(path, "cells.tsv.gz"), .scx_metadata(cells, if (is.null(schema)) NULL else schema$cells, "cell_id"))
  feature_meta <- .scx_resource(file.path(path, "features.tsv.gz"), .scx_metadata(features, if (is.null(schema)) NULL else schema$features, "feature_id"))
  payload <- list(manifest = manifest, cells = cells[[1]], features = features[[1]], meta = meta, feature_meta = feature_meta,
                  matrices = list(), embeddings = list())
  shape <- c(nrow(cells), nrow(features))
  for (i in seq_len(nrow(manifest))) {
    row <- manifest[i, ]; file <- .scx_safe_path(path, row$path)
    if (row$role %in% c("counts", "normalized", "cell_graph")) {
      x <- .scx_matrix(file, if (row$role == "cell_graph") rep(shape[1], 2L) else shape)
      if (row$role == "counts" && any(.scx_values(x) < 0)) stop(file, ": counts must be nonnegative")
      if (retain || identical(row$role, retain_expression)) payload$matrices[[row$path]] <- x
      rm(x)
    } else if (row$role == "embedding") {
      table <- .scx_table(file)
      .scx_resource(file, {
        if (ncol(table) < 2L || names(table)[1] != "cell_id" || !identical(table[[1]], cells[[1]])) stop("invalid embedding cell axis")
        values <- matrix(as.numeric(as.matrix(table[-1])), nrow = shape[1])
        if (any(!is.finite(values))) stop("embedding must be numeric and finite")
        if (retain) payload$embeddings[[row$path]] <- values
      })
    }
  }
  if (retain || !is.null(retain_expression)) payload else list(cells = shape[1], features = shape[2], files = nrow(manifest))
}

#' Validate the R-readable SCX profile without Python.
#' @export
validate_scx <- function(path) .scx_inspect(path)
