# Resolve nonexistent parent segments before creating anything or checking containment.
.scx_tsv_target <- function(output) {
  if (!is.character(output) || length(output) != 1L || is.na(output) || !nzchar(output)) stop("output must be one file path")
  output <- path.expand(output)
  parent <- dirname(output); suffix <- basename(output)
  while (!dir.exists(parent)) {
    next_parent <- dirname(parent)
    if (identical(next_parent, parent)) stop("cannot resolve output parent")
    suffix <- file.path(basename(parent), suffix)
    parent <- next_parent
  }
  target <- normalizePath(parent, winslash="/", mustWork=TRUE)
  if (.Platform$OS.type == "windows") suffix <- chartr("\\", "/", suffix)
  for (part in strsplit(suffix, "/", fixed=TRUE)[[1]]) {
    if (part == "..") target <- dirname(target)
    else if (part != "." && nzchar(part)) target <- file.path(target, part)
  }
  target
}

.scx_tsv_select <- function(ids, requested, label) {
  if (is.null(requested)) return(seq_along(ids))
  if (!is.character(requested) || !length(requested) || anyNA(requested) || anyDuplicated(requested))
    stop(label, " must be nonempty unique character IDs")
  selected <- match(requested, ids)
  if (anyNA(selected)) stop("unknown ", label, " IDs: ", paste(head(requested[is.na(selected)], 5L), collapse=", "))
  selected
}

.scx_tsv_field <- function(x) ifelse(grepl('["\t\r\n]', x), paste0('"', gsub('"', '""', x, fixed=TRUE), '"'), x)

#' Export an SCX expression matrix as CELL x FEATURE TSV or TSV.gz.
#' @export
to_tsv <- function(path, output, matrix="counts", cells=NULL, features=NULL) {
  if (!is.character(matrix) || length(matrix) != 1L || is.na(matrix) || !matrix %in% c("counts","normalized"))
    stop("matrix must be counts or normalized")
  root <- normalizePath(path, winslash="/", mustWork=TRUE)
  target <- .scx_tsv_target(output)
  root_key <- root; target_key <- target
  if (.Platform$OS.type == "windows") { root_key <- tolower(root); target_key <- tolower(target) }
  if (target_key == root_key || startsWith(target_key, paste0(root_key,"/"))) stop("TSV output must be outside the input SCX directory")
  if (.scx_exists(target)) stop("output already exists; choose a new path")
  data <- .scx_inspect(root, retain_expression=matrix)
  entry <- data$manifest[data$manifest$role == matrix, , drop=FALSE]
  if (!nrow(entry)) stop("SCX has no normalized matrix")
  values <- data$matrices[[entry$path]]
  data$matrices <- NULL
  rows <- .scx_tsv_select(data$cells, cells, "cells")
  columns <- .scx_tsv_select(data$features, features, "features")
  target <- .scx_output_path(target)
  stage <- tempfile(".scx-tsv-", tmpdir=dirname(target))
  con <- NULL; input <- NULL; destination <- NULL
  on.exit({
    if (!is.null(con)) close(con)
    if (!is.null(input)) close(input)
    if (!is.null(destination)) close(destination)
    unlink(stage)
  }, add=TRUE)
  con <- if (endsWith(target,".gz")) gzfile(stage,"wt",encoding="UTF-8") else file(stage,"wt",encoding="UTF-8")
  withCallingHandlers({
    # CSR gives direct access to each row's stored entries. Keep the source
    # sparse; the only full-width output buffer is a vector of text fields.
    values <- methods::as(values,"RsparseMatrix")
    column_map <- integer(ncol(values))
    column_map[columns] <- seq_along(columns)
    writeLines(paste(.scx_tsv_field(c("cell_id",data$features[columns])), collapse="\t"), con, useBytes=TRUE)
    for (row in rows) {
      fields <- rep.int("0",length(columns))
      from <- values@p[row]; to <- values@p[row+1]
      if (from < to) {
        entries <- seq.int(from+1L,to)
        positions <- column_map[values@j[entries]+1L]
        keep <- positions > 0L
        # Include explicit stored zeros as well: sprintf retains signed zero.
        fields[positions[keep]] <- sprintf("%.17g",values@x[entries[keep]])
      }
      writeLines(paste(c(.scx_tsv_field(data$cells[row]), fields),collapse="\t"), con, useBytes=TRUE)
    }
    close(con); con <- NULL
  }, warning=function(w) stop(conditionMessage(w),call.=FALSE))
  lock <- .scx_lock_path(target)
  lock_con <- tryCatch(suppressWarnings(file(lock,"wx")), error=function(e) stop("output is locked: ",lock))
  reserved <- FALSE; published <- FALSE
  on.exit(if (!reserved || published) unlink(lock), add=TRUE)
  tryCatch(writeLines(.scx_json(list(pid=Sys.getpid(), output=target, mode="tsv-copy")),lock_con), finally=close(lock_con))
  tryCatch({
    # Exclusive binary creation refuses existing files, links and racing writers.
    destination <- suppressWarnings(file(target,"wxb"))
    reserved <- TRUE
    input <- file(stage,"rb")
    repeat {
      bytes <- readBin(input,"raw",1024L*1024L)
      if (!length(bytes)) break
      writeBin(bytes,destination)
    }
    close(destination); destination <- NULL
    close(input); input <- NULL
    if (digest::digest(stage,file=TRUE,algo="sha256") != digest::digest(target,file=TRUE,algo="sha256")) stop("TSV copy checksum mismatch")
    published <- TRUE
  }, error=function(e) {
    if (reserved) stop("Incomplete TSV retained at ",target,"; lock retained at ",lock,". ",conditionMessage(e),call.=FALSE)
    stop("output already exists or cannot be created: ",target,". ",conditionMessage(e),call.=FALSE)
  })
  invisible(target)
}
