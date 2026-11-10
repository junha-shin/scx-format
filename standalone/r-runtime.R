# Pure R backend for the generated standalone script, not the installed package.
# Deliberately refuses replacement: base R has no portable no-replace dir rename.
.scx_native <- function(name, ...) {
  if (name != "scx_check_gzip") stop("Native helper unavailable in standalone mode")
  path <- list(...)[[1]]
  .scx_check_gzip_r(path)
}

.scx_symlink <- function(path) {
  link <- suppressWarnings(Sys.readlink(path))
  length(link) == 1L && !is.na(link) && nzchar(link)
}

.scx_check_gzip_r <- function(path) {
  # Verify a single gzip member's CRC and length with existing digest + base R.
  # Stream to scratch disk rather than retaining the decompressed matrix in RAM.
  input <- file(path, "rb")
  scratch <- tempfile("scx-gzip-")
  output <- NULL; decoded <- NULL
  on.exit({
    if (!is.null(decoded)) close(decoded)
    if (!is.null(output)) close(output)
    close(input)
    unlink(scratch)
  }, add = TRUE)
  size <- file.info(path)$size
  if (is.na(size) || size < 18) stop("truncated gzip")
  header <- readBin(input, "raw", 3L)
  if (!identical(header, as.raw(c(31, 139, 8)))) stop("invalid gzip header")
  seek(input, size - 8, origin = "start")
  trailer <- readBin(input, "raw", 8L)
  expected_crc <- paste(sprintf("%02x", as.integer(rev(trailer[1:4]))), collapse = "")
  expected_size <- sum(as.integer(trailer[5:8]) * 256^(0:3))
  output <- file(scratch, "wb")
  decoded <- gzfile(path, "rb")
  total <- 0
  withCallingHandlers(repeat {
    chunk <- readBin(decoded, "raw", 1024L * 1024L)
    if (!length(chunk)) break
    total <- total + length(chunk)
    writeBin(chunk, output)
  }, warning = function(w) stop(conditionMessage(w), call. = FALSE))
  close(output)
  # Don't call isOpen() on a closed, invalid connection in the exit handler.
  output <- NULL
  if (total %% 2^32 != expected_size || digest::digest(scratch, file = TRUE, algo = "crc32") != expected_crc)
    stop("gzip CRC/length mismatch; corrupt, truncated or unsupported concatenated gzip")
  invisible(TRUE)
}

.scx_copy <- function(source, target) {
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  if (!isTRUE(file.copy(source, target, overwrite = FALSE, copy.mode = TRUE))) stop("copy failed: ", target)
  if (digest::digest(source, file = TRUE, algo = "sha256") != digest::digest(target, file = TRUE, algo = "sha256"))
    stop("copied file checksum mismatch: ", target)
}

.scx_publish <- function(stage, path, overwrite) {
  if (overwrite) stop("Standalone R refuses overwrite=TRUE; choose a new output directory. The installed package supports staged replacement.")
  path <- .scx_output_path(path)
  lock <- .scx_lock_path(path)
  con <- tryCatch(suppressWarnings(file(lock, "wx", encoding = "UTF-8")), error = function(e)
    stop("output is locked or lock cannot be created: ", lock, call. = FALSE))
  reserved <- FALSE; published <- FALSE
  on.exit({
    if (!reserved || published) unlink(lock)
    # If interrupted during copy, keep the incomplete directory and lock.
    # The upstream writer cleans its private staging files; no old output is touched.
  }, add = TRUE)
  tryCatch(writeLines(.scx_json(list(pid = Sys.getpid(), output = path, mode = "standalone-copy")), con),
           finally = close(con))
  payload <- file.path(stage, "payload")
  validate_scx(payload)
  .scx_destination(path, FALSE)
  # mkdir reserves the leaf atomically and cannot replace any existing entry.
  if (!dir.create(path, showWarnings = FALSE)) stop("output already exists or cannot be reserved: ", path)
  reserved <- TRUE
  tryCatch({
    manifest <- .scx_manifest(payload)
    for (relative in manifest$path) .scx_copy(file.path(payload, relative), file.path(path, relative))
    # No reader can accept the new directory until the manifest appears last.
    pending <- file.path(path, ".scx-manifest-pending")
    .scx_copy(file.path(payload, "manifest.tsv"), pending)
    if (!file.rename(pending, file.path(path, "manifest.tsv"))) stop("manifest publication failed")
    published <- TRUE
  }, error = function(e) stop("Incomplete output retained at ", path, "; lock retained at ", lock,
       ". Original inputs are unchanged. Choose a new path; see standalone recovery guidance. ",
       conditionMessage(e), call. = FALSE))
  invisible(path)
}

# Reject unsupported replacement before doing any data export work.
.scx_write_impl <- write_scx
write_scx <- function(object, path, counts = "counts", normalized = "data",
                      analysis = FALSE, overwrite = FALSE,
                      assay = SeuratObject::DefaultAssay(object)) {
  if (!identical(overwrite, FALSE)) stop("Standalone R requires overwrite=FALSE; choose a new output directory")
  .scx_write_impl(object, path, counts = counts, normalized = normalized,
                  analysis = analysis, overwrite = FALSE, assay = assay)
}
