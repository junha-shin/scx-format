.scx_native <- function(name, ...) {
  if (!is.loaded(name, PACKAGE = "scxSeurat")) {
    if (!requireNamespace("scxSeurat", quietly = TRUE) || !is.loaded(name, PACKAGE = "scxSeurat"))
      stop("Install scxSeurat >= 0.1.4 with its native helpers; Python is not required", call. = FALSE)
  }
  .Call(name, ..., PACKAGE = "scxSeurat")
}

.scx_exists <- function(path) file.exists(path) || dir.exists(path) || .scx_symlink(path)
.scx_symlink <- function(path) { x <- Sys.readlink(path); !is.na(x) && nzchar(x) }
.scx_rename_new <- function(source, target) .scx_native("scx_rename_new", source, target)

.scx_output_path <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  file.path(normalizePath(dirname(path), winslash = "/", mustWork = TRUE), basename(path))
}

.scx_lock_path <- function(path) {
  # Same lock identity as Python 0.1.0a5 (including Windows case/separators).
  key <- if (.Platform$OS.type == "windows") tolower(chartr("/", "\\", path)) else path
  file.path(dirname(path), paste0(".scx-lock-", substr(digest::digest(charToRaw(enc2utf8(key)),
             algo = "sha256", serialize = FALSE), 1L, 20L)))
}

.scx_destination <- function(path, overwrite) {
  if (.scx_symlink(path)) stop("output symlinks are refused: ", path)
  if (.scx_exists(path)) {
    if (!overwrite) stop("output already exists: ", path, "; choose a new path or use overwrite = TRUE")
    if (!dir.exists(path)) stop("output must be an SCX directory: ", path)
    .scx_manifest(path)  # Refuse unrelated directories, even with opt-in.
  }
}

.scx_publish <- function(stage, path, overwrite) {
  path <- .scx_output_path(path)
  lock <- .scx_lock_path(path)
  con <- tryCatch(suppressWarnings(file(lock, "wx", encoding = "UTF-8")), error = function(e)
    stop("output is locked or lock cannot be created: ", path, "; inspect ", lock,
         " and docs/safe-io.md before recovery", call. = FALSE))
  backup <- file.path(stage, "previous")
  payload <- file.path(stage, "payload")
  # stage is an outer private directory; completed dataset is payload.
  published <- FALSE; preserve <- FALSE
  on.exit({
    if (dir.exists(backup) && !published) preserve <- TRUE
    if (!preserve) {
      unlink(stage, recursive = TRUE)
      unlink(lock)
    }
  }, add = TRUE)
  tryCatch({
    writeLines(.scx_json(list(pid = Sys.getpid(), output = path, stage = stage, backup = backup)), con)
    flush(con)
  }, finally = close(con))
  .scx_destination(path, overwrite)
  validate_scx(payload)
  .scx_destination(path, overwrite)
  tryCatch({
    if (.scx_exists(path)) .scx_rename_new(path, backup)
    .scx_rename_new(payload, path)
    published <- TRUE
  }, error = function(e) {
    if (dir.exists(backup)) {
      tryCatch(.scx_rename_new(backup, path), error = function(recovery) {
        preserve <<- TRUE
        stop("could not restore previous output; backup retained at ", backup,
             "; lock retained at ", lock, ". See docs/safe-io.md. ", conditionMessage(recovery), call. = FALSE)
      })
    }
    stop(e)
  }, interrupt = function(e) {
    # Suspend further interrupts during rollback; a hard kill retains artifacts.
    suspendInterrupts(if (dir.exists(backup)) .scx_rename_new(backup, path))
    stop(e)
  })
  invisible(path)
}
