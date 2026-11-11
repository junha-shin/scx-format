# Build installable artifacts and test the same offline path used by install.R.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
output <- args[1]; dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output, winslash = "/", mustWork = TRUE)
r <- file.path(R.home("bin"), "R")
run <- function(arguments) {
  status <- system2(r, arguments)
  if (status != 0L) stop("R package build/install failed")
}
run(c("CMD", "build", "adapters/seurat"))
source <- "scxSeurat_0.1.6.tar.gz"
before <- list.files(pattern = "^scxSeurat_.*\\.(zip|tgz|tar\\.gz)$")
run(c("CMD", "INSTALL", "--build", source))
binary <- setdiff(list.files(pattern = "^scxSeurat_.*\\.(zip|tgz|tar\\.gz)$"), before)
stopifnot(length(binary) == 1L)
Sys.setenv(SCX_INSTALL_DEFINE_ONLY = "1")
sys.source("install.R", envir = environment())
asset <- .scx_binary_asset()
stopifnot(!is.null(asset))
stopifnot(file.copy(source, file.path(output, source), overwrite = TRUE),
          file.copy(binary, file.path(output, asset), overwrite = TRUE))
files <- c(source, asset)
sums <- vapply(files, function(name) paste0(digest::digest(file.path(output, name), algo = "sha256", file = TRUE), "  ", name), character(1))
writeLines(sums, file.path(output, "SHA256SUMS"))
lib <- tempfile("scx-library-"); dir.create(lib)
result <- install_scx(lib = lib, assets_dir = output)
stopifnot(!length(result$dependencies_installed))
Sys.setenv(R_LIBS_USER = paste(c(lib, .libPaths()), collapse = .Platform$path.sep))
status <- system2(file.path(R.home("bin"), "Rscript"),
  c("adapters/seurat/scripts/r-only.R", shQuote(file.path(output, "installed-check"))))
stopifnot(status == 0L)
writeLines(jsonlite::toJSON(list(passed = TRUE, asset = asset, R = R.version.string,
  platform = R.version$platform, os = as.list(Sys.info()[c("sysname", "release", "version", "machine")]), dependencies_installed = as.list(result$dependencies_installed)), auto_unbox = TRUE),
  file.path(output, "installation.json"))
