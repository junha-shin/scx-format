# User entry point: source("https://raw.githubusercontent.com/junha-shin/scx-format/v0.1.0-alpha.6/install.R")
# Define-only mode is for offline tests/builds, not needed by users.
.scx_binary_asset <- function(system = Sys.info()[["sysname"]], arch = R.version$arch,
                               r_minor = paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = "."),
                               os_release = if (file.exists("/etc/os-release")) readLines("/etc/os-release", warn = FALSE) else character()) {
  if (r_minor != "4.5") return(NULL)
  arch <- if (arch %in% c("aarch64", "arm64")) "arm64" else arch
  if (system == "Windows" && arch == "x86_64") return("scxSeurat_0.1.5-windows-x86_64.zip")
  if (system == "Darwin" && arch %in% c("arm64", "x86_64"))
    return(paste0("scxSeurat_0.1.5-macos-", arch, ".tgz"))
  if (system == "Linux" && arch == "x86_64" && any(grepl('^ID="?ubuntu"?$', os_release))) {
    number <- sub('^VERSION_ID="?([^" ]+)"?$', "\\1", grep("^VERSION_ID=", os_release, value = TRUE))
    if (length(number) == 1L && !is.na(suppressWarnings(as.numeric(number))) && as.numeric(number) >= 22.04)
      return("scxSeurat_0.1.5-ubuntu22-x86_64.tar.gz")
  }
  NULL
}

.scx_verify_download <- function(path, sums) {
  entries <- strsplit(sums, "  ", fixed = TRUE)
  match <- Filter(function(x) length(x) == 2L && identical(x[2], basename(path)), entries)
  if (length(match) != 1L || !grepl("^[a-fA-F0-9]{64}$", match[[1]][1])) stop("No unique SHA-256 entry for ", basename(path))
  actual <- digest::digest(path, algo = "sha256", file = TRUE)
  if (!identical(tolower(match[[1]][1]), actual)) stop("Checksum mismatch; package was not installed")
  invisible(TRUE)
}

install_scx <- function(lib = NULL, source = FALSE,
                        repos = "https://cloud.r-project.org", assets_dir = NULL) {
  if ("scxSeurat" %in% loadedNamespaces()) stop("Restart R before installing or updating scxSeurat")
  if (!is.logical(source) || length(source) != 1L || is.na(source)) stop("source must be TRUE or FALSE")
  asset <- if (source) "scxSeurat_0.1.5.tar.gz" else .scx_binary_asset()
  if (is.null(asset)) stop("No binary for this R/OS combination. With R build tools and zlib installed, call install_scx(source=TRUE). See docs/install.md")
  if (is.null(lib)) {
    user <- strsplit(Sys.getenv("R_LIBS_USER"), .Platform$path.sep, fixed = TRUE)[[1]]
    lib <- if (length(user) && nzchar(user[1])) path.expand(user[1]) else tools::R_user_dir("scxSeurat", "data")
  }
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  lib <- normalizePath(lib, winslash = "/", mustWork = TRUE)
  if (file.access(lib, 2L)) stop("R library is not writable: ", lib)
  .libPaths(c(lib, .libPaths()))
  required <- c("SeuratObject", "Matrix", "jsonlite", "digest")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if ("SeuratObject" %in% setdiff(required, missing) && packageVersion("SeuratObject") < "5.0.0")
    stop("SCX requires SeuratObject >=5; existing packages were not upgraded")
  if (length(missing)) utils::install.packages(missing, lib = lib, repos = repos)
  if (!all(vapply(required, requireNamespace, logical(1), quietly = TRUE)) || packageVersion("SeuratObject") < "5.0.0")
    stop("Required R dependencies are unavailable; install a compatible SeuratObject >=5 environment")
  temp <- tempfile("scx-install-"); dir.create(temp)
  on.exit(unlink(temp, recursive = TRUE), add = TRUE)
  obtain <- function(name) {
    destination <- file.path(temp, name)
    if (is.null(assets_dir)) {
      utils::download.file(paste0("https://github.com/junha-shin/scx-format/releases/download/v0.1.0-alpha.6/", name),
                           destination, mode = "wb", quiet = TRUE)
    } else {
      if (!file.copy(file.path(assets_dir, name), destination)) stop("Missing local asset: ", name)
    }
    destination
  }
  sums <- readLines(obtain("SHA256SUMS"), warn = FALSE)
  available <- vapply(strsplit(sums, "  ", fixed = TRUE), function(x) if (length(x) == 2L) x[2] else "", character(1))
  if (!asset %in% available) stop("This release has no matching binary. With build tools, call install_scx(source=TRUE)")
  path <- obtain(asset)
  .scx_verify_download(path, sums)
  # R's binary installer derives the package directory from the archive name.
  # Keep platform labels in release assets, but use the standard name locally.
  extension <- if (endsWith(asset, ".zip")) ".zip" else if (endsWith(asset, ".tgz")) ".tgz" else ".tar.gz"
  install_path <- file.path(temp, paste0("scxSeurat_0.1.5", extension))
  if (!identical(path, install_path) && !file.rename(path, install_path)) stop("Cannot prepare verified package for installation")
  type <- if (source || Sys.info()[["sysname"]] == "Linux") "source" else if (.Platform$OS.type == "windows") "win.binary" else "mac.binary"
  utils::install.packages(install_path, lib = lib, repos = NULL, type = type)
  description <- file.path(lib, "scxSeurat", "DESCRIPTION")
  if (!file.exists(description) || read.dcf(description)[1, "Version"] != "0.1.5") stop("scxSeurat installation did not complete")
  message("Installed scxSeurat 0.1.5. Next: library(scxSeurat)")
  invisible(list(library = lib, asset = asset, dependencies_installed = missing))
}

if (!identical(Sys.getenv("SCX_INSTALL_DEFINE_ONLY"), "1")) install_scx()
