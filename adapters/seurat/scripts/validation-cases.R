library(scxSeurat)
args <- commandArgs(TRUE)
mode <- args[1]
if (mode == "cases") {
  root <- args[2]
  expected <- jsonlite::fromJSON(file.path(root, "expected.json"), simplifyVector = FALSE)
  results <- lapply(names(expected), function(name) {
    error <- NULL
    valid <- tryCatch({ suppressWarnings(validate_scx(file.path(root, name))); TRUE }, error = function(e) {
      error <<- conditionMessage(e); FALSE
    })
    stopifnot(identical(valid, expected[[name]]$R))
    list(name = name, accepted = valid, error = error)
  })
  writeLines(jsonlite::toJSON(list(passed = TRUE, cases = results), auto_unbox = TRUE, null = "null"), file.path(root, "R-results.json"))
} else if (mode == "locked") {
  object <- read_scx_seurat(args[2])
  error <- tryCatch({suppressWarnings(write_scx_seurat(object, args[3])); "not blocked"}, error = conditionMessage)
  stopifnot(grepl("locked", error, fixed = TRUE))
} else if (mode == "hold") {
  path <- getFromNamespace(".scx_output_path", "scxSeurat")(args[2])
  lock <- getFromNamespace(".scx_lock_path", "scxSeurat")(path)
  con <- file(lock, "wx"); writeLines("R validation test owns this lock", con); close(con)
  file.create(args[3])
  start <- proc.time()[[3]]
  while (!file.exists(args[4]) && proc.time()[[3]] - start < 30) Sys.sleep(0.05)
  unlink(lock)
} else stop("unknown validation mode")
