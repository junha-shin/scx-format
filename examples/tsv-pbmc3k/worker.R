# Each invocation measures exactly one no-install TSV export in a fresh process.
args <- commandArgs(TRUE)
stopifnot(length(args)==5L)
source(args[1])
stopifnot(!"scxSeurat" %in% loadedNamespaces())
config <- jsonlite::fromJSON(args[4],simplifyVector=TRUE)
start <- proc.time()[["elapsed"]]
to_tsv(args[2],args[3],matrix=config$matrix,cells=config$cells,features=config$features)
elapsed <- proc.time()[["elapsed"]]-start
versions <- vapply(c("SeuratObject","Matrix","jsonlite","digest"),function(x) as.character(packageVersion(x)),character(1))
writeLines(jsonlite::toJSON(list(operation_seconds=elapsed,R=R.version.string,versions=as.list(versions),
  gzip_compression_level=eval(formals(gzfile)$compression)),
  auto_unbox=TRUE,digits=17),args[5])
