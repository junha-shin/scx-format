# Diagnostic only: instrumentation and Rprof overhead exclude this run from timings.
# Rscript profile.R standalone/scx.R input.scx output.tsv profile-prefix
args <- commandArgs(TRUE)
stopifnot(length(args)==4L, !file.exists(args[3]))
source(args[1])
env <- environment(to_tsv)
timings <- new.env(parent=emptyenv())
timings$validation <- 0; timings$row_extract_and_dense <- 0
timings$number_format <- 0; timings$write_lines <- 0
inspect <- env$.scx_inspect
env$.scx_inspect <- function(...) {
  start <- proc.time()[["elapsed"]]
  on.exit(timings$validation <- timings$validation + proc.time()[["elapsed"]]-start)
  inspect(...)
}
env$as.matrix <- function(x,...) {
  start <- proc.time()[["elapsed"]]
  on.exit(timings$row_extract_and_dense <- timings$row_extract_and_dense + proc.time()[["elapsed"]]-start)
  base::as.matrix(x,...)
}
env$sprintf <- function(fmt,...) {
  args <- list(...)
  start <- proc.time()[["elapsed"]]
  value <- do.call(base::sprintf,c(list(fmt),args))
  if (identical(fmt,"%.17g")) timings$number_format <- timings$number_format + proc.time()[["elapsed"]]-start
  value
}
env$writeLines <- function(text,con=stdout(),sep="\n",useBytes=FALSE) {
  force(text)
  start <- proc.time()[["elapsed"]]
  base::writeLines(text,con,sep=sep,useBytes=useBytes)
  timings$write_lines <- timings$write_lines + proc.time()[["elapsed"]]-start
}
profile <- paste0(args[4],".Rprof")
Rprof(profile,interval=0.01)
start <- proc.time()[["elapsed"]]
tryCatch(to_tsv(args[2],args[3]),finally=Rprof(NULL))
elapsed <- proc.time()[["elapsed"]]-start
summary <- summaryRprof(profile)
saveRDS(summary,paste0(args[4],"-summary.rds"))
write.csv(head(summary$by.total,30L),paste0(args[4],"-total.csv"))
write.csv(head(summary$by.self,30L),paste0(args[4],"-self.csv"))
result <- list(diagnostic_only=TRUE, API_seconds=elapsed,
  instrumented_seconds=as.list(timings),
  note="Diagnostic wrappers/Rprof add overhead. Inclusive samples and wrapper timings are not additive. Row extraction is forced inside as.matrix; number-format arguments are forced before its timer.")
writeLines(jsonlite::toJSON(result,auto_unbox=TRUE,digits=17,pretty=TRUE),paste0(args[4],".json"))
cat(jsonlite::toJSON(result,auto_unbox=TRUE,digits=17),"\n")
