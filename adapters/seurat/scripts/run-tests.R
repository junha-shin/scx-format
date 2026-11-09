# Run from repository root: Rscript adapters/seurat/scripts/run-tests.R
source("adapters/seurat/R/scx.R")
testthat::test_dir("adapters/seurat/tests/testthat", stop_on_failure = TRUE)
