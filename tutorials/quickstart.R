# Save scx.R beside this script; use your existing Seurat environment and object.
# https://github.com/junha-shin/scx-format/releases/download/v0.1.0/scx.R
# object is your existing Seurat object (counts/data in the active assay).
source("scx.R")
write_scx(object, "sample.scx")
restored <- read_scx("sample.scx")
to_tsv("sample.scx", "counts.tsv")

# Counts-only object: write_scx(object, "sample.scx", normalized=NULL)
# Use fresh output names when running again. More help: detailed.R
