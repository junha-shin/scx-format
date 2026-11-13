# SCX quickstart — use your existing Seurat environment.
# Save scx.R beside this script, then run from that folder:
# https://github.com/junha-shin/scx-format/releases/download/v0.1.0-alpha.7/scx.R
source("scx.R")

# Try small data, or use your own Seurat object instead.
counts <- Matrix::Matrix(matrix(c(1,0,3,0,2,1), nrow=3), sparse=TRUE)
dimnames(counts) <- list(c("gene1","gene2","gene3"), c("cell1","cell2"))
object <- SeuratObject::CreateSeuratObject(counts=counts)
SeuratObject::LayerData(object, layer="data") <- log1p(counts)

# Save and read SCX.
write_scx(object, "example.scx")
restored <- read_scx("example.scx")
restored

# Export a plain CELL x FEATURE expression table.
to_tsv("example.scx", "counts.tsv")

# Existing outputs are protected: choose new names when running again.
# More options and verified downloads: detailed.R
