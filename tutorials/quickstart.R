# PBMC3k quickstart. Notebook-style walkthrough: quickstart.md

if (!file.exists("scx.R"))
  download.file("https://github.com/junha-shin/scx-format/releases/download/v0.1.0/scx.R", "scx.R", mode="wb", quiet=TRUE)
source("scx.R")
folders <- c("examples/seurat/data/hg19", "../examples/seurat/data/hg19",
             "pbmc3k/filtered_gene_bc_matrices/hg19")
folder <- folders[file.exists(file.path(folders,"matrix.mtx"))][1]
if (is.na(folder)) {
  archives <- c("examples/seurat/data/pbmc3k.tar.gz", "../examples/seurat/data/pbmc3k.tar.gz", "pbmc3k.tar.gz")
  archive <- archives[file.exists(archives)][1]
  if (is.na(archive)) {
    archive <- "pbmc3k.tar.gz"
    download.file("https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz", archive, mode="wb", quiet=TRUE)
  }
  dir.create("pbmc3k", showWarnings=FALSE)
  untar(archive, exdir="pbmc3k")
  folder <- folders[3]
}
counts <- Matrix::readMM(file.path(folder,"matrix.mtx"))
genes <- read.delim(file.path(folder,"genes.tsv"), header=FALSE)
cells <- readLines(file.path(folder,"barcodes.tsv"))
dimnames(counts) <- list(genes[[1]], cells)
object <- SeuratObject::CreateSeuratObject(counts=counts)
cat("Input (features x cells):", dim(object), "\n")
head(genes)


run <- tempfile("scx-pbmc3k-", tmpdir=".")
dir.create(run)
scx_dir <- file.path(run,"pbmc3k.scx")
write_scx(object, scx_dir, normalized=NULL)
restored <- read_scx(scx_dir)
cat("Restored Seurat (features x cells):", dim(restored), "\n")
cat("SCX orientation (cells x features):", rev(dim(restored)), "\n")


features <- rownames(counts)[order(Matrix::rowSums(counts), decreasing=TRUE)[1:5]]
preview_path <- file.path(run,"preview.tsv")
to_tsv(scx_dir, preview_path, cells=head(cells,5), features=features)
preview <- read.delim(preview_path, row.names=1, check.names=FALSE)
cat("TSV preview (cells x features):", dim(preview), "\n")
head(preview)
