# SCX quickstart with PBMC3k (R)

Use your existing Seurat environment. This lesson loads real PBMC3k raw counts,
saves and reads SCX, and previews a small TSV. No SCX package, Python or extra
notebook packages are needed to run the code. The rendered Markdown includes
actual output; this Rmd can also be knitted if you already use R Markdown.

## 1. Load SCX and PBMC3k

Reuse the repository's downloaded data when present. Otherwise download the
10x archive (about 8 MB) once. Run from the repository root, tutorials/, or a
folder of your choice.

```r
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
```

```text
Input (features x cells): 32738 2700 
               V1           V2
1 ENSG00000243485   MIR1302-10
2 ENSG00000237613      FAM138A
3 ENSG00000186092        OR4F5
4 ENSG00000238009 RP11-34P13.7
5 ENSG00000239945 RP11-34P13.8
6 ENSG00000237683   AL627309.1
```

## 2. Save SCX and read it back

Seurat holds FEATURE x CELL; SCX handles the transpose to CELL x FEATURE.
These are raw counts, so normalized=NULL. A full rerun creates a fresh folder
and keeps previous outputs.

```r
run <- tempfile("scx-pbmc3k-", tmpdir=".")
dir.create(run)
scx_dir <- file.path(run,"pbmc3k.scx")
write_scx(object, scx_dir, normalized=NULL)
restored <- read_scx(scx_dir)
cat("Restored Seurat (features x cells):", dim(restored), "\n")
cat("SCX orientation (cells x features):", rev(dim(restored)), "\n")
```

```text
Restored Seurat (features x cells): 32738 2700 
SCX orientation (cells x features): 2700 32738 
```

## 3. Preview an expression TSV

Select five cells and the five genes with the largest total counts. Only this
small preview is dense; the complete expression matrix remains sparse.

```r
features <- rownames(counts)[order(Matrix::rowSums(counts), decreasing=TRUE)[1:5]]
preview_path <- file.path(run,"preview.tsv")
to_tsv(scx_dir, preview_path, cells=head(cells,5), features=features)
preview <- read.delim(preview_path, row.names=1, check.names=FALSE)
cat("TSV preview (cells x features):", dim(preview), "\n")
head(preview)
```

```text
TSV preview (cells x features): 5 5 
                 ENSG00000251562 ENSG00000205542 ENSG00000166710
AAACATACAACCAC-1              49              47              76
AAACATTGAGCTAC-1             142              62              75
AAACATTGATCAGC-1             171             117              69
AAACCGTGCTTCCG-1              11             114              41
AAACCGTGTATGCG-1              22              21              35
                 ENSG00000147403 ENSG00000167526
AAACATACAACCAC-1              34              29
AAACATTGAGCTAC-1              92              45
AAACATTGATCAGC-1              49              16
AAACCGTGCTTCCG-1              22              15
AAACCGTGTATGCG-1               8               4
```

Expected: Seurat 32738 x 2700, SCX 2700 x 32738, preview 5 x 5.
For your own data, replace section 1 with your existing Seurat object. Default
export uses counts/data in the active assay; normalized=NULL is for counts only.
Keep your original analysis object.

[Python notebook](quickstart.ipynb) | [Runnable R script](quickstart.R) |
[Detailed options](detailed.R) | [Common errors](../docs/first-use.md)

Data: 10x PBMC3k ([provenance](../examples/seurat/README.md)); upstream dataset
terms apply. [R Markdown source](quickstart.Rmd).
