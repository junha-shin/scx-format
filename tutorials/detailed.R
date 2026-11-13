# SCX detailed example - R / Seurat
# Open in RStudio and run each section, or source the whole file.
# It uses the same small data as detailed.ipynb. No SCX or Python install is needed.
# Full reruns make fresh output folders. Earlier results are kept.
# A repeated write step refuses existing output.

# 1. Check your existing Seurat environment ----
tutorial_needed <- c("SeuratObject", "Matrix", "jsonlite", "digest")
tutorial_missing <- tutorial_needed[!vapply(tutorial_needed, requireNamespace, logical(1), quietly=TRUE)]
if (length(tutorial_missing)) stop("Choose your Seurat environment; missing: ", paste(tutorial_missing,collapse=", "))
if (packageVersion("SeuratObject") < "5.0.0") stop("This lesson needs SeuratObject >=5.")
cat("Ready: existing SeuratObject / Matrix / jsonlite / digest. No packages installed.\n")

# 2. Download and check alpha.7 files ----
# First use needs internet. For offline use, copy alpha.7 scx.R and SHA256SUMS
# into scx-tutorial/files/ under your working folder.
# Cached files are checked before source(). A bad hash stops the example.
# Inspect the cache and get the exact release file again. It is not replaced for you.
tutorial_release <- "v0.1.0-alpha.7"
tutorial_base <- paste0("https://github.com/junha-shin/scx-format/releases/download/",tutorial_release,"/")
tutorial_hashes <- c(
  "scx.R"="b6da0a10d7a912a4ed636fd51da9fdf4a4b84ca99d3b1d6e86826122773babce",
  "SHA256SUMS"="f7c277c380d04ed1f2766b0c9726d6b32c69811cf51a7b426808efb9cadfe914")
tutorial_home <- file.path(getwd(),"scx-tutorial")
tutorial_files <- file.path(tutorial_home,"files")
dir.create(tutorial_files,recursive=TRUE,showWarnings=FALSE)
tutorial_verified_file <- function(name) {
  path <- file.path(tutorial_files,name)
  if (file.exists(path)) {
    if (digest::digest(path,file=TRUE,algo="sha256") != tutorial_hashes[[name]])
      stop("Checksum mismatch: ",name,". Cached file was not executed or overwritten; inspect and download the pinned release again.")
  } else {
    temporary <- tempfile(".download-",tmpdir=tutorial_files)
    on.exit(unlink(temporary),add=TRUE)
    for (attempt in 1:2) {
      result <- tryCatch({
        utils::download.file(paste0(tutorial_base,name),temporary,method="libcurl",mode="wb",quiet=TRUE)
        TRUE
      },error=function(e) FALSE)
      if (result) break
      if (attempt==2L) stop("Download failed. Copy the alpha.7 files to scx-tutorial/files/ for offline use.")
      Sys.sleep(0.5)
    }
    if (digest::digest(temporary,file=TRUE,algo="sha256") != tutorial_hashes[[name]]) stop("Downloaded checksum mismatch: ",name)
    if (!file.copy(temporary,path,overwrite=FALSE)) stop("A cache file appeared or cannot be created. Inspect it, then rerun this section.")
  }
  normalizePath(path,winslash="/",mustWork=TRUE)
}
tutorial_checksum <- tutorial_verified_file("SHA256SUMS")
tutorial_script <- tutorial_verified_file("scx.R")
stopifnot(paste0(tutorial_hashes[["scx.R"]],"  scx.R") %in% readLines(tutorial_checksum,warn=FALSE))
tutorial_output <- file.path(tutorial_home,"output")
dir.create(tutorial_output,showWarnings=FALSE)
run_dir <- tempfile("run-",tmpdir=tutorial_output)
if (!dir.create(run_dir)) stop("Could not create a fresh run directory.")
cat("Verified alpha.7 bytes. Created a fresh run-* folder under scx-tutorial/output/.\n")

# 3. Load the checked file ----
source(tutorial_script,local=environment())
stopifnot(SCX_VERSION==tutorial_release)
cat("Loaded ",SCX_VERSION," from the verified file.\n",sep="")

# 4. Make small data ----
# Seurat uses features x cells. SCX and TSV use cells x features.
# Only this tiny display is dense. Keep real datasets sparse.
# Normalized values are log1p(counts) for this example, not a biological analysis.
# Seurat also adds QC metadata.
demo_counts <- Matrix::Matrix(matrix(c(1,0,3,0, 0,2,0,4, 5,1,0,0),nrow=4,ncol=3),sparse=TRUE)
dimnames(demo_counts) <- list(c("geneA","geneB","geneC","geneD"),c("cell1","cell2","cell3"))
demo <- SeuratObject::CreateSeuratObject(SeuratObject::CreateAssay5Object(counts=demo_counts),project="SCXdemo")
SeuratObject::LayerData(demo,layer="data") <- log1p(demo_counts)
demo$sample_number <- c(1L,NA_integer_,3L)
demo$flag <- c(TRUE,NA,FALSE)
demo$group <- ordered(c("B","A","B"),levels=c("B","A","unused"))
demo[["RNA"]] <- SeuratObject::AddMetaData(demo[["RNA"]],
  metadata=data.frame(symbol=c("A","B","C","D"),row.names=rownames(demo_counts)))
cat("Expected: 3 cells x 4 features, displayed CELL x FEATURE\n")
print(as.matrix(Matrix::t(demo_counts)))
print(demo[[]][,c("sample_number","flag","group")])

# 5. Save, read, and check ----
# Defaults use counts/data in the active assay. SCX checks data automatically.
# This example has no embeddings or graphs. Use analysis=TRUE to include them.
scx_dir <- file.path(run_dir,"demo.scx")
write_scx(demo,scx_dir)
restored <- read_scx(scx_dir)
stopifnot(all(dim(restored)==c(4,3)))
stopifnot(max(abs(SeuratObject::LayerData(restored,layer="counts")-demo_counts))==0)
normalized_error <- max(abs(SeuratObject::LayerData(restored,layer="data")-log1p(demo_counts)))
stopifnot(normalized_error<=1e-12,identical(demo[[]],restored[[]]),
  identical(demo[["RNA"]][[]],restored[["RNA"]][[]]))
cat("PASS: counts, normalized expression, both axes and typed metadata preserved.\n")
print(restored[[]][,c("sample_number","flag","group")])
# Read metadata with object[[]]. Seurat's $ can drop the ordered-factor class.

# 6. View counts TSV and compressed normalized TSV ----
counts_tsv <- file.path(run_dir,"counts.tsv")
normalized_tsv <- file.path(run_dir,"normalized.tsv.gz")
to_tsv(scx_dir,counts_tsv)
to_tsv(scx_dir,normalized_tsv,matrix="normalized")
counts_table <- read.delim(counts_tsv,check.names=FALSE,stringsAsFactors=FALSE)
connection <- gzfile(normalized_tsv,"rt")
normalized_table <- read.delim(connection,check.names=FALSE,stringsAsFactors=FALSE)
close(connection)
stopifnot(identical(counts_table$cell_id,colnames(demo_counts)),
  identical(names(counts_table)[-1],rownames(demo_counts)),
  max(abs(as.matrix(counts_table[-1])-as.matrix(Matrix::t(demo_counts))))==0,
  max(abs(as.matrix(normalized_table[-1])-as.matrix(Matrix::t(log1p(demo_counts)))))<=1e-12)
print(counts_table)
print(normalized_table)

# 7. Select and order IDs ----
# Expect cell3,cell1 by geneB,geneA: [1,5] and [0,1].
selected_cells <- c("cell3","cell1")
selected_features <- c("geneB","geneA")
subset_tsv <- file.path(run_dir,"subset.tsv.gz")
to_tsv(scx_dir,subset_tsv,cells=selected_cells,features=selected_features)
connection <- gzfile(subset_tsv,"rt")
subset_table <- read.delim(connection,check.names=FALSE,stringsAsFactors=FALSE)
close(connection)
stopifnot(identical(subset_table$cell_id,selected_cells),
  identical(names(subset_table)[-1],selected_features),
  max(abs(as.matrix(subset_table[-1])-matrix(c(1,0,5,1),nrow=2)))==0)
print(subset_table)

# 8. Check existing output protection ----
previous_hash <- digest::digest(counts_tsv,file=TRUE,algo="sha256")
refusal <- tryCatch({to_tsv(scx_dir,counts_tsv);NULL},error=conditionMessage)
stopifnot(!is.null(refusal),grepl("already exists",refusal,fixed=TRUE),
  digest::digest(counts_tsv,file=TRUE,algo="sha256")==previous_hash)
cat("Expected refusal: counts.tsv already exists; its bytes are preserved.\n")
# The R file refuses existing SCX folders and overwrite=TRUE.
# Run section 2 or the whole example for a new folder. Old results are kept.

# 9. Try your own RDS (leave FALSE on the first run) ----
# Set the real path and raw-count layer. Normalized values are optional.
# Check exact layer names. Join split layers with JoinLayers when needed.
# SCX does not choose how to combine layers.
use_own_data <- FALSE
input_rds <- "my-object.rds"
counts_layer <- "counts"
normalized_layer <- "data"  # set NULL for counts-only
if (use_own_data) {
  own <- readRDS(input_rds)
  if (!inherits(own,"Seurat")) stop("The RDS must contain a Seurat object.")
  assay <- SeuratObject::DefaultAssay(own)
  print(SeuratObject::Layers(own[[assay]],search=NA))
  own_scx <- file.path(run_dir,"own.scx")
  write_scx(own,own_scx,assay=assay,counts=counts_layer,normalized=normalized_layer)
  own_restored <- read_scx(own_scx)
  to_tsv(own_scx,file.path(run_dir,"own-preview.tsv.gz"),
    cells=head(colnames(own_restored),10),features=head(rownames(own_restored),10))
} else cat("Own-data example skipped. Enable it only when your path/layers are ready.\n")

# Results and next steps ----
tutorial_result <- list(passed=TRUE,SCX_VERSION=SCX_VERSION,shape=c(3L,4L),
  counts_equal=TRUE,normalized_max_absolute_error=normalized_error,
  typed_metadata_equal=TRUE,TSV_equal=TRUE,subset_order_equal=TRUE,existing_output_preserved=TRUE)
cat("SUCCESS: all demo checks passed. No SCX package was installed.\nFiles in this run:\n")
print(list.files(run_dir,recursive=TRUE))
# Results stay in scx-tutorial/output/run-*. The first TSV column is cell_id.
# Use subsets and gzip for large TSVs. Sparse input checks still need memory.
# TSV does not keep all metadata. Keep the original analysis object.
# More: https://github.com/junha-shin/scx-format/blob/main/docs/tsv-export.md
