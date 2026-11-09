# Rscript examples/seurat/prepare.R <10x hg19 directory> <new output directory>
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L, !dir.exists(args[2]))
source("adapters/seurat/R/scx.R")
set.seed(42)
counts <- Seurat::Read10X(args[1], gene.column = 1L)
object <- SeuratObject::CreateSeuratObject(counts, min.cells = 3L, project = "PBMC3k")
object <- Seurat::NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
object <- Seurat::FindVariableFeatures(object, nfeatures = 1000L, verbose = FALSE)
object <- Seurat::ScaleData(object, features = SeuratObject::VariableFeatures(object), verbose = FALSE)
object <- Seurat::RunPCA(object, npcs = 30L, seed.use = 42L, verbose = FALSE)
object <- Seurat::FindNeighbors(object, dims = 1:30, k.param = 15L, verbose = FALSE)
object <- Seurat::RunUMAP(object, dims = 1:30, n.neighbors = 15L, seed.use = 42L, verbose = FALSE)
total <- object$nCount_RNA
object$library_size <- ordered(ifelse(total >= median(total), "high", "low"), levels = c("low", "high", "unused"))
object$library_size[1] <- NA
object$nullable_integer <- as.integer(total)
object$nullable_integer[seq(1, ncol(object), 10)] <- NA_integer_
object$nullable_logical <- total > 0
object$nullable_logical[seq(1, ncol(object), 10)] <- NA
object$text <- rep(c("", "null", NA_character_, "literal"), length.out = ncol(object))
dir.create(args[2], recursive = TRUE)
saveRDS(object, file.path(args[2], "source.rds"))
# Independent reference dump in Seurat's native FEATURE x CELL orientation.
Matrix::writeMM(SeuratObject::LayerData(object, layer = "counts"), file.path(args[2], "reference_counts.mtx"))
Matrix::writeMM(SeuratObject::LayerData(object, layer = "data"), file.path(args[2], "reference_normalized.mtx"))
meta_reference <- function(frame) lapply(names(frame), function(name) {
  x <- frame[[name]]
  list(name = name, r_type = typeof(x), r_class = as.list(class(x)),
       values = as.list(if (is.factor(x)) as.character(x) else x),
       levels = if (is.factor(x)) as.list(levels(x)) else NULL,
       ordered = is.ordered(x))
})
reference <- list(cells = as.list(colnames(object)), features = as.list(rownames(object)),
                  obs = meta_reference(object[[]]), var = meta_reference(object[["RNA"]][[]]),
                  embeddings = lapply(SeuratObject::Reductions(object), function(name) list(name = name, values = SeuratObject::Embeddings(object[[name]]))),
                  graphs = as.list(SeuratObject::Graphs(object)),
                  versions = as.list(setNames(vapply(c("Seurat", "SeuratObject", "Matrix", "jsonlite"), function(p) as.character(packageVersion(p)), character(1)), c("Seurat", "SeuratObject", "Matrix", "jsonlite"))),
                  R = R.version.string)
for (name in SeuratObject::Graphs(object)) Matrix::writeMM(object[[name]], file.path(args[2], paste0("reference_graph_", name, ".mtx")))
writeLines(.scx_json(reference), file.path(args[2], "reference.json"), useBytes = TRUE)
writeLines(sub("[[:space:]]+$", "", capture.output(sessionInfo())), file.path(args[2], "R-session.txt"))
write.table(installed.packages()[, c("Package", "Version")], file.path(args[2], "R-packages.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write_scx_seurat(object, file.path(args[2], "seurat.scx"))
