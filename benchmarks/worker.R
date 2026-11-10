args <- commandArgs(trailingOnly = TRUE)
if (args[1] == "generate") {
  output <- args[2]; cells <- as.integer(args[3]); genes <- as.integer(args[4])
  draws <- as.integer(args[5]); set.seed(as.integer(args[6]))
  counts <- Matrix::sparseMatrix(i = sample.int(genes, cells * draws, replace = TRUE),
                                 j = rep(seq_len(cells), each = draws),
                                 x = as.numeric(sample.int(19L, cells * draws, replace = TRUE)),
                                 dims = c(genes, cells))
  dimnames(counts) <- list(paste0("gene", seq_len(genes)), paste0("cell", seq_len(cells)))
  object <- SeuratObject::CreateSeuratObject(SeuratObject::CreateAssay5Object(counts = counts))
  SeuratObject::LayerData(object, layer = "data") <- log1p(counts)
  object$batch <- ordered(paste0("b", seq_len(cells) %% 8), levels = paste0("b", 0:8))
  object$qc <- rep(c(NA, rep(TRUE, 12)), length.out = cells)
  object$number <- seq_len(cells); object$number[seq(1, cells, 17)] <- NA_integer_
  object$text <- rep(c("", NA_character_, "null", "literal"), length.out = cells)
  object[["RNA"]] <- SeuratObject::AddMetaData(object[["RNA"]], data.frame(feature_number = seq_len(genes), row.names = rownames(counts)))
  for (name in c("pca", "umap")) {
    dimensions <- if (name == "pca") 30L else 2L
    x <- matrix(rnorm(cells * dimensions), cells, dimensions,
                dimnames = list(colnames(counts), paste0(toupper(name), "_", seq_len(dimensions))))
    object[[name]] <- SeuratObject::CreateDimReducObject(x, assay = "RNA", key = paste0(toupper(name), "_"))
  }
  k <- min(15L, cells)
  graph <- Matrix::sparseMatrix(i = rep(seq_len(cells), each = k), j = sample.int(cells, cells*k, replace = TRUE),
                                x = runif(cells*k), dims = c(cells, cells), dimnames = list(colnames(counts), colnames(counts)))
  graph <- SeuratObject::as.Graph(graph); SeuratObject::DefaultAssay(graph) <- "RNA"
  object[["RNA_nn"]] <- graph
  saveRDS(object, output, compress = FALSE)
  info <- list(cells = cells, genes = genes, nnz = length(counts@x),
               R = R.version.string, SeuratObject = as.character(packageVersion("SeuratObject")),
               Matrix = as.character(packageVersion("Matrix")), jsonlite = as.character(packageVersion("jsonlite")))
  writeLines(jsonlite::toJSON(info, auto_unbox = TRUE), paste0(output, ".json"))
} else {
  source(args[2])
  if (args[1] == "export") {
    if ("python" %in% names(formals(write_scx_seurat))) {
      write_scx_seurat(readRDS(args[3]), args[4], python = args[5])
    } else {
      # Keep historical adapters usable for before/after benchmark reproduction.
      write_scx_seurat(readRDS(args[3]), args[4])
    }
  } else if (args[1] == "import") {
    saveRDS(read_scx_seurat(args[3], python = args[5]), args[4], compress = FALSE)
  } else stop("unknown benchmark operation")
}
