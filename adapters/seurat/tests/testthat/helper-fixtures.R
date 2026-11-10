tiny_seurat <- function(assay_class = "Assay5") {
  counts <- Matrix::Matrix(matrix(c(1, 0, 3, 0, 2, 4), 3, 2), sparse = TRUE)
  dimnames(counts) <- list(c("g1", "g2", "g3"), c("c1", "c2"))
  assay <- if (assay_class == "Assay") SeuratObject::CreateAssayObject(counts = counts) else SeuratObject::CreateAssay5Object(counts = counts)
  object <- SeuratObject::CreateSeuratObject(assay)
  SeuratObject::LayerData(object, layer = "data") <- log1p(counts)
  object$integer <- c(1L, NA_integer_)
  object$double <- c(0.12345678901234567, 1.234567890123456e-30)
  object$logical <- c(TRUE, NA)
  object$category <- ordered(c("b", NA), levels = c("b", "a", "unused"))
  object$text <- c("", "null")
  object[["RNA"]] <- SeuratObject::AddMetaData(object[["RNA"]], data.frame(score = c(1.2, NA, 3.4), row.names = rownames(counts)))
  embedding <- matrix(c(1.23456789, 2, 3, 4), 2, dimnames = list(colnames(counts), c("PC_1", "PC_2")))
  object[["pca"]] <- SeuratObject::CreateDimReducObject(embedding, assay = "RNA", key = "PC_")
  graph <- Matrix::Matrix(matrix(c(0, 2, 1, 0), 2), sparse = TRUE)
  dimnames(graph) <- list(colnames(counts), colnames(counts))
  graph <- SeuratObject::as.Graph(graph)
  SeuratObject::DefaultAssay(graph) <- "RNA"
  object[["RNA_nn"]] <- graph
  object
}
