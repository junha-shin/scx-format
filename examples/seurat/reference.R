# Independent native-Seurat reference snapshot; no SCX metadata encoder used.
save_reference <- function(object, output) {
  saveRDS(object, file.path(output, "source.rds"))
  Matrix::writeMM(SeuratObject::LayerData(object, layer = "counts"), file.path(output, "reference_counts.mtx"))
  Matrix::writeMM(SeuratObject::LayerData(object, layer = "data"), file.path(output, "reference_normalized.mtx"))
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
                    graphs = as.list(SeuratObject::Graphs(object)))
  for (name in SeuratObject::Graphs(object)) Matrix::writeMM(object[[name]], file.path(output, paste0("reference_graph_", name, ".mtx")))
  writeLines(.scx_json(reference), file.path(output, "reference.json"), useBytes = TRUE)
}
