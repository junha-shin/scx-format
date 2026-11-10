#' Write the basic SCX profile, optionally including analysis coordinates/graphs.
#' @export
write_scx <- function(object, path, counts = "counts", normalized = "data",
                      analysis = FALSE, overwrite = FALSE,
                      assay = SeuratObject::DefaultAssay(object)) {
  if (!is.logical(analysis) || length(analysis) != 1L || is.na(analysis)) stop("analysis must be TRUE or FALSE")
  omitted <- character()
  if (!analysis && (length(SeuratObject::Reductions(object)) || length(SeuratObject::Graphs(object))))
    omitted <- c(omitted, "embeddings/graphs (include with analysis=TRUE)")
  if (length(SeuratObject::Assays(object)) > 1L) omitted <- c(omitted, "other assays")
  layers <- SeuratObject::Layers(object[[assay]], search = NA)
  if (length(setdiff(layers, c(counts, normalized)))) omitted <- c(omitted, "other expression layers")
  loadings <- any(vapply(SeuratObject::Reductions(object), function(name)
    length(SeuratObject::Loadings(object[[name]])) > 0L, logical(1)))
  if (length(object@commands) || length(object@misc) || loadings)
    omitted <- c(omitted, "analysis internals outside the SCX profile")
  result <- write_scx_seurat(object, path, assay = assay, counts_layer = counts,
    normalized_layer = if (identical(counts, normalized)) NULL else normalized,
    reductions = if (analysis) SeuratObject::Reductions(object) else character(),
    graphs = if (analysis) SeuratObject::Graphs(object) else character(),
    overwrite = overwrite, warn_scope = FALSE)
  if (length(omitted)) warning("SCX omitted: ", paste(omitted, collapse = "; "), call. = FALSE)
  invisible(result)
}

#' Validate and read an SCX directory into Seurat.
#' @export
read_scx <- function(path, assay = "RNA", assay_class = c("Assay5", "Assay"))
  read_scx_seurat(path, assay = assay, assay_class = match.arg(assay_class))
