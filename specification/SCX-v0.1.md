# SCX v0.1-alpha specification

Status: draft. Scope: one scRNA-seq feature space.

## Layout and axes

An SCX dataset is a directory. It MUST contain `manifest.tsv`, `cells.tsv.gz`, `features.tsv.gz`, and `matrices/counts.mtx.gz`. All paths in the manifest are relative to the dataset root. The expression matrix MUST be CELL × FEATURE; the cell and feature table row orders define its axes. A cell graph is CELL × CELL, and an embedding is CELL × COMPONENT.

The cell and feature files MUST be UTF-8, tab-separated, gzip-compressed tables with a header. Their first columns MUST be `cell_id` and `feature_id`, respectively. IDs MUST be nonempty and unique. Additional columns are optional metadata; empty fields denote missing values. A metadata value containing tab or newline is unsupported in this alpha version.

The manifest MUST be UTF-8 TSV with the header `path\trole\tname`. Each file other than the manifest MUST have exactly one row. Required rows have roles `cells`, `features`, and `counts`, and names `cells`, `features`, and `counts`. Optional roles are `normalized`, `embedding`, and `cell_graph`; names MUST be unique within a role. Paths MUST be safe relative paths, and all referenced files MUST exist. Unlisted files have no SCX meaning.

Matrix files MUST be gzip-compressed Matrix Market coordinate files, numeric real or integer, with one-based coordinates. Duplicate coordinates are summed by Matrix Market readers. Counts MUST have nonnegative finite values. Normalized matrices have the same shape as counts; cell graphs are square with one row and column per cell. Embeddings are gzipped UTF-8 TSVs, with `cell_id` as the first column and one or more numeric component columns. Their row IDs MUST match `cells.tsv.gz` in order.

Recommended locations are `matrices/normalized.mtx.gz`, `embeddings/<name>.tsv.gz`, and `graphs/<name>.mtx.gz`. Names should use letters, digits, `_`, `-`, and `.`.

## AnnData mapping

The adapter takes counts from `layers['counts']`; `X` maps to normalized. Cell and feature tables map to `obs` and `var`, embeddings to `obsm`, and cell graphs to `obsp`. On import, `X` is normalized if available, otherwise counts. Other layers, `varm`, `varp`, `uns`, and pandas dtype fidelity are not covered by this alpha profile.

## Planned

Seurat support will transpose assay matrices from FEATURE × CELL to the canonical CELL × FEATURE orientation. Multi-assay, multiome, spatial data, and arbitrary object serialization are outside v0.1-alpha.
