# SCX v0.1-alpha specification

Status: draft. Scope: one scRNA-seq feature space.

The manifest MUST begin with `# scx_version=0.1-alpha` on its own line, followed by the TSV header. Earlier unversioned draft manifests must be regenerated. Matrix Market symmetry MUST be `general`. Table headers MUST be nonempty and unique, and every record MUST have the same number of fields.

Writers MUST refuse existing output destinations by default. Replacement MUST
require explicit opt-in, stage and validate the new dataset before modifying the
old one, and retain the previous dataset until publication succeeds. Failed
replacement MUST restore the old dataset or retain a recoverable backup and
report its path. Implementations MUST document interruption and concurrent-reader
limits; directory replacement need not be an atomic reader snapshot. See the
reference implementation's [safe I/O contract](../docs/safe-io.md).
New writers MUST emit the [metadata schema 1.0 extension](metadata-schema-v1.md), preserving numeric/boolean/categorical types and missing values. Legacy datasets without this extension load as text with a warning. Unsupported analysis objects omitted by an adapter MUST produce a warning.

## Layout and axes

An SCX dataset is a directory. It MUST contain `manifest.tsv`, `cells.tsv.gz`, `features.tsv.gz`, and `matrices/counts.mtx.gz`. All paths in the manifest are relative to the dataset root. The expression matrix MUST be CELL × FEATURE; the cell and feature table row orders define its axes. A cell graph is CELL × CELL, and an embedding is CELL × COMPONENT.

The cell and feature files MUST be UTF-8, tab-separated, gzip-compressed tables with a header. Their first columns MUST be `cell_id` and `feature_id`, respectively. IDs MUST be nonempty and unique. With a metadata schema, other fields contain JSON scalars, including explicit `null` for missing values. Without a schema, fields are legacy plain text and empty values cannot distinguish missing from empty strings. Decoded TSV fields MUST NOT contain literal tabs or newlines; JSON escaping permits these inside typed string values.

The manifest MUST be UTF-8 TSV with the header `path\trole\tname` after its version line. Each file other than the manifest MUST have exactly one row. Required rows have roles `cells`, `features`, and `counts`, and names `cells`, `features`, and `counts`. Additional roles are `normalized`, `embedding`, `cell_graph` and `metadata_schema`; names MUST be unique within a role. Paths MUST be safe relative paths, and all referenced files MUST exist. Unlisted files have no SCX meaning.

Matrix files MUST be gzip-compressed Matrix Market coordinate files, numeric real or integer, with one-based coordinates. Duplicate coordinates are summed by Matrix Market readers. Counts MUST have nonnegative finite values. Normalized matrices have the same shape as counts; cell graphs are square with one row and column per cell. Embeddings are gzipped UTF-8 TSVs, with `cell_id` as the first column and one or more numeric component columns. Their row IDs MUST match `cells.tsv.gz` in order.

Recommended locations are `matrices/normalized.mtx.gz`, `embeddings/<name>.tsv.gz`, and `graphs/<name>.mtx.gz`. Names should use letters, digits, `_`, `-`, and `.`.

## AnnData mapping

The adapter takes counts from `layers['counts']`; `X` maps to normalized. Typed cell and feature tables map to `obs` and `var`, embeddings to `obsm`, and cell graphs to `obsp`. On import, `X` is normalized if available, otherwise counts. Other layers, `varm`, `varp`, `uns`, and metadata types outside schema 1.0 are not covered by this alpha profile.

## Planned

The Seurat MVP adapter transposes assay matrices from FEATURE × CELL to the canonical CELL × FEATURE orientation and uses metadata schema 1.0 with explicit native R type mappings. Multi-assay, multiome, spatial data, and arbitrary object serialization are outside v0.1-alpha.
