# First use in your existing environment

Download scx.py or scx.R from v0.1.0 into your working folder.
Try the [Python notebook](../tutorials/quickstart.ipynb) or [R guide](../tutorials/quickstart.md).
They load PBMC3k, show matrix sizes, and export a small TSV.
Then replace the loading step with your own object.
The detailed examples use fixed alpha.7 files so their checks can be repeated.

| Environment | Existing dependencies | Tested scope |
|---|---|---|
| Scanpy / AnnData | Python >=3.10, NumPy >=1.24, pandas >=2, SciPy >=1.11.4, AnnData >=0.11.4 | AnnData 0.11.4, 0.12.3, 0.13.4; Linux, Windows, macOS standalone checks |
| Seurat | R, SeuratObject >=5, Matrix, jsonlite, digest | R 4.5.2, existing Assay and Assay5; Linux, Windows, macOS standalone checks |

These packages must already be present. Older environments may not meet the requirements.
SCX does not install or upgrade packages. R does not need Python or native SCX helpers.

## Counts choices

Python reads raw counts from `layers['counts']` and normalized values from `X`.
If `X` holds raw counts, use `write_scx(adata, "sample.scx", counts="X", normalized=None)`.
Use `counts="raw_counts"` to select a layer with that name.
Use `normalized=None` for counts only. SCX does not choose `raw.X` for you.

R uses counts/data in the active assay. Select other exact layer names when needed.
Use `normalized=NULL` for counts only. Join split Assay5 layers yourself when needed.

## Common errors

| Problem | Recovery |
|---|---|
| Cannot import scx or source scx.R | Save the release file in the current working folder; restart the Python kernel if another scx was loaded |
| Missing dependency | Choose your existing analysis environment with the dependencies above |
| Missing counts | Inspect adata.layers or SeuratObject::Layers(object[[assay]], search=NA), then select the actual raw-count source |
| Missing normalized layer | Select its actual name, or set normalized=None (R: NULL) for counts-only export |
| Input does not exist | Check the working folder and full input path; read_scx expects an SCX directory |
| Output already exists | Choose a new output name; existing results are protected |
| Corrupt or incomplete SCX | Restore/regenerate the named file from the original analysis object; see safe-io.md |

## Verification

Standalone CI installs analysis dependencies but not SCX.
It tests Assay, Assay5, and both conversion directions on three operating systems.
The checks block installed SCX imports and R subprocess calls.
PBMC3k checks compare full sparse counts, axis order, data hashes, and TSV previews.
Other tests cover bad input and existing output. Detailed examples also check
metadata, offline reruns, damaged files, and subset order.
The lessons and converters do not install dependencies.

Run `python scripts/check_tutorials.py --output .compat-results/first-use`
in a prepared test environment. Use a new output path.
The reports cover quickstart and detailed examples.
See [type limits](../compatibility/README.md) and [real-data tests](../examples/pbmc3k/README.md).
