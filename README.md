# SCX - Single Cell eXchange

Move single-cell data between Seurat and AnnData.
SCX stores each expression matrix as **cells x features**.

## Try PBMC3k

Use your existing Scanpy or Seurat environment. No SCX package install is needed.

- [Python notebook](tutorials/quickstart.ipynb)
- [R guide with results](tutorials/quickstart.md) | [R script](tutorials/quickstart.R)

The examples reuse local PBMC3k data or download it once.
They save and read 2700 cells x 32738 genes, then show a 5 x 5 TSV table.

## Use your own data

Download [scx.py](https://github.com/junha-shin/scx-format/releases/download/v0.1.0/scx.py)
or [scx.R](https://github.com/junha-shin/scx-format/releases/download/v0.1.0/scx.R)
into your working folder.
Check the file against [SHA256SUMS](https://github.com/junha-shin/scx-format/releases/download/v0.1.0/SHA256SUMS)
before you load it. See [setup and requirements](docs/install.md).

Python:

```python
import scx
scx.write_scx(adata, "sample.scx")
restored = scx.read_scx("sample.scx")
```

R:

```r
source("scx.R")
write_scx(seurat_object, "sample.scx")
restored <- read_scx("sample.scx")
```

Python reads counts from `layers['counts']` and normalized values from `X`.
R reads `counts` and `data` from the selected assay.
Counts are required. [Choose another source or export counts only](docs/first-use.md).

## What SCX keeps

- Counts and optional normalized values.
- Cell and feature metadata, including supported data types and missing values.
- Embeddings and cell graphs when you set `analysis=True` (R: `analysis=TRUE`).

SCX supports one assay. Other assays, loadings, spatial data, multiome data,
and other analysis objects are outside its scope. Keep your original object.

## Export TSV

Use `scx.to_tsv("sample.scx", "counts.tsv")` in Python,
or `to_tsv("sample.scx", "counts.tsv")` in R.
You can also export `.tsv.gz` or select cells and genes.
[See TSV options and memory needs](docs/tsv-export.md).

## Save safely

SCX checks data when it reads and writes. Existing outputs are protected.
The R file writes only to new folders.
[See overwrite options and recovery](docs/safe-io.md).

## More help

- [Common errors](docs/first-use.md) | [Detailed examples](tutorials/README.md#detailed-walkthrough)
- [Supported versions and limits](compatibility/README.md)
- [Format](specification/SCX-v0.1.md) | [Metadata types](specification/metadata-schema-v1.md)
- [Tests with real data](examples/pbmc3k/README.md) | [R tests](examples/seurat/README.md) | [Speed and memory](benchmarks/README.md)
- [Release checks](docs/release-verification.md) | [Changes](CHANGELOG.md)

## Development

```sh
pip install -e '.[test]'
pytest
```

After changing the adapters, run `python scripts/build_standalone.py`.
CI tests Python and R on Linux, Windows and macOS without installing SCX.
See the [R development guide](docs/r-only.md).
