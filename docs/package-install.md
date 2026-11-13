# Optional package installation

Package installation is optional. Use your existing Python or R environment.
Supported binary packages do not need a checkout, a new environment, or build tools.
They reuse compatible dependencies.

## Python / scverse

```sh
pip install "https://github.com/junha-shin/scx-format/releases/download/v0.1.0/scx_format-0.1.0-py3-none-any.whl"
```

Run this in your AnnData environment. In a notebook, you can use `%pip install`.
Requirements: Python >=3.10, NumPy >=1.24, pandas >=2, and SciPy >=1.11.4.
Object and H5AD conversion also need AnnData >=0.11.4.
If needed, install it with `pip install "anndata>=0.11.4"`.
Do not skip dependency checks with `--no-deps`.

## R / Seurat

```r
source("https://raw.githubusercontent.com/junha-shin/scx-format/v0.1.0/install.R")
library(scxSeurat)
```

The installer picks a binary, checks SHA-256, and installs it in your user library.
It installs only missing dependencies: SeuratObject >=5, Matrix, jsonlite, and digest.
It refuses SeuratObject <5 rather than upgrading it.
Restart R before updating a loaded scxSeurat package.

| Environment | Prebuilt coverage |
|---|---|
| R 4.5.x, Windows x86_64 | Binary ZIP |
| R 4.5.x, macOS arm64 | Binary TGZ; compatible macOS runtime required |
| R 4.5.x, Ubuntu >=22.04 x86_64 | Binary tarball built on Ubuntu 22.04 |
| Other R / OS / architecture | Source build below |

The macOS binary needs the build runner's macOS version or newer.
Check the release record first. Other Linux distributions and older macOS versions are not covered.
Missing CRAN dependencies may need build tools if no suitable binary is available.

To build from source:

```r
Sys.setenv(SCX_INSTALL_DEFINE_ONLY = "1")
source("https://raw.githubusercontent.com/junha-shin/scx-format/v0.1.0/install.R")
Sys.unsetenv("SCX_INSTALL_DEFINE_ONLY")
install_scx(source = TRUE)
```

Source builds need matching Rtools on Windows, a compiler and zlib headers on Linux,
or command-line developer tools on macOS. They do not need Python.

## Offline and verification

Download the wheel or matching R archive and SHA256SUMS from the
[release](https://github.com/junha-shin/scx-format/releases/tag/v0.1.0).
Check the Python wheel's hash before installing it.
For R, save install.R from the tag and define it as shown above.
Then call `install_scx(assets_dir="downloaded-assets", lib="your-library")`.
For offline use, dependencies must already be available.

Check `scx.__version__` or `packageVersion("scxSeurat")`.
See the [verification guide](release-verification.md) for a round-trip test.
Reads and writes still check the whole dataset.

These steps use stable 0.1.0 packages.
For use without an SCX install, see [the one-file tools](install.md).
