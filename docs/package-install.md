# Optional package installation

Install SCX in the R or Python environment you already use. A checkout, new
environment and development tools are unnecessary for supported binary packages.
Existing compatible dependencies are reused.

## Python / scverse

```sh
pip install "https://github.com/junha-shin/scx-format/releases/download/v0.1.0-alpha.7/scx_format-0.1.0a7-py3-none-any.whl"
```

Run this in your AnnData environment; notebooks can use `%pip install`.
Requirements: Python >=3.10, NumPy >=1.24, pandas >=2, SciPy >=1.11.4.
Object/H5AD conversion needs AnnData >=0.11.4. If absent, install it with
`pip install "anndata>=0.11.4"`. Do not bypass dependencies with `--no-deps`.

## R / Seurat

```r
source("https://raw.githubusercontent.com/junha-shin/scx-format/v0.1.0-alpha.7/install.R")
library(scxSeurat)
```

The installer selects a binary, verifies SHA-256 and installs into your user
library. Only missing dependencies are installed: SeuratObject >=5, Matrix,
jsonlite and digest. Existing SeuratObject <5 is refused rather than upgraded.
Restart R before updating a loaded scxSeurat.

| Environment | Prebuilt coverage |
|---|---|
| R 4.5.x, Windows x86_64 | Binary ZIP |
| R 4.5.x, macOS arm64 | Binary TGZ; compatible macOS runtime required |
| R 4.5.x, Ubuntu >=22.04 x86_64 | Binary tarball built on Ubuntu 22.04 |
| Other R / OS / architecture | Source build below |

The macOS binary requires the build runner's macOS version or newer; consult
the release verification record. Other Linux distributions and older macOS
versions are not covered by the binary promise. Missing CRAN dependencies may
need build tools when your mirror has no suitable binary.

For an explicit source build:

```r
Sys.setenv(SCX_INSTALL_DEFINE_ONLY = "1")
source("https://raw.githubusercontent.com/junha-shin/scx-format/v0.1.0-alpha.7/install.R")
Sys.unsetenv("SCX_INSTALL_DEFINE_ONLY")
install_scx(source = TRUE)
```

Source builds need matching Rtools on Windows, compiler/zlib headers on Linux,
or command-line developer tools on macOS. Python is unnecessary.

## Offline and verification

Download the wheel or matching R archive and `SHA256SUMS` from the
[release](https://github.com/junha-shin/scx-format/releases/tag/v0.1.0-alpha.7).
Verify Python's checksum before installing the local wheel. For R, save
`install.R` from the tag, define it as above, then call
`install_scx(assets_dir="downloaded-assets", lib="your-library")`.
Dependencies must already be available offline.

Check `scx.__version__` or `packageVersion("scxSeurat")`. See the
[verification guide](release-verification.md) for an installed-package round trip.
Every dataset read/write still performs full validation.

These instructions describe the immutable alpha.7 packages. For no-install usage,
use the [standalone scripts](install.md) from current main.
