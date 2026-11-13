# Verify 0.1.0 release files

This guide covers v0.1.0 files, Python package 0.1.0, and R package 0.1.7.
The format is 0.1-alpha; the metadata schema is 1.0.
See [setup](install.md) or the downloadable
[USAGE.md](https://github.com/junha-shin/scx-format/releases/download/v0.1.0/USAGE.md).
Old tags and files are unchanged.

## Check downloaded bytes first

Download scx.py or scx.R and SHA256SUMS from the
[release](https://github.com/junha-shin/scx-format/releases/tag/v0.1.0).
Compute SHA-256 and compare it with that file's entry:

```sh
sha256sum scx.py             # Linux
shasum -a 256 scx.R          # macOS
```

```powershell
Get-FileHash .\scx.py -Algorithm SHA256
Get-FileHash .\scx.R -Algorithm SHA256
```

Use `sha256sum -c SHA256SUMS` only when all listed files are present.
For one script, compare its hash with its own entry; optional packages are not needed.
Do not run a file whose hash differs. Keep the script and checksum with your analysis.
Release files use UTF-8 and LF. Editing them or changing line endings changes the hash.

After checking the hash, load the file and check its release ID:

```python
import scx
assert scx.SCX_VERSION == "v0.1.0"
```

```r
source("scx.R")
stopifnot(SCX_VERSION == "v0.1.0")
```

The version check does not check file integrity. Check the hash first.
These checks do not install SCX or dependencies.

## Repeat the staged-file tests

Use a checkout of the tag and your existing analysis environments.
Put checked release files in downloaded/, then run the CI checks:

```sh
python scripts/check_standalone.py --script downloaded/scx.py --output .compat-results/public-python
```

Set `SCX_SCRIPT=downloaded/scx.R`, then run:

```sh
Rscript scripts/check_standalone.R .compat-results/public-r .compat-results/public-python
python scripts/check_standalone_bridge.py .compat-results/public-r downloaded/scx.py
```

In PowerShell, use `$env:SCX_SCRIPT="downloaded/scx.R"`.
Choose new output folders. These tests block installed SCX loads and R subprocesses.
They check Assay/Assay5, metadata, axes, directed graphs, write profiles,
TSV/gzip/subsets, sparse edge cases, damaged files, existing output, and interrupted copies.
Successful runs write JSON reports.
`--require-uninstalled` also checks that the Python package is absent.
Normal users do not need to uninstall a package to use the one-file tools.

## Development and CI

```sh
python scripts/build_standalone.py --check
python scripts/stage_standalone_release.py dist/standalone
python scripts/stage_standalone_release.py dist/standalone --verify
python -m pytest -q
Rscript adapters/seurat/scripts/run-tests.R
```

Staging requires a new folder. It checks version agreement and writes UTF-8/LF files and hashes.
Bad, duplicate, or unsafe checksum entries fail.
CI tests these files without SCX installed on three systems, rechecks their hashes,
and saves files and reports. Main and tag CI must pass before a release is published.

Other jobs test dependency profiles, including SciPy 1.11.4, Seurat bridges,
R-only use, optional packages, and benchmark smoke cases.
See [optional package setup](package-install.md).
Run `python examples/quickstart/verify.py` for the installed-package example;
`--python-only` skips R.
PBMC3k and scaling reports are separate measurements with their own source hashes.
They are not new timing measurements for this release.

## Maintainer publication

1. Update `release.json`, Python/R package versions, installer, changelog and
   guides together. Regenerate the one-file scripts and run the tests.
2. Commit/push; verify all main CI jobs. Tag that exact commit and verify tag CI.
3. Use the canonical script artifact from Linux; check that staged script and
   usage hashes agree across Linux/Windows/macOS. Compare to tagged raw files.
4. Create a GitHub prerelease draft. Add scripts, USAGE, optional package assets,
   one combined `SHA256SUMS` and a sanitized `verification.json` with source,
   platform and main/tag CI evidence. Keep the no-install assets primary.
5. Verify uploaded sizes/hashes and tag target before publication. Download the
   public assets, verify hashes and repeat the no-install checks using those
   downloaded files. Attach public-install validation evidence if separate.
6. Do not replace published scripts, packages or tags. Any correction needs a
   new version. No PyPI/CRAN publication is performed by this process.
