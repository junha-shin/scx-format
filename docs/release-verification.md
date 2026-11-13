# Verify 0.1.0 release files

Standalone release **v0.1.0**, optional Python **0.1.0** / R **0.1.7**.
Format **0.1-alpha**, metadata schema **1.0**. Older tags/assets are unchanged.
The default installation-free instructions are in [setup](install.md) and the
downloadable [USAGE.md](https://github.com/junha-shin/scx-format/releases/download/v0.1.0/USAGE.md).

## Check downloaded bytes first

Download `scx.py` or `scx.R` and `SHA256SUMS` from the
[versioned release](https://github.com/junha-shin/scx-format/releases/tag/v0.1.0).
Compute SHA-256 and compare with the entry for that exact filename:

```sh
sha256sum scx.py             # Linux
shasum -a 256 scx.R          # macOS
```

```powershell
Get-FileHash .\scx.py -Algorithm SHA256
Get-FileHash .\scx.R -Algorithm SHA256
```

For an entire downloaded asset set, Linux `sha256sum -c SHA256SUMS` checks all
listed files. Downloading optional packages is unnecessary for checking a single
script; compare that file's line individually. A mismatch must be resolved before
executing code. Assets use canonical UTF-8/LF; editing or converting newlines
changes the checksum. Save the release file and checksum with the workflow.

After verification, import/source and check the release ID:

```python
import scx
assert scx.SCX_VERSION == "v0.1.0"
```

```r
source("scx.R")
stopifnot(SCX_VERSION == "v0.1.0")
```

This checks the loaded file version, not its integrity: hash checking comes
first. Neither test installs SCX or missing dependencies.

## Repeat the staged-file tests

From the tagged checkout in existing Python/R analysis environments, place the
verified files in `downloaded/` and run the same no-install checks used by CI:

```sh
python scripts/check_standalone.py --script downloaded/scx.py --output .compat-results/public-python
```

Set `SCX_SCRIPT=downloaded/scx.R` in your shell, then:

```sh
Rscript scripts/check_standalone.R .compat-results/public-r .compat-results/public-python
python scripts/check_standalone_bridge.py .compat-results/public-r downloaded/scx.py
```

For PowerShell, set the variable with `$env:SCX_SCRIPT="downloaded/scx.R"`.
Outputs must be fresh directories. The tests forbid installed SCX imports/loads
and R subprocesses, exercise Assay/Assay5, typed metadata, nonsquare axes,
asymmetric graphs, basic/analysis profiles, TSV/gzip/subsets, sparse zero/extreme
cases, checksum/corruption failures, overwrite protection and interrupted copies.
JSON reports are emitted on success. `--require-uninstalled` on the Python test
additionally proves the distribution is absent; normal users need not uninstall
anything because imports of it are already blocked by the test.

## Development and CI

```sh
python scripts/build_standalone.py --check
python scripts/stage_standalone_release.py dist/standalone
python scripts/stage_standalone_release.py dist/standalone --verify
python -m pytest -q
Rscript adapters/seurat/scripts/run-tests.R
```

Staging refuses an existing directory, checks package/release version agreement,
normalizes generated scripts/guide to LF, and writes SHA-256 entries. Tampered,
duplicate or unsafe entries fail validation. Three-platform CI stages these
exact files, tests them with SCX absent, verifies hashes again and retains them
with the reports. Main and tag CI must both pass before publishing.

Other CI jobs retain Python dependency profiles (including SciPy 1.11.4),
Seurat bridges, native R-only checks, three-platform optional package builds and
benchmark smoke tests. Optional packages use [their separate install guide](package-install.md).
The installed-package quickstart is `python examples/quickstart/verify.py`;
`--python-only` skips R. Real PBMC3k and scaling evidence remain separate recorded
measurements with their own source hashes, not new timings for this release.

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
