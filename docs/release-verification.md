# Verify the alpha release

The steps below are pinned to alpha.5 / R 0.1.3. Current main adds the
[R-only 0.1.4 validation procedure](r-only.md); published alpha.5 assets are
unchanged. Current CI also builds 0.1.4 and tests its native R-only paths.

Tag **v0.1.0-alpha.5** maps to Python **0.1.0a5** and R **0.1.3**.
Use the [unified installation guide](getting-started.md) to obtain the tagged
checkout, install the released wheel and R source package, and select the same
Python interpreter from R. No PyPI/CRAN publication is needed.

## External-user checks

The release includes `SHA256SUMS`. Verify downloaded package files before
installing, comparing the computed digest with the entry for that exact name:

```sh
# Linux
sha256sum scx_format-0.1.0a5-py3-none-any.whl
# macOS
shasum -a 256 scx_format-0.1.0a5-py3-none-any.whl
```

```powershell
Get-FileHash scx_format-0.1.0a5-py3-none-any.whl -Algorithm SHA256
```

After installing, from the tagged checkout run:

```sh
python examples/quickstart/verify.py --python-only --output quickstart-output-python
python examples/quickstart/verify.py --output quickstart-output-bridge
```

For a custom R installation append `--rscript "/path/to/Rscript"`. Both commands
require fresh outputs. The first verifies Python only; the second verifies the
installed R adapter and both cross-language directions. A successful run exits 0
and writes `report.json` with `passed: true`, versions, shapes and numeric errors.
`commands.json` records each public CLI/R subprocess and its output; it is written
even when a later command fails. No success report is emitted for a failed run.
Paths in command logs describe the local installation; inspect before sharing.

The tiny fixture is 3 cells × 4 features, with counts, log1p-normalized values,
one embedding, an asymmetric cell graph and nullable/ordered metadata. It
compares numeric values with tolerance 1e-12 and supported metadata exactly.
It checks transport, not biological results. It also does not replace fault
tests, compatibility tests or larger-data validation.

## Full repository tests

In a dedicated development environment, choose one profile and run:

```sh
python -m pip install -e ".[test]" -r compatibility/anndata-0.13.4.txt
python -m pytest -q
Rscript adapters/seurat/scripts/run-tests.R
Rscript adapters/seurat/scripts/compat-fixtures.R .r-tests/release-fixtures
python scripts/check_compatibility.py --fixtures .r-tests/release-fixtures --output .compat-results/release-bridge
```

R tests require testthat, SeuratObject, Matrix and jsonlite; install the R adapter
as described in the guide for the installed-package quickstart. To reproduce
other supported profiles, use their constraint files in **separate** Python
environments. AnnData 0.11.4 has separate Python 3.10 and 3.12 constraints; the
R bridge uses Python 3.12. See [compatibility](../compatibility/README.md).

CI runs Python compatibility tests, three R bridge environments, Windows safety
tests and Python/R benchmark smoke tests. The newest AnnData R job additionally
builds Python wheel/sdist and the R source package, installs the built packages,
and runs the cross-language quickstart without sourcing adapter code.
Reports and built packages are retained as CI artifacts. Tag pushes run the same
checks. CI does not automatically publish to GitHub, PyPI or CRAN.

Real-data validation and performance measurements are separate procedures:

- [Scanpy PBMC3k](../examples/pbmc3k/README.md) and
  [Seurat PBMC3k](../examples/seurat/README.md): pinned upstream checksums,
  provenance, additional analysis dependencies and recorded comparisons.
- [Scaling benchmarks](../benchmarks/README.md): workload, hardware/memory
  methods, source versions and all measured results. Historical speed figures
  predate the latest output-safety gates; they are not this release's timing.

## Maintainer build and release procedure

1. Finish changes on main and run the tests above. Update `CHANGELOG.md`, release
   notes and supported version mapping together. Use a clean checkout to avoid
   accidental development data in distributions.
2. Build Python packages with `python -m pip install build` then
   `python -m build`. Build the R source with `R CMD build adapters/seurat` and
   place the resulting `scxSeurat_0.1.3.tar.gz` in `dist/`.
3. In a fresh Python environment install the **wheel** plus the selected profile,
   and install the R tarball into a fresh R library. From outside the repository,
   verify that `scx.__file__` and `find.package('scxSeurat')` point to installed
   package libraries. Run the tagged-checkout quickstart using those environments.
   Repeat the Python check with the sdist in a separate fresh environment.
4. Record package SHA-256 hashes, full runtime dependency pins and quickstart reports. Do not include downloaded
   datasets or private absolute-path command logs in release assets.
5. Commit and push; wait for all CI checks. Create the annotated
   `v0.1.0-alpha.5` tag at the verified commit, push it and verify tag CI. Create a
   GitHub **prerelease** using the checked-in notes, then upload wheel, Python
   sdist, R source tarball, checksum list and sanitized verification report.
6. Check that asset names, sizes and SHA-256 digests match local packages and
   that the release/tag points to the tested commit. Download and exercise the
   public assets as a final install check. Do not replace an already published
   package silently; use a new version for any correction.
