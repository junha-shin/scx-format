"""Execute both public lessons twice, then test rejection of a damaged cache.

Requires the existing analysis dependencies, R, and tutorials/requirements-test.txt.
No user kernel registration or SCX package installation is performed.
"""
from __future__ import annotations

import argparse
import copy
import json
import os
from pathlib import Path
import subprocess
import sys

import nbformat
from nbclient import NotebookClient

ROOT = Path(__file__).resolve().parents[1]


def python_lesson(output: Path):
    work = output / "python"
    work.mkdir()
    kernel_root = output / "jupyter"
    kernel = kernel_root / "kernels" / "scx-tutorial-test"
    kernel.mkdir(parents=True)
    (kernel / "kernel.json").write_text(json.dumps({
        "argv": [sys.executable, "-m", "ipykernel_launcher", "-f", "{connection_file}"],
        "display_name": "Private tutorial test", "language": "python",
    }), encoding="utf-8")
    os.environ["JUPYTER_PATH"] = str(kernel_root) + os.pathsep + os.environ.get("JUPYTER_PATH", "")
    os.environ["JUPYTER_RUNTIME_DIR"] = str(output / "runtime")
    os.environ["IPYTHONDIR"] = str(output / "ipython")
    notebook = nbformat.read(ROOT / "tutorials/detailed.ipynb", as_version=4)
    nbformat.validate(notebook)
    client = NotebookClient(notebook, timeout=180, kernel_name="scx-tutorial-test",
                            resources={"metadata": {"path": str(work)}})

    def probe(source):
        # A test-only cell; never included in the published notebook.
        cell = nbformat.v4.new_code_cell(source)
        notebook.cells.append(cell)
        try:
            client.execute_cell(cell, len(notebook.cells) - 1)
        finally:
            notebook.cells.pop()

    with client.setup_kernel():
        probe("""
import sys, importlib.abc
class NoInstalledSCX(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, *args):
        if fullname == 'scx' or fullname.startswith('scx.'):
            raise AssertionError('An installed SCX package must not be used')
sys.meta_path.insert(0, NoInstalledSCX())
""")
        for index, cell in enumerate(notebook.cells):
            client.execute_cell(cell, index)
        executed = copy.deepcopy(notebook)
        probe("""
first_run = RUN_DIR
first_hash = hashlib.sha256(counts_tsv.read_bytes()).hexdigest()
first_cache_time = SCRIPT_PATH.stat().st_mtime_ns
def forbid_download(*args, **kwargs):
    raise AssertionError('Offline rerun attempted a download')
urllib.request.urlopen = forbid_download
""")
        for index, cell in enumerate(notebook.cells):
            client.execute_cell(cell, index)
        probe(r"""
assert RUN_DIR != first_run
assert hashlib.sha256((first_run / 'counts.tsv').read_bytes()).hexdigest() == first_hash
assert SCRIPT_PATH.stat().st_mtime_ns == first_cache_time
original = SCRIPT_PATH.read_bytes()
try:
    SCRIPT_PATH.write_bytes(original + b'\n# damaged cache\n')
    try:
        verified_file('scx.py')
    except RuntimeError as error:
        assert 'Checksum mismatch' in str(error)
    else:
        raise AssertionError('Damaged cache was accepted')
finally:
    SCRIPT_PATH.write_bytes(original)
tutorial_result.update(offline_rerun=True, prior_outputs_preserved=True, damaged_cache_rejected=True)
import json
Path('../python-report.json').write_text(json.dumps(tutorial_result, indent=2), encoding='utf-8')
""")
    for cell in executed.cells:
        cell.metadata.pop("execution", None)
    nbformat.write(executed, output / "executed.ipynb")
    return json.loads((output / "python-report.json").read_text(encoding="utf-8"))


def r_lesson(output: Path, rscript: str):
    work = output / "r"
    work.mkdir()
    wrapper = output / "check-r.R"
    wrapper.write_text('''args <- commandArgs(TRUE)
env <- new.env(parent=.GlobalEnv)
# Fail if the lesson tries to use the optional SCX package or Python.
trace(base::loadNamespace, tracer=quote({
  if (identical(package,"scxSeurat")) stop("Installed SCX package forbidden")
}), print=FALSE)
trace(base::system, tracer=quote(stop("Subprocess forbidden")), print=FALSE)
trace(base::system2, tracer=quote(stop("Subprocess forbidden")), print=FALSE)
sys.source(args[1],envir=env)
first_run <- env$run_dir
first_hash <- digest::digest(file.path(first_run,"counts.tsv"),file=TRUE,algo="sha256")
first_time <- file.info(env$tutorial_script)$mtime
trace(utils::download.file, tracer=quote(stop("Offline rerun attempted download")),print=FALSE)
sys.source(args[1],envir=env)
stopifnot(env$run_dir!=first_run,
  digest::digest(file.path(first_run,"counts.tsv"),file=TRUE,algo="sha256")==first_hash,
  identical(file.info(env$tutorial_script)$mtime,first_time))
original <- readBin(env$tutorial_script,"raw",n=file.info(env$tutorial_script)$size)
writeBin(c(original,charToRaw("\\n# damaged cache\\n")),env$tutorial_script)
refusal <- tryCatch({env$tutorial_verified_file("scx.R");NULL},error=conditionMessage)
writeBin(original,env$tutorial_script)
stopifnot(!is.null(refusal),grepl("Checksum mismatch",refusal,fixed=TRUE))
report <- c(env$tutorial_result,list(offline_rerun=TRUE,prior_outputs_preserved=TRUE,damaged_cache_rejected=TRUE))
jsonlite::write_json(report,args[2],auto_unbox=TRUE,pretty=TRUE)
''', encoding="utf-8")
    result = subprocess.run([rscript, "--vanilla", str(wrapper),
                             str(ROOT / "tutorials/detailed.R"), str(output / "r-report.json")],
                            cwd=work, capture_output=True, text=True, encoding="utf-8", errors="replace")
    (output / "r-console.txt").write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(result.stdout + result.stderr)
    return json.loads((output / "r-report.json").read_text(encoding="utf-8"))


def quick_lessons(output: Path, rscript: str):
    import shutil
    import pandas as pd
    for language, filename in (("python", "scx.py"), ("r", "scx.R")):
        work = output / ("quick-" + language)
        work.mkdir()
        shutil.copyfile(ROOT / "standalone" / filename, work / filename)
    work = output / "quick-python"
    notebook = nbformat.read(ROOT / "tutorials/quickstart.ipynb", as_version=4)
    nbformat.validate(notebook)
    notebook.cells.insert(0, nbformat.v4.new_code_cell("""
import anndata as ad
import numpy as np
from scipy import sparse
counts = sparse.csr_matrix([[1,0,3],[0,2,1]])
adata = ad.AnnData(X=np.log1p(counts.astype(float)))
adata.layers['counts'] = counts
"""))
    notebook.cells.append(nbformat.v4.new_code_cell("""
from pathlib import Path
import hashlib
assert Path(scx.__file__).resolve() == Path('scx.py').resolve()
np.testing.assert_array_equal(restored.layers['counts'].toarray(), counts.toarray())
np.testing.assert_allclose(restored.X.toarray(), adata.X.toarray(), rtol=0, atol=1e-12)
missing = adata.copy()
del missing.layers['counts']
try:
    scx.write_scx(missing, 'missing-counts.scx')
except ValueError as error:
    assert 'counts' in str(error) and 'normalized=None' in str(error)
else:
    raise AssertionError('Missing counts accepted')
assert not Path('missing-counts.scx').exists()
try:
    scx.read_scx('absent.scx')
except (OSError, ValueError):
    pass
else:
    raise AssertionError('Missing input accepted')
original = hashlib.sha256(Path('counts.tsv').read_bytes()).hexdigest()
try:
    scx.to_tsv('sample.scx', 'counts.tsv')
except FileExistsError:
    pass
else:
    raise AssertionError('Expected overwrite refusal')
assert hashlib.sha256(Path('counts.tsv').read_bytes()).hexdigest() == original
"""))
    NotebookClient(notebook, timeout=180, kernel_name="scx-tutorial-test",
                   resources={"metadata": {"path": str(work)}}).execute()
    notebook.cells.pop()
    notebook.cells.pop(0)
    for cell in notebook.cells:
        cell.metadata.pop("execution", None)
    nbformat.write(notebook, output / "quickstart-executed.ipynb")
    wrapper = output / "check-quick.R"
    wrapper.write_text('''args <- commandArgs(TRUE)
source("scx.R")
counts <- Matrix::Matrix(matrix(c(1,0,3,0,2,1),nrow=3),sparse=TRUE)
dimnames(counts) <- list(c("gene1","gene2","gene3"),c("cell1","cell2"))
object <- SeuratObject::CreateSeuratObject(counts=counts)
SeuratObject::LayerData(object,layer="data") <- log1p(counts)
source(args[1])
stopifnot(max(abs(SeuratObject::LayerData(restored,layer="counts")-counts))==0,
  max(abs(SeuratObject::LayerData(restored,layer="data")-log1p(counts)))<=1e-12)
bad <- object
SeuratObject::DefaultAssay(bad) <- "RNA"
SeuratObject::LayerData(bad,layer="counts") <- NULL
missing <- tryCatch({write_scx(bad,"missing-counts.scx");NULL},error=conditionMessage)
stopifnot(!is.null(missing),grepl("counts",missing),!dir.exists("missing-counts.scx"))
absent <- tryCatch({read_scx("absent.scx");NULL},error=conditionMessage)
stopifnot(!is.null(absent))
previous <- digest::digest("counts.tsv",file=TRUE,algo="sha256")
refusal <- tryCatch({to_tsv("sample.scx","counts.tsv");NULL},error=conditionMessage)
stopifnot(!is.null(refusal),grepl("already exists",refusal,fixed=TRUE),
  previous==digest::digest("counts.tsv",file=TRUE,algo="sha256"))
''', encoding="utf-8")
    result = subprocess.run([rscript, "--vanilla", str(wrapper), str(ROOT / "tutorials/quickstart.R")],
                            cwd=output / "quick-r", capture_output=True, text=True,
                            encoding="utf-8", errors="replace")
    (output / "quick-r-console.txt").write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(result.stdout + result.stderr)
    expected = [[1, 0, 3], [0, 2, 1]]
    for language in ("python", "r"):
        table = pd.read_csv(output / ("quick-" + language) / "counts.tsv", sep="\t", index_col=0)
        assert table.values.tolist() == expected
    return {"passed": True, "counts_and_normalized_roundtrip": True,
            "TSV_values_equal": True, "overwrite_protection": True,
            "missing_counts_rejected": True, "missing_input_rejected": True,
            "SCX_VERSION": json.loads((ROOT / "release.json").read_text(encoding="utf-8"))["tag"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="New output directory")
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--save-executed", type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    python = python_lesson(output)
    r = r_lesson(output, args.rscript)
    assert python["shape"] == r["shape"] == [3, 4]
    assert python["SCX_VERSION"] == r["SCX_VERSION"] == "v0.1.0-alpha.7"
    for language in (python, r):
        assert all(language[key] is True for key in (
            "passed", "counts_equal", "typed_metadata_equal", "TSV_equal", "subset_order_equal",
            "existing_output_preserved", "offline_rerun", "prior_outputs_preserved", "damaged_cache_rejected"))
    # Compare actual text exports across both languages, including compression.
    import pandas as pd
    p_run = sorted((output / "python/scx-tutorial/output").iterdir())[0]
    r_run = sorted((output / "r/scx-tutorial/output").iterdir())[0]
    for name in ("counts.tsv", "normalized.tsv.gz", "subset.tsv.gz"):
        p = pd.read_csv(p_run / name, sep="\t", index_col=0)
        r_frame = pd.read_csv(r_run / name, sep="\t", index_col=0)
        pd.testing.assert_frame_equal(p, r_frame, check_exact=False, atol=1e-12, rtol=0)
    report = {"passed": True, "cross_language_TSV_equal": True, "python": python, "R": r, "quickstart": quick_lessons(output, args.rscript)}
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if args.save_executed:
        args.save_executed.write_bytes((output / "executed.ipynb").read_bytes())
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
