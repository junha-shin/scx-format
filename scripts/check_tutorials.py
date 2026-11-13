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
    """Run the actual PBMC3k lessons, preserving notebook/Markdown output."""
    import shutil
    import re
    import hashlib
    import pandas as pd
    for language, filename in (("python", "scx.py"), ("r", "scx.R")):
        work = output / ("quick-" + language)
        work.mkdir()
        shutil.copyfile(ROOT / "standalone" / filename, work / filename)
    # Reuse previously downloaded real data, without committing raw datasets.
    for source, destination in (
        (ROOT / "examples/pbmc3k/data/pbmc3k_raw.h5ad", output / "quick-python/pbmc3k_raw.h5ad"),
        (ROOT / "examples/seurat/data/pbmc3k.tar.gz", output / "quick-r/pbmc3k.tar.gz"),
    ):
        if source.exists():
            shutil.copyfile(source, destination)
    notebook = nbformat.read(ROOT / "tutorials/quickstart.ipynb", as_version=4)
    nbformat.validate(notebook)
    notebook.cells.append(nbformat.v4.new_code_cell("""
assert adata.shape == restored.shape == (2700, 32738)
assert (restored.layers['counts'] != adata.X).nnz == 0
assert restored.obs_names.tolist() == adata.obs_names.tolist()
assert restored.var_names.tolist() == adata.var_names.tolist()
np.testing.assert_array_equal(preview.values, adata[:5, features].X.toarray())
assert preview.index.tolist() == adata.obs_names[:5].tolist()
assert preview.columns.tolist() == features
assert not (run / 'unexpected.tsv').exists()
"""))
    NotebookClient(notebook, timeout=300, kernel_name="scx-tutorial-test",
                   resources={"metadata": {"path": str(output / "quick-python")}}).execute()
    notebook.cells.pop()
    for cell in notebook.cells:
        cell.metadata.pop("execution", None)
    nbformat.write(notebook, output / "quickstart-executed.ipynb")
    rmd = (ROOT / "tutorials/quickstart.Rmd").read_text(encoding="utf-8")
    blocks = re.findall(r"```\{r[^\n]*\}\n(.*?)```", rmd, re.S)
    assert blocks
    script = "# PBMC3k quickstart. Notebook-style walkthrough: quickstart.md\n\n" + "\n\n".join(blocks)
    assert (ROOT / "tutorials/quickstart.R").read_text(encoding="utf-8") == script
    wrapper = output / "check-quick.R"
    lines = ['results <- list()', 'env <- new.env(parent=.GlobalEnv)']
    for i, block in enumerate(blocks, 1):
        lines.append(f'results[[{i}]] <- capture.output(eval(parse(text={json.dumps(block)}),envir=env))')
    lines.extend([
        'stopifnot(all(dim(env$restored)==c(32738,2700)))',
        'stopifnot(max(abs(SeuratObject::LayerData(env$restored,layer="counts")-env$counts))==0)',
        'stopifnot(identical(colnames(env$restored),colnames(env$counts)),identical(rownames(env$restored),rownames(env$counts)))',
        'stopifnot(max(abs(as.matrix(env$preview)-as.matrix(Matrix::t(env$counts[env$features,head(env$cells,5),drop=FALSE]))))==0)',
        'jsonlite::write_json(results,"../quick-r-output.json",auto_unbox=FALSE)',
    ])
    wrapper.write_text("\n".join(lines)+"\n", encoding="utf-8")
    result = subprocess.run([rscript, "--vanilla", str(wrapper)], cwd=output / "quick-r",
                            capture_output=True, text=True, encoding="utf-8", errors="replace")
    (output / "quick-r-console.txt").write_text(result.stdout + result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(result.stdout + result.stderr)
    results = iter(json.loads((output / "quick-r-output.json").read_text(encoding="utf-8")))
    def rendered(match):
        return "```r\n" + match.group(1) + "```\n\n```text\n" + "\n".join(next(results)) + "\n```"
    markdown = re.sub(r"```\{r[^\n]*\}\n(.*?)```", rendered, rmd, flags=re.S)
    markdown = re.sub(r"\A---\n.*?\n---\n", "# SCX quickstart with PBMC3k (R)\n", markdown, count=1, flags=re.S)
    (output / "quickstart-rendered.md").write_text(markdown, encoding="utf-8")
    p_run = next((output / "quick-python").glob("scx-pbmc3k-*"))
    r_run = next((output / "quick-r").glob("scx-pbmc3k-*"))
    p_table = pd.read_csv(p_run / "preview.tsv", sep="\t", index_col=0)
    r_table = pd.read_csv(r_run / "preview.tsv", sep="\t", index_col=0)
    pd.testing.assert_frame_equal(p_table, r_table, check_dtype=False, check_exact=True)
    for path, expected in (
        (output / "quick-python/pbmc3k_raw.h5ad", "89a96f1beaa2dd83a687666d3f19a4513ac27a2a2d12581fcd77afed7ea653a1"),
        (output / "quick-r/pbmc3k.tar.gz", "847d6ebd9a1ec9a768f2be7e40ca42cbfe75ebeb6d76a4c24167041699dc28b5"),
    ):
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected
    return {"passed": True, "dataset": "PBMC3k raw counts", "shape": [2700,32738],
            "preview_shape": list(p_table.shape), "full_counts_equal": True,
            "axis_order_equal": True, "cross_language_preview_equal": True,
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
