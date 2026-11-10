"""Exercise the downloaded Python file with installed SCX imports forbidden."""
import argparse
import csv
import gzip
import importlib.abc
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import warnings


class NoInstalledSCX(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "scx" or fullname.startswith("scx."):
            raise AssertionError("Standalone script attempted an installed SCX import")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--require-uninstalled", action="store_true")
    args = parser.parse_args()
    if args.require_uninstalled:
        try:
            importlib.metadata.distribution("scx-format")
        except importlib.metadata.PackageNotFoundError:
            pass
        else:
            raise AssertionError("CI environment must not have scx-format installed")
    installed = sorted((d.metadata["Name"], d.version) for d in importlib.metadata.distributions())
    sys.meta_path.insert(0, NoInstalledSCX())
    script = Path(__file__).resolve().parents[1] / "standalone/scx.py"
    spec = importlib.util.spec_from_file_location("scx_script", script)
    scx = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scx)
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    a = ad.AnnData(sparse.csr_matrix([[0., 1.2], [2.3, 0.], [3.4, 1.]]),
        obs=pd.DataFrame(index=["c1", "c2", "c3"]), var=pd.DataFrame(index=["g1", "g2"]))
    a.layers["counts"] = sparse.csr_matrix([[0, 2], [3, 0], [4, 1]])
    a.obs["number"] = pd.array([1, None, 3], dtype="Int32")
    a.obs["flag"] = pd.array([True, False, None], dtype="boolean")
    a.obs["group"] = pd.Categorical(["b", None, "a"], categories=["b", "a", "unused"], ordered=True)
    a.var["score"] = [0.25, float("nan")]
    a.obsm["pca"] = np.array([[1., 2.], [3., 4.], [5., 6.]])
    a.obsp["nn"] = sparse.csr_matrix([[0, 2, 0], [1, 0, 3], [0, 0, 0]])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        scx.write_scx(a, out / "basic.scx")
    basic = scx.read_scx(out / "basic.scx")
    assert not basic.obsm and not basic.obsp
    pd.testing.assert_frame_equal(a.obs, basic.obs)
    pd.testing.assert_frame_equal(a.var, basic.var)
    np.testing.assert_array_equal(a.layers["counts"].toarray(), basic.layers["counts"].toarray())
    scx.write_scx(a, out / "full.scx", analysis=True)
    full = scx.read_scx(out / "full.scx")
    np.testing.assert_array_equal(a.obsm["pca"], full.obsm["pca"])
    np.testing.assert_array_equal(a.obsp["nn"].toarray(), full.obsp["nn"].toarray())
    scx.to_tsv(out / "full.scx", out / "counts.tsv")
    scx.to_tsv(out / "full.scx", out / "normalized.tsv.gz", matrix="normalized",
               cells=['c3','c1'], features=['g2','g1'])
    with gzip.open(out / "normalized.tsv.gz", 'rt', encoding='utf-8') as handle:
        table = list(csv.reader(handle, delimiter='\t'))
    assert table[0] == ['cell_id','g2','g1'] and [r[0] for r in table[1:]] == ['c3','c1']
    np.testing.assert_array_equal(np.array([r[1:] for r in table[1:]], dtype=float), [[1.,3.4],[1.2,0.]])
    try:
        scx.to_tsv(out / "full.scx", out / "counts.tsv")
    except FileExistsError:
        pass
    else:
        raise AssertionError("TSV was overwritten")
    before = (out / "full.scx/matrices/counts.mtx.gz").read_bytes()
    try:
        scx.write_scx(a, out / "full.scx")
    except FileExistsError:
        pass
    else:
        raise AssertionError("Existing output was overwritten")
    bad = a.copy(); bad.layers["counts"].data[0] = -1
    try:
        scx.write_scx(bad, out / "full.scx", overwrite=True)
    except ValueError:
        pass
    else:
        raise AssertionError("Negative counts accepted")
    assert (out / "full.scx/matrices/counts.mtx.gz").read_bytes() == before
    a.write_h5ad(out / "input.h5ad")
    def cli(*argv):
        subprocess.run([sys.executable, str(script), *map(str, argv)], check=True, cwd=out)
    cli("from-h5ad", out / "input.h5ad", out / "cli.scx", "--analysis")
    cli("validate", out / "cli.scx")
    cli("to-h5ad", out / "cli.scx", out / "restored.h5ad")
    pd.testing.assert_frame_equal(a.obs, ad.read_h5ad(out / "restored.h5ad").obs)
    path = out / "full.scx/matrices/counts.mtx.gz"
    path.write_bytes(before[:-8])
    try:
        scx.read_scx(out / "full.scx")
    except (ValueError, OSError, EOFError):
        pass
    else:
        raise AssertionError("Truncated gzip accepted")
    assert "scx" not in sys.modules
    assert installed == sorted((d.metadata["Name"], d.version) for d in importlib.metadata.distributions())
    report = {"passed": True, "installed_scx_imports": 0, "package_installations": 0,
              "shape": [3, 2], "typed_metadata": True, "analysis": True,
              "overwrite_protection": True, "corrupt_gzip_rejected": True, "CLI": True, "TSV": True}
    (out / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
