"""Validate installed packages through the public CLI and a tiny R/Python bridge.

No downloads, analysis pipeline or repository source imports are needed. This is
a transport smoke test; PBMC3k examples provide separate real-data evidence.
"""
import argparse
from importlib.metadata import version
import json
from pathlib import Path
import platform
import subprocess
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scx.h5ad import write_h5ad
from scx import read_scx, write_scx


def sample():
    counts = sparse.csr_matrix([[1., 0, 2, 0], [0, 3, 0, 1], [4, 0, 0, 2]])
    obs = pd.DataFrame(index=[f'cell{i}' for i in range(1, 4)])
    obs['batch'] = pd.Categorical(['b', 'a', None], categories=['b', 'a', 'unused'], ordered=True)
    obs['passing'] = pd.array([True, None, False], dtype='boolean')
    obs['number'] = pd.array([1, None, 3], dtype='Int32')
    obs['note'] = pd.array(['', None, 'null'], dtype='string')
    var = pd.DataFrame({'score': [0.5, np.nan, 1.5, 2.5]}, index=[f'gene{i}' for i in range(1, 5)])
    a = ad.AnnData(X=counts.log1p(), obs=obs, var=var)
    a.layers['counts'] = counts
    a.obsm['pca'] = np.array([[1., 2], [3, 4], [5, 6]])
    a.obsp['RNA_nn'] = sparse.csr_matrix([[0., 1, 0], [0, 0, 2], [3, 0, 0]])
    return a


def compare(expected, actual):
    assert expected.obs_names.equals(actual.obs_names)
    assert expected.var_names.equals(actual.var_names)
    errors = {}
    for name, left, right in [('counts', expected.layers['counts'], actual.layers['counts']),
                               ('normalized', expected.X, actual.X),
                               *[(f'embedding:{k}', expected.obsm[k], actual.obsm[k]) for k in expected.obsm],
                               *[(f'graph:{k}', expected.obsp[k], actual.obsp[k]) for k in expected.obsp]]:
        assert left.shape == right.shape, name
        delta = sparse.csr_matrix(left).astype(float) - sparse.csr_matrix(right).astype(float)
        error = float(abs(delta.data).max()) if delta.nnz else 0.0
        assert error <= 1e-12, (name, error)
        errors[name] = error
    assert set(actual.obsm) == set(expected.obsm) and set(actual.obsp) == set(expected.obsp)
    for axis in ('obs', 'var'):
        left, right = getattr(expected, axis).copy(), getattr(actual, axis)
        for name in left:
            series = left[name]
            # H5AD can infer string storage for textual category labels;
            # values, order, unused levels and missing codes must still match.
            if isinstance(series.dtype, pd.CategoricalDtype):
                left[name] = pd.Categorical.from_codes(series.cat.codes,
                    categories=pd.Index(series.cat.categories, dtype=right[name].cat.categories.dtype),
                    ordered=series.cat.ordered)
        pd.testing.assert_frame_equal(left, right, check_exact=True)
    return {'shape': list(actual.shape), 'max_absolute_errors': errors, 'metadata_and_axes': 'passed'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('quickstart-output'))
    parser.add_argument('--python-only', action='store_true', help='skip the two R/Python bridge cases')
    parser.add_argument('--rscript', default='Rscript')
    args = parser.parse_args()
    if args.output.exists():
        parser.error('output already exists; choose a fresh --output directory')
    output = args.output.resolve()
    output.mkdir(parents=True)
    log = []

    def run(command):
        result = subprocess.run(list(map(str, command)), capture_output=True, text=True)
        log.append({'command': list(map(str, command)), 'returncode': result.returncode,
                    'stdout': result.stdout, 'stderr': result.stderr})
        if result.returncode:
            raise RuntimeError(result.stderr or result.stdout)

    def scx(*command):
        run([sys.executable, '-m', 'scx.cli', *command])

    checks = {}
    try:
        source = output/'input.h5ad'
        write_h5ad(sample(), source)
        before = ad.read_h5ad(source)
        write_scx(before, output/'basic.scx')
        basic_expected = before.copy()
        basic_expected.obsm.clear(); basic_expected.obsp.clear()
        checks['basic_object_api'] = compare(basic_expected, read_scx(output/'basic.scx'))
        scx('from-h5ad', source, output/'python.scx', '--analysis')
        scx('validate', output/'python.scx')
        scx('to-h5ad', output/'python.scx', output/'python.h5ad')
        checks['python_h5ad_scx_h5ad'] = compare(before, ad.read_h5ad(output/'python.h5ad'))
        if not args.python_only:
            script = Path(__file__).with_name('seurat.R').resolve()
            run([args.rscript, str(script), 'convert', output/'python.scx', output/'via-r.scx', sys.executable])
            scx('to-h5ad', output/'via-r.scx', output/'via-r.h5ad')
            checks['anndata_seurat_anndata'] = compare(before, ad.read_h5ad(output/'via-r.h5ad'))
            run([args.rscript, str(script), 'create', output/'source.rds', output/'r.scx', sys.executable])
            scx('to-h5ad', output/'r.scx', output/'r.h5ad')
            checks['seurat_anndata'] = compare(before, ad.read_h5ad(output/'r.h5ad'))
            scx('from-h5ad', output/'r.h5ad', output/'r-back.scx', '--analysis')
            run([args.rscript, str(script), 'verify', output/'source.rds', output/'r-back.scx', sys.executable])
            checks['seurat_anndata_seurat'] = 'passed'
        report = {'passed': True, 'python_only': args.python_only, 'python': platform.python_version(),
                  'platform': platform.system(),
                  'versions': {name: version(name) for name in ('scx-format', 'anndata', 'pandas', 'numpy', 'scipy')},
                  'checks': checks}
        if not args.python_only:
            report['R'] = json.loads((output/'R-versions.json').read_text())
        (output/'report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
        print(json.dumps(report, indent=2))
    finally:
        (output/'commands.json').write_text(json.dumps(log, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
