"""Reproduce a real-data SCX round trip; raw data and outputs stay local."""
from __future__ import annotations

import argparse
import hashlib
from importlib.metadata import distributions, version
import json
from pathlib import Path
import platform
import subprocess
import sys
import urllib.request

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

URL = 'https://falexwolf.de/data/pbmc3k_raw.h5ad'
SHA256 = '89a96f1beaa2dd83a687666d3f19a4513ac27a2a2d12581fcd77afed7ea653a1'


def compare_matrix(left, right):
    left, right = sparse.csr_matrix(left), sparse.csr_matrix(right)
    assert left.shape == right.shape
    delta = left.astype(np.float64) - right.astype(np.float64)
    error = float(abs(delta).max()) if delta.nnz else 0.0
    assert error <= 1e-6, f'max absolute error {error}'
    return {'shape': list(left.shape), 'nnz_before': left.nnz, 'nnz_after': right.nnz,
            'max_absolute_error': error, 'passed': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('examples/pbmc3k/output'))
    parser.add_argument('--data', type=Path, default=Path('examples/pbmc3k/data/pbmc3k_raw.h5ad'))
    args = parser.parse_args()
    if args.output.exists():
        parser.error('output exists; choose a new --output directory')
    args.data.parent.mkdir(parents=True, exist_ok=True)
    if not args.data.exists():
        urllib.request.urlretrieve(URL, args.data)
    digest = hashlib.sha256(args.data.read_bytes()).hexdigest()
    if digest != SHA256:
        raise ValueError('Source checksum mismatch; refusing changed input')
    args.output.mkdir(parents=True)
    requirements = sorted({f'{d.metadata["Name"]}=={d.version}' for d in distributions()
                           if d.metadata['Name'].lower() not in {'scx-format'}})
    (args.output/'requirements.txt').write_text('\n'.join(requirements) + '\n', encoding='utf-8')
    print('Preparing PBMC3k analysis...', flush=True)
    a = ad.read_h5ad(args.data)
    original_shape = list(a.shape)
    sc.pp.filter_genes(a, min_cells=3)
    a.layers['counts'] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=10000)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=1000, flavor='seurat')
    sc.pp.pca(a, n_comps=30, random_state=0)
    sc.pp.neighbors(a, n_neighbors=15, n_pcs=30, random_state=0)
    sc.tl.umap(a, random_state=0)
    print('Analysis complete; validating CLI round trip...', flush=True)
    # A derived QC category exercises common categorical metadata, not cell labels.
    totals = np.asarray(a.layers['counts'].sum(axis=1)).ravel()
    a.obs['library_size_group'] = pd.Categorical(np.where(totals >= np.median(totals), 'high', 'low'))
    a.obs['total_counts'] = totals
    a.obs['qc_pass'] = totals > 0
    a.obs['nullable_count'] = pd.array([int(x) if i % 10 else None for i, x in enumerate(totals)], dtype='Int64')
    a.obs['nullable_flag'] = pd.array([True if i % 10 else None for i in range(a.n_obs)], dtype='boolean')
    a.obs['library_size_group'] = a.obs['library_size_group'].cat.add_categories(['unused']).cat.as_ordered()
    a.obs.loc[a.obs_names[0], 'library_size_group'] = np.nan
    source, target, folder = args.output/'input.h5ad', args.output/'restored.h5ad', args.output/'dataset.scx'
    a.write_h5ad(source)
    # Actual command entry points, rather than only in-memory adapter calls.
    cli = [sys.executable, '-m', 'scx.cli']
    logs = []
    for command in [('from-h5ad', source, folder), ('validate', folder), ('to-h5ad', folder, target)]:
        result = subprocess.run(cli + [str(x) for x in command], capture_output=True, text=True)
        logs.append({'command': command[0], 'returncode': result.returncode,
                     'stdout': result.stdout, 'stderr': result.stderr.replace(str(Path(__file__).resolve().parents[2]), '<repo>')})
        if result.returncode:
            raise RuntimeError(result.stderr)
    before, after = ad.read_h5ad(source), ad.read_h5ad(target)
    checks = {}
    assert before.obs_names.equals(after.obs_names)
    assert before.var_names.equals(after.var_names)
    checks['cell_order'] = checks['feature_order'] = True
    checks['counts'] = compare_matrix(before.layers['counts'], after.layers['counts'])
    checks['normalized'] = compare_matrix(before.X, after.X)
    for key in before.obsm:
        checks[key] = compare_matrix(before.obsm[key], after.obsm[key])
    for key in before.obsp:
        checks[key] = compare_matrix(before.obsp[key], after.obsp[key])
    pd.testing.assert_frame_equal(before.obs, after.obs, check_exact=True)
    pd.testing.assert_frame_equal(before.var, after.var, check_exact=True)
    checks['metadata_values_dtypes_categories_and_missing'] = True
    report = {'source_url': URL, 'source_sha256': digest, 'original_shape': original_shape,
              'tested_shape': list(before.shape), 'seed': 0, 'tolerance': 1e-6,
              'python': platform.python_version(), 'platform': platform.platform(),
              'versions': {p: version(p) for p in ['scx-format', 'scanpy', 'anndata', 'numpy', 'pandas', 'scipy', 'scikit-learn', 'umap-learn', 'numba', 'h5py']},
              'checks': checks,
              'metadata_dtypes': {'before_obs': before.obs.dtypes.astype(str).to_dict(),
                                  'after_obs': after.obs.dtypes.astype(str).to_dict()},
              'known_losses': {'uns': list(before.uns), 'varm': list(before.varm)},
              'commands': logs, 'passed': True}
    (args.output/'report.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(checks, indent=2))


if __name__ == '__main__':
    main()
