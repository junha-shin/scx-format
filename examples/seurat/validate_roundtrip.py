"""Run Seurat PBMC3k analysis and compare independent R references with AnnData."""
import argparse
import hashlib
from importlib.metadata import distributions, version
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tarfile
import urllib.request

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from scx.anndata import from_anndata, to_anndata
from scx.validate import validate

URL = 'https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz'
SHA256 = '847d6ebd9a1ec9a768f2be7e40ca42cbfe75ebeb6d76a4c24167041699dc28b5'


def compare(left, right, tolerance=1e-12):
    left, right = sparse.csr_matrix(left), sparse.csr_matrix(right)
    assert left.shape == right.shape
    delta = left.astype(float) - right.astype(float)
    error = float(abs(delta).max()) if delta.nnz else 0.0
    assert error <= tolerance, error
    return {'shape': list(left.shape), 'max_absolute_error': error, 'passed': True}


def expected_metadata(columns, ids, h5ad=False):
    result = pd.DataFrame(index=pd.Index(ids))
    for column in columns:
        values = column['values']
        if column['levels'] is not None:
            result[column['name']] = pd.Categorical(values, categories=pd.Index(column['levels'], dtype=None if h5ad else object), ordered=column['ordered'])
        else:
            nullable = any(x is None for x in values)
            dtype = {'integer': 'Int32' if nullable else 'int32', 'double': 'float64',
                     'logical': 'boolean' if nullable else 'bool', 'character': 'string' if h5ad and nullable else 'object'}[column['r_type']]
            result[column['name']] = pd.Series(values, dtype=dtype).array
    return result


def validate_export(output):
    reference = json.loads((output/'reference.json').read_text(encoding='utf-8'))
    validate(output/'seurat.scx')
    in_memory = to_anndata(output/'seurat.scx')
    pd.testing.assert_frame_equal(in_memory.obs, expected_metadata(reference['obs'], reference['cells']), check_exact=True)
    pd.testing.assert_frame_equal(in_memory.var, expected_metadata(reference['var'], reference['features']), check_exact=True)
    target = output/'restored.h5ad'
    subprocess.run([sys.executable, '-m', 'scx.cli', 'to-h5ad', str(output/'seurat.scx'), str(target)], check=True)
    restored = ad.read_h5ad(target)
    assert restored.obs_names.tolist() == reference['cells']
    assert restored.var_names.tolist() == reference['features']
    checks = {'cell_order': True, 'feature_order': True}
    checks['counts'] = compare(mmread(output/'reference_counts.mtx', spmatrix=True).T, restored.layers['counts'], tolerance=0)
    checks['normalized'] = compare(mmread(output/'reference_normalized.mtx', spmatrix=True).T, restored.X)
    for embedding in reference['embeddings']:
        checks[embedding['name']] = compare(np.array(embedding['values']), restored.obsm[embedding['name']])
    for name in reference['graphs']:
        checks[name] = compare(mmread(output/f'reference_graph_{name}.mtx', spmatrix=True), restored.obsp[name])
    # AnnData/H5AD may infer the string storage dtype for category labels;
    # numeric widths, category values/order and nullability remain strict.
    pd.testing.assert_frame_equal(restored.obs, expected_metadata(reference['obs'], reference['cells'], h5ad=True), check_exact=True)
    pd.testing.assert_frame_equal(restored.var, expected_metadata(reference['var'], reference['features'], h5ad=True), check_exact=True)
    checks['metadata_values_dtypes_categories_missing'] = True
    # Also exercise the reverse path back into Seurat from AnnData-exported SCX.
    from_anndata(restored, output/'ann.scx')
    return reference, checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rscript', default='Rscript')
    parser.add_argument('--data', type=Path, default=Path('examples/seurat/data/pbmc3k.tar.gz'))
    parser.add_argument('--output', type=Path, default=Path('examples/seurat/output'))
    args = parser.parse_args()
    if args.output.exists():
        parser.error('choose a fresh output directory')
    args.data.parent.mkdir(parents=True, exist_ok=True)
    if not args.data.exists():
        urllib.request.urlretrieve(URL, args.data)
    digest = hashlib.sha256(args.data.read_bytes()).hexdigest()
    assert digest == SHA256, 'source checksum mismatch'
    raw = args.data.parent/'hg19'
    raw.mkdir(exist_ok=True)
    with tarfile.open(args.data) as archive:
        for name in ('matrix.mtx', 'genes.tsv', 'barcodes.tsv'):
            with archive.extractfile(f'filtered_gene_bc_matrices/hg19/{name}') as source, (raw/name).open('wb') as target:
                shutil.copyfileobj(source, target)
    subprocess.run([args.rscript, '--vanilla', 'examples/seurat/prepare.R', str(raw.resolve()), str(args.output.resolve())], check=True)
    reference, checks = validate_export(args.output)
    subprocess.run([args.rscript, '--vanilla', 'examples/seurat/check_reverse.R', str(args.output.resolve()), sys.executable], check=True)
    checks['AnnData_SCX_Seurat_reverse'] = True
    report = {'source_url': URL, 'source_sha256': digest, 'tested_shape': checks['counts']['shape'],
              'R': reference['R'], 'R_versions': reference['versions'], 'python': platform.python_version(),
              'python_versions': {p: version(p) for p in ('scx-format', 'anndata', 'numpy', 'pandas', 'scipy')},
              'seed': 42, 'checks': checks, 'passed': True,
              'limitations': ['single assay', 'R native dtype mapping', 'no identities, loadings or analysis history']}
    pins = sorted({f'{d.metadata["Name"]}=={d.version}' for d in distributions()
                   if d.metadata['Name'].lower() != 'scx-format'})
    (args.output/'python-requirements.txt').write_text('\n'.join(pins)+'\n', encoding='utf-8')
    (args.output/'report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
