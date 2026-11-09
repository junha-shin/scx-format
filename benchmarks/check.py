"""Check numeric values and typed metadata after an untimed benchmark run."""
import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def same_matrix(left, right):
    assert left.shape == right.shape
    delta = left.astype(np.float64) - right.astype(np.float64)
    if hasattr(delta, 'data') and hasattr(delta, 'nnz'):
        error = float(abs(delta).max()) if delta.nnz else 0.0
    else:
        error = float(np.max(np.abs(delta))) if delta.size else 0.0
    assert error <= 1e-12, error
    return error


def check(output, inputs):
    report = json.loads((output/'report.json').read_text())
    assert not report['profiled'], 'profiled runs are diagnostic, not timing evidence'
    results = []
    for record in report['records']:
        assert record['wall_seconds'] > 0 and record['peak_bytes'] > 0 and record['output_bytes'] > 0
        if record['command'] != 'to-h5ad':
            continue
        n, repeat = record['cells'], record['repeat']
        before = ad.read_h5ad(inputs/f"{n}-{report['genes']}-{report['draws_per_cell']}-{report['seed']}.h5ad")
        after = ad.read_h5ad(output/f'{n}-{repeat}.h5ad')
        pd.testing.assert_frame_equal(before.obs, after.obs, check_exact=True)
        pd.testing.assert_frame_equal(before.var, after.var, check_exact=True)
        errors = {'counts': same_matrix(before.layers['counts'], after.layers['counts']), 'normalized': same_matrix(before.X, after.X)}
        for key in before.obsm:
            errors[key] = same_matrix(before.obsm[key], after.obsm[key])
        for key in before.obsp:
            errors[key] = same_matrix(before.obsp[key], after.obsp[key])
        results.append({'cells': n, 'repeat': repeat, 'errors': errors, 'typed_metadata': True})
        del before, after
    assert len(results) == len(report['cells']) * report['repeats']
    result = {'passed': True, 'checks': results}
    (output/'verification.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--inputs', type=Path, default=Path('benchmarks/work/inputs'))
    args = parser.parse_args()
    print(json.dumps(check(args.output, args.inputs), indent=2))
