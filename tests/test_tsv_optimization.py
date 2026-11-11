"""Before/after summaries must not compare different data or changed outputs."""
import copy
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('tsv_comparison',ROOT/'examples/tsv-pbmc3k/compare_optimization.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_same_input_comparison_and_ratio():
    before = json.loads((ROOT/'examples/tsv-pbmc3k/results.json').read_text(encoding='utf-8'))
    after = copy.deepcopy(before)
    for row in after['summary']: row['median_wall_seconds'] /= 2
    result = module.compare(before,after)
    assert result['passed'] and all(x['speedup']==2 for x in result['rows'])


@pytest.mark.parametrize('mutation,pattern',[
    ('input','input_files'),('runtime','runtime'),('values','numeric error'),('bytes','bytes changed')])
def test_incompatible_or_incorrect_runs_are_refused(mutation,pattern):
    before = json.loads((ROOT/'examples/tsv-pbmc3k/results.json').read_text(encoding='utf-8'))
    after = copy.deepcopy(before)
    if mutation=='input': next(iter(after['input_files'].values()))['sha256'] = 'wrong'
    elif mutation=='runtime':
        for row in after['records']:
            if row['engine']=='R': row['R'] = 'different runtime'
    elif mutation=='values': after['checks'][0]['max_absolute_errors']['R_vs_source'] = 1
    else:
        for row in after['records']:
            if row['engine']=='R' and row['case']=='counts-full-tsv': row['output_sha256'] = 'wrong'
    with pytest.raises(ValueError,match=pattern): module.compare(before,after)
