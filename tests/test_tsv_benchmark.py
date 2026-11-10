"""The benchmark must detect matching-but-wrong exports as well as misaligned axes."""
import importlib.util
from pathlib import Path

import pytest
from scipy import sparse

spec = importlib.util.spec_from_file_location('tsv_benchmark',Path(__file__).resolve().parents[1]/'examples/tsv-pbmc3k/benchmark.py')
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)


def test_reference_comparison_catches_matching_wrong_values(tmp_path):
    a,b = tmp_path/'a.tsv',tmp_path/'b.tsv'
    a.write_text('cell_id\tg1\tg2\nc1\t0\t2\n'); b.write_bytes(a.read_bytes())
    result = benchmark.compare_pair(a,b,['c1'],['g1','g2'],sparse.csr_matrix([[0,2]]),0)
    assert result['numeric_fields_checked']==2
    with pytest.raises(ValueError,match='numeric error'):
        benchmark.compare_pair(a,b,['c1'],['g1','g2'],sparse.csr_matrix([[1,2]]),0)


@pytest.mark.parametrize('text',['cell_id\tg2\tg1\nc1\t0\t2\n','cell_id\tg1\tg2\nc2\t0\t2\n',
    'cell_id\tg1\tg2\n','cell_id\tg1\tg2\nc1\t0\t2\nc2\t0\t0\n',
    'cell_id\tg1\tg2\nc1\t0\tnan\n'])
def test_reference_comparison_rejects_bad_structure(tmp_path,text):
    a,b = tmp_path/'a.tsv',tmp_path/'b.tsv'
    a.write_text(text); b.write_bytes(a.read_bytes())
    with pytest.raises(ValueError):
        benchmark.compare_pair(a,b,['c1'],['g1','g2'],sparse.csr_matrix([[0,2]]),0)
