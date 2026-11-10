"""Cross-check R-native acceptance and shared locks against Python SCX."""
import argparse
import gzip
import json
from pathlib import Path
import shutil
import subprocess
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scx.anndata import from_anndata
from scx.core import read_table, write_table
from scx.publish import staged_output
from scx.validate import validate


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--rscript', default='Rscript')
args = parser.parse_args()
root = args.output.resolve()
if root.exists():
    parser.error('choose a fresh output directory')
root.mkdir(parents=True)
a = ad.AnnData(sparse.csr_matrix([[1., 0, 2], [0, 3, 0]]),
               obs=pd.DataFrame(index=['c1', 'c2']), var=pd.DataFrame(index=['g1', 'g2', 'g3']))
a.layers['counts'] = a.X.copy()
a.obs['number'] = pd.array([1, None], dtype='Int32')
a.obs['flag'] = pd.array([True, None], dtype='boolean')
a.obs['group'] = pd.Categorical(['a', None], categories=['a', 'unused'], ordered=True)
a.obsm['pca'] = np.array([[1., 2], [3, 4]])
a.obsp['graph'] = sparse.csr_matrix([[0., 1], [2, 0]])
from_anndata(a, root/'valid')
expected = {'valid': {'Python': True, 'R': True}}

def case(name, edit, py=False, r=False):
    output = root/name
    shutil.copytree(root/'valid', output)
    edit(output)
    expected[name] = {'Python': py, 'R': r}

def matrix(output, text):
    with gzip.open(output/'matrices/counts.mtx.gz', 'wt') as f:
        f.write(text)

def schema(output, edit):
    path = output/'metadata.schema.json'
    value = json.loads(path.read_text())
    edit(value)
    path.write_text(json.dumps(value))

def metadata(output, column, token):
    path = output/'cells.tsv.gz'
    table = read_table(path)
    table.loc[0, column] = token
    write_table(path, table)

case('truncated_gzip', lambda p: (p/'matrices/counts.mtx.gz').write_bytes((p/'matrices/counts.mtx.gz').read_bytes()[:-5]))
case('bad_crc', lambda p: (p/'cells.tsv.gz').write_bytes((p/'cells.tsv.gz').read_bytes()[:-8] + b'\0'*8))
case('missing', lambda p: (p/'features.tsv.gz').unlink())
case('short_mtx', lambda p: matrix(p, '%%MatrixMarket matrix coordinate real general\n2 3 2\n1 1 1\n'))
case('index', lambda p: matrix(p, '%%MatrixMarket matrix coordinate real general\n2 3 1\n3 1 1\n'))
case('shape', lambda p: matrix(p, '%%MatrixMarket matrix coordinate real general\n3 2 1\n1 1 1\n'))
case('nonfinite', lambda p: matrix(p, '%%MatrixMarket matrix coordinate real general\n2 3 1\n1 1 NaN\n'))
case('negative', lambda p: matrix(p, '%%MatrixMarket matrix coordinate real general\n2 3 1\n1 1 -1\n'))
case('symmetry', lambda p: matrix(p, '%%MatrixMarket matrix coordinate real symmetric\n2 3 1\n1 1 1\n'))
case('duplicate_coords', lambda p: matrix(p, '%%MatrixMarket matrix coordinate real general\n2 3 2\n1 1 1\n1 1 2\n'), True, True)
case('integer_cancellation_R_limit', lambda p: matrix(p, '%%MatrixMarket matrix coordinate integer general\n2 3 3\n1 1 9007199254740991\n1 1 2\n1 1 -9007199254740991\n'), True, False)
case('bool_number', lambda p: metadata(p, 'flag', '1'))
case('integer_overflow', lambda p: metadata(p, 'number', '2147483648'))
case('fraction_integer', lambda p: metadata(p, 'number', '1.5'))
case('category_code', lambda p: metadata(p, 'group', '99'))
case('fractional_category_code', lambda p: metadata(p, 'group', '0.0'))
def float32_case(p):
    schema(p, lambda s: s['cells']['columns'][0].update(dtype='float32'))
    metadata(p, 'number', '0.1')
case('lossy_float32', float32_case)
case('bad_json', lambda p: metadata(p, 'number', '{'))
case('bad_schema', lambda p: (p/'metadata.schema.json').write_text('[]'))
case('schema_null', lambda p: schema(p, lambda s: s.update(cells=None)))
case('manifest_width', lambda p: (p/'manifest.tsv').write_text((p/'manifest.tsv').read_text()+'bad\tembedding\tx\textra\n'))
case('unsafe_path', lambda p: (p/'manifest.tsv').write_text((p/'manifest.tsv').read_text()+'../outside\tembedding\tx\n'))
case('wide_integer_R_limit', lambda p: schema(p, lambda s: s['cells']['columns'][0].update(dtype='Int64')), True, False)
case('numeric_category_R_limit', lambda p: schema(p, lambda s: s['cells']['columns'][2].update(categories_dtype='int32', categories=[1, 2])), True, False)
case('legacy_text', lambda p: (p/'manifest.tsv').write_text(''.join(line for line in (p/'manifest.tsv').read_text().splitlines(keepends=True) if '\tmetadata_schema\t' not in line)), True, True)
for name, result in expected.items():
    try:
        validate(root/name)
        accepted = True
    except (ValueError, OSError):
        accepted = False
    assert accepted == result['Python'], (name, accepted)

script = str(Path('adapters/seurat/scripts/validation-cases.R').resolve())
(root/'expected.json').write_text(json.dumps(expected), encoding='utf-8')
def run(*command):
    subprocess.run([args.rscript, script, *map(str, command)], check=True)
run('cases', root)

# Both languages must recognize each other's actual on-disk lock identity.
target = root/'shared-output'
with staged_output(target, directory=True):
    run('locked', root/'valid', target)
target.rmdir()
ready, release = root/'ready', root/'release'
worker = subprocess.Popen([args.rscript, script, 'hold', str(target), str(ready), str(release)])
try:
    deadline = time.monotonic() + 30
    while not ready.exists():
        assert worker.poll() is None, 'R lock holder failed'
        assert time.monotonic() < deadline, 'R lock holder timed out'
        time.sleep(0.05)
    try:
        from_anndata(a, target)
    except ValueError as exc:
        assert 'locked' in str(exc)
    else:
        raise AssertionError('Python ignored R writer lock')
finally:
    release.touch()
    worker.wait(timeout=30)
assert worker.returncode == 0
report = {'passed': True, 'cases': expected, 'shared_locks_both_directions': True,
          'R_subset_limits': 'Wide integer metadata and numeric/boolean category levels remain rejected by native R.'}
(root/'report.json').write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
