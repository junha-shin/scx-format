"""Compare Assay/Assay5 fixtures through the active AnnData environment."""
import argparse
import hashlib
from importlib.metadata import version
import importlib.util
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--fixtures', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--rscript', default='Rscript')
parser.add_argument('--dataset', default='small-fixture')
args = parser.parse_args()
if args.output.exists():
    parser.error('output exists; choose a fresh directory')
spec = importlib.util.spec_from_file_location('seurat_example', 'examples/seurat/validate_roundtrip.py')
example = importlib.util.module_from_spec(spec)
spec.loader.exec_module(example)
args.output.mkdir(parents=True)
results = {}
for assay in ('Assay', 'Assay5'):
    output = args.output/assay
    shutil.copytree(args.fixtures/assay, output)
    _, checks = example.validate_export(output)
    subprocess.run([args.rscript, '--vanilla', 'examples/seurat/check_reverse.R', str(output), sys.executable], check=True)
    checks['reverse_assay_class_and_data'] = True
    checks['source_snapshot_sha256'] = hashlib.sha256((output/'source.rds').read_bytes()).hexdigest()
    results[assay] = checks
report = {'dataset': args.dataset, 'python': platform.python_version(), 'platform': platform.system(),
          'versions': {p: version(p) for p in ('scx-format', 'anndata', 'pandas', 'numpy', 'scipy', 'h5py', 'zarr')},
          'R': json.loads((args.fixtures/'R-versions.json').read_text()),
          'cases': results, 'passed': True}
(args.output/'report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
print(json.dumps(report, indent=2))
