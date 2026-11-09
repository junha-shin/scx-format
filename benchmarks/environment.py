"""Record hardware and package pins without embedding local installation paths."""
import argparse
from importlib.metadata import distributions
import json
from pathlib import Path
import platform
import psutil

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
versions = {d.metadata['Name']: d.version for d in distributions()}
data = {'platform': platform.platform(), 'python': platform.python_version(), 'processor': platform.processor(),
        'physical_cores': psutil.cpu_count(logical=False), 'logical_cpus': psutil.cpu_count(),
        'physical_memory_bytes': psutil.virtual_memory().total, 'versions': dict(sorted(versions.items())),
        'filesystem': 'local workspace; disk type and OS cache state not characterized'}
(args.output/'environment.json').write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8')
pins = sorted(f'{name}=={value}' for name, value in versions.items() if name.lower() != 'scx-format')
(args.output/'python-requirements.txt').write_text('\n'.join(pins)+'\n', encoding='utf-8')
