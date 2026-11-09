"""Compare the independent R fixture and actual H5AD; produce reverse SCX."""
import importlib.util
from pathlib import Path
import sys

spec = importlib.util.spec_from_file_location('seurat_example', 'examples/seurat/validate_roundtrip.py')
example = importlib.util.module_from_spec(spec)
spec.loader.exec_module(example)
_, checks = example.validate_export(Path(sys.argv[1]))
print(checks)
