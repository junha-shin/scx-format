"""Check generated scripts stay current and exercise the standalone contract."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def test_standalone_bundles_are_current():
    subprocess.run([sys.executable, str(ROOT / "scripts/build_standalone.py"), "--check"], check=True)


def test_standalone_python_contract(tmp_path):
    subprocess.run([sys.executable, str(ROOT / "scripts/check_standalone.py"),
                    "--output", str(tmp_path / "standalone")], check=True)
