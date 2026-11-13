"""Reject non-ASCII README text and invalid UTF-8 in tracked README files."""
from pathlib import Path
import subprocess

root = Path(__file__).resolve().parents[1]
result = subprocess.run(["git", "-c", f"safe.directory={root.as_posix()}", "ls-files", "-z"],
                        cwd=root, check=True, capture_output=True)
files = [root / p.decode("utf-8") for p in result.stdout.split(b"\0") if p and Path(p.decode("utf-8")).name == "README.md"]
for path in files:
    text = path.read_text(encoding="utf-8")
    if not text.isascii():
        raise SystemExit(f"{path.relative_to(root)}: README must contain ASCII text only; see CONTRIBUTING.md")
print(f"README encoding check passed ({len(files)} files)")
