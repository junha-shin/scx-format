"""Prepare canonical LF standalone release files and SHA-256 checksums."""
import argparse
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def verify(folder):
    entries = (folder/'SHA256SUMS').read_text(encoding='utf-8').splitlines()
    seen = set()
    parsed = []
    if not entries:
        raise ValueError('Empty checksum list')
    for line in entries:
        match = re.fullmatch(r'([0-9a-f]{64})  ([A-Za-z0-9_.-]+)',line)
        if not match or match[2] in seen:
            raise ValueError('Invalid or duplicate checksum entry')
        seen.add(match[2])
        parsed.append((match[1],match[2]))
    for checksum,name in parsed:
        if hashlib.sha256((folder/name).read_bytes()).hexdigest() != checksum:
            raise ValueError('Checksum mismatch: '+name)
    return seen


def stage(folder):
    release = json.loads((ROOT/'release.json').read_text(encoding='utf-8'))
    if f'__version__ = "{release["python_version"]}"' not in (ROOT/'src/scx/__init__.py').read_text():
        raise ValueError('Python package/release version mismatch')
    if f'version = "{release["python_version"]}"' not in (ROOT/'pyproject.toml').read_text():
        raise ValueError('Python build/release version mismatch')
    if f'Version: {release["R_version"]}\n' not in (ROOT/'adapters/seurat/DESCRIPTION').read_text():
        raise ValueError('R package/release version mismatch')
    folder.mkdir(parents=True,exist_ok=False)
    files = {name:(ROOT/'standalone'/name).read_text(encoding='utf-8') for name in ('scx.py','scx.R')}
    files['USAGE.md'] = (ROOT/'docs/standalone-usage.md').read_text(encoding='utf-8')
    for name, content in files.items():
        if release['tag'] not in content:
            raise ValueError('Missing release identifier: '+name)
        (folder/name).write_bytes(content.encode('utf-8'))
    lines = [hashlib.sha256((folder/name).read_bytes()).hexdigest()+'  '+name for name in sorted(files)]
    (folder/'SHA256SUMS').write_bytes(('\n'.join(lines)+'\n').encode('ascii'))
    verify(folder)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder',type=Path)
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    verify(args.folder) if args.verify else stage(args.folder)
