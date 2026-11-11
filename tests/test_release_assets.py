"""Release staging must be reproducible, versioned, and reject damaged downloads."""
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('stage_release',ROOT/'scripts/stage_standalone_release.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_staged_release_is_reproducible_and_detects_tampering(tmp_path):
    a,b = tmp_path/'a',tmp_path/'b'
    module.stage(a); module.stage(b)
    assert module.verify(a)=={'scx.py','scx.R','USAGE.md'}
    for name in ('scx.py','scx.R','USAGE.md','SHA256SUMS'):
        assert (a/name).read_bytes()==(b/name).read_bytes()
        assert b'\r\n' not in (a/name).read_bytes()
    with pytest.raises(FileExistsError): module.stage(a)
    (a/'scx.R').write_bytes((a/'scx.R').read_bytes()+b'\n# altered\n')
    with pytest.raises(ValueError,match='Checksum mismatch'): module.verify(a)


@pytest.mark.parametrize('entry',['../scx.R','scx.R\n','/scx.R'])
def test_unsafe_checksum_entries_are_rejected(tmp_path,entry):
    (tmp_path/'SHA256SUMS').write_text('a'*64+'  '+entry+'\n')
    with pytest.raises(ValueError): module.verify(tmp_path)
