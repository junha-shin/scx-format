"""Fault injection verifies existing data survive failed publication."""
import errno
import importlib
import json
import subprocess
import sys

from anndata import read_h5ad
import numpy as np
import pytest

from scx.anndata import from_anndata, to_anndata
from scx.cli import main
from scx.h5ad import write_h5ad
from scx.publish import staged_output
from test_stability import sample


def snapshot(root):
    return {str(p.relative_to(root)): p.read_bytes() for p in root.rglob('*') if p.is_file()}


def assert_clean(parent):
    assert not list(parent.glob('.scx-*'))


def test_replace_dataset_and_refuse_default(tmp_path):
    root = tmp_path/'out.scx'
    original, changed = sample(), sample()
    from_anndata(original, root)
    before = snapshot(root)
    changed.layers['counts'][0, 0] = 42
    with pytest.raises(FileExistsError):
        from_anndata(changed, root)
    assert snapshot(root) == before
    from_anndata(changed, root, overwrite=True)
    np.testing.assert_array_equal(to_anndata(root).layers['counts'].toarray(), changed.layers['counts'].toarray())
    assert_clean(tmp_path)


def test_replacement_requires_boolean_opt_in(tmp_path):
    root = tmp_path/'out.scx'
    from_anndata(sample(), root)
    before = snapshot(root)
    with pytest.raises(ValueError, match='overwrite must be True or False'):
        from_anndata(sample(), root, overwrite='false')
    assert snapshot(root) == before
    assert_clean(tmp_path)


@pytest.mark.parametrize('existing', [False, True])
@pytest.mark.parametrize('failure', [OSError(errno.ENOSPC, 'disk full'), PermissionError('write denied'), KeyboardInterrupt()])
def test_interrupted_scx_write_preserves_previous(tmp_path, monkeypatch, existing, failure):
    root = tmp_path/'out.scx'
    if existing:
        from_anndata(sample(), root)
    before = snapshot(root) if existing else None
    adapter = importlib.import_module('scx.anndata')
    def fail(path, matrix):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b'partial matrix')
        raise failure
    monkeypatch.setattr(adapter, 'write_matrix', fail)
    with pytest.raises(type(failure)):
        from_anndata(sample(), root, overwrite=existing)
    assert snapshot(root) == before if existing else not root.exists()
    assert_clean(tmp_path)


@pytest.mark.parametrize('failure', [PermissionError('publication denied'), KeyboardInterrupt()])
def test_directory_publish_failure_rolls_back(tmp_path, monkeypatch, failure):
    root = tmp_path/'out.scx'
    from_anndata(sample(), root)
    before = snapshot(root)
    module = importlib.import_module('scx.publish')
    rename = module._rename_new
    def fail(source, target):
        if source.name == 'payload':
            raise failure
        return rename(source, target)
    monkeypatch.setattr(module, '_rename_new', fail)
    with pytest.raises(type(failure)):
        from_anndata(sample(), root, overwrite=True)
    assert snapshot(root) == before
    assert_clean(tmp_path)


def test_failed_rollback_keeps_backup_and_lock(tmp_path, monkeypatch):
    root = tmp_path/'out.scx'
    from_anndata(sample(), root)
    before = snapshot(root)
    module = importlib.import_module('scx.publish')
    rename = module._rename_new
    def fail(source, target):
        if source.name in {'payload', 'previous'}:
            raise PermissionError('destination unavailable')
        return rename(source, target)
    monkeypatch.setattr(module, '_rename_new', fail)
    with pytest.raises(RuntimeError, match='backup retained'):
        from_anndata(sample(), root, overwrite=True)
    stages = list(tmp_path.glob('.scx-stage-*'))
    assert len(stages) == 1
    assert snapshot(stages[0]/'previous') == before
    assert len(list(tmp_path.glob('.scx-lock-*'))) == 1
    with pytest.raises(ValueError, match='output is locked'):
        from_anndata(sample(), root)


def test_racing_directory_is_never_replaced(tmp_path, monkeypatch):
    root = tmp_path/'out.scx'
    module = importlib.import_module('scx.publish')
    rename = module._rename_new
    def raced(source, target):
        if source.name == 'payload':
            target.mkdir()  # Even an empty directory must survive a race.
        return rename(source, target)
    monkeypatch.setattr(module, '_rename_new', raced)
    with pytest.raises(FileExistsError):
        from_anndata(sample(), root)
    assert root.is_dir() and not list(root.iterdir())
    assert_clean(tmp_path)


def test_file_race_and_writer_lock(tmp_path):
    target = tmp_path/'out.h5ad'
    with pytest.raises(FileExistsError):
        with staged_output(target) as temporary:
            temporary.write_bytes(b'new')
            with pytest.raises(ValueError, match='output is locked'):
                with staged_output(target):
                    pytest.fail('a second writer acquired the same destination')
            target.write_bytes(b'other writer')
    assert target.read_bytes() == b'other writer'
    assert_clean(tmp_path)


def test_refuses_unrelated_directory(tmp_path):
    unrelated = tmp_path/'unrelated'
    unrelated.mkdir()
    (unrelated/'important.txt').write_text('keep me')
    with pytest.raises(ValueError, match='manifest.tsv'):
        from_anndata(sample(), unrelated, overwrite=True)
    assert (unrelated/'important.txt').read_text() == 'keep me'


@pytest.mark.parametrize('dangling', [False, True])
def test_refuses_output_symlinks(tmp_path, dangling):
    unrelated = tmp_path/'unrelated'
    if not dangling:
        unrelated.mkdir()
    link = tmp_path/'link'
    try:
        link.symlink_to(unrelated, target_is_directory=True)
    except OSError:
        pytest.skip('creating symlinks requires local Windows privileges')
    with pytest.raises(ValueError, match='symlinks are refused'):
        from_anndata(sample(), link, overwrite=True)


def test_h5ad_replace_and_failed_writer(tmp_path, monkeypatch):
    target = tmp_path/'out.h5ad'
    write_h5ad(sample(), target)
    before = target.read_bytes()
    with pytest.raises(FileExistsError):
        write_h5ad(sample(), target)
    assert target.read_bytes() == before
    changed = sample()
    changed.X[0, 0] = 24
    write_h5ad(changed, target, overwrite=True)
    assert read_h5ad(target).X[0, 0] == 24
    before = target.read_bytes()
    module = importlib.import_module('scx.h5ad')
    def fail(adata, path):
        path.write_bytes(b'partial hdf5')
        raise OSError(errno.ENOSPC, 'disk full')
    monkeypatch.setattr(module, '_write_h5ad', fail)
    with pytest.raises(OSError):
        write_h5ad(sample(), target, overwrite=True)
    assert target.read_bytes() == before
    assert_clean(tmp_path)


def test_h5ad_reopen_verification_preserves_previous(tmp_path, monkeypatch):
    target = tmp_path/'out.h5ad'
    write_h5ad(sample(), target)
    before = target.read_bytes()
    module = importlib.import_module('scx.h5ad')
    monkeypatch.setattr(module, '_write_h5ad', lambda a, p: p.write_bytes(b'not an hdf5 file'))
    with pytest.raises(OSError):
        write_h5ad(sample(), target, overwrite=True)
    assert target.read_bytes() == before
    assert_clean(tmp_path)


def test_cli_overwrite_and_input_dataset_protection(tmp_path, capsys):
    source, root, target = tmp_path/'in.h5ad', tmp_path/'out.scx', tmp_path/'out.h5ad'
    sample().write_h5ad(source)
    main(['from-h5ad', str(source), str(root)])
    main(['from-h5ad', str(source), str(root), '--overwrite'])
    main(['to-h5ad', str(root), str(target)])
    before = target.read_bytes()
    with pytest.raises(SystemExit) as exc:
        main(['to-h5ad', str(root), str(target)])
    assert exc.value.code == 1
    assert '--overwrite' in capsys.readouterr().err
    assert target.read_bytes() == before
    main(['to-h5ad', str(root), str(target), '--overwrite'])
    before = snapshot(root)
    with pytest.raises(SystemExit):
        main(['to-h5ad', str(root), str(root/'matrices/counts.mtx.gz'), '--overwrite'])
    assert 'outside the input SCX' in capsys.readouterr().err
    assert snapshot(root) == before


@pytest.mark.parametrize('during_commit', [False, True])
def test_hard_process_exit_leaves_recoverable_artifacts(tmp_path, during_commit):
    root = tmp_path/'out.scx'
    from_anndata(sample(), root)
    before = snapshot(root)
    # A hard exit bypasses finally blocks. Verify both pre-commit and the
    # directory-replacement gap instead of claiming automatic crash recovery.
    code = '''
import os, sys
from pathlib import Path
import scx.publish as publish
root = Path(sys.argv[1])
with publish.staged_output(root, directory=True, overwrite=True) as payload:
    (payload/'partial').write_text('interrupted')
    if sys.argv[2] == 'True':
        publish._rename_new(root, payload.parent/'previous')
    os._exit(9)
'''
    result = subprocess.run([sys.executable, '-c', code, str(root), str(during_commit)], capture_output=True, timeout=30)
    assert result.returncode == 9, result.stderr.decode()
    locks = list(tmp_path.glob('.scx-lock-*'))
    assert len(locks) == 1
    info = json.loads(locks[0].read_text())
    assert (tmp_path/info['stage']).exists()
    if during_commit:
        assert not root.exists()
        assert snapshot(tmp_path/info['backup']) == before
    else:
        assert snapshot(root) == before
    with pytest.raises(ValueError, match='output is locked'):
        from_anndata(sample(), root, overwrite=True)
