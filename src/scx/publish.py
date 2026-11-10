"""Stage and publish outputs on the destination filesystem, with rollback."""
from contextlib import contextmanager
import argparse
import ctypes
import errno
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import warnings

from .errors import message


def _exists(path):
    return os.path.lexists(path)


def _rename_new(source, target):
    """Rename a directory without replacing even an empty racing destination."""
    if os.name == 'nt':
        os.rename(source, target)  # Windows rename never replaces a destination.
        return
    libc = ctypes.CDLL(None, use_errno=True)
    if sys.platform.startswith('linux') and hasattr(libc, 'renameat2'):
        function = libc.renameat2
        function.argtypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint]
        result = function(-100, os.fsencode(source), -100, os.fsencode(target), 1)
    elif sys.platform == 'darwin' and hasattr(libc, 'renamex_np'):
        function = libc.renamex_np
        function.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_uint]
        result = function(os.fsencode(source), os.fsencode(target), 4)
    else:
        raise OSError(errno.ENOTSUP, 'filesystem/platform lacks safe no-replace directory rename', str(target))
    if result:
        code = ctypes.get_errno()
        raise OSError(code, os.strerror(code), str(target))


def _check_destination(path, directory, overwrite):
    if path.is_symlink():
        raise ValueError(f"output symlinks are refused: {path}")
    if _exists(path):
        if not overwrite:
            raise FileExistsError(errno.EEXIST, 'output already exists', str(path))
        if path.is_dir() != directory or (not directory and not path.is_file()):
            raise ValueError(f"output has the wrong type: {path}")
        if directory:
            # Never recursively replace an unrelated folder with --overwrite.
            from .core import read_manifest
            read_manifest(path)


@contextmanager
def staged_output(destination, *, directory=False, overwrite=False):
    """Yield a private payload path; publish only after the caller succeeds.

    Directory replacement has a brief name gap; it is rollback-safe, not a
    reader-atomic swap. Hard termination/power loss require manual recovery.
    """
    if type(overwrite) is not bool:
        raise ValueError('overwrite must be True or False')
    destination = Path(destination).absolute()
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Resolve the parent, never the leaf: output symlinks must be refused.
    destination = destination.parent.resolve() / destination.name
    key = hashlib.sha256(os.fsencode(os.path.normcase(str(destination)))).hexdigest()[:20]
    lock = destination.parent / f'.scx-lock-{key}'
    try:
        handle = lock.open('x', encoding='utf-8')
    except FileExistsError as exc:
        raise ValueError(f"output is locked: {destination}; inspect {lock} and follow docs/safe-io.md before recovery") from exc
    stage = None
    preserve = False
    published = False
    try:
        with handle:
            _check_destination(destination, directory, overwrite)
            stage = Path(tempfile.mkdtemp(prefix='.scx-stage-', dir=destination.parent))
            handle.write(json.dumps({'pid': os.getpid(), 'output': str(destination),
                                     'stage': str(stage), 'backup': str(stage/'previous')}))
            handle.flush()
            os.fsync(handle.fileno())
        payload, backup = stage/'payload', stage/'previous'
        if directory:
            payload.mkdir()
        yield payload
        _check_destination(destination, directory, overwrite)
        if not directory:
            if overwrite:
                os.replace(payload, destination)
            else:
                # Atomic create-if-absent, on the same filesystem as the stage.
                os.link(payload, destination)
        else:
            try:
                if _exists(destination):
                    _rename_new(destination, backup)
                _rename_new(payload, destination)
            except BaseException as exc:
                if backup.exists():
                    try:
                        _rename_new(backup, destination)
                    except OSError as recovery:
                        preserve = True
                        raise RuntimeError(f"could not restore previous output; backup retained at {backup}; "
                                           f"lock retained at {lock}. Follow docs/safe-io.md. "
                                           f"Publication failed: {exc}; restoration failed: {recovery}") from exc
                raise
        published = True
    finally:
        if stage is not None and (stage/'previous').exists() and not published:
            preserve = True
        if not preserve:
            try:
                if stage is not None:
                    shutil.rmtree(stage)
            except OSError as exc:
                warnings.warn(f"output cleanup incomplete at {stage}: {exc}; see docs/safe-io.md", UserWarning)
            finally:
                lock.unlink(missing_ok=True)


def main(argv=None):
    """Internal validated publication bridge for the R exporter."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage')
    parser.add_argument('output')
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args(argv)
    try:
        from .validate import validate
        source = Path(args.stage).resolve()
        target = Path(args.output).absolute()
        target = target.parent.resolve() / target.name
        if source == target or source in target.parents or target in source.parents:
            raise ValueError('stage and output must be separate sibling directories')
        if source.parent != target.parent.resolve():
            raise ValueError('stage and output must be on the same parent filesystem')
        with staged_output(target, directory=True, overwrite=args.overwrite) as payload:
            validate(source)
            payload.rmdir()
            _rename_new(source, payload)
    except (ValueError, OSError, RuntimeError) as exc:
        parser.exit(1, f"scx: {message(exc)}\n")
    except KeyboardInterrupt:
        parser.exit(130, 'scx: interrupted. See docs/safe-io.md for recovery.\n')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
