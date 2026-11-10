"""File context and actionable messages for expected transport failures."""
from contextlib import contextmanager
import csv
import errno
import zlib


class SCXError(ValueError):
    """Invalid SCX input, including a damaged transport file."""


@contextmanager
def resource(path):
    try:
        yield
    except SCXError:
        raise
    except (ValueError, TypeError, OverflowError, EOFError, UnicodeError, csv.Error, OSError, zlib.error) as exc:
        if isinstance(exc, OSError) and not isinstance(exc, (FileNotFoundError,)) and exc.errno is not None:
            raise
        raise SCXError(f"{path}: {exc or type(exc).__name__}") from exc


def message(exc):
    if isinstance(exc, FileExistsError):
        hint = "Choose a new output path or pass --overwrite to replace an existing output."
    elif isinstance(exc, FileNotFoundError):
        hint = "Check the input path and restore any missing source files."
    elif isinstance(exc, PermissionError):
        hint = "Check output permissions and close programs holding the destination open."
    elif isinstance(exc, OSError) and exc.errno == errno.ENOSPC:
        hint = "Free disk space; replacement needs space for the new and previous outputs."
    elif isinstance(exc, ImportError):
        hint = "Install scx-format[anndata] in this Python environment."
    elif isinstance(exc, OSError) and exc.errno == errno.ENOTSUP:
        hint = "Use a supported local filesystem with safe rename and hard-link support."
    elif isinstance(exc, SCXError):
        hint = "Restore or regenerate the named input file and retry. For SCX directories, run scx validate."
    else:
        hint = "Check the input and output paths and the supported SCX profile."
    return f"{exc}\n{hint}"
