"""SCX v0.1-alpha reference implementation."""

from .api import read_scx, write_scx, validate_scx

__version__ = "0.1.0a7"
__all__ = ["read_scx", "write_scx", "validate_scx", "to_tsv", "__version__"]


def to_tsv(path, output, *, matrix="counts", cells=None, features=None):
    """Export a selected expression matrix as protected CELL x FEATURE TSV[.gz]."""
    from .tsv import to_tsv as export
    return export(path, output, matrix=matrix, cells=cells, features=features)
