"""SCX v0.1-alpha reference implementation."""

from .api import read_scx, write_scx, validate_scx

__version__ = "0.1.0a6"
__all__ = ["read_scx", "write_scx", "validate_scx", "__version__"]
