"""Small object API, with optional analysis coordinates/graphs."""
import warnings


def write_scx(object, path, *, counts="counts", normalized="X", analysis=False, overwrite=False):
    """Write AnnData safely. Basic = expression + typed obs/var.

    counts selects layers[counts], or explicitly X. normalized selects X, a
    layer, or None. X is never inferred to be raw counts. analysis=True adds
    supported obsm/obsp. The old from_anndata adapter retains its full defaults.
    """
    if type(analysis) is not bool:
        raise ValueError("analysis must be True or False")
    from .anndata import from_anndata
    omitted = []
    if not analysis and (len(object.obsm) or len(object.obsp)):
        omitted.append("embeddings/graphs (include with analysis=True)")
    if (set(object.layers) - {counts, normalized, None}) or len(object.uns) or len(object.varm) or len(object.varp) or object.raw is not None:
        omitted.append("extra layers/analysis internals outside the SCX profile")
    from_anndata(object, path, counts=counts,
                 normalized=None if counts == normalized else normalized,
                 analysis=analysis, overwrite=overwrite, warn_scope=False)
    if omitted:
        warnings.warn("SCX omitted: " + "; ".join(omitted), UserWarning, stacklevel=2)


def read_scx(path):
    """Validate and read every supported resource present in an SCX directory."""
    from .anndata import to_anndata
    return to_anndata(path)


def validate_scx(path):
    """Optional diagnostic; read_scx/write_scx already validate automatically."""
    from .validate import validate
    return validate(path)
