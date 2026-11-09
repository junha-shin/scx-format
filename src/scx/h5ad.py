"""H5AD transport compatibility across supported AnnData versions."""
from contextlib import nullcontext
import warnings

import anndata as ad
import pandas as pd


def write_h5ad(adata, path):
    """Retain missing text with an explicit H5AD nullable-string mask."""
    changed = []
    for axis in ('obs', 'var'):
        frame = getattr(adata, axis)
        for name in frame:
            series = frame[name]
            if (series.dtype == object or isinstance(series.dtype, pd.StringDtype)) and series.isna().any():
                if not changed:
                    adata = adata.copy()
                getattr(adata, axis)[name] = series.astype(pd.StringDtype(storage='python'))
                changed.append(f'{axis}.{name}')
    if changed:
        warnings.warn('H5AD stores missing text as nullable string dtype: ' + ', '.join(changed), UserWarning)
    # AnnData 0.11/0.12 requires explicit opt-in for nullable string encoding.
    context = ad.settings.override(allow_write_nullable_strings=True) if hasattr(ad.settings, 'allow_write_nullable_strings') else nullcontext()
    with context:
        adata.write_h5ad(path, convert_strings_to_categoricals=False)
