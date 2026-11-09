"""Versioned, typed metadata encoded as JSON scalars inside TSV fields."""
import json
import re
import math

import numpy as np
import pandas as pd

from .core import read_table, write_table

NUMERIC = re.compile(r'^(?:u?int(?:8|16|32|64)|float(?:32|64)|U?Int(?:8|16|32|64)|Float(?:32|64))$')


def scalar(value):
    if pd.isna(value):
        return None
    return value.item() if isinstance(value, np.generic) else value


def descriptor(series):
    dtype = series.dtype
    if isinstance(dtype, pd.CategoricalDtype):
        return {'dtype': 'category', 'ordered': dtype.ordered,
                'categories_dtype': str(series.cat.categories.dtype),
                'categories': [scalar(x) for x in series.cat.categories]}
    name = str(dtype)
    if name not in {'object', 'str', 'string', 'bool', 'boolean'} and not NUMERIC.fullmatch(name):
        raise ValueError(f'unsupported metadata dtype: {name}')
    if name in {'object', 'str', 'string'} and any(not isinstance(x, str) for x in series.dropna()):
        raise ValueError('object metadata must contain only strings or missing values')
    return {'dtype': name}


def write_metadata(path, frame, id_column):
    if not frame.columns.is_unique or any(not isinstance(c, str) or not c for c in frame.columns):
        raise ValueError('metadata column names must be unique nonempty strings')
    table = pd.DataFrame({id_column: frame.index.astype(str)})
    columns = []
    for name in frame:
        if name == id_column:
            raise ValueError(f'reserved metadata column: {name}')
        series = frame[name]
        spec = {'name': name, **descriptor(series)}
        columns.append(spec)
        values = series.cat.codes if spec['dtype'] == 'category' else series
        encoded = [None if spec['dtype'] == 'category' and x == -1 else scalar(x) for x in values]
        try:
            table[name] = [json.dumps(x, ensure_ascii=False, allow_nan=False) for x in encoded]
        except (TypeError, ValueError) as exc:
            raise ValueError(f'unsupported metadata value in {name}') from exc
    write_table(path, table)
    return {'index_name': frame.index.name, 'columns': columns}


def decode(values, dtype):
    if dtype in {'object', 'str', 'string'}:
        if any(x is not None and not isinstance(x, str) for x in values):
            raise ValueError('string metadata contains non-string value')
    elif dtype in {'bool', 'boolean'}:
        if any(x is not None and type(x) is not bool for x in values):
            raise ValueError('boolean metadata must contain true/false/null')
    elif NUMERIC.fullmatch(dtype):
        integer = 'int' in dtype.lower()
        if any(x is not None and (type(x) not in ({int} if integer else {int, float})) for x in values):
            raise ValueError('invalid numeric metadata value')
        if any(type(x) is float and not math.isfinite(x) for x in values):
            raise ValueError('nonfinite numeric metadata value')
    else:
        raise ValueError(f'unsupported metadata dtype: {dtype}')
    if dtype in {'bool'} or dtype.startswith(('int', 'uint')):
        if None in values:
            raise ValueError('nonnullable dtype contains null')
    result = pd.Series(values, dtype=dtype)
    # Detect silent integer wraparound and lossy numeric casts.
    for old, new in zip(values, result):
        if old is not None and scalar(new) != old:
            raise ValueError('metadata value cannot be represented by declared dtype')
    return result


def read_metadata(path, spec, id_column):
    table = read_table(path)
    if not isinstance(spec, dict) or set(spec) != {'index_name', 'columns'}:
        raise ValueError('invalid metadata table schema')
    if spec['index_name'] is not None and not isinstance(spec['index_name'], str):
        raise ValueError('index name must be text or null')
    if not isinstance(spec['columns'], list):
        raise ValueError('schema columns must be a list')
    names = [c['name'] for c in spec['columns']]
    if list(table.columns) != [id_column, *names]:
        raise ValueError('metadata schema columns do not match TSV')
    result = pd.DataFrame(index=pd.Index(table[id_column].tolist(), name=spec['index_name']))
    for column in spec['columns']:
        name, dtype = column['name'], column['dtype']
        try:
            values = [json.loads(x, parse_constant=lambda x: (_ for _ in ()).throw(ValueError('nonfinite JSON'))) for x in table[name]]
            if dtype == 'category':
                if set(column) != {'name', 'dtype', 'categories', 'categories_dtype', 'ordered'} or type(column['ordered']) is not bool:
                    raise ValueError('invalid categorical schema')
                if not isinstance(column['categories'], list):
                    raise ValueError('categories must be a list')
                categories = decode(column['categories'], column['categories_dtype'])
                if categories.isna().any() or categories.duplicated().any():
                    raise ValueError('categories must be nonmissing and unique')
                if any(x is not None and (type(x) is not int or not 0 <= x < len(categories)) for x in values):
                    raise ValueError('categorical code out of range')
                series = pd.Series(pd.Categorical.from_codes([-1 if x is None else x for x in values], categories=categories, ordered=column['ordered']))
            else:
                if set(column) != {'name', 'dtype'}:
                    raise ValueError('invalid column schema')
                series = decode(values, dtype)
            result[name] = series.array
        except (ValueError, TypeError, OverflowError) as exc:
            raise ValueError(f'invalid typed metadata column: {name}: {exc}') from exc
    return result


def read_schema(path):
    try:
        schema = json.loads(path.read_text(encoding='utf-8'))
        if set(schema) != {'schema_version', 'cells', 'features'} or schema['schema_version'] != '1.0':
            raise ValueError('unsupported metadata schema version')
        return schema
    except (KeyError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError('invalid metadata schema') from exc
