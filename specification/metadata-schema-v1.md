# SCX metadata schema 1.0

This optional extension to the versioned SCX v0.1-alpha manifest restores typed
cell/feature metadata. New reference writers MUST emit it. A manifest row MUST
be `metadata.schema.json`, role `metadata_schema`, name `metadata`. Older readers
that do not recognize this extension reject its role rather than misread data.
New readers accept legacy manifests without this row as plain text, with a
warning: legacy files cannot recover missing dtype information.

## Sidecar

UTF-8 JSON, with exactly these top-level keys:

```json
{
  "schema_version": "1.0",
  "cells": {
    "index_name": null,
    "columns": [
      {"name": "total_counts", "dtype": "float32"},
      {"name": "qc_pass", "dtype": "boolean"},
      {"name": "group", "dtype": "category", "ordered": true,
       "categories_dtype": "object", "categories": ["low", "high", "unused"]}
    ]
  },
  "features": {"index_name": null, "columns": []}
}
```

Each table descriptor MUST contain `index_name` (string or null) and `columns`
(ordered list). Column names MUST exactly match the TSV header after the ID
column, in order, with no duplicates. IDs remain literal, nonempty, unique text.
Index names are preserved, but index IDs remain strings as required by SCX.

## Types and values

All metadata TSV fields after the ID MUST contain one JSON scalar. A JSON string
is quoted inside the TSV field using normal TSV/CSV quoting. Decode the TSV
field first, then parse its JSON. This permits escaped tabs/newlines in string
values while the decoded TSV field itself contains no literal control separator.
No Python/R serialization or executable dtype expression is allowed.

| dtype | JSON scalar | Restored type |
|---|---|---|
| `int8/16/32/64`, `uint8/16/32/64` | integer, no null | fixed-width integer |
| `Int8/16/32/64`, `UInt8/16/32/64` | integer or null | nullable integer |
| `float32/64` | finite number or null | floating point; null becomes NaN |
| `Float32/64` | finite number or null | nullable floating point |
| `bool` | true/false, no null | boolean |
| `boolean` | true/false/null | nullable boolean |
| `object`, `str`, `string` | JSON string or null | text with declared logical dtype |
| `category` | zero-based integer category code or null | categorical |

For example, empty text is JSON `""`, literal text null is JSON `"null"`, and a
missing value is JSON `null`. These MUST remain distinct. Nulls preserve missing
positions, not Python-specific missing sentinels or NaN payload bits. String
storage backends such as Python vs Arrow are not part of the portable schema.

The `category` descriptor also MUST contain `ordered` (boolean),
`categories_dtype` (a supported noncategory dtype above), and `categories`
(JSON scalar list). Categories MUST be unique and nonmissing. Order and unused
categories MUST be retained even when `ordered` is false. Nonmissing codes MUST
be integers in range. Numeric/boolean categories preserve their value type.

Readers MUST reject unknown schema versions/dtypes, mismatched columns, invalid
JSON, wrong scalar types, integer overflow, lossy numeric casts, and invalid
category codes. Writers MUST reject unsupported mixed object columns, datetime,
arbitrary objects, complex values and infinities. Floating values are formatted
with sufficient precision to recover the declared numeric value.

## AnnData mapping and limits

The schema applies to `obs` and `var` only. Adapters restore values, dtype widths,
nullability, categorical categories/order and missing masks. The SCX CLI disables
H5AD's automatic string-to-categorical conversion. For cross-version H5AD
compatibility it writes missing text using nullable-string encoding, with a
warning: object/str storage may become pandas `string` without losing values or
missing masks. AnnData's writing option is scoped to that write. Direct external
H5AD writes may still use AnnData defaults; see the compatibility tables. `uns`, `varm`, other
unsupported objects and metadata column-index names are outside this extension.
