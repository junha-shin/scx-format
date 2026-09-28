"""SCX command line interface."""
from __future__ import annotations

import argparse
from .validate import validate


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="scx")
    commands = parser.add_subparsers(dest="command", required=True)
    check = commands.add_parser("validate", help="validate an SCX directory")
    check.add_argument("directory")
    export = commands.add_parser("from-h5ad", help="convert AnnData to SCX")
    export.add_argument("input")
    export.add_argument("output")
    imported = commands.add_parser("to-h5ad", help="convert SCX to AnnData")
    imported.add_argument("input")
    imported.add_argument("output")
    args = parser.parse_args(argv)
    try:
        if args.command == "validate":
            result = validate(args.directory)
            print(f"Valid SCX: {result['cells']} cells, {result['features']} features")
        elif args.command == "from-h5ad":
            import anndata as ad
            from .anndata import from_anndata
            from_anndata(ad.read_h5ad(args.input), args.output)
        else:
            from .anndata import to_anndata
            to_anndata(args.input).write_h5ad(args.output)
    except (ValueError, OSError, ImportError) as exc:
        parser.exit(1, f"scx: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
