"""SCX command line interface."""
from __future__ import annotations

import argparse
from pathlib import Path
from .validate import validate
from .errors import message, resource


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="scx")
    commands = parser.add_subparsers(dest="command", required=True)
    check = commands.add_parser("validate", help="validate an SCX directory")
    check.add_argument("directory")
    export = commands.add_parser("from-h5ad", help="convert AnnData to SCX")
    export.add_argument("input")
    export.add_argument("output")
    export.add_argument("--overwrite", action="store_true", help="replace an existing SCX output after validation")
    export.add_argument("--analysis", action="store_true", help="also export embeddings and cell graphs")
    export.add_argument("--counts", default="counts", help="counts layer name, or X only when X contains raw counts")
    export.add_argument("--normalized", default="X", help="normalized layer name, X, or none")
    imported = commands.add_parser("to-h5ad", help="convert SCX to AnnData")
    imported.add_argument("input")
    imported.add_argument("output")
    imported.add_argument("--overwrite", action="store_true", help="replace an existing H5AD output after writing")
    args = parser.parse_args(argv)
    try:
        if args.command == "validate":
            result = validate(args.directory)
            print(f"Valid SCX: {result['cells']} cells, {result['features']} features")
        elif args.command == "from-h5ad":
            import anndata as ad
            from .api import write_scx
            with resource(args.input):
                adata = ad.read_h5ad(args.input)
            write_scx(adata, args.output, overwrite=args.overwrite, analysis=args.analysis,
                      counts=args.counts, normalized=None if args.normalized == "none" else args.normalized)
        else:
            from .anndata import to_anndata
            from .h5ad import write_h5ad
            output = Path(args.output)
            if output.resolve().is_relative_to(Path(args.input).resolve()):
                raise ValueError("H5AD output must be outside the input SCX directory")
            write_h5ad(to_anndata(args.input), output, overwrite=args.overwrite)
    except (ValueError, OSError, ImportError, EOFError, RuntimeError) as exc:
        parser.exit(1, f"scx: {message(exc)}\n")
    except KeyboardInterrupt:
        parser.exit(130, "scx: interrupted. Inspect the output and see docs/safe-io.md for recovery.\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
