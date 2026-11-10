# Safe writes, replacement and recovery

## Standalone scripts

The no-install Python script uses the reference Python safety implementation.
The R source file supports **new outputs only**, not replacement/rollback.
It reserves a new directory with `dir.create`, verifies each copied staged file
with SHA-256, and publishes `manifest.tsv` last. Before that point readers reject
the incomplete directory. Copy failure/interruption retains the incomplete
folder and shared lock; inspect after the writer stops and use a fresh path.
R checks single-member gzip CRC and length using scratch disk. These are basic
safety guarantees, not the native package's no-replace directory rename or its
full gzip profile. See [the standalone safety table](install.md#shared-default-and-safety).
The replacement/native-helper details below describe the installed packages.


Current main's **R adapter 0.1.4** now implements validation and publication
without Python, using compiled R helpers. See [R-only operation](r-only.md) for
its installation, validation subset and tests. It shares Python's lock identity,
backup layout and no-replace publication contract. The version-specific Python
requirement below describes the original R 0.1.3 / alpha.5 release.

Available in Python 0.1.0a5 and scxSeurat 0.1.3. The on-disk format and CELL ×
FEATURE orientation are unchanged. This applies to `from_anndata`, `write_h5ad`,
the conversion CLI and `write_scx_seurat`. Core `write_table`, `write_matrix`
and `write_manifest` are low-level building blocks; callers must stage them.

## Normal use

Existing outputs are refused unless replacement is explicit:

```sh
scx from-h5ad input.h5ad result.scx --overwrite
scx to-h5ad result.scx result.h5ad --overwrite
```

```python
from scx.anndata import from_anndata
from scx.h5ad import write_h5ad
from_anndata(adata, "result.scx", overwrite=True)
write_h5ad(adata, "result.h5ad", overwrite=True)
```

```r
write_scx_seurat(object, "result.scx", overwrite = TRUE,
                 python = "/path/to/python")
```

Seurat export now requires Python with **scx-format >= 0.1.0a5** installed,
including when no overwrite is requested. This enables full Python validation
and the shared publication code. Import retains its existing Python requirement.
Pass the interpreter explicitly if it is not on PATH.

The destination type must match: SCX is a directory, H5AD a regular file.
Output symlinks, including dangling links, are refused. Directory replacement
requires a readable v0.1-alpha manifest with valid required entries; unrelated
folders are not recursively replaced. Other file types are refused. H5AD output
from the CLI must be outside its input SCX directory. `--overwrite` replaces the
whole output, including unlisted files in an SCX directory; save personal files
elsewhere. It does not merge datasets.

## Publication and failures

1. An exclusive `.scx-lock-<hash>` next to the output prevents cooperating SCX
   writers from publishing to the same resolved destination concurrently.
2. Files are written under a private sibling `.scx-stage-<random>/payload` on
   the destination filesystem. The existing output stays in place during writing
   and validation. SCX runs the complete validator; H5AD uses AnnData's backed
   reader to check readability, shape and cell/feature IDs after writing.
   This H5AD check may read metadata, layers, embeddings and graphs again; it
   is not a byte checksum or a second independent comparison of all values.
3. New outputs use create-if-absent operations, including a no-replace directory
   rename. A racing destination, even an empty directory, is preserved. Existing
   H5AD files use a same-filesystem atomic replacement. Existing SCX directories
   are moved to `stage/previous`, then the validated payload takes their name.
4. If directory publication raises an exception, the old directory is restored.
   If restoration also fails, `previous` and the lock are retained and the error
   gives their paths. Successfully published outputs and restored old outputs
   remain intact; normal staging files and locks are removed.

An SCX directory replacement has a brief interval when its final name is absent.
Pause readers during replacement, or use a new versioned destination. Locks only
coordinate SCX writers; unrelated programs can still modify paths. There is no
general concurrent-reader snapshot or distributed/network-filesystem guarantee.

Staging requires free space for the new output while keeping the old one. Moving
the old directory to `previous` does not copy it. Extra H5AD verification adds
I/O and can increase memory. Historical benchmarks under `benchmarks/results`
describe 0.1.0a4/0.1.2; they are not measurements of this additional safety gate.

Regular Python exceptions and Ctrl-C clean staging or roll back. SIGKILL,
`os._exit`, ending the process in Task Manager, host failure and power loss bypass
cleanup. SCX does not promise durable transactions after power loss: file/directory
data are not all fsynced, and stale locks are never removed automatically.

## Manual recovery after a hard interruption

1. Stop writers and confirm that the PID in the lock is no longer running this
   conversion (PIDs can be reused). Read its JSON: `output`, `stage` and `backup`
   are absolute paths. An empty lock can mean interruption before staging began.
2. If the output exists, validate it (`scx validate` for a directory; reopen an
   H5AD with AnnData). Keep any `previous` backup until the chosen output is
   confirmed. A `payload` may be incomplete; never assume it is a finished output.
3. If output is absent and `stage/previous` exists, validate that backup and move
   it back to the exact output path only when that path is still absent. Do not
   overwrite an independently created output. If neither exists, preserve any
   artifacts needed for diagnosis and rerun from the original source to a new path.
4. Once the selected output is verified, remove only the identified abandoned
   stage and its lock. Never bulk-delete `.scx-*` paths while a writer is active.

R may also leave its outer `.scx-<random>` staging directory after a hard exit
before Python publication begins. It contains no previous output; preserve the
original and inspect/validate the stage before deciding to remove it.

## Errors and checks

| Error | Meaning and action |
|---|---|
| output already exists | Use a new path or explicitly opt in to replacement. |
| output is locked | Another writer or abandoned run owns this destination; follow recovery above. |
| missing file / invalid manifest | Restore the named file or regenerate SCX; do not edit axes independently. |
| gzip EOF, CRC or decompression failure | Compressed transport is damaged; restore/regenerate the named file. |
| truncated MTX / invalid index / shape mismatch | Regenerate the matrix with its canonical cell and feature tables. |
| invalid typed metadata | Check the named TSV column and `metadata.schema.json` together. |
| permission denied | Check permissions and close software holding destination files open. |
| no space left | Free space for staging; the previous output is retained. |
| rollback failed | Keep `previous` and the lock; recover manually using the error's paths. |

Expected CLI failures exit 1 with a file/reason and an action; Ctrl-C exits 130.
Parser/usage errors exit 2. Python APIs continue to raise exceptions. The validator
checks the entire gzip stream including its trailer. It does not authenticate
inputs or detect every syntactically valid numerical alteration without checksums.

CI covers fault-injected disk-full, permission, interruption, publish/rollback
failure, writer races, actual subprocess hard exit, corrupt gzip/UTF-8/TSV/MTX/
manifest/schema, and preservation of previous files. Linux and Windows exercise
their native publication paths. Directory rename uses Linux `renameat2` or Windows
rename; macOS has an unverified `renamex_np` implementation. Unsupported systems
or filesystems fail rather than falling back to a replacing rename. File creation
requires same-filesystem hard-link support. Local-filesystem tests do not certify
NFS, SMB, cloud sync folders or power-loss recovery.

## Alpha stabilization validation

Local Windows 11 / Python 3.12.14 / AnnData 0.13.4 verification for this change:
84 Python cases passed; two symlink cases skipped because the account cannot
create Windows symlinks. R 4.5.2 / SeuratObject 5.4.0 passed all 45 assertions.
The full Python suite, including the symlink cases, also runs on Linux CI.

```sh
pip install -e '.[test]' -r compatibility/anndata-0.13.4.txt
pytest -q
Rscript adapters/seurat/scripts/run-tests.R
```

Previously generated real-data inputs from the linked examples were reused;
biological analysis was not rerun. New SCX and H5AD outputs were created and then
explicitly replaced using this version's writers:

| Source | Cells × features | Replacement / round trip |
|---|---:|---|
| [Scanpy PBMC3k typed input](../examples/pbmc3k/README.md) | 2,700 × 13,714 | SCX and H5AD replacement passed; supported matrix errors 0 |
| [Seurat PBMC3k saved source object](../examples/seurat/README.md) | 2,700 × 13,714 | R export/replacement and SCX/H5AD replacement passed; supported matrix errors 0 |

Counts, normalized values, every embedding and cell graph, axis order and
metadata were compared. SCX metadata matched exact dtypes. H5AD retains its
documented transport behavior: nullable text becomes `string`, and AnnData may
infer `str` storage for textual category labels while preserving values, order
and unused categories. See [dtype support](../compatibility/README.md#dtype-support).
