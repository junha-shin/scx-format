# Safe writes, replacement and recovery

## Standalone scripts

The Python file uses the Python safety code.
The R file writes new outputs only. It cannot replace or roll back a folder.
It reserves a new folder, checks copied files with SHA-256, then writes manifest.tsv last.
Readers reject the folder before that point.
A failed or interrupted copy leaves a partial folder and lock. Inspect them only after the writer stops.
R uses scratch disk to check single-member gzip CRC and length.
The native package supports more gzip forms and native no-replace folder rename.
See [safety differences](install.md#shared-default-and-safety).
The replacement steps below apply to Python and the optional native R package.


The current R package checks and publishes data without Python, using native helpers.
It shares Python's lock identity, backup layout, and no-replace behavior.
See [R-only use](r-only.md). The older Python requirement applied to R 0.1.3 in alpha.5.

The safety work started in Python 0.1.0a5 and scxSeurat 0.1.3.
The format and cells x features layout did not change.
It covers from_anndata, write_h5ad, the CLI, and write_scx_seurat.
Low-level write_table, write_matrix, and write_manifest do not stage output for callers.

## Normal use

Existing output is refused unless you request replacement:

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

The python argument above is kept for old calls and ignored by the current R package.
Only historical R 0.1.3 required Python with scx-format >=0.1.0a5 for validated reads and writes.
Current R use does not need a Python interpreter.

SCX output must be a folder; H5AD output must be a regular file.
Symlinks, including broken ones, and other file types are refused.
Folder replacement needs a readable 0.1-alpha manifest with valid required entries.
It cannot replace an unrelated folder.
CLI H5AD output must be outside the input SCX folder.
Overwrite replaces the whole output, including unlisted files. It does not merge datasets.
Keep personal files elsewhere.

## Publication and failures

1. A .scx-lock-<hash> file beside the output gives one SCX writer control of that path.
2. New data are written in a private sibling .scx-stage-<random>/payload folder.
   The old output stays in place while the new data are written and checked.
   SCX runs the full validator. H5AD is reopened with AnnData's backed reader
   to check readability, shape, and cell/feature IDs.
   This may reread metadata, layers, embeddings, and graphs.
   It is not a checksum or an independent check of every value.
3. New output is created only if its path is free. A racing target, even an
   empty folder, is kept. Existing H5AD uses same-filesystem atomic replacement.
   Existing SCX is moved to stage/previous; the checked payload takes its name.
4. If folder publication fails, SCX restores the old folder.
   If restore also fails, it keeps previous and the lock and reports their paths.
   On success or successful rollback, normal staging files and locks are removed.

Replacing an SCX folder leaves a brief gap when its final name is absent.
Pause readers or write to a new versioned path.
Locks coordinate SCX writers only. Other programs can still change files.
There is no concurrent-reader snapshot or network-filesystem guarantee.

Keep enough free space for the new data and old output.
Moving the old folder to previous does not copy it.
H5AD read-back adds I/O and may use more memory.
Historical benchmark results describe 0.1.0a4/0.1.2 and do not measure this safety step.

Normal Python exceptions and Ctrl-C clean up staging or roll back.
SIGKILL, os._exit, Task Manager termination, system failure, and power loss bypass cleanup.
Not all files and folders are fsynced. There is no durable power-loss transaction.
Stale locks are never removed automatically.

## Manual recovery after a hard interruption

1. Stop writers. Check that the lock's PID is no longer running this conversion;
   PIDs can be reused. Read output, stage, and backup from the lock JSON.
   They are absolute paths. An empty lock may mean staging had not started.
2. If output exists, check it with scx validate or reopen H5AD with AnnData.
   Keep previous until you have checked the chosen output.
   A payload may be incomplete; do not assume it is ready.
3. If output is absent and stage/previous exists, check the backup first.
   Move it back only if the output path is still absent.
   Do not overwrite a file made by another program.
   If neither exists, keep useful evidence and rerun from the source to a new path.
4. After checking the chosen output, remove only its abandoned stage and lock.
   Never bulk-delete .scx-* paths while a writer is active.

A hard R exit may also leave an outer .scx-<random> staging folder.
In the old Python-based path, this could happen before publication began.
It contains no previous output. Keep the source and inspect the stage before removing it.

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

Expected CLI failures exit 1 with the file, reason, and an action.
Ctrl-C exits 130; usage errors exit 2. Python APIs raise exceptions.
Gzip checks include the full stream and trailer.
Validation does not authenticate input or detect every valid-looking value change without hashes.

CI injects disk-full, permission, interruption, publish, rollback, and race failures.
It also tests hard exits and damaged gzip, UTF-8, TSV, MTX, manifest, and schema files.
Previous output must be kept. Linux and Windows exercise native write paths.
macOS has a renamex_np path; the earlier stabilization report did not verify it.
Unsupported systems fail instead of using an unsafe replacing rename.
File creation needs same-filesystem hard links.
Local tests do not certify NFS, SMB, cloud sync, or power-loss recovery.

## Alpha stabilization validation

The original stabilization check used Windows 11, Python 3.12.14, and AnnData 0.13.4.
It passed 84 Python cases; two symlink cases were skipped because the account lacked permission.
R 4.5.2 / SeuratObject 5.4.0 passed 45 assertions.
These are historical counts. The full Python suite also runs on Linux CI.

```sh
pip install -e '.[test]' -r compatibility/anndata-0.13.4.txt
pytest -q
Rscript adapters/seurat/scripts/run-tests.R
```

The checks reused saved real-data inputs. They did not rerun biological analysis.
They created, then replaced, SCX and H5AD outputs with that version's writers:

| Source | Cells × features | Replacement / round trip |
|---|---:|---|
| [Scanpy PBMC3k typed input](../examples/pbmc3k/README.md) | 2,700 × 13,714 | SCX and H5AD replacement passed; supported matrix errors 0 |
| [Seurat PBMC3k saved source object](../examples/seurat/README.md) | 2,700 × 13,714 | R export/replacement and SCX/H5AD replacement passed; supported matrix errors 0 |

The checks compared counts, normalized values, embeddings, graphs, axes, and metadata.
SCX kept exact metadata types.
H5AD may turn nullable text into string and infer different storage for text category labels.
Values, order, missing values, and unused categories remain intact.
See [type support](../compatibility/README.md#dtype-support).
