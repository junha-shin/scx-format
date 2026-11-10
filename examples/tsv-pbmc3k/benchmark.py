"""Reproduce PBMC3k TSV equivalence and isolated-process performance measurements."""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tarfile
import time
from importlib.metadata import version

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = Path(__file__).resolve()
URL = 'https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz'
ARCHIVE_SHA256 = '847d6ebd9a1ec9a768f2be7e40ca42cbfe75ebeb6d76a4c24167041699dc28b5'


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            value.update(chunk)
    return value.hexdigest()


def api():
    sys.path.insert(0, str(ROOT/'standalone'))
    import scx
    assert Path(scx.__file__).resolve() == ROOT/'standalone/scx.py'
    return scx


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False)+'\n', encoding='utf-8')


def prepare(args):
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse
    from scipy.io import mmread
    from inspect import signature
    if digest(args.archive) != ARCHIVE_SHA256:
        raise ValueError('PBMC3k archive SHA-256 mismatch')
    # Read only explicitly named members; never extract archive paths to disk.
    with tarfile.open(args.archive, 'r:gz') as archive:
        prefix = 'filtered_gene_bc_matrices/hg19/'
        options = {'spmatrix':True} if 'spmatrix' in signature(mmread).parameters else {}
        counts = sparse.csr_matrix(mmread(archive.extractfile(prefix+'matrix.mtx'), **options).T)
        genes = [line.decode('utf-8').rstrip('\n').split('\t') for line in archive.extractfile(prefix+'genes.tsv')]
        cells = [line.decode('utf-8').strip() for line in archive.extractfile(prefix+'barcodes.tsv')]
    keep = np.asarray((counts != 0).sum(axis=0)).ravel() >= 3
    counts = counts[:, keep].tocsr()
    genes = [gene for gene, selected in zip(genes, keep) if selected]
    totals = np.asarray(counts.sum(axis=1)).ravel()
    normalized = (sparse.diags(10000.0/totals) @ counts.astype(np.float64)).tocsr()
    normalized.data = np.log1p(normalized.data)
    data = ad.AnnData(normalized, obs=pd.DataFrame(index=cells),
        var=pd.DataFrame({'symbol': [g[1] for g in genes]}, index=[g[0] for g in genes]))
    data.layers['counts'] = counts
    api().write_scx(data, args.output)
    save(args.output.parent/(args.output.name+'-provenance.json'), {
        'dataset': '10x PBMC3k filtered hg19', 'url': URL, 'archive_sha256': ARCHIVE_SHA256,
        'original_shape': [len(cells), len(keep)], 'shape': list(data.shape), 'counts_nnz': counts.nnz,
        'feature_filter': 'expressed (nonzero counts) in at least 3 cells; Ensembl feature IDs',
        'normalized': 'log1p(counts / filtered-cell total * 10000), float64',
        'analysis_objects': 'none; expression and axis metadata only',
    })
    print(json.dumps({'prepared': str(args.output), 'shape': list(data.shape), 'nnz': counts.nnz}))


def worker(args):
    scx = api()
    # Warm code imports only, never input data, before the API timer.
    import importlib
    importlib.import_module('_scx_standalone.tsv')
    config = json.loads(args.config.read_text(encoding='utf-8'))
    start = time.perf_counter()
    scx.to_tsv(args.input, args.output, matrix=config['matrix'], cells=config['cells'], features=config['features'])
    from inspect import signature
    save(args.metrics, {'operation_seconds': time.perf_counter()-start,
        'gzip_compression_level':signature(gzip.open).parameters['compresslevel'].default,
        'python': platform.python_version(), 'versions': {p: version(p) for p in ('numpy','pandas','scipy')}})


def execute(command, log, env):
    import psutil
    peak = 0; tree_peak = 0
    start = time.perf_counter()
    with log.open('w', encoding='utf-8') as stream:
        process = subprocess.Popen(command, stdout=stream, stderr=stream, env=env)
        monitored = psutil.Process(process.pid)
        while process.poll() is None:
            try:
                members = [monitored,*monitored.children(recursive=True)]
                rss = 0
                for member in members:
                    try:
                        memory = member.memory_info()
                        rss += memory.rss
                        # venv Python and Rscript may be tiny launchers on Windows.
                        # The largest process's OS high-water mark captures the real
                        # worker, rather than incorrectly reporting only its launcher.
                        if sys.platform == 'win32': peak = max(peak,memory.peak_wset)
                    except (psutil.NoSuchProcess,psutil.AccessDenied):
                        pass
                tree_peak = max(tree_peak,rss)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
            time.sleep(0.01)
    duration = time.perf_counter()-start
    if process.returncode:
        raise RuntimeError(f'Worker failed; see {log}\n{log.read_text(encoding="utf-8", errors="replace")[-4000:]}')
    if sys.platform != 'win32': peak = tree_peak
    if peak <= 0:
        raise RuntimeError('Worker memory could not be measured')
    return duration, peak, tree_peak


def open_tsv(path):
    return (gzip.open if str(path).endswith('.gz') else open)(path, 'rt', encoding='utf-8', newline='')


def compare_pair(left, right, cell_ids, feature_ids, expected, tolerance):
    """Check every field/zero against the source, without retaining dense TSVs."""
    import numpy as np
    errors = {'python_vs_source': 0.0, 'R_vs_source': 0.0, 'python_vs_R': 0.0}
    with open_tsv(left) as a, open_tsv(right) as b:
        readers = [csv.reader(a, delimiter='\t'), csv.reader(b, delimiter='\t')]
        for reader in readers:
            if next(reader, None) != ['cell_id', *feature_ids]:
                raise ValueError('TSV header/feature order mismatch')
        for row, cell_id in enumerate(cell_ids):
            values = []
            for reader in readers:
                fields = next(reader, None)
                if fields is None or len(fields) != len(feature_ids)+1 or fields[0] != cell_id:
                    raise ValueError('TSV cell order/width/truncation mismatch')
                numbers = np.asarray(fields[1:], dtype=float)
                if not np.isfinite(numbers).all():
                    raise ValueError('TSV nonfinite values')
                values.append(numbers)
            reference = expected[row:row+1].toarray().ravel()
            for key, x, y in [('python_vs_source', values[0], reference),
                              ('R_vs_source', values[1], reference), ('python_vs_R', *values)]:
                error = float(np.max(np.abs(x-y)))
                errors[key] = max(errors[key], error)
                if error > tolerance:
                    raise ValueError(f'{key}: numeric error {error} exceeds {tolerance}')
        if any(next(reader, None) is not None for reader in readers):
            raise ValueError('Extra TSV rows')
    return {'passed': True, 'shape': [len(cell_ids),len(feature_ids)], 'tolerance': tolerance,
            'max_absolute_errors': errors, 'numeric_fields_checked': len(cell_ids)*len(feature_ids),
            'all_ids_and_zero_entries_checked': True}


def run(args):
    scx = api()
    import psutil
    import numpy as np
    args.output.mkdir(parents=True, exist_ok=False)
    data = scx.read_scx(args.input)  # Reference retained only in the unmeasured controller.
    cells, features = list(data.obs_names), list(data.var_names)
    rows = list(range(0,len(cells), max(1,len(cells)//100)))[:100][::-1]
    columns = list(range(0,len(features), max(1,len(features)//200)))[:200][::-1]
    cases = []
    for matrix in ('counts','normalized'):
        for gz in (False,True):
            cases.append({'name': matrix+('-full-gzip' if gz else '-full-tsv'), 'matrix': matrix,
                          'gzip': gz, 'cells': None, 'features': None})
        cases.append({'name': matrix+'-subset-gzip','matrix':matrix,'gzip':True,
                      'cells':[cells[i] for i in rows], 'features':[features[i] for i in columns]})
    env = os.environ.copy()
    for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
        env[name] = '1'
    records, checks = [], []
    for repeat in range(args.repeats):
        for case in cases:
            config = args.output/(case['name']+'.json')
            if not config.exists():
                save(config,case)
            outputs = {}
            # Reverse engine order on alternate repeats to reduce fixed order bias.
            for engine in (('python','R') if repeat%2 == 0 else ('R','python')):
                stem = f'{case["name"]}-{repeat}-{engine}'
                output = args.output/(stem+('.tsv.gz' if case['gzip'] else '.tsv'))
                metrics = args.output/(stem+'-metrics.json')
                if engine == 'python':
                    command = [sys.executable,str(SCRIPT),'worker','--input',str(args.input.resolve()),
                               '--output',str(output.resolve()),'--config',str(config.resolve()),'--metrics',str(metrics.resolve())]
                else:
                    command = [args.rscript,'--vanilla',str(SCRIPT.with_name('worker.R')),
                               str(ROOT/'standalone/scx.R'),str(args.input.resolve()),str(output.resolve()),str(config.resolve()),str(metrics.resolve())]
                wall, peak, tree_peak = execute(command,args.output/(stem+'.log'),env)
                record = {'case':case['name'],'repeat':repeat,'engine':engine,'wall_seconds':wall,
                          **json.loads(metrics.read_text(encoding='utf-8')), 'peak_bytes':peak,
                          'sampled_process_tree_peak_bytes':tree_peak,
                          'output_bytes':output.stat().st_size,'output_sha256':digest(output)}
                records.append(record); outputs[engine] = output
                print(json.dumps(record),flush=True)
            matrix = data.layers['counts'] if case['matrix']=='counts' else data.X
            if case['cells'] is not None:
                selected = matrix[rows,:][:,columns].tocsr()
            else:
                selected = matrix
            check = compare_pair(outputs['python'],outputs['R'],case['cells'] or cells,case['features'] or features,
                                 selected,0 if case['matrix']=='counts' else 1e-12)
            checks.append({'case':case['name'],'repeat':repeat,**check})
    from statistics import median
    summary = []
    for case in cases:
        for engine in ('python','R'):
            group = [r for r in records if r['case']==case['name'] and r['engine']==engine]
            summary.append({'case':case['name'],'engine':engine,
                'median_wall_seconds':median(r['wall_seconds'] for r in group),
                'median_operation_seconds':median(r['operation_seconds'] for r in group),
                'max_peak_bytes':max(r['peak_bytes'] for r in group),
                'output_bytes':group[0]['output_bytes']})
    provenance = args.input.parent/(args.input.name+'-provenance.json')
    report = {'passed':True,'dataset':json.loads(provenance.read_text()) if provenance.exists() else {'label':'synthetic smoke'},
        'shape':list(data.shape),'counts_nnz':data.layers['counts'].nnz, 'repeats':args.repeats,
        'platform':platform.platform(),'machine':platform.machine(),'CPU':platform.processor(),
        'logical_CPUs':psutil.cpu_count(),'physical_memory_bytes':psutil.virtual_memory().total,
        'memory_method':'largest worker/descendant process OS peak working set (Windows); sampled process-tree RSS also recorded' if sys.platform=='win32' else 'process-tree RSS sampled every 10 ms; observed lower bound',
        'controller_versions':{'Python':platform.python_version(), **{p:version(p) for p in ('psutil','anndata','numpy','pandas','scipy','h5py')}},
        'timing':'wall includes process startup, imports/source, validation, export, compression and publication; operation starts after code imports/source; excludes preparation/comparison',
        'thread_limits':{'BLAS_OpenMP_requested':1,'SciPy_Matrix_Market':'automatic'},
        'input_files':{str(p.relative_to(args.input)).replace('\\','/'):{'bytes':p.stat().st_size,'sha256':digest(p)} for p in sorted(args.input.rglob('*')) if p.is_file()},
        'implementation_sha256':{p:digest(ROOT/p) for p in ['standalone/scx.py','standalone/scx.R','examples/tsv-pbmc3k/benchmark.py','examples/tsv-pbmc3k/worker.R']},
        'cases':cases,'records':records,'summary':summary,'checks':checks}
    save(args.output/'report.json',report)
    lines = ['| Case | Engine | Wall s (median) | API s (median) | Peak MiB (max) | File MiB |',
             '|---|---|---:|---:|---:|---:|']
    for row in summary:
        lines.append(f'| {row["case"]} | {row["engine"]} | {row["median_wall_seconds"]:.3f} | {row["median_operation_seconds"]:.3f} | {row["max_peak_bytes"]/2**20:.1f} | {row["output_bytes"]/2**20:.3f} |')
    (args.output/'summary.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps({'passed':True,'comparisons':len(checks),'report':str(args.output/'report.json')}),flush=True)


def smoke(args):
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse
    args.output.mkdir(parents=True,exist_ok=False)
    data = ad.AnnData(sparse.csr_matrix(np.log1p([[0,2,3],[4,0,1],[1,5,0],[0,0,1]])),
                      obs=pd.DataFrame(index=['c1','c2','c3','c4']),var=pd.DataFrame(index=['g1','g2','g3']))
    data.layers['counts'] = sparse.csr_matrix([[0,2,3],[4,0,1],[1,5,0],[0,0,1]])
    source = args.output/'input.scx'
    api().write_scx(data,source)
    args.input = source; args.output = args.output/'run'; args.repeats = 1
    run(args)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest='mode',required=True)
    p = modes.add_parser('prepare')
    p.add_argument('--archive',type=Path,required=True); p.add_argument('--output',type=Path,required=True)
    for name in ('worker','run','smoke'):
        p = modes.add_parser(name); p.add_argument('--output',type=Path,required=True)
        if name != 'smoke': p.add_argument('--input',type=Path,required=True)
        if name == 'worker':
            p.add_argument('--config',type=Path,required=True); p.add_argument('--metrics',type=Path,required=True)
        else:
            p.add_argument('--rscript',default='Rscript')
            if name == 'run': p.add_argument('--repeats',type=int,default=2)
    args = parser.parse_args()
    if hasattr(args,'repeats') and args.repeats < 1: parser.error('--repeats must be positive')
    globals()[args.mode](args)


if __name__ == '__main__':
    main()
