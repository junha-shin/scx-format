"""Seeded sparse workload; isolated CLI workers with OS peak-memory metrics."""
import argparse
import json
import hashlib
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
from importlib.metadata import version


def peak_bytes():
    if sys.platform == 'win32':
        import psutil
        return int(psutil.Process().memory_info().peak_wset)
    import resource
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(peak if sys.platform == 'darwin' else peak * 1024)


def worker(args):
    sys.path.insert(0, str(args.implementation.resolve()))
    from scx.cli import main
    import scx
    from threadpoolctl import threadpool_limits
    threadpool_limits(limits=1)
    assert Path(scx.__file__).resolve().is_relative_to(args.implementation.resolve())
    if args.profile:
        import cProfile
        profiler = cProfile.Profile()
        profiler.enable()
    start = time.perf_counter()
    main(args.command)
    duration = time.perf_counter() - start
    peak = peak_bytes()
    if args.profile:
        profiler.disable()
        profiler.dump_stats(str(args.profile))
    args.metrics.write_text(json.dumps({'operation_seconds': duration, 'peak_bytes': peak, 'engine_version': scx.__version__}), encoding='utf-8')


def make_input(path, cells, genes, per_cell, seed):
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse
    rng = np.random.default_rng(seed)
    # Independent random positions; duplicates are summed. No dense cells x genes allocation.
    indptr = np.arange(cells + 1, dtype=np.int64) * per_cell
    indices = rng.integers(0, genes, cells * per_cell, dtype=np.int32)
    values = rng.integers(1, 20, cells * per_cell).astype(np.float32)
    counts = sparse.csr_matrix((values, indices, indptr), shape=(cells, genes))
    counts.sum_duplicates(); counts.sort_indices()
    normalized = counts.copy()
    normalized.data = np.log1p(normalized.data)
    obs = pd.DataFrame(index=[f'cell_{i}' for i in range(cells)])
    obs['batch'] = pd.Categorical([f'b{i % 8}' for i in range(cells)], categories=[f'b{i}' for i in range(9)], ordered=True)
    obs['qc'] = pd.array([None if i % 13 == 0 else True for i in range(cells)], dtype='boolean')
    obs['total'] = np.asarray(counts.sum(axis=1)).ravel()
    obs['number'] = pd.array([None if i % 17 == 0 else i for i in range(cells)], dtype='Int32')
    a = ad.AnnData(normalized, obs=obs, var=pd.DataFrame(index=[f'gene_{i}' for i in range(genes)]))
    a.layers['counts'] = counts
    a.obsm['X_pca'] = rng.normal(size=(cells, 30)).astype(np.float32)
    a.obsm['X_umap'] = rng.normal(size=(cells, 2)).astype(np.float32)
    k = min(15, cells)
    a.obsp['neighbors'] = sparse.csr_matrix((rng.uniform(size=cells*k), rng.integers(0, cells, cells*k), np.arange(cells+1)*k), shape=(cells, cells))
    a.write_h5ad(path, convert_strings_to_categoricals=False)
    return int(counts.nnz)


def size(path):
    return path.stat().st_size if path.is_file() else sum(p.stat().st_size for p in path.rglob('*') if p.is_file())


def source_hash(path):
    digest = hashlib.sha256()
    for file in sorted(path.rglob('*.py')):
        digest.update(str(file.relative_to(path)).replace('\\', '/').encode())
        digest.update(file.read_bytes())
    return digest.hexdigest()


def run(args):
    args.output.mkdir(parents=True, exist_ok=False)
    inputs = args.inputs.resolve(); inputs.mkdir(parents=True, exist_ok=True)
    records = []
    for cells in args.cells:
        source = inputs/f'{cells}-{args.genes}-{args.per_cell}-{args.seed}.h5ad'
        # Inputs are generated once and reused byte-for-byte in baseline/improved runs.
        if not source.exists():
            make_input(source, cells, args.genes, args.per_cell, args.seed)
        digest = hashlib.sha256(source.read_bytes()).hexdigest()
        for repeat in range(args.repeats):
            folder = args.output/f'{cells}-{repeat}.scx'
            target = args.output/f'{cells}-{repeat}.h5ad'
            for command in [('from-h5ad', source, folder), ('validate', folder), ('to-h5ad', folder, target)]:
                metrics = args.output/f'{cells}-{repeat}-{command[0]}.json'
                cmd = [sys.executable, str(Path(__file__).resolve()), 'worker', '--implementation', str(args.implementation.resolve()), '--metrics', str(metrics), '--command', *map(str, command)]
                if args.profile and cells == args.cells[0] and repeat == 0:
                    cmd[cmd.index('--command'):cmd.index('--command')] = ['--profile', str(args.output/f'{command[0]}.prof')]
                start = time.perf_counter()
                result = subprocess.run(cmd, capture_output=True, text=True)
                elapsed = time.perf_counter() - start
                if result.returncode:
                    raise RuntimeError(result.stderr)
                record = {'cells': cells, 'repeat': repeat, 'command': command[0], 'wall_seconds': elapsed, **json.loads(metrics.read_text()),
                          'output_bytes': size(target if command[0] == 'to-h5ad' else folder), 'input_h5ad_bytes': size(source)}
                record['input_sha256'] = digest
                records.append(record)
                print(json.dumps(record), flush=True)
    report = {'label': args.label, 'seed': args.seed, 'genes': args.genes, 'draws_per_cell': args.per_cell,
              'cells': args.cells, 'repeats': args.repeats, 'profiled': args.profile,
              'thread_limits': {'blas_openmp_requested': 1, 'matrix_market_parallelism': 'default/automatic (not explicitly limited)'},
              'platform': platform.platform(), 'python': platform.python_version(),
              'implementation_sha256': source_hash(args.implementation),
              'memory_method': 'OS peak RSS / Windows peak working set of isolated CLI worker',
              'versions': {p: version(p) for p in ('scx-format', 'anndata', 'numpy', 'pandas', 'scipy', 'h5py')},
              'records': records}
    (args.output/'report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')


def r_run(args):
    import psutil
    args.output.mkdir(parents=True, exist_ok=False)
    args.inputs.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env['PYTHONPATH'] = os.pathsep.join([str(args.implementation.resolve()), str(Path('benchmarks/thread-limit').resolve())])
    env['SCX_BENCHMARK_THREADS'] = '1'
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        env[name] = '1'
    records = []
    info = None
    for cells in args.cells:
        source = args.inputs/f'r-{cells}-{args.genes}-{args.per_cell}-{args.seed}.rds'
        if not source.exists():
            subprocess.run([args.rscript, '--vanilla', 'benchmarks/worker.R', 'generate', str(source), str(cells), str(args.genes), str(args.per_cell), str(args.seed)], check=True, env=env)
        info = json.loads(Path(str(source)+'.json').read_text())
        digest = hashlib.sha256(source.read_bytes()).hexdigest()
        for repeat in range(args.repeats):
            folder = args.output/f'{cells}-{repeat}.scx'
            target = args.output/f'{cells}-{repeat}.rds'
            for mode, input_path, output_path in [('export', source, folder), ('import', folder, target)]:
                command = [args.rscript, '--vanilla', 'benchmarks/worker.R', mode, str(args.adapter.resolve()), str(input_path), str(output_path), sys.executable]
                log = args.output/f'{cells}-{repeat}-{mode}.log'
                peak = 0
                start = time.perf_counter()
                with log.open('w', encoding='utf-8') as stream:
                    process = subprocess.Popen(command, stdout=stream, stderr=stream, env=env)
                    parent = psutil.Process(process.pid)
                    while process.poll() is None:
                        try:
                            members = [parent, *parent.children(recursive=True)]
                        except psutil.NoSuchProcess:
                            members = []
                        rss = 0
                        for member in members:
                            try:
                                rss += member.memory_info().rss
                            except psutil.NoSuchProcess:
                                pass
                        peak = max(peak, rss)
                        time.sleep(0.01)
                elapsed = time.perf_counter() - start
                if process.returncode:
                    raise RuntimeError(log.read_text())
                assert peak > 0
                record = {'cells': cells, 'repeat': repeat, 'command': f'seurat-{mode}', 'wall_seconds': elapsed,
                          'peak_bytes': peak, 'output_bytes': size(output_path), 'input_rds_bytes': size(source), 'input_sha256': digest}
                records.append(record)
                print(json.dumps(record), flush=True)
            subprocess.run([args.rscript, '--vanilla', 'benchmarks/verify-r.R', str(source), str(target)], check=True, env=env)
    report = {'label': args.label, 'seed': args.seed, 'genes': args.genes, 'draws_per_cell': args.per_cell,
              'cells': args.cells, 'repeats': args.repeats, 'platform': platform.platform(),
              'thread_limits': {'blas_openmp_requested': 1, 'matrix_market_parallelism': 'default/automatic (not explicitly limited)'},
              'memory_method': 'sampled process-tree RSS at 10 ms (includes Python validator child)',
              'adapter_sha256': hashlib.sha256(args.adapter.read_bytes()).hexdigest(),
              'python_implementation_sha256': source_hash(args.implementation),
              'R': info, 'records': records, 'round_trip_passed': True}
    (args.output/'report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='mode', required=True)
    w = sub.add_parser('worker')
    w.add_argument('--implementation', type=Path, required=True)
    w.add_argument('--metrics', type=Path, required=True)
    w.add_argument('--profile', type=Path)
    w.add_argument('--command', nargs=argparse.REMAINDER, required=True)
    for mode in ('run', 'r-run'):
        r = sub.add_parser(mode)
        r.add_argument('--implementation', type=Path, default=Path('src'))
        r.add_argument('--output', type=Path, required=True)
        r.add_argument('--inputs', type=Path, default=Path('benchmarks/work/inputs'))
        r.add_argument('--label', default='current')
        r.add_argument('--cells', type=int, nargs='+', default=[1000, 10000, 50000])
        r.add_argument('--genes', type=int, default=2000)
        r.add_argument('--per-cell', type=int, default=100)
        r.add_argument('--seed', type=int, default=42)
        r.add_argument('--repeats', type=int, default=2 if mode == 'run' else 1)
        if mode == 'run':
            r.add_argument('--profile', action='store_true')
        else:
            r.add_argument('--rscript', default='Rscript')
            r.add_argument('--adapter', type=Path, default=Path('adapters/seurat/R/scx.R'))
    args = parser.parse_args()
    if args.mode == 'worker':
        worker(args)
    else:
        if min(args.cells) < 1 or args.genes < 1 or args.per_cell < 1 or args.repeats < 1:
            parser.error('dimensions and repeats must be positive')
        (r_run if args.mode == 'r-run' else run)(args)


if __name__ == '__main__':
    main()
