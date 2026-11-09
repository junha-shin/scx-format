"""Benchmark-only BLAS/OpenMP limit for the Python validator spawned by R.

Matrix Market's lazy native thread pool retains its default/automatic setting.
"""
import os

if os.environ.get('SCX_BENCHMARK_THREADS') == '1':
    import scipy.io
    from threadpoolctl import threadpool_limits
    _benchmark_limits = threadpool_limits(limits=1)
