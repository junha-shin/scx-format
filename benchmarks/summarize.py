"""Generate comparison tables from raw measured runs; never time profiled runs."""
import argparse
import csv
import json
from pathlib import Path
import statistics


def load(path):
    report = json.loads(path.read_text())
    if report.get('profiled'):
        raise ValueError('profiled runs cannot be used as timing evidence')
    return report


def groups(report):
    result = {}
    for row in report['records']:
        result.setdefault((row['cells'], row['command']), []).append(row)
    return result


def median(rows, field):
    return statistics.median(r[field] for r in rows)


def render(before, after, title):
    for key in ('seed', 'genes', 'draws_per_cell'):
        assert before[key] == after[key]
    b, a = groups(before), groups(after)
    lines = [f'## {title}', '', '| Cells | Operation | Before s | After s | Speedup | Before peak MiB | After peak MiB | Before output MiB | After output MiB |',
             '|---:|---|---:|---:|---:|---:|---:|---:|---:|']
    rows = []
    for key in sorted(a):
        current, old = a[key], b.get(key)
        if old:
            assert {r['input_sha256'] for r in old} == {r['input_sha256'] for r in current}, 'input changed'
        row = {'cells': key[0], 'command': key[1], 'before_seconds': median(old, 'wall_seconds') if old else None,
               'after_seconds': median(current, 'wall_seconds'),
               'before_peak_mib': max(r['peak_bytes'] for r in old)/2**20 if old else None,
               'after_peak_mib': max(r['peak_bytes'] for r in current)/2**20,
               'before_output_mib': median(old, 'output_bytes')/2**20 if old else None,
               'after_output_mib': median(current, 'output_bytes')/2**20}
        row['speedup'] = row['before_seconds']/row['after_seconds'] if old else None
        rows.append(row)
        def fmt(value):
            return 'not measured' if value is None else f'{value:.2f}'
        lines.append(f"| {key[0]:,} | {key[1]} | " + ' | '.join(fmt(row[k]) for k in ('before_seconds', 'after_seconds', 'speedup', 'before_peak_mib', 'after_peak_mib', 'before_output_mib', 'after_output_mib')) + ' |')
    return lines, rows


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--improved', type=Path, required=True)
    parser.add_argument('--r-baseline', type=Path, required=True)
    parser.add_argument('--r-improved', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    lines = ['# Measured scaling results', '', 'Wall time includes process startup, imports and I/O. Python: median of two runs; R: one exploratory run. Peaks: maximum observed per case. R memory is sampled process-tree RSS; Python memory is OS process peak. See README for methods and limits.', '']
    data = []
    for before, after, title in [(args.baseline, args.improved, 'Python CLI'), (args.r_baseline, args.r_improved, 'Seurat R')]:
        table, rows = render(load(before), load(after), title)
        lines.extend(table + [''])
        data.extend({'engine': title, **r} for r in rows)
    lines += ['R baseline 50,000-cell cases were not completed: the old finite check can allocate a dense 50,000 × 50,000 logical graph (~10 GB). No measured time, peak or speedup is inferred for those cases.', '',
              'SCX gzip level changed from 9 to 6 in Python. File size differences are measured, not assumed; H5AD and RDS outputs use their documented uncompressed defaults.', '']
    (args.output/'results.md').write_text('\n'.join(lines), encoding='utf-8')
    with (args.output/'summary.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(data[0]))
        writer.writeheader(); writer.writerows(data)
