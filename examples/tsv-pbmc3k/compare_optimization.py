"""Compare a validated baseline and optimized run; refuse mismatched inputs/methods."""
import argparse
import json
from pathlib import Path


def compare(before, after):
    for report in (before,after):
        if not report['checks'] or not report['records'] or not report['summary']:
            raise ValueError('Missing completed comparison/measurement records')
        if not report['passed'] or not all(x['passed'] for x in report['checks']):
            raise ValueError('Both runs must pass every source/language comparison')
        for check in report['checks']:
            if max(check['max_absolute_errors'].values()) > check['tolerance']:
                raise ValueError('Recorded numeric error exceeds tolerance')
    for key in ('input_files','shape','cases','thread_limits','memory_method','platform'):
        if before[key] != after[key]:
            raise ValueError(f'Cannot compare different {key}')
    def environment(report):
        return {(r['R'],tuple(sorted(r['versions'].items())),r['gzip_compression_level'])
                for r in report['records'] if r['engine']=='R'}
    if environment(before) != environment(after):
        raise ValueError('R runtime/dependencies/compression changed')
    rows = []
    for old in before['summary']:
        if old['engine'] != 'R': continue
        new = next(x for x in after['summary'] if x['engine']=='R' and x['case']==old['case'])
        if not old['case'].endswith('gzip'):
            old_hashes = {r['output_sha256'] for r in before['records'] if r['engine']=='R' and r['case']==old['case']}
            new_hashes = {r['output_sha256'] for r in after['records'] if r['engine']=='R' and r['case']==old['case']}
            if old_hashes != new_hashes: raise ValueError('Uncompressed R TSV bytes changed')
        rows.append({'case':old['case'],'before_wall_seconds':old['median_wall_seconds'],
                     'after_wall_seconds':new['median_wall_seconds'],
                     'speedup':old['median_wall_seconds']/new['median_wall_seconds'],
                     'before_peak_bytes':old['max_peak_bytes'],'after_peak_bytes':new['max_peak_bytes'],
                     'before_output_bytes':old['output_bytes'],'after_output_bytes':new['output_bytes']})
    return {'passed':True,'same_input_files':True,'same_R_environment_and_gzip':True,
            'plain_R_TSV_bytes_identical':True,'rows':rows,
            'baseline_repeats':before['repeats'],'optimized_repeats':after['repeats'],
            'note':'Historical baseline and a separate optimized session; OS cache was not reset. See runtime/source hashes in the original reports.'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('before',type=Path); p.add_argument('after',type=Path); p.add_argument('output',type=Path)
    args = p.parse_args()
    result = compare(json.loads(args.before.read_text(encoding='utf-8')),json.loads(args.after.read_text(encoding='utf-8')))
    with args.output.open('x',encoding='utf-8') as handle:
        handle.write(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__=='__main__': main()
