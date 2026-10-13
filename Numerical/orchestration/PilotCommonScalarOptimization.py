"""Read-only Stage B scalar coordinate pilot for the 16 T3 models.

Run dry: python -m Numerical.orchestration.PilotCommonScalarOptimization
Run: python -m Numerical.orchestration.PilotCommonScalarOptimization --run
Resume: python -m Numerical.orchestration.PilotCommonScalarOptimization --run --resume

Seven candidates: baseline plus +/- fractional changes in lambdaT3, mS1Sq, mS2Sq.
Yukawas and fermion masses stay at the latest common least-squares optimum.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from copy import deepcopy
from pathlib import Path

from Numerical.orchestration.RunCommonParameterComparisonV2 import (
    PROJECT_ROOT, SCENARIO, OBS, candidate_for_target, models_and_cases, load,
)
from Numerical.orchestration.OptimizeCommonYukawas import make_evaluator
from Numerical.fitting.OscillationFit import OscillationFitTarget

SOURCE = PROJECT_ROOT / 'Reports/output/full/common_yukawa_least_squares' / SCENARIO / 'optimized_common_input.json'
OUTPUT = PROJECT_ROOT / 'Reports/output/full/common_scalar_pilot' / SCENARIO
FIELDS = ('lambdaT3', 'mS1Sq', 'mS2Sq')


def atomic_json(path: Path, data):
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    os.replace(temp, path)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run', action='store_true')
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--fraction', type=float, default=0.05)
    args = ap.parse_args(argv)
    if not (0 < args.fraction < 0.5):
        ap.error('--fraction must be between 0 and 0.5')
    if not SOURCE.is_file():
        raise FileNotFoundError(SOURCE)
    common = load(SOURCE)
    ordinary = common['shared_ordinary']
    configs, cases = models_and_cases()
    names = sorted(configs)
    if len(names) != 16:
        raise RuntimeError('Expected exactly 16 models')
    source = deepcopy(load(configs[common['source_model']]))
    for key in ('y1', 'y2', 'lambdaT3', 'MF', 'mS1Sq', 'mS2Sq'):
        source['base_state']['ordinary'][key] = deepcopy(ordinary[key])
    candidates = [('baseline', None, 1.0)]
    for field in FIELDS:
        for suffix, scale in (('minus', 1-args.fraction), ('plus', 1+args.fraction)):
            candidates.append((field+'_'+suffix, field, scale))
    print('Source:', SOURCE, flush=True)
    print('Output:', OUTPUT, flush=True)
    print('Candidates:', len(candidates), 'models:', len(names), 'evaluations:', len(candidates)*len(names), flush=True)
    if not args.run:
        print('Dry run; add --run to evaluate.', flush=True)
        return 0
    if OUTPUT.exists() and not args.resume:
        raise FileExistsError('Existing output; use --resume: '+str(OUTPUT))
    if args.resume and not OUTPUT.exists():
        raise FileNotFoundError('Cannot resume missing output: '+str(OUTPUT))
    OUTPUT.mkdir(parents=True, exist_ok=True)
    manifest = {'source': str(SOURCE), 'fraction': args.fraction, 'models': names,
                'candidates': [c[0] for c in candidates]}
    manifest_path = OUTPUT / 'manifest.json'
    if manifest_path.exists():
        if load(manifest_path) != manifest:
            raise ValueError('Resume settings/source differ from saved manifest')
    else:
        atomic_json(manifest_path, manifest)
    target = OscillationFitTarget.from_json(PROJECT_ROOT / 'configs/nufit_6_1_ic24_no.json')
    if target.ordering != 'NO' or tuple(target.observable_names) != OBS:
        raise ValueError('Unexpected NuFIT target')
    payloads = {name: candidate_for_target(source, load(configs[name])) for name in names}
    evaluators = {name: make_evaluator(payloads[name], cases[name], target) for name in names}
    records_path = OUTPUT / 'evaluations.csv'
    columns = ['candidate','model','status','chi2','error',*OBS]
    rows = {}
    if records_path.exists():
        with records_path.open(newline='', encoding='utf-8') as handle:
            for row in csv.DictReader(handle):
                rows[(row['candidate'], row['model'])] = row
    def checkpoint():
        tmp = records_path.with_suffix('.csv.tmp')
        with tmp.open('w', newline='', encoding='utf-8') as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            for key in sorted(rows):
                writer.writerow(rows[key])
        os.replace(tmp, records_path)
    for label, field, scale in candidates:
        for name in names:
            key = (label, name)
            if key in rows and rows[key]['status'] == 'Success':
                print(label, name, 'Reused', flush=True)
                continue
            payload = payloads[name]
            original = deepcopy(payload['base_state']['ordinary'])
            try:
                if field is not None:
                    spec = payload['base_state']['ordinary'][field]
                    if field == 'lambdaT3':
                        spec['real'] = float(spec['real']) * scale
                    else:
                        payload['base_state']['ordinary'][field] = float(spec) * scale
                fit = evaluators[name]({})
                chi2 = float(fit.chi2)
                if not math.isfinite(chi2):
                    raise ValueError('Nonfinite chi2')
                row = {'candidate':label, 'model':name, 'status':'Success',
                       'chi2':chi2, 'error':'', **dict(zip(fit.observable_names, map(float,fit.prediction)))}
            except Exception as exc:
                row = {'candidate':label, 'model':name, 'status':'Failed',
                       'chi2':'', 'error':f'{type(exc).__name__}: {exc}',
                       **{k:'' for k in OBS}}
            finally:
                payload['base_state']['ordinary'].clear()
                payload['base_state']['ordinary'].update(original)
            rows[key] = row
            checkpoint()
            print(label, name, row['status'], 'chi2='+str(row['chi2']), flush=True)
    rankings = []
    for label, field, scale in candidates:
        subset = [rows[(label,n)] for n in names]
        successful = [r for r in subset if r['status'] == 'Success']
        rankings.append({'candidate':label, 'parameter':field, 'scale':scale,
                         'successful_models':len(successful),
                         'mean_chi2':sum(float(r['chi2']) for r in successful)/16 if len(successful)==16 else None})
    valid = [r for r in rankings if r['mean_chi2'] is not None]
    best = min(valid, key=lambda r:r['mean_chi2']) if valid else None
    atomic_json(OUTPUT/'summary.json', {'rankings':rankings, 'best':best,
        'note':'One-at-a-time scalar pilot; Yukawas and fermion masses fixed. Not a joint optimisation.'})
    if best is not None:
        chosen = deepcopy(common)
        field = best['parameter']
        if field is not None:
            if field == 'lambdaT3':
                chosen['shared_ordinary'][field]['real'] *= best['scale']
            else:
                chosen['shared_ordinary'][field] *= best['scale']
        chosen['average_chi2'] = best['mean_chi2']
        chosen['objective'] = 'best of seven one-at-a-time scalar pilot candidates'
        atomic_json(OUTPUT/'best_common_input.json', chosen)
    print('Best:', best, flush=True)
    print('Saved:', OUTPUT, flush=True)
    return 0 if len(valid)==len(candidates) else 2

if __name__ == '__main__':
    raise SystemExit(main())
