"""A(0) common-input T3 T1/T2/T3 sensitivity and original envelope-HTML rebuild.

Uses the canonical BenchmarkSensitivity.run_t3_sensitivity_scan and the original
Reports/output/display/common_input_model_envelopes_all.html as an HTML template.
Does not overwrite that original or recompute the original fitted benchmarks.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
from pathlib import Path
import re
import time

import numpy as np
from Numerical.fitting.BenchmarkSensitivity import run_t3_sensitivity_scan
from Numerical.orchestration import RunA0SharedBenchmark as a0
from Numerical.orchestration import RunSmallLambdaRefit as base
from Numerical.plotting.BuildFittedModelEnvelopes import OBS, LABELS, JUNO_FUTURE_6Y, short_name

ROOT = base.STUDY / 'a0_shared_sensitivity'
ORIGINAL = base.ROOT / 'Reports/output/display/common_input_model_envelopes_all.html'
DISPLAY = base.ROOT / 'Reports/output/display_a0_shared_fixed/common_input_model_envelopes_all.html'
DIRECTIONS = ('T1', 'T2', 'T3')


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temp.replace(path)


def settings(args):
    return {'gamma_min': args.gamma_min, 'gamma_max': args.gamma_max,
            'n_points': args.points, 'spacing': 'log'}


def run(args):
    if not (0 < args.gamma_min < 1 < args.gamma_max):
        raise ValueError('Require gamma_min < 1 < gamma_max with positive gamma_min')
    if args.points < 3 or args.points % 2 != 1:
        raise ValueError('--points must be odd and >= 3')
    shared, models = a0.load_inputs(args.artifact_map)
    if args.selected_only:
        accepted = set()
        comparison = base.load(a0.ROOT / 'summary.json')['models']
        accepted = {item['model'] for item in comparison
                    if item.get('status') == 'Success' and item.get('pass_strict_5p86')}
        models = [m for m in models if m[0]['model'] in accepted]
    if not models:
        raise ValueError('No selected models')
    config = settings(args)
    fingerprint = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()
    records = []
    for k, (row, payload, artifacts) in enumerate(models, 1):
        name = row['model']
        folder = ROOT / name
        result_file = folder / 'sensitivity_diagnostics.json'
        if args.resume and result_file.is_file():
            old = base.load(result_file)
            if old.get('study_sha256') == fingerprint and old.get('status') == 'Success':
                records.append({'model': name, 'status': 'Success', 'verified_chi2': old.get('centre_chi2'), 'source': 'Resumed'})
                print(f'[{k}/{len(models)}] SKIP {name}: checkpoint exists', flush=True)
                continue
        cfg = copy.deepcopy(payload)
        cfg['sensitivity'] = {'enabled': True, **config}
        started = time.monotonic()
        print(f'[{k}/{len(models)}] SCAN {name}: 3 directions x {args.points} gamma values', flush=True)
        try:
            result = run_t3_sensitivity_scan(
                base_payload=cfg, uv_rgbeta_path=artifacts[0],
                intermediate_rgbeta_path=artifacts[1], final_weinberg_path=artifacts[2],
                oscillation_target_path=base.TARGET)
            result['model'] = name
            result['study_sha256'] = fingerprint
            result['study_metadata'] = {'donor': a0.DONOR, 'fixed_lambda': float(cfg['base_state']['ordinary']['lambdaT3']['real']), **config}
            middle = []
            for direction in DIRECTIONS:
                points = result['directions'][direction]['points']
                centres = [p for p in points if abs(p['gamma'] - 1.0) < 1e-12 and p.get('status') == 'Success']
                if len(centres) != 1:
                    raise RuntimeError(f'{direction} missing successful gamma=1 centre')
                middle.append(centres[0])
            if not np.allclose([middle[0]['chi2']], [middle[1]['chi2'], middle[2]['chi2']], rtol=1e-5, atol=1e-6):
                raise RuntimeError('Direction centres have inconsistent chi2')
            result['centre_chi2'] = float(middle[0]['chi2'])
            write(result_file, result)
            record = {'model': name, 'status': 'Success', 'verified_chi2': result['centre_chi2'], 'failed_points': result['n_failed_points']}
            print(f'  DONE {name}: chi2(gamma=1)={result["centre_chi2"]:.6g}, failed gamma points={result["n_failed_points"]}, elapsed={time.monotonic()-started:.1f}s', flush=True)
        except Exception as exc:
            record = {'model': name, 'status': 'Failed', 'message': f'{type(exc).__name__}: {exc}'}
            print(f'  FAILED {name}: {record["message"]}', flush=True)
        records.append(record)
        write(ROOT / 'comparison_summary.json', records)
    write(ROOT / 'comparison_summary.json', records)
    print(f'Saved {sum(x["status"] == "Success" for x in records)}/{len(records)} model scans: {ROOT}', flush=True)


def collect(source):
    rows_file = source / 'comparison_summary.json'
    rows = base.load(rows_file)
    models, skipped, target = [], [], None
    for row in rows:
        name = row['model']
        if row['status'] != 'Success':
            skipped.append(row)
            continue
        path = source / name / 'sensitivity_diagnostics.json'
        if not path.is_file():
            skipped.append({'model':name,'status':'MissingSavedScan'})
            continue
        scan = base.load(path)
        dirs, centres = {}, []
        try:
            for key in DIRECTIONS:
                pts = [{'gamma':float(p['gamma']), **{o:float(p['prediction'][o]) for o in OBS}}
                       for p in scan['directions'][key]['points']
                       if p.get('status') == 'Success' and all(o in p.get('prediction',{}) for o in OBS)]
                pts.sort(key=lambda x:x['gamma'])
                centre = [p for p in pts if abs(p['gamma']-1)<1e-12]
                if len(centre)!=1 or len(pts)<3:
                    raise ValueError(f'{key} missing centre or insufficient successful points')
                centres.append(centre[0]); dirs[key]=pts
            for c in centres[1:]:
                if not np.allclose([c[o] for o in OBS], [centres[0][o] for o in OBS], rtol=1e-7, atol=1e-10):
                    raise ValueError('gamma=1 predictions mismatch')
            if target is None: target=scan['target']
            models.append({'model':name,'status':'Success', 'chi2':row.get('verified_chi2'),
                           'fit_stage':'A(0) shared UV Yukawas', 'label':short_name(name),
                           'centre':centres[0], 'directions':dirs})
        except (KeyError, ValueError, TypeError) as exc:
            skipped.append({'model':name, 'status':'InvalidSavedScan', 'message':str(exc)})
    return {'models':models,'skipped':skipped,'target':target or {}, 'labels':LABELS,
            'observables':list(OBS),'source':str(source),
            'future_sigma':JUNO_FUTURE_6Y,
            'future_source':'JUNO 6-year projection, CPC 46 (2022) 123001, Table 6; NuFIT centre retained',
            'envelope_definition':'2x empirical sample standard-deviation covariance ellipse of saved A(0) common-UV-input gamma sensitivity points; descriptive only'}


def build(args):
    template = args.template.resolve()
    if not template.is_file(): raise FileNotFoundError(f'Missing original common-input envelope HTML: {template}')
    source = args.source.resolve()
    data = collect(source)
    if not data['models']: raise RuntimeError('No valid A(0) gamma scans available')
    page = template.read_text(encoding='utf-8-sig')
    start_marker, end_marker = 'const D=', ';const $='
    start = page.find(start_marker)
    end = page.find(end_marker, start + len(start_marker))
    if start < 0 or end < 0: raise ValueError('Original HTML no longer has the verified data-block markers')
    encoded = json.dumps(data, ensure_ascii=False, allow_nan=False).replace('</', '<\\/')
    page = page[:start + len(start_marker)] + encoded + page[end:]
    page = page.replace('Common-input T3 prediction-space sensitivity', 'A(0) shared-input T3 prediction-space sensitivity')
    page = page.replace('common-input models plotted', 'A(0) shared-input models plotted')
    page = page.replace('about each common-input centre', 'about each A(0) shared-input centre')
    page = page.replace('common-input centre', 'A(0) shared-input centre')
    page = page.replace('Source: '+str(args.template), 'Source: '+str(source))
    dest = args.output.resolve()
    if dest == template: raise ValueError('Refusing to overwrite the original historical envelope')
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(page, encoding='utf-8')
    write(dest.parent / 'a0_sensitivity_display_manifest.json', {
        'template':str(template), 'template_sha256':hashlib.sha256(template.read_bytes()).hexdigest(),
        'source':str(source), 'model_count':len(data['models']), 'skipped':data['skipped'],
        'warning':'Envelopes are empirical gamma-scan covariances, not parameter confidence regions or statistical model errors.'})
    print(f'Wrote {dest}; models={len(data["models"])}, skipped={len(data["skipped"])}',flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('run','build','all'))
    parser.add_argument('--gamma-min',type=float,default=0.9)
    parser.add_argument('--gamma-max',type=float,default=1.1)
    parser.add_argument('--points',type=int,default=5,help='Odd number of gamma points; first pass defaults to five')
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--selected-only',action='store_true',help='Only five strict-pass models, otherwise all 16')
    parser.add_argument('--artifact-map',type=Path)
    parser.add_argument('--template',type=Path,default=ORIGINAL)
    parser.add_argument('--source',type=Path,default=ROOT)
    parser.add_argument('--output',type=Path,default=DISPLAY)
    args=parser.parse_args()
    if args.action in ('run','all'):run(args)
    if args.action in ('build','all'):build(args)

if __name__ == '__main__':main()
