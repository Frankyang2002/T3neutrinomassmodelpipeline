"""Summarise saved fermion-start T1/T2/T3 scans without rerunning RGEs.

Outputs one row per case/mode/model/direction/observable, with the signed and
absolute maximum departure from the saved gamma=1 prediction. Mass-squared
splittings are in eV^2; mixing-angle observables are dimensionless.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SOURCE = ROOT / 'output' / 'studies' / 'fermion_start_lambda_comparison'
OBSERVABLES = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12', 'sin2_theta13', 'sin2_theta23')
DIRECTIONS = ('T1', 'T2', 'T3')
CASES = ('5e-9', '5e-10')
MODES = ('optimised', 'shared_a_m2')
FIELDS = ('case', 'mode', 'model', 'direction', 'observable', 'centre', 'max_abs_change',
          'max_fractional_change', 'max_percent_change', 'gamma_at_max',
          'signed_change_at_max', 'min_prediction', 'max_prediction',
          'valid_points', 'failed_points', 'gamma_min', 'gamma_max', 'relative_status')


def finite_number(value):
    try:
        v = float(value)
        return v if math.isfinite(v) else None
    except (ValueError, TypeError, OverflowError):
        return None


def summarise_scan(scan, case, mode, model, relative_floor):
    if scan.get('status') != 'Success':
        raise ValueError('Scan status is not Success')
    rows = []
    for direction in DIRECTIONS:
        points = scan.get('directions', {}).get(direction, {}).get('points', [])
        valid, failed = [], 0
        for point in points:
            gamma = finite_number(point.get('gamma'))
            prediction = point.get('prediction')
            if (point.get('status') != 'Success' or gamma is None or gamma <= 0 or
                    not isinstance(prediction, dict)):
                failed += 1
                continue
            values = {name: finite_number(prediction.get(name)) for name in OBSERVABLES}
            if any(value is None for value in values.values()):
                failed += 1
                continue
            valid.append((gamma, values))
        centres = [values for gamma, values in valid if math.isclose(gamma, 1.0, rel_tol=0, abs_tol=1e-9)]
        if len(centres) != 1:
            raise ValueError(f'{direction}: expected exactly one valid gamma=1 centre, found {len(centres)}')
        centre = centres[0]
        if len(valid) < 3:
            raise ValueError(f'{direction}: fewer than three valid scan points')
        for name in OBSERVABLES:
            base = centre[name]
            gamma_max, extreme = max(valid, key=lambda p: abs(p[1][name] - base))
            diff = extreme[name] - base
            near_zero = abs(base) <= relative_floor
            fraction = None if near_zero else abs(diff) / abs(base)
            rows.append(dict(case=case, mode=mode, model=model, direction=direction,
                             observable=name, centre=base, max_abs_change=abs(diff),
                             max_fractional_change=fraction,
                             max_percent_change=None if fraction is None else 100*fraction,
                             gamma_at_max=gamma_max, signed_change_at_max=diff,
                             min_prediction=min(v[name] for _, v in valid),
                             max_prediction=max(v[name] for _, v in valid),
                             valid_points=len(valid), failed_points=failed,
                             gamma_min=min(g for g, _ in valid), gamma_max=max(g for g, _ in valid),
                             relative_status='undefined_near_zero_centre' if near_zero else 'defined'))
    return rows


def write_csv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, default=DEFAULT_SOURCE)
    p.add_argument('--output', type=Path, default=None)
    p.add_argument('--case', choices=(*CASES, 'all'), default='all')
    p.add_argument('--mode', choices=(*MODES, 'all'), default='all')
    p.add_argument('--relative-floor', type=float, default=1e-30,
                   help='Below this absolute baseline, fractional sensitivity is undefined.')
    args = p.parse_args(argv)
    if not math.isfinite(args.relative_floor) or args.relative_floor < 0:
        p.error('--relative-floor must be a finite nonnegative number')
    cases = CASES if args.case == 'all' else (args.case,)
    modes = MODES if args.mode == 'all' else (args.mode,)
    output = args.output or (args.source / 'sensitivity_metrics')
    rows, errors, counts = [], [], {}
    for case in cases:
        for mode in modes:
            folder = args.source / case / 'reports_existing' / mode / 'sensitivity'
            files = sorted(folder.glob('*/sensitivity_diagnostics.json'))
            counts[f'{case}/{mode}'] = len(files)
            if len(files) != 16:
                errors.append(f'{case}/{mode}: found {len(files)} scans, expected 16')
            for filepath in files:
                try:
                    scan = json.loads(filepath.read_text(encoding='utf-8-sig'))
                    rows.extend(summarise_scan(scan, case, mode, filepath.parent.name,
                                               args.relative_floor))
                except (ValueError, TypeError, KeyError, json.JSONDecodeError) as exc:
                    errors.append(f'{case}/{mode}/{filepath.parent.name}: {exc}')
    if not rows:
        raise SystemExit('No valid sensitivity data found at ' + str(args.source))
    write_csv(output / 'sensitivity_by_model_direction_observable.csv', rows, FIELDS)
    aggregate = []
    for case in cases:
        for mode in modes:
            for name in OBSERVABLES:
                matched = [r for r in rows if r['case'] == case and r['mode'] == mode
                           and r['observable'] == name and r['max_fractional_change'] is not None]
                if not matched:
                    continue
                worst = max(matched, key=lambda r: r['max_fractional_change'])
                aggregate.append({k: worst[k] for k in ('case', 'mode', 'observable', 'model',
                    'direction', 'max_fractional_change', 'max_percent_change',
                    'max_abs_change', 'gamma_at_max')})
    write_csv(output / 'largest_sensitivity_by_observable.csv', aggregate,
              ('case', 'mode', 'observable', 'model', 'direction', 'max_fractional_change',
               'max_percent_change', 'max_abs_change', 'gamma_at_max'))
    metadata = {'scan_counts': counts, 'rows': len(rows), 'errors': errors,
                'relative_floor': args.relative_floor,
                'definition': 'max_gamma |O(gamma)-O(1)| / |O(1)| (when baseline > floor)',
                'input': str(args.source), 'output': str(output)}
    (output / 'summary.json').write_text(json.dumps(metadata, indent=2)+'\n', encoding='utf-8')
    print('Saved', output / 'sensitivity_by_model_direction_observable.csv')
    print('Saved', output / 'largest_sensitivity_by_observable.csv')
    print('Scans:', counts, '| metric rows:', len(rows))
    for row in aggregate:
        print(f"{row['case']:6} {row['mode']:11} {row['observable']:22} max={row['max_percent_change']:.6g}% "
              f"{row['model']} {row['direction']} gamma={row['gamma_at_max']:.4g}")
    if errors:
        print('Warnings/errors:')
        for error in errors:
            print(' ', error)
        raise SystemExit(1)


if __name__ == '__main__':
    main()
