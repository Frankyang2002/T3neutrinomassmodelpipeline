"""Run three prefactor-preserving T3 scans and plot NuFIT correlation panels.

Usage:
  python -m Numerical.orchestration.PlotCommonYukawaCorrelations
  python -m Numerical.orchestration.PlotCommonYukawaCorrelations --run --points 7 --gamma-min 0.9 --gamma-max 1.1

Reads the latest common least-squares Yukawas; never edits fitted benchmarks.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from copy import deepcopy
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from Numerical.orchestration.RunCommonParameterComparisonV2 import (
    PROJECT_ROOT, SCENARIO, OBS, candidate_for_target, models_and_cases, load, evaluate,
)
from Numerical.fitting.BenchmarkSensitivity import scaled_config_payload
from Numerical.fitting.OscillationFit import OscillationFitTarget

INPUT = PROJECT_ROOT / 'Reports/output/full/common_yukawa_least_squares' / SCENARIO
OUTPUT = PROJECT_ROOT / 'Reports/output/full/common_yukawa_correlation_scans' / SCENARIO
COLORS = {'T1': '#2376c9', 'T2': '#ed713d', 'T3': '#269a57'}


def _write_json(path, obj):
    with path.open('x', encoding='utf-8') as f:
        json.dump(obj, f, indent=2, allow_nan=False)
        f.write('\n')


def _plot(model, rows, target, out):
    center = dict(zip(target.observable_names, target.central_values))
    sigma = dict(zip(target.observable_names, np.sqrt(np.diag(target.covariance))))
    panels = [('sin2_theta13', 'delta_m3l_sq_ev2', r'$\sin^2\theta_{13}$', r'$\Delta m_{31}^2\;[\mathrm{eV}^2]$'),
              ('sin2_theta23', 'delta_m3l_sq_ev2', r'$\sin^2\theta_{23}$', r'$\Delta m_{31}^2\;[\mathrm{eV}^2]$')]
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.8), constrained_layout=True)
    for ax, (x, y, xl, yl) in zip(axes, panels):
        for n, alpha in [(3, .11), (1, .24)]:
            ax.axvspan(center[x]-n*sigma[x], center[x]+n*sigma[x], color='#268ed0', alpha=alpha, lw=0)
            ax.axhspan(center[y]-n*sigma[y], center[y]+n*sigma[y], color='#268ed0', alpha=alpha, lw=0)
        ax.scatter([center[x]], [center[y]], marker='+', s=90, color='black', zorder=5, label='NuFIT centre')
        for direction, color in COLORS.items():
            group = sorted((r for r in rows if r['direction'] == direction and r['status'] == 'Success'), key=lambda r: r['gamma'])
            if not group:
                continue
            ax.plot([r[x] for r in group], [r[y] for r in group], '-o', ms=3, lw=1.5, color=color, label=direction)
            bm = min(group, key=lambda r: abs(r['gamma']-1))
            if abs(bm['gamma']-1) < 1e-10:
                ax.scatter([bm[x]], [bm[y]], s=60, marker='x', color=color, zorder=6)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.grid(alpha=.17)
        ax.ticklabel_format(axis='y', style='sci', scilimits=(-3,3))
    axes[0].legend(loc='best', fontsize=8)
    fig.suptitle(model.replace('T3_', 'T3 '), fontsize=11)
    fig.savefig(out / f'{model}_correlations.png', dpi=170)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', action='store_true')
    parser.add_argument('--resume', action='store_true', help='Resume existing partial scan without overwriting completed rows')
    parser.add_argument('--points', type=int, default=7)
    parser.add_argument('--gamma-min', type=float, default=.9)
    parser.add_argument('--gamma-max', type=float, default=1.1)
    args = parser.parse_args(argv)
    if args.points < 3 or args.points % 2 == 0 or not (0 < args.gamma_min < 1 < args.gamma_max):
        parser.error('Require odd --points >= 3 and 0 < gamma-min < 1 < gamma-max')
    configs, cases = models_and_cases()
    common_path = INPUT / 'optimized_common_input.json'
    if not common_path.is_file():
        raise FileNotFoundError('Latest least-squares common input missing: ' + str(common_path))
    common = load(common_path)
    ordinary = common['shared_ordinary']
    target = OscillationFitTarget.from_json(PROJECT_ROOT / 'configs/nufit_6_1_ic24_no.json')
    if tuple(target.observable_names) != OBS or target.ordering != 'NO':
        raise ValueError('Unexpected NuFIT target')
    print('Input:', common_path)
    print('Output:', OUTPUT)
    print('Models:', len(configs), 'directions: T1 T2 T3', 'points:', args.points)
    print('Estimated evaluations:', len(configs)*3*args.points)
    print('CP phase is not available from the existing five-observable fit interface; only two physical panels are produced.')
    if not args.run:
        print('Dry run only. Add --run to calculate trajectories.')
        return 0
    if OUTPUT.exists() and not args.resume:
        raise FileExistsError('Output exists; add --resume to continue: ' + str(OUTPUT))
    if args.resume and OUTPUT.exists() and (OUTPUT / 'scan_summary.json').exists():
        raise FileExistsError('Completed scan already exists: ' + str(OUTPUT))
    # Validate all model configs and inputs before creating outputs.
    source_name = common['source_model']
    source = load(configs[source_name])
    source = deepcopy(source)
    for key in ('y1','y2','lambdaT3','MF','mS1Sq','mS2Sq'):
        source['base_state']['ordinary'][key] = deepcopy(ordinary[key])
    payloads = {name: candidate_for_target(source, load(path)) for name, path in configs.items()}
    OUTPUT.mkdir(parents=True, exist_ok=True)
    gammas = np.linspace(args.gamma_min, args.gamma_max, args.points)
    csv_path = OUTPUT / 'scan_points.csv'
    fieldnames = ['model', 'direction', 'gamma', 'status', 'error', 'chi2', *OBS]
    all_rows = []
    completed = {}
    if args.resume and csv_path.exists():
        with csv_path.open(newline='', encoding='utf-8') as f:
            for old in csv.DictReader(f):
                if old['status'] == 'Success':
                    completed[(old['model'], old['direction'], round(float(old['gamma']), 12))] = old
        print('Recovered completed points:', len(completed), flush=True)
    def checkpoint():
        temporary = csv_path.with_suffix('.csv.tmp')
        with temporary.open('w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_rows)
        os.replace(temporary, csv_path)
    for name in sorted(payloads):
        model_rows = []
        for direction in ('T1','T2','T3'):
            for gamma in gammas:
                key = (name, direction, round(float(gamma), 12))
                if key in completed:
                    row = dict(completed[key])
                    row['gamma'] = float(gamma)
                    for k in ['chi2', *OBS]:
                        row[k] = float(row[k])
                    model_rows.append(row)
                    all_rows.append(row)
                    checkpoint()
                    print(f'{name} {direction} gamma={gamma:.4f}: Reused', flush=True)
                    continue
                row = {'model': name, 'direction': direction, 'gamma': float(gamma), 'status': 'Success', 'error': '', 'chi2': ''}
                row.update({key: '' for key in OBS})
                try:
                    candidate = scaled_config_payload(payloads[name], direction, float(gamma))
                    chi2, prediction = evaluate(candidate, cases[name], target)
                    row.update(chi2=chi2, **prediction)
                except Exception as exc:
                    row['status'] = 'Failed'
                    row['error'] = f'{type(exc).__name__}: {exc}'
                model_rows.append(row)
                all_rows.append(row)
                checkpoint()
                print(f'{name} {direction} gamma={gamma:.4f}: {row["status"]} chi2={row["chi2"]}', flush=True)
        _plot(name, model_rows, target, OUTPUT)
        checkpoint()
    if any(r['status'] != 'Success' for r in all_rows):
        print('Some evaluations failed; scan_points.csv saved for resumption.', flush=True)
        return 2
    _write_json(OUTPUT / 'scan_summary.json', {
        'input': str(common_path), 'model_count': len(payloads), 'directions': {
            'T1': 'lambdaT3 *= gamma; y1,y2 /= sqrt(gamma)',
            'T2': 'lambdaT3 *= gamma; y2 /= gamma',
            'T3': 'y1 *= gamma; y2 /= gamma',
        }, 'gamma_min': args.gamma_min, 'gamma_max': args.gamma_max,
        'points_per_direction': args.points, 'evaluations': len(all_rows),
        'successes': sum(r['status']=='Success' for r in all_rows),
        'nuFIT_source': target.source,
        'regions': 'Gaussian marginal 1sigma and 3sigma rectangular bands, not official joint contours',
        'cp_phase': 'Unavailable from existing five-observable pipeline; no delta_CP panel fabricated',
        'scalars': 'Fixed to latest least-squares common input, except lambdaT3 in T1/T2 scans',
        'representation_specific_couplings': 'Preserved from each model reference config',
    })
    print('Finished:', OUTPUT)
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
