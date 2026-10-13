"""Analyze fitted T3 trajectories as finite-difference shifts in NuFIT sigma units.

Run: python -m Numerical.plotting.NormalizedFittedSensitivity [--scenario smallY_smallL]
Reads existing sensitivity JSON; does not rerun or modify physics calculations.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OBS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12', 'sin2_theta13', 'sin2_theta23')
LABELS = (r'$\Delta m^2_{21}$', r'$\Delta m^2_{3\ell}$', r'$\sin^2\theta_{12}$', r'$\sin^2\theta_{13}$', r'$\sin^2\theta_{23}$')
MODEL_RE = re.compile(r'T3_dS1_(\d+)_dS2_(\d+)_dF_(\d+)_alpha_([mp])(\d+)')


def representation(model):
    m = MODEL_RE.fullmatch(model)
    if not m:
        raise ValueError(f'Unexpected model identifier: {model}')
    d1, d2, df, sign, alpha = m.groups()
    return f'({d1},{d2},{df}), {"-" if sign == "m" else "+"}{alpha}'


def analyze(input_dir: Path, output_dir: Path) -> None:
    files = sorted(input_dir.glob('T3_*/sensitivity_diagnostics.json'))
    if not files:
        raise FileNotFoundError(f'No sensitivity_diagnostics.json files in {input_dir}')
    rows = []
    sigmas_reference = None
    for path in files:
        data = json.loads(path.read_text(encoding='utf-8-sig'))
        if data.get('status') != 'Success':
            raise ValueError(f'Unsuccessful scan: {path}')
        target = data['target']
        names = target['observable_names']
        cov = np.asarray(target['covariance'], dtype=float)
        if cov.shape != (len(names), len(names)):
            raise ValueError(f'Invalid target covariance: {path}')
        sigma = dict(zip(names, np.sqrt(np.diag(cov))))
        if not all(np.isfinite(sigma[k]) and sigma[k] > 0 for k in OBS):
            raise ValueError(f'Invalid target sigmas: {path}')
        sigvec = np.array([sigma[k] for k in OBS])
        if sigmas_reference is None:
            sigmas_reference = sigvec
        elif not np.allclose(sigvec, sigmas_reference, rtol=1e-10, atol=0):
            raise ValueError(f'Inconsistent NuFIT covariance: {path}')
        for direction in ('T1', 'T2', 'T3'):
            pts = sorted(data['directions'][direction]['points'], key=lambda p: p['gamma'])
            gamma = np.asarray([float(p['gamma']) for p in pts])
            if len(pts) < 3 or not np.all(np.diff(gamma) > 0):
                raise ValueError(f'Invalid gamma grid: {path} {direction}')
            matrix = np.array([[float(p['prediction'][name]) for name in OBS] for p in pts])
            if not np.all(np.isfinite(matrix)):
                raise ValueError(f'Nonfinite predictions: {path} {direction}')
            if not np.any(np.isclose(gamma, 1., atol=1e-10, rtol=0)):
                raise ValueError(f'Gamma=1 baseline missing: {path} {direction}')
            base = matrix[np.argmin(abs(gamma-1.))]
            # Central finite-difference slope in units of sigma per unit gamma.
            slope = (matrix[-1]-matrix[0]) / (gamma[-1]-gamma[0]) / sigvec
            # Maximum displacement relative to the actual gamma=1 benchmark.
            maximum = np.max(np.abs((matrix-base)/sigvec), axis=0)
            for i, name in enumerate(OBS):
                rows.append(dict(model=path.parent.name, representation=representation(path.parent.name),
                                 direction=direction, observable=name,
                                 sigma=float(sigvec[i]), gamma_min=float(gamma[0]),
                                 gamma_max=float(gamma[-1]), baseline=float(base[i]),
                                 endpoint_difference=float(matrix[-1,i]-matrix[0,i]),
                                 slope_sigma_per_gamma=float(slope[i]),
                                 max_abs_shift_sigma=float(maximum[i])))
    output_dir.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with (output_dir/'normalized_sensitivity.csv').open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader(); w.writerows(rows)
    (output_dir/'normalized_sensitivity.json').write_text(json.dumps(rows, indent=2)+'\n', encoding='utf-8')
    models = sorted({r['model'] for r in rows})
    for direction in ('T1','T2','T3'):
        values = np.array([[next(r['slope_sigma_per_gamma'] for r in rows if r['model']==m and r['direction']==direction and r['observable']==obs) for obs in OBS] for m in models])
        fig, ax = plt.subplots(figsize=(10, max(5, .38*len(models)+1.5)))
        vmax = float(np.max(abs(values))) or 1.0
        im = ax.imshow(values, aspect='auto', cmap='coolwarm', vmin=-vmax, vmax=vmax)
        ax.set_xticks(range(len(OBS)), LABELS)
        ax.set_yticks(range(len(models)), [representation(m) for m in models], fontsize=8)
        ax.set_title(f'{direction}: normalized endpoint slope (NuFIT σ per unit γ)')
        fig.colorbar(im, ax=ax, label=r'$[O(\gamma_{max})-O(\gamma_{min})]/[(\gamma_{max}-\gamma_{min})\sigma]$')
        fig.tight_layout()
        fig.savefig(output_dir/f'{direction}_normalized_heatmap.png', dpi=180)
        plt.close(fig)
    # Class-level aggregation: no physical attribution to multiplicity is implied.
    group_rows = []
    for direction in ('T1','T2','T3'):
        for dims in sorted({tuple(MODEL_RE.fullmatch(m).group(i) for i in (1,2,3)) for m in models}):
            subset = [r for r in rows if r['direction']==direction and tuple(MODEL_RE.fullmatch(r['model']).group(i) for i in (1,2,3))==dims]
            for obs in OBS:
                x = [r for r in subset if r['observable']==obs]
                group_rows.append({'representation_dims':','.join(dims), 'direction':direction, 'observable':obs,
                                   'n_models':len(x), 'median_abs_slope_sigma_per_gamma':float(np.median([abs(r['slope_sigma_per_gamma']) for r in x])),
                                   'max_abs_shift_sigma':float(max(r['max_abs_shift_sigma'] for r in x))})
    with (output_dir/'class_summary.csv').open('w', newline='', encoding='utf-8') as f:
        w=csv.DictWriter(f, fieldnames=list(group_rows[0]));w.writeheader();w.writerows(group_rows)
    print(f'Analyzed {len(models)} models; {len(rows)} model/direction/observable measurements')
    print(f'Saved CSV, JSON and 3 heatmaps to {output_dir}')
    for direction in ('T1','T2','T3'):
        r = max((r for r in rows if r['direction']==direction),key=lambda r:r['max_abs_shift_sigma'])
        print(f'{direction}: largest displacement {r["max_abs_shift_sigma"]:.6g} sigma; {r["model"]} / {r["observable"]}')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scenario',default='smallY_smallL')
    parser.add_argument('--input-dir',type=Path)
    parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args()
    input_dir=args.input_dir or ROOT/'Reports/output/full/fitted_sensitivity'/args.scenario
    output_dir=args.output_dir or ROOT/'Reports/output/full/normalized_fitted_sensitivity'/args.scenario
    analyze(input_dir,output_dir)

if __name__=='__main__':
    main()
