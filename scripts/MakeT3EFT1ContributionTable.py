"""Read-only five-model EFT1 induced-vs-finite Weinberg contribution table.

Project path: scripts/MakeT3EFT1ContributionTable.py
Run: python -X utf8 -m scripts.MakeT3EFT1ContributionTable
Reads saved simultaneous-IVP cases; never runs RGEs or overwrites inputs.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
STUDY = ROOT / 'output' / 'studies' / 'general_scalar_extended'
SOURCES = {
    'A_-4': STUDY / 'experimental_simultaneous_eft1_seven_term' / 'simultaneous_A_m4_v2',
    'B_-3': STUDY / 'experimental_simultaneous_eft1_general_nonmixed' / 'simultaneous_nonmixed13_v1',
    'C_-3': STUDY / 'experimental_simultaneous_eft1_general_nonmixed' / 'simultaneous_nonmixed13_v1',
    'D_-2': STUDY / 'experimental_simultaneous_eft1_general_nonmixed' / 'simultaneous_nonmixed13_v1',
    'E_-4': STUDY / 'experimental_simultaneous_eft1_general_nonmixed' / 'simultaneous_nonmixed13_v1',
}
DEFAULT_OUT = STUDY / 'reports' / 'simultaneous_eft1_five_contributions_v1'
OBS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2')


def matrix(raw, label):
    if isinstance(raw, dict) and 'real' in raw and 'imag' in raw:
        z = np.asarray(raw['real'], dtype=float) + 1j * np.asarray(raw['imag'], dtype=float)
    else:
        z = np.asarray(raw, dtype=complex)
    if z.shape != (3, 3) or not np.isfinite(z).all():
        raise ValueError(f'{label}: expected finite complex 3x3 matrix')
    return z


def norm(a):
    return float(np.linalg.norm(a, 'fro'))


def pct(a, b):
    if b == 0:
        raise ValueError('Zero EFT1-off observable; relative percentage undefined')
    return float(100.0 * (a - b) / b)


def load(model):
    folder = SOURCES[model]
    path = folder / model.replace('-', 'm').replace('+', 'p') / 'cases.json'
    with path.open('r', encoding='utf-8-sig') as f:
        data = json.load(f)
    if data.get('status') != 'Computed' or data.get('model') != model:
        raise ValueError(f'{path}: unexpected model/status')
    cases = data.get('cases')
    if not isinstance(cases, dict) or len(cases) != 16:
        raise ValueError(f'{path}: expected 16 cases')
    for key in ('1111', '1001'):
        if key not in cases or not isinstance(cases[key].get('prediction'), dict):
            raise ValueError(f'{path}: missing {key} prediction')
    return path, cases['1111'], cases['1001']


def row_for(model):
    path, full, off = load(model)
    induced = matrix(full['c5_induced_wilson_log'], f'{model}/induced')
    finite = matrix(full['c5_joint_scalar_bubble'], f'{model}/full finite')
    frozen_finite = matrix(off['c5_joint_scalar_bubble'], f'{model}/frozen finite')
    c_full = matrix(full['c5_at_scalar_threshold'], f'{model}/full C5')
    c_off = matrix(off['c5_at_scalar_threshold'], f'{model}/frozen C5')
    n_induced, n_finite = norm(induced), norm(finite)
    if n_induced <= 0 or n_finite <= 0:
        raise ValueError(f'{model}: zero induced or finite norm; undefined fit')
    if norm(c_full - (induced + finite)) > 1e-9 * max(norm(c_full), 1e-30):
        raise ValueError(f'{model}: full C5 decomposition mismatch')
    if norm(c_off - frozen_finite) > 1e-9 * max(norm(c_off), 1e-30):
        raise ValueError(f'{model}: EFT1-off C5 is not the saved finite matching')
    k = np.vdot(finite, induced) / np.vdot(finite, finite)
    residual = norm(induced - k * finite) / n_induced
    ratio = n_induced / n_finite
    finite_drift = norm(finite - frozen_finite) / norm(frozen_finite)
    factor = abs(1 + k)**2
    values = {
        'model': model.replace('_', '(').replace('(', '(', 1) + ')',
        'source_cases_json': str(path.relative_to(ROOT)),
        'R_induced_over_finite': ratio,
        'k_real': float(k.real),
        'k_imag': float(k.imag),
        'epsilon_proportionality': residual,
        'finite_matching_drift_rel_full_vs_1001': finite_drift,
        'predicted_delta_m2_pct_proportional_full_finite': 100.0 * (factor - 1),
    }
    for name in OBS:
        full_val = float(full['prediction'][name])
        off_val = float(off['prediction'][name])
        if not (math.isfinite(full_val) and math.isfinite(off_val)):
            raise ValueError(f'{model}: nonfinite {name}')
        values[f'full_{name}'] = full_val
        values[f'eft1_off_{name}'] = off_val
        values[f'actual_{name}_pct'] = pct(full_val, off_val)
    return values


def fmt(x, digits=3):
    return f'{x:.{digits}g}'


def latex_num(x):
    return fmt(x, 4).replace('e+', r'\times 10^{').replace('e-', r'\times 10^{-') if False else fmt(x, 4)


def make_markdown(rows):
    lines = [
        '| Model | $R$ | $k$ | $\\epsilon$ | Finite drift | Predicted $\\delta(\\Delta m^2)$ | Actual solar | Actual atmospheric |',
        '|:--|--:|--:|--:|--:|--:|--:|--:|',
    ]
    for row in rows:
        k = f"{fmt(row['k_real'])}" + (f"{row['k_imag']:+.2g}i" if abs(row['k_imag']) > 1e-12 else '')
        lines.append('| ' + ' | '.join((
            row['model'], fmt(row['R_induced_over_finite']), k,
            fmt(row['epsilon_proportionality']),
            fmt(100 * row['finite_matching_drift_rel_full_vs_1001']) + r'%',
            fmt(row['predicted_delta_m2_pct_proportional_full_finite']) + r'%',
            fmt(row['actual_delta_m21_sq_ev2_pct']) + r'%',
            fmt(row['actual_delta_m3l_sq_ev2_pct']) + r'%',
        )) + ' |')
    lines.extend(['',
        '$R=\\|C_5^{\\mathrm{ind}}\\|_F/\\|C_5^{\\mathrm{fin,full}}\\|_F$; '
        '$k=\\mathrm{Tr}[(C_5^{\\mathrm{fin,full}})^\\dagger C_5^{\\mathrm{ind}}]/\\|C_5^{\\mathrm{fin,full}}\\|_F^2$; '
        '$\\epsilon=\\|C_5^{\\mathrm{ind}}-kC_5^{\\mathrm{fin,full}}\\|_F/\\|C_5^{\\mathrm{ind}}\\|_F$.',
        'Predicted correction is $100(|1+k|^2-1)\\%$. Actual corrections use `1111` versus `1001`.',
        'Finite drift (%) is $100\\|C_5^{\\mathrm{fin,full}}-C_5^{\\mathrm{fin,1001}}\\|_F/\\|C_5^{\\mathrm{fin,1001}}\\|_F$.',
        'The prediction is justified only for small proportionality residual and small finite drift; it neglects differences in SMEFT evolution.',
        'Results are exploratory and inherit the matching and scalar-potential approximations of the saved runs.',
    ])
    return '\n'.join(lines) + '\n'


def make_latex(rows):
    lines = [
        r'\begin{table}[htbp]',
        r'\centering',
        r'\small',
        r'\begin{tabular}{lrrrrrrr}',
        r'\hline',
        r'Model & $R$ & $\mathrm{Re}\,k$ & $\epsilon$ & drift (\%) & pred. (\%) & solar (\%) & atm. (\%) \\',
        r'\hline',
    ]
    for row in rows:
        parts = [row['model'].replace('_', ''), fmt(row['R_induced_over_finite']),
                 fmt(row['k_real']), fmt(row['epsilon_proportionality']),
                 fmt(100 * row['finite_matching_drift_rel_full_vs_1001']),
                 fmt(row['predicted_delta_m2_pct_proportional_full_finite']),
                 fmt(row['actual_delta_m21_sq_ev2_pct']),
                 fmt(row['actual_delta_m3l_sq_ev2_pct'])]
        lines.append(' & '.join(parts) + r' \\')
    lines += [r'\hline', r'\end{tabular}',
              r'\caption{Induced versus finite Weinberg coefficients at the common scalar threshold for five T3 benchmarks. '
              r'The predicted mass-squared shift assumes proportional coefficients; the actual solar and atmospheric shifts compare the simultaneous full (1111) and EFT1-frozen (1001) cases. '
              r'The drift column gives the relative change in finite matching between these cases. All results inherit the exploratory matching approximations.}',
              r'\label{tab:eft1-five-contributions}', r'\end{table}']
    return '\n'.join(lines) + '\n'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, default=DEFAULT_OUT)
    args = p.parse_args()
    rows = [row_for(model) for model in SOURCES]
    args.output.mkdir(parents=True, exist_ok=True)
    csv_path = args.output / 'eft1_five_contributions.csv'
    with csv_path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (args.output / 'eft1_five_contributions.md').write_text(make_markdown(rows), encoding='utf-8')
    (args.output / 'eft1_five_contributions.tex').write_text(make_latex(rows), encoding='utf-8')
    print(make_markdown(rows))
    print('Saved CSV, Markdown and LaTeX to:', args.output)


if __name__ == '__main__':
    main()
