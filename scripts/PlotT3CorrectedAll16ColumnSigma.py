"""Column-wise z-score heatmaps of corrected LLSS stage ablations for all 16 T3 models.

Project path: scripts/PlotT3CorrectedAll16ColumnSigma.py

Each column (stage x observable) is standardised ACROSS MODELS:
    d_mj = O_full,mj - O_stage_off,mj
    z_mj = (d_mj - mean_models(d_.j)) / std_models(d_.j)
with sample standard deviation (ddof=1). No experimental NuFIT errors are used.
A zero-variance column has z=0 for all available models, explicitly flagged in CSV.
Missing models remain blank when --allow-partial is supplied.
This is relative dispersion across models, not the absolute stage effect or
statistical significance. No physics or fits are recomputed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = ROOT / 'output/studies/general_scalar_extended/corrected_llss_eft1_freeze'
OUTPUT_ROOT = ROOT / 'Reports/output/thesis_qualitative/corrected_all16_column_sigma'
FORMAT = 't3-independent-llss-U_E_W_S-ablation-v1'
MODELS = (
    'A_-4', 'A_-2', 'A_+0', 'B_-3', 'B_-1', 'B_+1',
    'C_-3', 'C_-1', 'C_+1', 'D_-2', 'D_+0', 'D_+2',
    'E_-4', 'E_-2', 'E_+0', 'E_+2',
)
STAGES = (
    ('UV OFF', '0111'),
    ('EFT1 OFF (E+W)', '1001'),
    ('SMEFT OFF', '1110'),
)
MASS = (
    ('delta_m21_sq_ev2', r'$\Delta m^2_{21}$'),
    ('delta_m3l_sq_ev2', r'$\Delta m^2_{3\ell}$'),
)
MIXING = (
    ('sin2_theta12', r'$\sin^2\theta_{12}$'),
    ('sin2_theta13', r'$\sin^2\theta_{13}$'),
    ('sin2_theta23', r'$\sin^2\theta_{23}$'),
)


def read(path: Path) -> dict:
    return json.loads(path.read_text(encoding='utf-8-sig'))


def obs_value(cases: dict, case: str, obs: str, model: str) -> float:
    if case not in cases or not isinstance(cases[case].get('prediction'), dict):
        raise ValueError(f'{model}: missing prediction for case {case}')
    value = float(cases[case]['prediction'][obs])
    if not math.isfinite(value):
        raise ValueError(f'{model}: nonfinite {obs} for case {case}')
    return value


def standardize(raw: np.ndarray, *, kind: str, observable_defs: tuple) -> tuple[np.ndarray, list[dict]]:
    z = np.full_like(raw, np.nan)
    stats = []
    for j in range(raw.shape[1]):
        valid = np.isfinite(raw[:, j])
        n = int(np.sum(valid))
        if n < 2:
            raise ValueError(f'{kind} column {j}: need at least 2 computed models; found {n}')
        x = raw[valid, j]
        mean = float(np.mean(x))
        sigma = float(np.std(x, ddof=1))
        if not math.isfinite(sigma):
            raise ValueError(f'{kind} column {j}: invalid standard deviation')
        constant = (sigma == 0.0)
        z[valid, j] = 0.0 if constant else (x - mean) / sigma
        nobs = len(observable_defs)
        stage, key = STAGES[j // nobs]
        obs = observable_defs[j % nobs][0]
        stats.append({
            'kind': kind, 'stage_off': stage, 'case_off': key,
            'observable': obs, 'n_models': n, 'mean_raw_shift': mean,
            'sample_standard_deviation_raw_shift': sigma,
            'constant_column': constant,
            'raw_shift_units': 'eV^2' if kind == 'mass' else 'absolute sin^2(theta)',
        })
    return z, stats


def plot_matrix(z: np.ndarray, defs: tuple, *, kind: str, dest: Path) -> None:
    nobs = len(defs)
    labels = [f'{stage}\n{label}' for stage, _ in STAGES for _, label in defs]
    lim = max(2.5, float(np.nanmax(np.abs(z))))
    fig, ax = plt.subplots(figsize=(max(11.5, z.shape[1] * 1.37 + 3), 9.3))
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#dddddd')
    image = ax.imshow(np.ma.masked_invalid(z), aspect='auto', cmap=cmap,
                      norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim),
                      interpolation='nearest')
    ax.set_xticks(np.arange(z.shape[1]), labels, fontsize=9)
    ax.set_yticks(np.arange(len(MODELS)), MODELS, fontsize=9)
    title = 'Mass-squared differences' if kind == 'mass' else 'Mixing angles'
    ax.set_title(f'{title}: column-wise z-scores across T3 models\n'
                 r'$z=(\Delta O-\overline{\Delta O}_{\rm column})/s_{\rm column}$; '
                 'each stage-observable column standardised separately',
                 fontsize=11, pad=18)
    ax.set_ylabel('T3 model')
    ax.set_xlabel('Disabled running stage and observable')
    ax.set_xticks(np.arange(-.5, z.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(MODELS), 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=.9)
    ax.tick_params(which='minor', bottom=False, left=False)
    for y in (3, 6, 9, 12):
        ax.axhline(y - .5, color='black', linewidth=1.2)
    for x in range(nobs, z.shape[1], nobs):
        ax.axvline(x - .5, color='black', linewidth=1.5)
    for i in range(z.shape[0]):
        for j in range(z.shape[1]):
            v = z[i, j]
            if math.isfinite(v):
                ax.text(j, i, f'{v:+.2f}', ha='center', va='center', fontsize=7.5,
                        color='white' if abs(v) > .65 * lim else 'black')
    fig.colorbar(image, ax=ax, fraction=.03, pad=.02).set_label('Within-column z-score (sample SD units)')
    fig.tight_layout()
    fig.savefig(dest, dpi=240, bbox_inches='tight')
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-id', default='corrected_llss16_v4')
    parser.add_argument('--allow-partial', action='store_true')
    args = parser.parse_args()
    src = SOURCE_ROOT / args.run_id
    summary_path = src / 'summary.json'
    if not summary_path.is_file():
        raise FileNotFoundError(f'Missing corrected campaign summary: {summary_path}')
    summary = read(summary_path)
    if summary.get('format') != FORMAT:
        raise ValueError(f'Unexpected summary format: {summary.get("format")!r}')

    raw = {
        'mass': np.full((16, len(STAGES) * len(MASS)), np.nan),
        'mixing': np.full((16, len(STAGES) * len(MIXING)), np.nan),
    }
    missing = []
    values = {}
    for i, model in enumerate(MODELS):
        tag = model.replace('+', 'p').replace('-', 'm')
        case_path = src / tag / 'cases.json'
        if not case_path.is_file():
            missing.append(model)
            continue
        payload = read(case_path)
        if payload.get('model') != model or payload.get('status') != 'Computed':
            missing.append(model)
            continue
        cases = payload.get('cases', {})
        for stage_i, (stage, off_key) in enumerate(STAGES):
            for kind, observables in (('mass', MASS), ('mixing', MIXING)):
                for obs_i, (name, _) in enumerate(observables):
                    full = obs_value(cases, '1111', name, model)
                    off = obs_value(cases, off_key, name, model)
                    diff = full - off
                    j = stage_i * len(observables) + obs_i
                    raw[kind][i, j] = diff
                    values[(kind, i, j)] = (full, off)
    if missing and not args.allow_partial:
        raise ValueError(f'Incomplete saved campaign: {16-len(missing)}/16, missing {missing}; '
                         'rerun with --allow-partial to leave those rows blank')

    output = OUTPUT_ROOT / args.run_id
    output.mkdir(parents=True, exist_ok=True)
    all_stats = []
    all_rows = []
    for kind, defs in (('mass', MASS), ('mixing', MIXING)):
        z, stats = standardize(raw[kind], kind=kind, observable_defs=defs)
        all_stats.extend(stats)
        plot_matrix(z, defs, kind=kind, dest=output / f'{kind}_column_sigma_heatmap.png')
        for i, model in enumerate(MODELS):
            for j in range(z.shape[1]):
                if not np.isfinite(z[i, j]):
                    continue
                stage_i, obs_i = divmod(j, len(defs))
                stage, key = STAGES[stage_i]
                full, off = values[(kind, i, j)]
                st = stats[j]
                all_rows.append({
                    'model': model, 'kind': kind, 'stage_off': stage,
                    'case_full': '1111', 'case_off': key,
                    'observable': defs[obs_i][0],
                    'full_value': full, 'stage_off_value': off,
                    'raw_shift_full_minus_off': float(raw[kind][i, j]),
                    'column_mean': st['mean_raw_shift'],
                    'column_sample_std': st['sample_standard_deviation_raw_shift'],
                    'within_column_z': float(z[i, j]),
                    'constant_column': st['constant_column'],
                })
    with (output / 'column_sigma_values.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)
    with (output / 'column_statistics.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(all_stats[0]))
        writer.writeheader()
        writer.writerows(all_stats)
    manifest = {
        'source_run': args.run_id, 'summary_sha256': hashlib.sha256(summary_path.read_bytes()).hexdigest(),
        'models': list(MODELS), 'missing_models': missing,
        'full_case': '1111', 'stage_off_cases': dict(STAGES),
        'raw_difference': 'O_full - O_stage_off (NOT percent)',
        'standardisation': 'For EACH (stage,observable) column separately across available models: '
                           'z=(difference-mean_column)/sample_std_column; ddof=1',
        'zero_variance_policy': 'Set available z-scores to 0; flag column in CSV',
        'scope': 'EFT1 E+W OFF is LLSS approximation, not complete EFT1 freeze',
        'physics_recomputed': False,
        'note': 'A z-score of zero means the column average, NOT zero running effect. '
                'This is a between-model standardisation, NOT NuFIT significance.',
    }
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
    print(f'Computed models: {16-len(missing)}/16; missing: {missing}')
    print('Wrote:', output)
    for name in ('mass_column_sigma_heatmap.png', 'mixing_column_sigma_heatmap.png',
                 'column_sigma_values.csv', 'column_statistics.csv', 'manifest.json'):
        print(' ', output / name)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
