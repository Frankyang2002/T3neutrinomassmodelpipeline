"""All-16 corrected LLSS stage-freeze shifts in NuFIT 1-sigma units.

Destination: scripts/PlotT3CorrectedAll16SigmaHeatmaps.py

No physics is recomputed. Uses corrected_llss_eft1_freeze/<run-id> case
1111 versus UV=0111, combined EFT1 E+W=1001, and SMEFT=1110.

z = (observable_full - observable_stage_off) / sigma_observable.
The sigma values come from configs/nufit_6_1_ic24_no.json, the same
Gaussianized NuFIT 6.1 NO fit target as the original numerical pipeline.
This is a standardized stage shift, NOT a model's deviation from NuFIT best fit
or a full correlated Delta-chi-square.
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
OUTPUT_ROOT = ROOT / 'Reports/output/thesis_qualitative/corrected_all16_sigma_heatmaps'
DEFAULT_TARGET = ROOT / 'configs/nufit_6_1_ic24_no.json'
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
    if case not in cases or 'prediction' not in cases[case] or obs not in cases[case]['prediction']:
        raise ValueError(f'{model}: missing {obs} for case {case}')
    x = float(cases[case]['prediction'][obs])
    if not math.isfinite(x):
        raise ValueError(f'{model}: nonfinite {obs} for case {case}')
    return x


def nice(z: float) -> str:
    if not math.isfinite(z):
        return 'N/A'
    if abs(z) < 1e-3 and z != 0:
        return f'{z:+.1e}'
    if abs(z) < 10:
        return f'{z:+.2f}'
    if abs(z) < 1e4:
        return f'{z:+.1f}'
    return f'{z:+.1e}'


def plot_heatmap(data: np.ndarray, obs: tuple, title: str, dest: Path) -> None:
    ncol = data.shape[1]
    limdata = np.abs(data[np.isfinite(data)])
    lim = max(float(np.max(limdata)) if limdata.size else 0., 1e-12)
    fig, ax = plt.subplots(figsize=(max(11, ncol * 1.45 + 3), 9.5))
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#dddddd')
    im = ax.imshow(np.ma.masked_invalid(data), aspect='auto', cmap=cmap,
                   norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim))
    labels = [f'{stage}\n{label}' for stage, _ in STAGES for _, label in obs]
    ax.set_xticks(range(ncol), labels, fontsize=9)
    ax.set_yticks(range(len(MODELS)), MODELS, fontsize=9)
    ax.set_title(title + '\nFull minus frozen, divided by NuFIT 1$\\sigma$ uncertainty', fontsize=12, pad=15)
    ax.set_ylabel('T3 model')
    ax.set_xlabel('Disabled running stage and observable')
    ax.set_xticks(np.arange(-.5, ncol, 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(MODELS), 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=.8)
    ax.tick_params(which='minor', bottom=False, left=False)
    for y in (3, 6, 9, 12):
        ax.axhline(y - .5, color='black', linewidth=1.3)
    for x in range(len(obs), ncol, len(obs)):
        ax.axvline(x - .5, color='black', linewidth=1.5)
    for i in range(len(MODELS)):
        for j in range(ncol):
            val = data[i, j]
            ax.text(j, i, nice(val), ha='center', va='center', fontsize=7.3,
                    color=('white' if math.isfinite(val) and abs(val) > .72 * lim else 'black'))
    fig.colorbar(im, ax=ax, fraction=.029, pad=.02).set_label('Signed shift (NuFIT 1$\\sigma$ units)')
    fig.tight_layout()
    fig.savefig(dest, dpi=240, bbox_inches='tight')
    plt.close(fig)


def plot_grouped_rms(rows: list[dict], dest: Path) -> list[dict]:
    allobs = MASS + MIXING
    names = [v[0] for v in allobs]
    labels = [v[1] for v in allobs]
    summary = []
    values = np.full((len(allobs), len(STAGES)), np.nan)
    for i, (name, _) in enumerate(allobs):
        for j, (stage, _) in enumerate(STAGES):
            x = np.asarray([float(r['shift_sigma']) for r in rows
                            if r['observable'] == name and r['stage_off'] == stage], dtype=float)
            if x.size:
                value = float(np.sqrt(np.mean(x*x)))
                values[i, j] = value
                summary.append({'observable': name, 'stage_off': stage,
                                'n_models': int(x.size), 'rms_sigma_shift': value,
                                'mean_signed_sigma_shift': float(np.mean(x)),
                                'max_abs_sigma_shift': float(np.max(np.abs(x)))})
    fig, ax = plt.subplots(figsize=(12, 6))
    positions = np.arange(len(allobs))
    width = .24
    for j, (stage, _) in enumerate(STAGES):
        ax.bar(positions + (j-1)*width, values[:, j], width=width, label=stage)
    ax.set_xticks(positions, labels, fontsize=11)
    ax.set_ylabel(r'RMS $|\Delta O|/\sigma_O$ across models')
    ax.set_title('Stage-grouped RMS shifts across the 16 T3 models')
    ax.legend(frameon=False)
    ax.grid(axis='y', alpha=.25)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(dest, dpi=240, bbox_inches='tight')
    plt.close(fig)
    return summary


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-id', default='corrected_llss16_v4')
    ap.add_argument('--target', type=Path, default=DEFAULT_TARGET)
    ap.add_argument('--allow-partial', action='store_true')
    args = ap.parse_args()
    src = SOURCE_ROOT / args.run_id
    spath = src / 'summary.json'
    if not spath.is_file():
        raise FileNotFoundError(f'Missing corrected run summary: {spath}')
    summary = read(spath)
    if summary.get('format') != FORMAT:
        raise ValueError('Input is not corrected independent LLSS campaign')
    target = args.target.resolve()
    if not target.is_file():
        raise FileNotFoundError(f'Missing experimental uncertainty file: {target}')
    t = read(target)
    if t.get('ordering', '').upper() != 'NO':
        raise ValueError('The corrected study uses NO; require NO target uncertainties')
    sigma = t.get('one_sigma_errors')
    if not isinstance(sigma, dict):
        raise ValueError('Target lacks one_sigma_errors')
    for key, _ in MASS + MIXING:
        if key not in sigma or not math.isfinite(float(sigma[key])) or float(sigma[key]) <= 0:
            raise ValueError(f'Invalid positive 1-sigma error for {key}')
    matrices = {'mass': np.full((16, 6), np.nan), 'mixing': np.full((16, 9), np.nan)}
    rows, absent = [], []
    for i, model in enumerate(MODELS):
        tag = model.replace('+', 'p').replace('-', 'm')
        path = src / tag / 'cases.json'
        if not path.is_file():
            absent.append(model)
            continue
        result = read(path)
        if result.get('model') != model or result.get('status') != 'Computed':
            absent.append(model)
            continue
        cases = result.get('cases', {})
        for j, (stage, offkey) in enumerate(STAGES):
            for kind, observables in (('mass', MASS), ('mixing', MIXING)):
                for k, (obs, _) in enumerate(observables):
                    full = obs_value(cases, '1111', obs, model)
                    off = obs_value(cases, offkey, obs, model)
                    sig = float(sigma[obs])
                    shift = (full - off) / sig
                    matrices[kind][i, j*len(observables)+k] = shift
                    rows.append({'model': model, 'stage_off': stage,
                                 'case_full': '1111', 'case_off': offkey,
                                 'observable': obs, 'full_value': full,
                                 'stage_off_value': off, 'sigma_observable': sig,
                                 'shift_sigma': shift})
    if absent and not args.allow_partial:
        raise ValueError(f'Incomplete corrected run: {16-len(absent)}/16, missing {absent}')
    if not rows:
        raise ValueError('No computed cases')
    dest = OUTPUT_ROOT / args.run_id
    dest.mkdir(parents=True, exist_ok=True)
    plot_heatmap(matrices['mass'], MASS, 'Mass-squared differences: 16-model stage shifts',
                 dest / 'mass_squared_sigma_heatmap.png')
    plot_heatmap(matrices['mixing'], MIXING, 'Mixing angles: 16-model stage shifts',
                 dest / 'mixing_sigma_heatmap.png')
    grouped = plot_grouped_rms(rows, dest / 'grouped_rms_sigma_shifts.png')
    with (dest / 'sigma_heatmap_values.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (dest / 'grouped_rms_sigma_shifts.csv').open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(grouped[0]))
        writer.writeheader()
        writer.writerows(grouped)
    (dest / 'manifest.json').write_text(json.dumps({
        'source_run': args.run_id, 'source_format': summary['format'],
        'target': str(target), 'target_sha256': hashlib.sha256(target.read_bytes()).hexdigest(),
        'target_description': t.get('source', ''), 'sigma_values': {k: float(sigma[k]) for k, _ in MASS+MIXING},
        'models': list(MODELS), 'missing_models': absent,
        'full_case': '1111', 'stage_off_cases': dict(STAGES),
        'quantity': 'z=(O_full-O_frozen)/sigma_observable',
        'grouped_quantity': 'RMS(z) across computed models; not a statistical error bar',
        'eft1_scope': 'E+W off: renormalisable freeze plus remove LLSS log, not complete EFT1 Wilson evolution',
        'not_a_correlated_delta_chi_square': True,
        'physics_recomputed': False
    }, indent=2)+'\n', encoding='utf-8')
    print(f'Computed models: {16-len(absent)}/16; missing: {absent}')
    for name in ('mass_squared_sigma_heatmap.png', 'mixing_sigma_heatmap.png',
                 'grouped_rms_sigma_shifts.png', 'sigma_heatmap_values.csv'):
        print(dest / name)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
