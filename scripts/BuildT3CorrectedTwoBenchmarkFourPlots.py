"""Four original-style stage-ablation heatmaps for A(-4) and D(-2).

Project path: scripts/BuildT3CorrectedTwoBenchmarkFourPlots.py

Reads ONLY the saved corrected LLSS 16-model campaign. Produces two plots
(mass-squared splittings and mixing angles) per selected MODEL benchmark.
This is NOT the historical two-Yukawa-DONOR-across-16-models experiment.
No physics is rerun or altered. General seven-term A(-4) trial is NOT used.

Stage comparisons: UV 1111 vs 0111, EFT1 1111 vs 1001, SMEFT 1111 vs 1110.
Mass: 100*(full-off)/abs(off) in percent; mixing: full-off absolute.
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
INPUT_BASE = ROOT / 'output/studies/general_scalar_extended/corrected_llss_eft1_freeze'
OUTPUT_BASE = ROOT / 'Reports/output/thesis_qualitative/corrected_two_benchmark_four_plots'
SCHEMA = 't3-independent-llss-U_E_W_S-ablation-v1'
MODELS = ('A_-4', 'D_-2')
STAGES = (('UV', '0111'), ('EFT1 (E+W off)', '1001'), ('SMEFT', '1110'))
GROUPS = {
    'masssq': (
        ('delta_m21_sq_ev2', r'$\Delta m_{21}^2$'),
        ('delta_m3l_sq_ev2', r'$\Delta m_{3\ell}^2$'),
    ),
    'mixing': (
        ('sin2_theta12', r'$\sin^2\theta_{12}$'),
        ('sin2_theta13', r'$\sin^2\theta_{13}$'),
        ('sin2_theta23', r'$\sin^2\theta_{23}$'),
    ),
}

def read_json(path: Path) -> dict:
    if not path.is_file():
        raise FileNotFoundError(f'Missing corrected LLSS data: {path}')
    with path.open(encoding='utf-8-sig') as fp:
        obj = json.load(fp)
    if not isinstance(obj, dict):
        raise ValueError(f'Invalid JSON object: {path}')
    return obj


def measure(record: dict, case: str, field: str, model: str) -> float:
    if case not in record:
        raise ValueError(f'{model}: missing case {case}')
    pred = record[case].get('prediction')
    if not isinstance(pred, dict) or field not in pred:
        raise ValueError(f'{model}: {case} missing {field}')
    value = float(pred[field])
    if not math.isfinite(value):
        raise ValueError(f'{model}: {case} {field} nonfinite')
    return value


def data_rows(model: str, cases: dict) -> list[dict]:
    rows = []
    for stage, off in STAGES:
        for group, metrics in GROUPS.items():
            for field, _ in metrics:
                full = measure(cases, '1111', field, model)
                frozen = measure(cases, off, field, model)
                if group == 'masssq':
                    if frozen == 0:
                        raise ValueError(f'{model}: {stage} {field}: frozen denominator is zero')
                    delta = 100.0 * (full-frozen) / abs(frozen)
                    unit = 'percent'
                else:
                    delta = full-frozen
                    unit = 'absolute_sin2_difference'
                rows.append(dict(model=model, benchmark_type='model_specific_fitted_point',
                                 stage=stage, reference_case='1111', frozen_case=off,
                                 observable=field, full_value=full, off_value=frozen,
                                 change=delta, unit=unit))
    return rows


def render(rows: list[dict], model: str, group: str, destination: Path) -> None:
    metrics = GROUPS[group]
    columns = [(stage, metric, label) for stage, _ in STAGES for metric, label in metrics]
    lookup = {(row['stage'], row['observable']): row['change'] for row in rows}
    vals = np.asarray([[lookup[stage, metric] for stage, metric, _ in columns]], dtype=float)
    max_abs = max(float(np.max(abs(vals))), 1e-15)
    norm = TwoSlopeNorm(vmin=-max_abs, vcenter=0, vmax=max_abs)
    fig, ax = plt.subplots(figsize=(max(10, len(columns)*1.35), 3.8))
    im = ax.imshow(vals, cmap='RdBu_r', norm=norm, aspect='auto')
    ax.set_xticks(range(len(columns)), [f'{stage}\n{label}' for stage, _, label in columns], fontsize=9)
    ax.set_yticks([0], [model], fontsize=11)
    ax.set_xticks(np.arange(-0.5, len(columns), 1), minor=True)
    ax.set_yticks([-0.5, 0.5], minor=True)
    ax.grid(which='minor', color='black', linewidth=.6)
    ax.tick_params(which='minor', bottom=False, left=False)
    for x in range(len(metrics), len(columns), len(metrics)):
        ax.axvline(x-.5, color='black', linewidth=1.6)
    for j, v in enumerate(vals[0]):
        label = f'{v:+.3g}%' if group == 'masssq' else f'{v:+.4g}'
        ax.text(j, 0, label, ha='center', va='center', fontsize=9,
                color='white' if abs(v) > .60*max_abs else 'black')
    if group == 'masssq':
        ax.set_title(f'{model} | mass-squared differences | corrected LLSS stage ablations', pad=12)
        unit = r'$100\,(O_{\mathrm{full}}-O_{\mathrm{off}})/|O_{\mathrm{off}}|$ (%)'
    else:
        ax.set_title(f'{model} | mixing observables | corrected LLSS stage ablations', pad=12)
        unit = r'$\sin^2\theta_{\mathrm{full}}-\sin^2\theta_{\mathrm{off}}$ (absolute)'
    fig.colorbar(im, ax=ax, fraction=.035, pad=.03).set_label(unit)
    fig.text(.5, .02,
             'EFT1 off: renormalisable EFT1 frozen + LLSS log removed; finite scalar matching retained. Exploratory approximation.',
             ha='center', fontsize=8)
    fig.subplots_adjust(bottom=.26)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=260, bbox_inches='tight')
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-id', default='corrected_llss16_v4')
    ap.add_argument('--output', type=Path)
    args = ap.parse_args()
    source = INPUT_BASE / args.run_id
    summary_path = source / 'summary.json'
    summary = read_json(summary_path)
    if summary.get('format') != SCHEMA:
        raise ValueError('Refusing non-corrected or one-scale-anchored results: unexpected schema')
    out = (args.output if args.output is not None else OUTPUT_BASE/args.run_id).resolve()
    all_rows = []
    fingerprints = {'summary.json': hashlib.sha256(summary_path.read_bytes()).hexdigest()}
    planned = []
    for model in MODELS:
        name = model.replace('+', 'p').replace('-', 'm')
        path = source / name / 'cases.json'
        record = read_json(path)
        if record.get('status') != 'Computed' or record.get('model') != model:
            raise ValueError(f'{model}: corrected calculation not computed')
        cases = record.get('cases')
        if not isinstance(cases, dict) or len(cases) != 16:
            raise ValueError(f'{model}: expected sixteen stored cases')
        fingerprints[f'{name}/cases.json'] = hashlib.sha256(path.read_bytes()).hexdigest()
        rows = data_rows(model, cases)
        all_rows.extend(rows)
        for group in ('masssq', 'mixing'):
            planned.append((rows, model, group, out/model/f'oscillation_{group}_heatmap.png'))
    out.mkdir(parents=True, exist_ok=True)
    with (out/'comparison_values.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(all_rows[0]));writer.writeheader();writer.writerows(all_rows)
    for rows, model, group, path in planned:
        render(rows, model, group, path)
        print('Figure:', path)
    meta = {
        'source_run': args.run_id, 'source_root': str(source), 'sha256': fingerprints,
        'benchmark_models': list(MODELS),
        'benchmark_definition': 'two model-specific saved fitted configurations, NOT two donors across sixteen models',
        'comparisons': {stage: ['1111', off] for stage, off in STAGES},
        'mass_metric': '100*(full-off)/abs(off)',
        'mixing_metric': 'full-off (absolute sin-squared difference)',
        'approximations': ['LLSS scalar-pair leading log and finite bubble only',
                           'no independent complete EFT1 Wilson mixing',
                           'no certified general threshold matching'],
        'general_seven_term_trial_included': False,
        'recomputed_physics': False,
    }
    (out/'manifest.json').write_text(json.dumps(meta, indent=2)+'\n', encoding='utf-8')
    print('Four plots created:', out)
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
