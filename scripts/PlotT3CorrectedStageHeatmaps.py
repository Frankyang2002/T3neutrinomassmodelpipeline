"""Plot three corrected LLSS stage comparisons (UV, combined EFT1, SMEFT).

Project path: scripts/PlotT3CorrectedStageHeatmaps.py

Requires a completed or partial run made by RunT3CorrectedEFT1Freeze.py.
It never recomputes physics, rewrites saved results, fills missing models with
zero, or interprets the combined EFT1-off comparison as certified full Wilson
transport. The EFT1 column compares 1111 with 1001 (E=0 and W=0).

Example (PowerShell, from project root):
  python -X utf8 -m scripts.PlotT3CorrectedStageHeatmaps --run-id corrected_llss16_v4
  python -X utf8 -m scripts.PlotT3CorrectedStageHeatmaps --run-id corrected_llss_A4_diagnostic_v1 --allow-partial
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INPUT_BASE = ROOT / 'output/studies/general_scalar_extended/corrected_llss_eft1_freeze'
OUTPUT_BASE = ROOT / 'Reports/output/thesis_qualitative/corrected_three_stage_heatmaps'
MODELS = tuple(f'{family}_{a:+d}' for family, values in {
    'A': (-4, -2, 0), 'B': (-3, -1, 1), 'C': (-3, -1, 1),
    'D': (-2, 0, 2), 'E': (-4, -2, 0, 2),
}.items() for a in values)
STAGES = (
    ('UV', '0111'),
    ('EFT1 (E+W off)', '1001'),
    ('SMEFT', '1110'),
)
OBSERVABLES = (
    ('delta_m21_sq_ev2', r'$\Delta m^2_{21}$', 'mass'),
    ('delta_m3l_sq_ev2', r'$\Delta m^2_{3\ell}$', 'mass'),
    ('sin2_theta12', r'$\sin^2\theta_{12}$', 'mixing'),
    ('sin2_theta13', r'$\sin^2\theta_{13}$', 'mixing'),
    ('sin2_theta23', r'$\sin^2\theta_{23}$', 'mixing'),
)
REFERENCE = '1111'
EXPECTED_FORMAT = 't3-independent-llss-U_E_W_S-ablation-v1'


def load(path: Path) -> dict:
    with path.open('r', encoding='utf-8-sig') as stream:
        return json.load(stream)


def model_tag(model: str) -> str:
    return model.replace('+', 'p').replace('-', 'm')


def prediction(cases: dict, case: str, field: str, model: str) -> float:
    if case not in cases:
        raise ValueError(f'{model}: required case {case} is missing')
    entry = cases[case]
    if not isinstance(entry, dict) or not isinstance(entry.get('prediction'), dict):
        raise ValueError(f'{model}: case {case} lacks a prediction')
    value = float(entry['prediction'][field])
    if not math.isfinite(value):
        raise ValueError(f'{model}: non-finite {field} in case {case}')
    return value


def colour_limit(data: np.ndarray) -> float:
    valid = np.abs(data[np.isfinite(data)])
    return float(max(np.max(valid), 1e-15)) if valid.size else 1.0


def pretty_value(v: float, kind: str) -> str:
    if not math.isfinite(v):
        return 'N/A'
    if kind == 'mass':
        return f'{v:+.2f}%' if abs(v) >= 0.01 else f'{v:+.1e}%'
    if abs(v) < 0.0001:
        return f'{v:+.1e}'
    return f'{v:+.4f}'


def draw_matrix(ax, matrix: np.ndarray, models: tuple[str, ...], heading: str,
                unit: str, kind: str):
    import matplotlib.colors as colors
    from matplotlib.cm import ScalarMappable
    nrow, ncol = matrix.shape
    limit = colour_limit(matrix)
    # Linear scaling preserves the meaning of zero and signs. Different observables
    # use independent scales; numeric labels give exact values.
    norm = colors.Normalize(vmin=-limit, vmax=limit)
    cmap = __import__('matplotlib').colormaps['RdBu_r'].copy()
    cmap.set_bad('#dedede')
    ax.imshow(np.ma.masked_invalid(matrix), aspect='auto', cmap=cmap, norm=norm,
              interpolation='nearest')
    ax.set_yticks(np.arange(nrow), labels=[m.replace('_', '(') + ')' for m in models])
    ax.set_xticks(np.arange(ncol), labels=[s for s, _ in STAGES])
    ax.tick_params(axis='x', labelrotation=25, labelsize=8)
    ax.tick_params(axis='y', labelsize=8)
    ax.set_title(heading + f'\n{unit}', fontsize=12)
    ax.set_xticks(np.arange(-.5, ncol, 1), minor=True)
    ax.set_yticks(np.arange(-.5, nrow, 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=.9)
    ax.tick_params(which='minor', bottom=False, left=False)
    for i in range(nrow):
        for j in range(ncol):
            v = matrix[i, j]
            label = pretty_value(v, kind)
            foreground = 'white' if math.isfinite(v) and abs(v) >= limit * .73 else 'black'
            ax.text(j, i, label, ha='center', va='center', fontsize=7.0,
                    color=foreground)
    return ScalarMappable(norm=norm, cmap=cmap)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id', required=True, help='Corrected LLSS campaign ID')
    p.add_argument('--allow-partial', action='store_true',
                   help='Plot available models only, showing missing ones as N/A')
    p.add_argument('--output-dir', type=Path, default=None)
    args = p.parse_args(argv)
    if not args.run_id or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-' for c in args.run_id):
        p.error('Invalid --run-id')
    folder = INPUT_BASE / args.run_id
    summary_path = folder / 'summary.json'
    if not summary_path.is_file():
        raise FileNotFoundError(f'No corrected LLSS summary: {summary_path}')
    summary = load(summary_path)
    if summary.get('format') != EXPECTED_FORMAT:
        raise ValueError('Refusing non-corrected / one-scale-anchored source data: '
                         + repr(summary.get('format')))
    if summary.get('one_scale_C5_used_as_W_on_reference') is not False:
        raise ValueError('Source metadata does not certify the corrected W construction')
    out = args.output_dir if args.output_dir else OUTPUT_BASE / args.run_id
    matrices = {name: np.full((len(MODELS), len(STAGES)), np.nan, dtype=float)
                for name, _, _ in OBSERVABLES}
    rows = []
    absent = []
    for i, model in enumerate(MODELS):
        path = folder / model_tag(model) / 'cases.json'
        if not path.is_file():
            absent.append(model)
            continue
        record = load(path)
        if record.get('model') != model or record.get('status') != 'Computed':
            raise ValueError(f'{model}: cases.json is not a computed model record')
        if record.get('matching_approximation') != 'independent_LLSS_plus_bubble':
            raise ValueError(f'{model}: cases lack corrected matching provenance')
        cases = record.get('cases')
        if not isinstance(cases, dict):
            raise ValueError(f'{model}: cases field is invalid')
        for j, (stage, off_case) in enumerate(STAGES):
            for name, label, kind in OBSERVABLES:
                full = prediction(cases, REFERENCE, name, model)
                frozen = prediction(cases, off_case, name, model)
                absolute = full - frozen
                if kind == 'mass':
                    if abs(frozen) < 1e-30:
                        raise ValueError(f'{model}: near-zero frozen mass-squared denominator for {name}')
                    plotted = 100.0 * absolute / abs(frozen)
                    unit = '% of |frozen|'
                else:
                    plotted = absolute
                    unit = 'absolute difference'
                matrices[name][i, j] = plotted
                rows.append({'model': model, 'stage': stage, 'full_case': REFERENCE,
                             'frozen_case': off_case, 'observable': name,
                             'full_value': f'{full:.17g}', 'frozen_value': f'{frozen:.17g}',
                             'full_minus_frozen': f'{absolute:.17g}',
                             'heatmap_value': f'{plotted:.17g}', 'heatmap_unit': unit})
    if absent and not args.allow_partial:
        raise ValueError(f'{len(absent)} models lack computed corrected results: '
                         + ', '.join(absent) + '. Rerun with --allow-partial '
                         'to produce a visibly incomplete figure.')
    if len(absent) == len(MODELS):
        raise ValueError('No computed models found')
    out.mkdir(parents=True, exist_ok=True)
    with (out / 'heatmap_values.csv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'savefig.dpi': 200})
    for name, label, kind in OBSERVABLES:
        fig, ax = plt.subplots(figsize=(11, 9.0), constrained_layout=True)
        unit = ('100 × (full − off) / |off| [%]' if kind == 'mass'
                else 'full − off [absolute sin²θ units]')
        draw_matrix(ax, matrices[name], MODELS, label, unit, kind)
        fig.suptitle('Corrected LLSS stage ablations: ' + args.run_id,
                     fontsize=13, fontweight='bold')
        fig.savefig(out / f'{name}.png', bbox_inches='tight')
        plt.close(fig)
    fig, axs = plt.subplots(1, 5, figsize=(25, 10), constrained_layout=True)
    for ax, (name, label, kind) in zip(axs, OBSERVABLES):
        unit = '% versus frozen' if kind == 'mass' else 'absolute Δsin²θ'
        draw_matrix(ax, matrices[name], MODELS, label, unit, kind)
    fig.suptitle('Corrected LLSS: UV off (0111), combined EFT1 off (1001), SMEFT off (1110)\n'
                 'All compared with full (1111); missing models are N/A', fontsize=14, fontweight='bold')
    fig.savefig(out / 'combined_5observable_heatmaps.png', bbox_inches='tight')
    plt.close(fig)
    manifest = {
        'source_run_id': args.run_id,
        'source_summary': str(summary_path),
        'source_format': summary['format'],
        'baseline_case': REFERENCE,
        'comparison_cases': {stage: case for stage, case in STAGES},
        'mass_definition': '100 * (full - frozen) / abs(frozen)',
        'mixing_definition': 'full sin2theta - frozen sin2theta (absolute)',
        'missing_models': absent,
        'computed_models': len(MODELS) - len(absent),
        'interpretation': 'Exploratory LLSS+scalar-bubble, not full physical one-loop matching',
        'EFT1_off_scope': 'E=0 freezes renormalisable EFT1 state; W=0 omits conditional LLSS-induced leading log; finite scalar threshold bubble retained',
        'EFT1_off_caveat': 'Combined 1001 is an exploratory proxy, not independent freezing of every EFT1 Wilson coefficient or complete physical matching',
        'W_scope': 'W is not a separate plotted stage; it is set to zero together with E in the EFT1-off column',
        'heatmap_values_csv': str(out / 'heatmap_values.csv'),
    }
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    print(f'Plotted {len(MODELS)-len(absent)}/16 corrected LLSS models; missing {len(absent)}')
    if absent:
        print('N/A models:', ', '.join(absent))
    print('Images and exact plotted data:', out)
    for model in ('A_-4',):
        if model in absent:
            continue
        idx = MODELS.index(model)
        jj = [k for k, (stage, _) in enumerate(STAGES) if stage.startswith('EFT1')][0]
        print('A_-4 combined EFT1-off (1111 - 1001) mixing-angle shifts:',
              {name: float(matrices[name][idx, jj]) for name, _, kind in OBSERVABLES
               if kind == 'mixing'})
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
