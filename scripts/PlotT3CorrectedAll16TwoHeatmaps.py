"""Two original-style 16-model heatmaps from corrected LLSS ablation cases.

Destination: scripts/PlotT3CorrectedAll16TwoHeatmaps.py

This script reads saved cases only. Full reference is 1111; stage-off cases
are UV=0111, EFT1 E+W=1001, SMEFT=1110. EFT1 is the existing LLSS-only
approximation, not a complete Wilson-RGE stage freeze.
"""
from __future__ import annotations

import argparse
import csv
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
PLOT_ROOT = ROOT / 'Reports/output/thesis_qualitative/corrected_all16_two_heatmaps'
MODELS = (
    'A_-4', 'A_-2', 'A_+0', 'B_-3', 'B_-1', 'B_+1',
    'C_-3', 'C_-1', 'C_+1', 'D_-2', 'D_+0', 'D_+2',
    'E_-4', 'E_-2', 'E_+0', 'E_+2',
)
STAGES = (('UV OFF', '0111'), ('EFT1 OFF (E+W)', '1001'), ('SMEFT OFF', '1110'))
MASS = (
    ('delta_m21_sq_ev2', r'$\Delta m^2_{21}$'),
    ('delta_m3l_sq_ev2', r'$\Delta m^2_{3\ell}$'),
)
MIXING = (
    ('sin2_theta12', r'$\sin^2\theta_{12}$'),
    ('sin2_theta13', r'$\sin^2\theta_{13}$'),
    ('sin2_theta23', r'$\sin^2\theta_{23}$'),
)
FORMAT = 't3-independent-llss-U_E_W_S-ablation-v1'


def read(path: Path) -> dict:
    return json.loads(path.read_text(encoding='utf-8-sig'))


def prediction(cases: dict, key: str, obs: str, model: str) -> float:
    if key not in cases or 'prediction' not in cases[key]:
        raise ValueError(f'{model}: missing {key} prediction')
    result = float(cases[key]['prediction'][obs])
    if not math.isfinite(result):
        raise ValueError(f'{model}: nonfinite {obs} for {key}')
    return result


def graph(data: np.ndarray, kinds: tuple, title: str, unit: str, path: Path):
    fig, ax = plt.subplots(figsize=(max(10, data.shape[1]*1.30+3.2), 9.4))
    nonzero = np.abs(data[np.isfinite(data)])
    lim = max(float(np.max(nonzero)) if nonzero.size else 0., 1e-14)
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#dedede')
    im = ax.imshow(np.ma.masked_invalid(data), aspect='auto', cmap=cmap,
                   norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim))
    labels = [f'{stage}\n{label}' for stage, _ in STAGES for _, label in kinds]
    ax.set_xticks(np.arange(data.shape[1]), labels, fontsize=9)
    ax.set_yticks(np.arange(len(MODELS)), MODELS, fontsize=9)
    ax.set_title(title + '\n' + unit, fontsize=13, pad=14)
    ax.set_ylabel('T3 model')
    ax.set_xlabel('Disabled running stage / observable')
    ax.set_xticks(np.arange(-.5, data.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(MODELS), 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=.7)
    ax.tick_params(which='minor', bottom=False, left=False)
    for y in (3,6,9,12):
        ax.axhline(y-.5, color='black', linewidth=1.3)
    for x in (len(kinds), 2*len(kinds)):
        ax.axvline(x-.5, color='black', linewidth=1.5)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data[i,j]
            if not math.isfinite(val):
                label = 'N/A'
            elif unit.startswith('Signed percent'):
                label = f'{val:+.2f}' if abs(val) >= 0.01 else f'{val:+.1e}'
            else:
                label = f'{val:+.4f}' if abs(val) >= .0001 else f'{val:+.1e}'
            ax.text(j,i,label,ha='center',va='center',fontsize=7.3,
                    color='white' if math.isfinite(val) and abs(val) > .72*lim else 'black')
    fig.colorbar(im, ax=ax, fraction=.027, pad=.02).set_label(unit)
    fig.tight_layout()
    fig.savefig(path, dpi=240, bbox_inches='tight')
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-id', default='corrected_llss16_v4')
    parser.add_argument('--allow-partial', action='store_true')
    args = parser.parse_args()
    src = SOURCE_ROOT / args.run_id
    summary_path = src / 'summary.json'
    if not summary_path.is_file():
        raise FileNotFoundError(f'Missing corrected run: {summary_path}')
    summary = read(summary_path)
    if summary.get('format') != FORMAT:
        raise ValueError(f'Not a corrected LLSS dataset: {summary.get("format")!r}')
    matrices = {'mass': np.full((16,6),np.nan), 'mixing':np.full((16,9),np.nan)}
    rows = []
    absent = []
    for i, model in enumerate(MODELS):
        tag = model.replace('+','p').replace('-','m')
        file = src / tag / 'cases.json'
        if not file.is_file():
            absent.append(model)
            continue
        record = read(file)
        if record.get('model') != model or record.get('status') != 'Computed':
            absent.append(model)
            continue
        cases = record['cases']
        for j,(stage, off_key) in enumerate(STAGES):
            for kind, observations in (('mass',MASS),('mixing',MIXING)):
                for k,(obs,_) in enumerate(observations):
                    full = prediction(cases,'1111',obs,model)
                    frozen = prediction(cases,off_key,obs,model)
                    if kind=='mass':
                        if frozen == 0:
                            raise ValueError(f'{model} {obs}: zero frozen denominator')
                        difference = 100*(full-frozen)/abs(frozen)
                    else:
                        difference = full-frozen
                    matrices[kind][i,j*len(observations)+k]=difference
                    rows.append({'model':model,'stage_off':stage,'case_full':'1111',
                                 'case_off':off_key,'observable':obs,
                                 'full_value':full,'stage_off_value':frozen,
                                 'change':difference,
                                 'change_unit':'percent' if kind=='mass' else 'absolute_sin2'})
    if absent and not args.allow_partial:
        raise ValueError(f'Incomplete dataset: {len(MODELS)-len(absent)}/16 computed; missing {absent}. '
                         'Use --allow-partial only for diagnostic N/A cells.')
    if not rows:
        raise ValueError('No computed cases available')
    dest = PLOT_ROOT / args.run_id
    dest.mkdir(parents=True,exist_ok=True)
    graph(matrices['mass'], MASS, 'Mass-squared differences: all 16 T3 models',
          'Signed percent change (%)', dest/'mass_squared_differences.png')
    graph(matrices['mixing'], MIXING, 'Mixing observables: all 16 T3 models',
          r'Absolute $\Delta\sin^2\theta$', dest/'mixing_angles.png')
    with (dest/'heatmap_values.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    (dest/'manifest.json').write_text(json.dumps({
        'source_run':args.run_id,'source_format':summary['format'],
        'models':list(MODELS),'missing_models':absent,
        'reference_case':'1111','stage_off_cases':dict(STAGES),
        'mass_comparison':'100*(full-off)/abs(off)',
        'mixing_comparison':'full-off (absolute sin^2 angle)',
        'eft1_scope':'E+W off: renormalizable freeze + omit LLSS leading log; not full Wilson running',
        'physics_recomputed':False
    },indent=2)+'\n',encoding='utf-8')
    print(f'Models: {16-len(absent)}/16; missing: {absent}')
    print('Mass:',dest/'mass_squared_differences.png')
    print('Mixing:',dest/'mixing_angles.png')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
