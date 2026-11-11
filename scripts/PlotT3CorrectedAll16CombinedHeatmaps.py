"""Four all-16 corrected LLSS heatmaps: original-scale changes and column z-scores.

Project path: scripts/PlotT3CorrectedAll16CombinedHeatmaps.py

Read-only physics: load saved U/E/W/S cases; no fits or RGEs are rerun.
Figures: mass signed percentage changes (symmetric-log colour scale), mixing absolute sin^2 changes,
         mass columnwise z-scores, mixing columnwise z-scores.
Each column z-score uses RAW full-minus-off differences across available models.
EFT1 off on plots means E=0,W=0 (case 1001) in the LLSS approximation;
this is NOT the full physical EFT1 Wilson freeze. See manifest.
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
from matplotlib.colors import TwoSlopeNorm, SymLogNorm
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INPUT_ROOT = ROOT / 'output/studies/general_scalar_extended/corrected_llss_eft1_freeze'
OUTPUT_ROOT = ROOT / 'Reports/output/thesis_qualitative/corrected_all16_percentage_and_column_sigma'
EXPECTED_FORMAT = 't3-independent-llss-U_E_W_S-ablation-v1'
MODELS = (
    'A_-4', 'A_-2', 'A_+0', 'B_-3', 'B_-1', 'B_+1',
    'C_-3', 'C_-1', 'C_+1', 'D_-2', 'D_+0', 'D_+2',
    'E_-4', 'E_-2', 'E_+0', 'E_+2',
)
STAGES = (('UV off', '0111'), ('EFT1 off', '1001'), ('SMEFT off', '1110'))
MASS = (('delta_m21_sq_ev2', r'$\Delta m^2_{21}$'),
        ('delta_m3l_sq_ev2', r'$\Delta m^2_{3\ell}$'))
MIXING = (('sin2_theta12', r'$\sin^2\theta_{12}$'),
          ('sin2_theta13', r'$\sin^2\theta_{13}$'),
          ('sin2_theta23', r'$\sin^2\theta_{23}$'))


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding='utf-8-sig'))


def get_obs(cases: dict, key: str, name: str, model: str) -> float:
    if key not in cases or not isinstance(cases[key].get('prediction'), dict):
        raise ValueError(f'{model}: missing case {key} prediction')
    value = float(cases[key]['prediction'][name])
    if not math.isfinite(value):
        raise ValueError(f'{model}: nonfinite {name} for case {key}')
    return value


def column_standardize(raw: np.ndarray) -> tuple[np.ndarray, list[dict]]:
    z = np.full_like(raw, np.nan)
    stats = []
    for j in range(raw.shape[1]):
        good = np.isfinite(raw[:, j])
        n = int(good.sum())
        if n < 2:
            raise ValueError(f'Column {j} has only {n} models; need >= 2')
        x = raw[good, j]
        mean = float(x.mean())
        std = float(x.std(ddof=1))
        constant = std == 0.0
        z[good, j] = 0.0 if constant else (x - mean) / std
        stats.append({'n_models': n, 'raw_mean': mean, 'raw_sample_std': std,
                      'constant_column': constant})
    return z, stats


def compact_value(v: float, mode: str) -> str:
    """Short annotations that remain legible at thesis text width.

    Full precision is retained in heatmap_values.csv.
    """
    if mode == 'z':
        return f'{v:+.1f}'
    if mode == 'milli':
        if abs(v) < 0.005: return '0'
        return f'{v:+.1f}' if abs(v) < 100 else f'{v:+.0f}'
    av = abs(v)
    if av < 5e-12:
        return '0'
    if mode == 'percent':
        if av >= 1000 or av < 0.01:
            return f'{v:+.1e}'
        if av >= 100:
            return f'{v:+.0f}'
        if av >= 10:
            return f'{v:+.1f}'
        return f'{v:+.2f}'
    if mode == 'absolute':
        if av >= 0.1: return f'{v:+.2f}'
        if av >= 0.01: return f'{v:+.3f}'
        if av >= 0.001: return f'{v:+.3f}'
        return f'{v:+.1e}'
    if av >= 0.01:
        return f'{v:+.3f}'
    return f'{v:+.1e}'


def render(a: np.ndarray, observables: tuple, *, title: str, subtitle: str,
           color_label: str, mode: str, dest: Path,
           figure_width_in: float = 9.6, thesis_width_in: float = 6.3,
           thesis_body_pt: float = 11.0) -> None:
    """Produce graphics whose fonts match thesis body text AFTER LaTeX scaling.

    Example: a 9.6-inch source graphic included at 6.3 inches is resized by
    6.3/9.6; a source font of 16.8 pt becomes 11 pt in the PDF.
    """
    if thesis_width_in <= 0 or figure_width_in <= 0 or thesis_body_pt <= 0:
        raise ValueError('Plot width and thesis font size must be positive')
    scale = thesis_width_in / figure_width_in
    font = thesis_body_pt / scale
    nobs = len(observables)
    # Preserve the actual sin^2(theta) observable labels for all three stages.
    labels = [label for _stage, _ in STAGES for _, label in observables]
    finite = np.abs(a[np.isfinite(a)])
    limit = max(float(finite.max()) if finite.size else 0.0, 1e-13)
    if mode == 'z':
        limit = max(limit, 2.5)
    norm = (SymLogNorm(linthresh=1.0, linscale=1.0, vmin=-limit,
                       vmax=limit, base=10)
            if mode == 'percent' and limit > 1.0
            else TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit))

    # Keep the entire diagram compact: LaTeX scales image AND text together.
    # Thesis figure captions already explain the physics; duplicating a long
    # title/subtitle inside the image would reduce the useful plot area.
    fig, ax = plt.subplots(figsize=(figure_width_in, 9.2))
    fig.subplots_adjust(left=.105, right=.79, top=.91, bottom=.10)
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('#dddddd')
    im = ax.imshow(np.ma.masked_invalid(a), aspect='auto', interpolation='nearest',
                   cmap=cmap, norm=norm)
    ax.set_xticks(np.arange(a.shape[1]), labels, fontsize=font*.74 if nobs == 3 else font*.90)
    # Show each stage ONCE, centered over its observable columns.
    for stage_idx, (stage, _code) in enumerate(STAGES):
        xpos = (stage_idx*nobs + (nobs-1)/2) / (a.shape[1]-1)
        # imshow extends from -0.5 to ncols-0.5; thus use an axes fraction.
        xpos = (stage_idx+.5)/len(STAGES)
        ax.text(xpos, 1.023, stage, transform=ax.transAxes, ha='center',
                va='bottom', fontsize=font, fontweight='bold', clip_on=False)
    ax.set_yticks(np.arange(len(MODELS)), MODELS, fontsize=font)
    ax.tick_params(axis='both', which='major', length=0, pad=6)
    ax.set_xticks(np.arange(-.5, a.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(MODELS), 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=.9)
    ax.tick_params(which='minor', bottom=False, left=False)
    for y in (3, 6, 9, 12):
        ax.axhline(y-.5, color='black', linewidth=1.4)
    for x in (nobs, nobs*2):
        ax.axvline(x-.5, color='black', linewidth=1.8)
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            v = a[i, j]
            if not math.isfinite(v):
                continue
            label = compact_value(v, mode)
            rgba = cmap(norm(v))
            luminance = .2126*rgba[0] + .7152*rgba[1] + .0722*rgba[2]
            ax.text(j, i, label, ha='center', va='center', fontsize=font*.91,
                    fontweight='medium',
                    color='white' if luminance < .48 else 'black')

    # Put the colourbar into an explicitly reserved margin, avoiding tight_layout
    # changing the effective font scale. Its label is a short unit string.
    cax = fig.add_axes([.815, .20, .026, .65])
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(color_label, fontsize=font*.90, labelpad=9)
    cbar.ax.tick_params(labelsize=font*.76)
    # Avoid Mathematplotlib auto-generating very long tick labels for huge values.
    if mode == 'percent':
        from matplotlib.ticker import FuncFormatter
        cbar.ax.yaxis.set_major_formatter(FuncFormatter(
            lambda x, pos: ('0' if abs(x)<1e-12 else f'{x:.0e}'
                            if abs(x)>=1000 or abs(x)<.01 else f'{x:g}')))
    fig.savefig(dest, dpi=280, facecolor='white')
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run-id', default='corrected_llss16_v4')
    ap.add_argument('--figure-width-in', type=float, default=9.6,
                    help='Unscaled source figure width in inches')
    ap.add_argument('--thesis-width-in', type=float, default=6.3,
                    help='Width of includegraphics in thesis, in inches')
    ap.add_argument('--thesis-body-pt', type=float, default=11.0,
                    help='Body text font size to match after thesis resizing')
    ap.add_argument('--allow-partial', action='store_true',
                    help='Leave genuinely missing model rows blank instead of failing')
    args = ap.parse_args()
    base = INPUT_ROOT / args.run_id
    summary_path = base / 'summary.json'
    if not summary_path.is_file():
        raise FileNotFoundError(f'Missing summary: {summary_path}')
    summary = load(summary_path)
    if summary.get('format') != EXPECTED_FORMAT:
        raise ValueError(f'Unexpected corrected campaign format: {summary.get("format")!r}')

    raw = {kind: np.full((16, 3*len(defs)), np.nan)
           for kind, defs in (('mass', MASS), ('mixing', MIXING))}
    normal = {kind: np.full_like(arr, np.nan) for kind, arr in raw.items()}
    records = []
    absent = []
    for i, model in enumerate(MODELS):
        case_path = base / model.replace('+', 'p').replace('-', 'm') / 'cases.json'
        if not case_path.is_file():
            absent.append(model)
            continue
        result = load(case_path)
        if result.get('model') != model or result.get('status') != 'Computed':
            absent.append(model)
            continue
        cases = result.get('cases', {})
        for stage_i, (stage, off_case) in enumerate(STAGES):
            for kind, defs in (('mass', MASS), ('mixing', MIXING)):
                for obs_i, (obs, _) in enumerate(defs):
                    full = get_obs(cases, '1111', obs, model)
                    off = get_obs(cases, off_case, obs, model)
                    d = full-off
                    # Preserve prior corrected report convention: mass %,
                    # mixing absolute delta(sin^2 theta), not mixing %.
                    if kind == 'mass':
                        if off == 0.0:
                            raise ValueError(f'{model}/{stage}/{obs}: zero percent denominator')
                        displayed = 100.0*d/abs(off)
                        unit = 'signed percent'
                    else:
                        displayed = d
                        unit = 'absolute sin2 shift'
                    j = stage_i*len(defs)+obs_i
                    raw[kind][i,j] = d
                    normal[kind][i,j] = displayed
                    records.append({'model':model, 'kind':kind, 'stage_off':stage,
                                    'full_case':'1111', 'off_case':off_case,
                                    'observable':obs, 'full_value':full,
                                    'off_value':off, 'raw_full_minus_off':d,
                                    'normal_display_value':displayed,
                                    'normal_display_units':unit,
                                    'matrix_row':i, 'matrix_column':j})
    if absent and not args.allow_partial:
        raise ValueError(f'Incomplete corrected run: {16-len(absent)}/16; missing {absent}. '
                         'Pass --allow-partial to show missing cells as blank.')
    if not records:
        raise ValueError('No computed models')

    out = OUTPUT_ROOT / args.run_id
    out.mkdir(parents=True, exist_ok=True)
    filenames = []
    column_stats = []
    for kind, defs in (('mass', MASS), ('mixing', MIXING)):
        z, stats = column_standardize(raw[kind])
        for j, st in enumerate(stats):
            stage_index, obs_index = divmod(j, len(defs))
            column_stats.append({'kind':kind,'stage_off':STAGES[stage_index][0],
                                 'off_case':STAGES[stage_index][1],
                                 'observable':defs[obs_index][0],**st})
        for record in records:
            if record['kind']==kind:
                j = record['matrix_column']
                record['column_mean_raw'] = stats[j]['raw_mean']
                record['column_sample_std_raw'] = stats[j]['raw_sample_std']
                record['column_z_score'] = float(z[record['matrix_row'],j])
        if kind == 'mass':
            normal_name, z_name = 'mass_percentage_heatmap.png','mass_column_sigma_heatmap.png'
            title='Mass-squared differences — all 16 T3 models'
            subtitle='Signed percentage changes: 100 (full - stage off) / |stage off|; symmetric-log colours (linear within ±1%)'
            color_label='Signed change (%)'
            normal_mode='percent'
        else:
            normal_name, z_name = 'mixing_changes_heatmap.png','mixing_column_sigma_heatmap.png'
            title='Mixing observables — all 16 T3 models'
            subtitle=r'Absolute shift: $\sin^2\theta_{\rm full}-\sin^2\theta_{\rm stage\ off}$'
            color_label=r'$\Delta\sin^2\theta$'
            normal_mode='absolute'
        # Display raw differences in sin^2(theta); no 10^-3 unit conversion.
        plot_normal = normal[kind]
        render(plot_normal, defs, title=title, subtitle=subtitle,
               color_label=color_label, mode=normal_mode, dest=out/normal_name,
               figure_width_in=args.figure_width_in, thesis_width_in=args.thesis_width_in,
               thesis_body_pt=args.thesis_body_pt)
        render(z, defs, title=title+' (column-standardised)',
               subtitle=r'$z=(\Delta O-\overline{\Delta O}_{\rm column})/s_{\rm column}$',
               color_label='Within-column z-score (sample SD units)', mode='z',
               dest=out/z_name, figure_width_in=args.figure_width_in,
               thesis_width_in=args.thesis_width_in, thesis_body_pt=args.thesis_body_pt)
        filenames.extend((normal_name,z_name))

    with (out/'heatmap_values.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0]))
        writer.writeheader();writer.writerows(records)
    with (out/'column_statistics.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(column_stats[0]))
        writer.writeheader();writer.writerows(column_stats)
    (out/'manifest.json').write_text(json.dumps({
        'source_run':args.run_id,
        'summary_sha256':hashlib.sha256(summary_path.read_bytes()).hexdigest(),
        'models':list(MODELS),'missing_models':absent,
        'full_case':'1111','stage_off_cases':dict(STAGES),
        'display_labels': [s for s,_ in STAGES],
        'normal_mass':'100*(full-off)/abs(off), signed percentage',
        'mass_colour_normalization':'SymLogNorm(base=10, linthresh=1 percentage point, linscale=1) if maximum abs exceeds 1%; otherwise linear',
        'mixing_colour_normalization':'linear diverging colour map',
        'column_sigma_colour_normalization':'linear diverging colour map',
        'normal_mixing':'full-off, absolute difference in sin^2 theta (not percentage)',
        'column_sigma':'z=(raw difference-column mean)/sample std (ddof=1), per stage/observable column across available models',
        'column_zero_meaning':'average across models, not no running effect',
        'constant_column_policy':'z=0 if sample standard deviation is identically zero',
        'EFT1_off_definition':'case 1001: both renormalisable EFT1 parameter evolution and LLSS logarithm disabled; existing LLSS approximation, not complete EFT1 operator freeze',
        'SMEFT_off_definition':'case 1110: SM+Weinberg evolution disabled',
        'general_seven_term_trial_included':False,
        'physics_recomputed':False,'figures':filenames,
        'thesis_presentation':{'figure_width_in':args.figure_width_in,
             'includegraphics_width_in':args.thesis_width_in,
             'target_body_font_pt':args.thesis_body_pt,
             'internal_titles_omitted':'Thesis caption supplies the figure title',
             'mixing_figure_unit':'Raw delta(sin^2 theta), no multiplier',
             'cell_labels':'Compact annotations; full precision in CSV'},
    },indent=2)+'\n',encoding='utf-8')
    print(f'Models: {16-len(absent)}/16; missing: {absent}')
    print(f'Output: {out}')
    for name in filenames:
        print(' ', out/name)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
