"""Plot existing fitted T3 gamma sensitivities, without recomputing physics.

Run from the repository root:
    python -m Numerical.plotting.PlotFittedSensitivityFigures

The plotted quantity is the change under the prescribed UV gamma deformation
relative to gamma=1. It is NOT a running-vs-no-running quantum correction.
Non-successful points are masked, so curves never bridge failed gamma points.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / 'Reports' / 'output' / 'full' / 'fitted_sensitivity' / 'smallY_smallL'
OBS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12', 'sin2_theta13', 'sin2_theta23')
LABELS = {
    'delta_m21_sq_ev2': r'$\Delta m^2_{21}$',
    'delta_m3l_sq_ev2': r'$\Delta m^2_{3\ell}$',
    'sin2_theta12': r'$\sin^2\theta_{12}$',
    'sin2_theta13': r'$\sin^2\theta_{13}$',
    'sin2_theta23': r'$\sin^2\theta_{23}$',
}
DIRECTIONS = ('T1', 'T2', 'T3')
STYLES = {'T1': ('#286f9b', '-'), 'T2': ('#b75c39', '--'), 'T3': ('#41836d', ':')}
RE_MODEL = re.compile(r'^T3_dS1_(\d+)_dS2_(\d+)_dF_(\d+)_alpha_([mp])(\d+)$')


def model_name(key: str) -> str:
    match = RE_MODEL.fullmatch(key)
    if not match:
        return key
    a, b, c, sign, magnitude = match.groups()
    alpha = ('-' if sign == 'm' else '+') + magnitude
    return f'({a},{b},{c}), alpha={alpha}'


def valid_prediction(p: dict, observable: str) -> float:
    if p.get('status', 'Success') != 'Success':
        return math.nan
    try:
        value = float(p['prediction'][observable])
        return value if math.isfinite(value) else math.nan
    except (TypeError, ValueError, KeyError):
        return math.nan


def baseline(points: list[dict], observable: str) -> float:
    matches = [p for p in points if math.isclose(float(p['gamma']), 1.0, rel_tol=0, abs_tol=1e-12)]
    if len(matches) != 1:
        raise ValueError(f'Expected exactly one gamma=1 point, found {len(matches)}')
    val = valid_prediction(matches[0], observable)
    if not math.isfinite(val):
        raise ValueError(f'Invalid gamma=1 value for {observable}')
    return val


def series(points: list[dict], observable: str) -> tuple[np.ndarray, np.ndarray, int]:
    xs = np.asarray([float(p['gamma']) for p in points], dtype=float)
    if np.any(~np.isfinite(xs)) or np.any(np.diff(xs) <= 0):
        raise ValueError('Gamma samples must be finite and strictly increasing')
    ys = np.asarray([valid_prediction(p, observable) for p in points])
    ref = baseline(points, observable)
    if abs(ref) < 1e-40:
        raise ValueError(f'Cannot calculate a percentage change relative to zero {observable}')
    deviations = 100.0 * (ys / ref - 1.0)
    failures = int(np.count_nonzero(~np.isfinite(ys)))
    return xs, deviations, failures


def load_models(folder: Path) -> tuple[list[dict], list[str]]:
    items, warnings = [], []
    for path in sorted(folder.glob('*/sensitivity_diagnostics.json')):
        try:
            obj = json.loads(path.read_text(encoding='utf-8-sig'))
            for direction in DIRECTIONS:
                if direction not in obj.get('directions', {}):
                    raise ValueError(f'Missing {direction}')
                points = obj['directions'][direction]['points']
                for obs in OBS:
                    baseline(points, obs)
            items.append({'model': path.parent.name, 'json': path, 'data': obj})
        except (OSError, KeyError, ValueError, TypeError) as exc:
            warnings.append(f'{path}: {type(exc).__name__}: {exc}')
    return items, warnings


def style_axes(ax, all_x: list[np.ndarray]) -> None:
    ax.axhline(0, color='0.5', alpha=.7, lw=.8)
    ax.axvline(1, color='0.6', alpha=.5, lw=.8)
    ax.grid(alpha=.17)
    if all_x:
        x = np.concatenate(all_x)
        if x.min() > 0 and x.max() / x.min() > 10:
            ax.set_xscale('log')
    ax.set_xlabel(r'UV deformation $\gamma$')


def draw_one(ax, data: dict, observable: str, show_legend: bool = True) -> None:
    xx = []
    for direction in DIRECTIONS:
        points = data['directions'][direction]['points']
        xs, ys, _ = series(points, observable)
        color, line = STYLES[direction]
        # NaN gaps ensure failed points are not joined across.
        ax.plot(xs, ys, linestyle=line, marker='o', markersize=2.8,
                color=color, lw=1.4, label=direction)
        xx.append(xs)
    style_axes(ax, xx)
    ax.set_ylabel(r'$100\,[O(\gamma)/O(1)-1]$ (%)')
    if show_legend:
        ax.legend(fontsize=8, loc='best')


def per_model_plots(item: dict, target: Path) -> None:
    key = item['model']
    data = item['data']
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    for ax, observable in zip(axes, OBS[:2]):
        draw_one(ax, data, observable)
        ax.set_title(LABELS[observable])
    fig.suptitle(model_name(key) + ' — fitted UV sensitivity (not an RG-only correction)')
    fig.savefig(target / f'{key}_mass_sensitivity.png', dpi=170)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.1), constrained_layout=True)
    for ax, observable in zip(axes, OBS[2:]):
        draw_one(ax, data, observable)
        ax.set_title(LABELS[observable])
    fig.suptitle(model_name(key) + ' — fitted UV sensitivity (not an RG-only correction)')
    fig.savefig(target / f'{key}_angle_sensitivity.png', dpi=170)
    plt.close(fig)


def overview(items: list[dict], target: Path, observable: str) -> None:
    cols = 4
    rows = math.ceil(len(items) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(17, rows * 3.05), constrained_layout=True,
                             squeeze=False)
    for ax, item in zip(axes.flat, items):
        draw_one(ax, item['data'], observable, show_legend=False)
        ax.set_title(model_name(item['model']), fontsize=9)
        ax.set_ylabel('Change (%)')
    for ax in list(axes.flat)[len(items):]:
        ax.axis('off')
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside upper center', ncol=3)
    fig.suptitle(f'{LABELS[observable]}: relative variation from fitted $\\gamma=1$ baseline', y=1.02)
    fig.savefig(target / f'overview_{observable}.png', dpi=180, bbox_inches='tight')
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=DEFAULT_INPUT)
    parser.add_argument('--output', type=Path, default=None)
    args = parser.parse_args()
    source = args.input.resolve()
    dest = (args.output or source / 'thesis_figures').resolve()
    items, warnings = load_models(source)
    if not items:
        raise SystemExit(f'No usable sensitivity_diagnostics.json files found under {source}')
    dest.mkdir(parents=True, exist_ok=True)
    summary = []
    for item in items:
        per_model_plots(item, dest)
        for direction in DIRECTIONS:
            points = item['data']['directions'][direction]['points']
            n_failed = sum(p.get('status', 'Success') != 'Success' for p in points)
            for observable in OBS:
                xs, ys, _ = series(points, observable)
                finite = np.isfinite(ys)
                summary.append({'model': item['model'], 'representation': model_name(item['model']),
                    'direction': direction, 'observable': observable, 'n_points': len(points),
                    'n_failed_points': n_failed, 'n_valid_points': int(finite.sum()),
                    'gamma_min_valid': float(xs[finite].min()),
                    'gamma_max_valid': float(xs[finite].max()),
                    'max_abs_pct_change_from_gamma1': float(np.abs(ys[finite]).max()),
                    'interpretation': 'UV gamma sensitivity; not isolated RG correction'})
    for observable in OBS:
        overview(items, dest, observable)
    with (dest / 'sensitivity_plot_summary.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    print(f'Plotted {len(items)} models into: {dest}')
    print(f'Created {2*len(items)+len(OBS)} PNG figures plus sensitivity_plot_summary.csv')
    for warning in warnings:
        print('SKIPPED:', warning)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
