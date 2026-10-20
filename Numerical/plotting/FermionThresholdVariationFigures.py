"""Presentation-only figures regenerated from saved heavy-fermion scan JSON.

Project path: Numerical/plotting/FermionThresholdVariationFigures.py

The horizontal coordinate is the fermion/scalar threshold ratio mu_F/mu_S:
M_F and mu_F are varied together, while mu_S remains fixed. It is not
an RGE evolution coordinate. No physics calculation or saved data is changed.

Historical (unchanged) source directories:
  Reports/output/full/staged_rg_comparison/<scenario>/scan.json
  Reports/output/full/weinberg_running_comparison/<scenario>/scan.json
  Reports/output/full/fermion_scale/<scenario>/scan.json
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from Numerical.orchestration.ModelNames import shorten

ROOT = Path(__file__).resolve().parents[2]
OBSERVABLES = (
    'delta_m21_sq_ev2', 'delta_m3l_sq_ev2',
    'sin2_theta12', 'sin2_theta13', 'sin2_theta23',
)
OBS_LABELS = {
    'delta_m21_sq_ev2': r'$\Delta m_{21}^2$',
    'delta_m3l_sq_ev2': r'$\Delta m_{3\ell}^2$',
    'sin2_theta12': r'$\sin^2\theta_{12}$',
    'sin2_theta13': r'$\sin^2\theta_{13}$',
    'sin2_theta23': r'$\sin^2\theta_{23}$',
}
# Changing the ordering changes the attribution. Preserve the saved definitions.
INTERVALS = (
    ('UV', '000', '100'),
    ('EFT1', '100', '110'),
    ('Weinberg', '110', '111'),
)


def model_legend(raw_model: str) -> str:
    """Canonical ModelNames mapping, with the T3A-m2 presentation spelling."""
    return shorten(raw_model).replace('T3-', 'T3', 1)


def load_rows(source: Path) -> list[dict]:
    source = Path(source)
    data = json.loads(source.read_text(encoding='utf-8-sig'))
    if not isinstance(data, list) or not all(isinstance(row, dict) for row in data):
        raise ValueError(f'Expected saved list of scan records: {source}')
    return data


def _plt():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def _ratio(row: dict) -> float | None:
    """Model-level MissingArtifacts/MissingFit rows have no scan coordinate.

    Those records must not enter the sort. By contrast, unsuccessful numeric
    points *with* a coordinate remain present and split plotted line segments.
    """
    try:
        ratio = float(row['mass_ratio'])
    except (KeyError, ValueError, TypeError, OverflowError):
        return None
    return ratio if math.isfinite(ratio) and ratio > 0 else None


def _segments(rows: list[dict], getter):
    """Finite successful samples only; never connect across failed scan points."""
    ordered = sorted(((ratio, row) for row in rows
                      if (ratio := _ratio(row)) is not None), key=lambda p: p[0])
    sections, current = [], []
    for x, row in ordered:
        value = getter(row) if row.get('status') == 'Success' else None
        try:
            y = float(value) if value is not None else math.nan
        except (TypeError, ValueError, OverflowError):
            y = math.nan
        if math.isfinite(y):
            current.append((x, y))
        elif current:
            sections.append(current)
            current = []
    if current:
        sections.append(current)
    return sections


def _model_groups(rows: list[dict]):
    names = sorted({row['model'] for row in rows if isinstance(row.get('model'), str)},
                   key=model_legend)
    return {model: [row for row in rows if row.get('model') == model]
            for model in names}


def _require_plottable(rows: list[dict]) -> None:
    if not any(r.get('status') == 'Success' and _ratio(r) is not None for r in rows):
        raise ValueError('No successful saved numerical points with a valid mass_ratio')


def _ratio_axis(ax) -> None:
    ax.set_xscale('log')
    ax.set_xlabel(r'$\mu_F/\mu_S$  ($M_F$ and $\mu_F$ varied together; $\mu_S$ fixed)')
    ax.grid(alpha=0.22)


def _legend(fig, ax, *, external: bool = True) -> None:
    handles, labels = ax.get_legend_handles_labels()
    # Retain order but never print the same model twice because of scan gaps.
    unique = dict(zip(labels, handles))
    if not unique:
        return
    if external:
        fig.legend(unique.values(), unique.keys(), loc='outside right center',
                   fontsize=8, ncol=1, title='T3 models')
    else:
        ax.legend(unique.values(), unique.keys(), fontsize=8, ncol=2,
                  loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0,
                  title='T3 models')


def stage_switch_off_plots(rows: list[dict], output: Path) -> int:
    """Use saved incremental percentages: 000→100, 100→110, 110→111.

    Each displayed value compares a final low-energy observable after enabling
    one additional RGE segment with the corresponding frozen-segment result.
    Matching is part of both calculations. It is not an energy derivative.
    """
    _require_plottable(rows)
    plt = _plt()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    by_model = _model_groups(rows)
    count = 0
    for obs in OBSERVABLES:
        fig, axes = plt.subplots(3, 1, figsize=(12.8, 10.9), sharex=True,
                                 constrained_layout=True)
        for ax, (stage, before, after) in zip(axes, INTERVALS):
            for model, points in by_model.items():
                sections = _segments(points, lambda row: (
                    (row.get('increment_percent') or {}).get(stage) or {}).get(obs))
                for i, section in enumerate(sections):
                    ax.plot([x for x, _ in section], [y for _, y in section],
                            '-o', markersize=3, lw=1.35,
                            label=model_legend(model) if i == 0 else None)
            ax.axhline(0, lw=0.7, color='black')
            ax.set_title(f'{stage}: {before} → {after}  (matching retained)',
                         fontsize=10, loc='left')
            ax.set_ylabel(f'{stage} running enabled\nvs frozen (%)')
            ax.grid(alpha=0.22)
        _ratio_axis(axes[-1])
        _legend(fig, axes[0])
        fig.suptitle('Fermion-threshold variation — RGE running versus frozen segments\n'
                     + f'{OBS_LABELS[obs]} at the final low-energy endpoint', fontsize=12.4)
        fig.savefig(output / f'fermion_threshold_running_effect_{obs}.png', dpi=180)
        plt.close(fig)
        count += 1
    print(f'Wrote {count} fermion-threshold running/frozen figures to {output}')
    return count


def weinberg_running_plots(rows: list[dict], output: Path) -> int:
    """Compare saved final-EFT C5 evolved/frozen low-energy observables."""
    _require_plottable(rows)
    plt = _plt()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    by_model = _model_groups(rows)
    for obs in OBSERVABLES:
        fig, ax = plt.subplots(figsize=(11.8, 6.0), constrained_layout=True)
        for model, points in by_model.items():
            sections = _segments(points, lambda row: (
                row.get('weinberg_running_percent') or {}).get(obs))
            for i, section in enumerate(sections):
                ax.plot([x for x, _ in section], [y for _, y in section],
                        '-o', markersize=3, lw=1.35,
                        label=model_legend(model) if i == 0 else None)
        ax.axhline(0, color='black', lw=0.7)
        _ratio_axis(ax)
        ax.set_ylabel('Weinberg running enabled vs frozen (%)')
        ax.set_title('Fermion-threshold variation — Weinberg running versus frozen $C_5$\n'
                     f'Final low-energy {OBS_LABELS[obs]} (UV/EFT1 running and matching included)')
        _legend(fig, ax, external=False)
        fig.savefig(output / f'weinberg_running_{obs}.png', dpi=180,
                    bbox_inches='tight')
        plt.close(fig)
    print(f'Wrote {len(OBSERVABLES)} Weinberg figures to {output}')
    return len(OBSERVABLES)


def fitted_mass_plots(rows: list[dict], output: Path) -> int:
    """Saved mass-scale effects relative to factor=1, NOT an RGE switch-off."""
    _require_plottable(rows)
    plt = _plt()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    by_model = _model_groups(rows)
    for obs in OBSERVABLES:
        fig, ax = plt.subplots(figsize=(11.8, 6.0), constrained_layout=True)
        for model, points in by_model.items():
            sections = _segments(points, lambda row: (
                row.get('relative_percent') or {}).get(obs))
            for i, section in enumerate(sections):
                ax.plot([x for x, _ in section], [y for _, y in section],
                        '-o', markersize=3, lw=1.35,
                        label=model_legend(model) if i == 0 else None)
        ax.axhline(0, color='black', lw=0.7)
        _ratio_axis(ax)
        ax.set_ylabel('Change versus same-model factor=1 baseline (%)')
        ax.set_title(f'Fitted fermion-threshold variation: final {OBS_LABELS[obs]}\n'
                     'Mass/threshold variation (not an RGE switch-off comparison)')
        _legend(fig, ax, external=False)
        fig.savefig(output / f'fermion_scale_{obs}.png', dpi=180,
                    bbox_inches='tight')
        plt.close(fig)
    print(f'Wrote {len(OBSERVABLES)} fitted-mass figures to {output}')
    return len(OBSERVABLES)


def rebuild(kind: str, scenario: str = 'smallY_smallL') -> None:
    full = ROOT / 'Reports' / 'output' / 'full'
    display = ROOT / 'Reports' / 'output' / 'display'
    if kind == 'fermion-threshold':
        stage_switch_off_plots(load_rows(full / 'staged_rg_comparison' / scenario / 'scan.json'),
                               display / 'fermion_threshold_variation' / scenario)
    elif kind == 'weinberg':
        weinberg_running_plots(load_rows(full / 'weinberg_running_comparison' / scenario / 'scan.json'),
                               display / 'weinberg_running' / scenario)
    elif kind == 'fitted-mass':
        folder = full / 'fermion_scale' / scenario
        fitted_mass_plots(load_rows(folder / 'scan.json'), folder)
    else:
        raise ValueError(f'Unknown presentation study: {kind}')


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', choices=('fermion-threshold', 'weinberg', 'fitted-mass'),
                        action='append', help='Default: all three saved-data plot groups')
    parser.add_argument('--scenario', default='smallY_smallL')
    args = parser.parse_args(argv)
    for study in dict.fromkeys(args.study or ('fermion-threshold', 'weinberg', 'fitted-mass')):
        rebuild(study, args.scenario)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
