"""Presentation-only plots for saved fitted heavy-fermion threshold scans.

Project path: Numerical/plotting/FermionThresholdVariationFigures.py

No RGEs, fitting, matching, or source JSON/CSV are recalculated here.
Saved raw sources (historical folder IDs retained for compatibility):
  Reports/output/full/staged_rg_comparison/<scenario>/scan.json
  Reports/output/full/weinberg_running_comparison/<scenario>/scan.json
  Reports/output/full/fermion_scale/<scenario>/scan.json

The staged percentages are endpoint differences between stage-on and
corresponding segment-frozen calculations, NOT O(mu) trajectories.
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
INTERVALS = (
    ('UV', '000', '100', 'UV segment: 000 → 100'),
    ('EFT1', '100', '110', 'EFT1 segment: 100 → 110'),
    ('Weinberg', '110', '111', 'Weinberg segment: 110 → 111'),
)


def model_legend(raw_model: str) -> str:
    """Use requested T3A-m2 visual convention; preserve physical model IDs."""
    canonical = shorten(raw_model)  # verified project mapping in ModelNames.py
    return canonical.replace('T3-', 'T3', 1)


def load_rows(source: Path) -> list[dict]:
    source = Path(source)
    data = json.loads(source.read_text(encoding='utf-8-sig'))
    if not isinstance(data, list):
        raise ValueError(f'Expected saved list of scan rows: {source}')
    return data


def _plt():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    return plt


def _segments(rows: list[dict], getter):
    """Do not connect curves across missing / failed saved scan points."""
    sections, section = [], []
    for row in sorted(rows, key=lambda r: float(r['mass_ratio'])):
        value = getter(row)
        try:
            ok = row.get('status') == 'Success' and value is not None and math.isfinite(float(value))
        except (ValueError, TypeError):
            ok = False
        if ok:
            section.append((float(row['mass_ratio']), float(value)))
        elif section:
            sections.append(section)
            section = []
    if section:
        sections.append(section)
    return sections


def _model_groups(rows: list[dict]):
    return {
        model: [row for row in rows if row.get('model') == model]
        for model in sorted({row['model'] for row in rows if row.get('model')})
    }


def _ratio_axis(ax) -> None:
    ax.set_xscale('log')
    ax.set_xlabel(r'$\mu_F/\mu_S$  (vary $M_F$ and $\mu_F$ together; $\mu_S$ fixed)')
    ax.grid(alpha=0.22)


def stage_switch_off_plots(rows: list[dict], output: Path) -> int:
    """Render saved 000→100→110→111 running vs frozen comparisons."""
    plt = _plt()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    count = 0
    by_model = _model_groups(rows)
    for obs in OBSERVABLES:
        fig, axes = plt.subplots(3, 1, figsize=(12.2, 10.2), sharex=True, constrained_layout=True)
        for ax, (key, before, after, label) in zip(axes, INTERVALS):
            for model, model_rows in by_model.items():
                def getter(r):
                    return (r.get('increment_percent') or {}).get(key, {}).get(obs)
                sections = _segments(model_rows, getter)
                for i, section in enumerate(sections):
                    ax.plot([x for x, _ in section], [y for _, y in section],
                            '-o', ms=3.0, lw=1.35,
                            label=model_legend(model) if i == 0 else None)
            ax.axhline(0, lw=0.7, color='black')
            ax.set_title(label, loc='left', fontsize=10)
            ax.set_ylabel('Running on vs\nsegment frozen (%)')
            ax.grid(alpha=0.22)
        _ratio_axis(axes[-1])
        # Legend outside axes: legible even with all 16 fitted representations.
        handles, labels = axes[0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc='outside right center', ncol=1,
                       fontsize=8, title='T3 models')
        fig.suptitle(f'Fermion-threshold variation: {OBS_LABELS[obs]}\n'
                     'Final low-energy observable; RGE segment enabled vs same segment frozen',
                     fontsize=13)
        fig.savefig(output / f'fermion_threshold_running_effect_{obs}.png', dpi=180)
        plt.close(fig)
        count += 1
    print(f'Wrote {count} fermion-threshold RGE switch-off figures to {output}')
    return count


def weinberg_running_plots(rows: list[dict], output: Path) -> int:
    """Isolate the final C5 running for scans which also vary mu_F and M_F."""
    plt = _plt()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    by_model = _model_groups(rows)
    count = 0
    for obs in OBSERVABLES:
        fig, ax = plt.subplots(figsize=(11.8, 6.0), constrained_layout=True)
        for model, model_rows in by_model.items():
            sections = _segments(model_rows, lambda r: (r.get('weinberg_running_percent') or {}).get(obs))
            for i, section in enumerate(sections):
                ax.plot([x for x, _ in section], [y for _, y in section],
                        '-o', markersize=3, lw=1.35,
                        label=model_legend(model) if i == 0 else None)
        ax.axhline(0, color='black', lw=0.7)
        _ratio_axis(ax)
        ax.set_ylabel('Evolved vs frozen $C_5$ (%)')
        ax.set_title('Fermion-threshold variation: final Weinberg-EFT running\n'
                     f'Effect on {OBS_LABELS[obs]} at the low-energy endpoint')
        if ax.get_legend_handles_labels()[0]:
            ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0.,
                      fontsize=8, title='T3 models')
        fig.savefig(output / f'weinberg_running_{obs}.png', dpi=180, bbox_inches='tight')
        plt.close(fig)
        count += 1
    print(f'Wrote {count} Weinberg frozen-C5 comparison figures to {output}')
    return count


def fitted_mass_plots(rows: list[dict], output: Path) -> int:
    """Re-label the separate fitted heavy-mass scan (baseline gamma=1)."""
    plt = _plt()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    by_model = _model_groups(rows)
    count = 0
    for obs in OBSERVABLES:
        fig, ax = plt.subplots(figsize=(11.8, 6.0), constrained_layout=True)
        for model, model_rows in by_model.items():
            sections = _segments(model_rows, lambda r: (r.get('relative_percent') or {}).get(obs))
            for i, section in enumerate(sections):
                ax.plot([x for x, _ in section], [y for _, y in section],
                        '-o', markersize=3, lw=1.35,
                        label=model_legend(model) if i == 0 else None)
        ax.axhline(0, color='black', lw=0.7)
        _ratio_axis(ax)
        ax.set_ylabel('Change vs same-model baseline (%)')
        ax.set_title(f'Fitted fermion-threshold variation: {OBS_LABELS[obs]}\n'
                     'Mass/threshold variation, NOT an RGE switch-off comparison')
        if ax.get_legend_handles_labels()[0]:
            ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0.,
                      fontsize=8, title='T3 models')
        fig.savefig(output / f'fermion_scale_{obs}.png', dpi=180, bbox_inches='tight')
        plt.close(fig)
        count += 1
    print(f'Wrote {count} fitted fermion-scale figures to {output}')
    return count


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
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--study', choices=('fermion-threshold', 'weinberg', 'fitted-mass'),
                    action='append', help='Default: all three saved-data plot groups')
    ap.add_argument('--scenario', default='smallY_smallL',
                    help='Historical saved-data folder ID, not the four-scenario study')
    args = ap.parse_args(argv)
    for study in dict.fromkeys(args.study or ('fermion-threshold', 'weinberg', 'fitted-mass')):
        rebuild(study, args.scenario)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
