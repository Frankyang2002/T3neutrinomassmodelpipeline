"""Publishable saved Zee sensitivity figures: 2x3 focus and 4x4 full sets.

Project path: scripts/PlotT3SplitZeeSensitivity.py

Default produces three figures from the existing general-T3 T1/T2/T3
sensitivity scan (no refitting, RGEs, matching, or new scan points):
  - 6 panels: solar top row, atmospheric bottom row for B_-1,C_-1,A_-4
  - 16 panels: solar plane for all models
  - 16 panels: atmospheric plane for all models

Both PNG and PDF are written by the repository's original _save_figure helper.
The saved source uses conditional one-scale C5 and is not the corrected
sequential LLSS EFT1-freeze calculation.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from scripts import PlotT3GeneralZeeOvals as original

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / 'Reports/output/thesis_general_5_3/zee_solar_atmospheric_split_v1'
FOCUS_DEFAULT = ('B_-1', 'C_-1', 'A_-4')
PLANES = {
    'solar': ('sin2_theta12', 'delta_m21_sq_ev2'),
    'atmospheric': ('sin2_theta23', 'delta_m3l_sq_ev2'),
}


def parse_models(value: str) -> tuple[str, ...]:
    group = tuple(x.strip() for x in value.split(',') if x.strip())
    if not group or len(set(group)) != len(group):
        raise ValueError('Model selection empty or has duplicates')
    unknown = sorted(set(group) - set(original.MODEL_ORDER))
    if unknown:
        raise ValueError(f'Unknown models: {unknown}')
    return group


def compact_legend(handles):
    """Keep actual handles from the original plotting code, shorten only text."""
    replacements = {
        'NuFIT 6.1 ~3σ (joint Gaussian)': 'NuFIT 6.1 ~3σ',
        'NuFIT 6.1 ~1σ (joint Gaussian)': 'NuFIT 6.1 ~1σ',
        'JUNO 2026 1σ (measured, joint Gaussian)': 'JUNO 2026 1σ',
        'JUNO 2026 1σ (measured)': 'JUNO 2026 1σ',
        'JUNO 6-y 1σ (projection)': 'JUNO 6-y 1σ',
        'Model response oval (descriptive)': 'Model response oval',
    }
    for handle in handles or []:
        if hasattr(handle, 'get_label') and hasattr(handle, 'set_label'):
            label = handle.get_label()
            if label in replacements:
                handle.set_label(replacements[label])
    return handles


def add_legend(fig, handles, ncols: int) -> None:
    if not handles:
        return
    fig.legend(
        handles=compact_legend(handles),
        loc='upper center', bbox_to_anchor=(0.5, 0.98),
        ncol=min(len(handles), 7),
        fontsize=8.7 if ncols == 3 else 8.2,
        frameon=False, columnspacing=1.2, handlelength=1.7,
    )


def draw_panel(ax, model, data, nufit, *, plane, view, juno_mode):
    x, y = PLANES[plane]
    return original._plot_panel(
        ax, model, data, nufit, x, y,
        show_model_oval=True,
        juno_mode=juno_mode,
        zoom=view,
    )


def save_figure(fig, destination: Path) -> None:
    original._save_figure(fig, destination, dpi=200)
    import matplotlib.pyplot as plt
    plt.close(fig)


def draw_focus(plots, nufit, models, out, *, view, juno_mode):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    n = len(models)
    fig, axes = plt.subplots(2, n, figsize=(5.6*n + 0.5, 10.0), squeeze=False)
    fig.set_facecolor('white')
    handles = None
    for row, plane in enumerate(('solar', 'atmospheric')):
        for col, model in enumerate(models):
            handles = draw_panel(axes[row, col], model, plots[model], nufit,
                                 plane=plane, view=view, juno_mode=juno_mode)
    # Row identifiers are kept outside the panel area, not in their titles.
    fig.text(0.017, 0.684, 'Solar', fontsize=12, weight='semibold',
             ha='center', va='center', rotation=90)
    fig.text(0.017, 0.300, 'Atmospheric', fontsize=12, weight='semibold',
             ha='center', va='center', rotation=90)
    fig.subplots_adjust(left=0.080, right=0.985, bottom=0.075,
                        top=0.885, wspace=0.34, hspace=0.36)
    add_legend(fig, handles, n)
    save_figure(fig, out)


def draw_all16(plots, nufit, out, *, plane, view, juno_mode):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(4, 4, figsize=(18.4, 16.5), squeeze=False)
    fig.set_facecolor('white')
    handles = None
    for ax, model in zip(axes.flat, original.MODEL_ORDER):
        handles = draw_panel(ax, model, plots[model], nufit,
                             plane=plane, view=view, juno_mode=juno_mode)
    fig.subplots_adjust(left=0.068, right=0.982, bottom=0.061,
                        top=0.916, wspace=0.45, hspace=0.51)
    add_legend(fig, handles, 4)
    save_figure(fig, out)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=original.DEFAULT_SOURCE,
                        help='Directory with saved points.csv and summary.json')
    parser.add_argument('--nufit', type=Path, default=original.DEFAULT_TARGET)
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--focus-models', default=','.join(FOCUS_DEFAULT),
                        help='Models shown in the 2-row joint figure (default: B_-1,C_-1,A_-4)')
    parser.add_argument('--view', choices=('experiment', 'trajectory'), default='experiment')
    parser.add_argument('--juno-mode', choices=('measured', 'future', 'both'), default='measured')
    parser.add_argument('--figures', choices=('all', 'focus', 'all16'), default='all',
                        help='Which figures to create')
    args = parser.parse_args(argv)

    try:
        focus_models = parse_models(args.focus_models)
    except ValueError as exc:
        parser.error(str(exc))
    for path in (args.source/'points.csv', args.source/'summary.json', args.nufit):
        if not path.is_file():
            parser.error(f'Missing input: {path}')

    # Validate all 16 only if needed; default produces all three figures.
    needed = (list(original.MODEL_ORDER) if args.figures != 'focus'
              else list(focus_models))
    plots, _metadata = original.load_points(
        args.source/'points.csv', args.source/'summary.json', models=needed)
    nufit = original.load_nufit(args.nufit)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {}
    suffix = args.view

    if args.figures in ('all', 'focus'):
        tag = '_'.join(x.replace('-', 'm').replace('+', 'p') for x in focus_models)
        out = args.output_dir/f'zee_focus_solar_atmospheric_{tag}_{suffix}.png'
        draw_focus(plots, nufit, focus_models, out,
                   view=args.view, juno_mode=args.juno_mode)
        outputs['focus_2row'] = {'models': list(focus_models),
                                'png': str(out), 'pdf': str(out.with_suffix('.pdf'))}
        print('Saved combined solar + atmospheric:', out, flush=True)

    if args.figures in ('all', 'all16'):
        for plane in ('solar', 'atmospheric'):
            out = args.output_dir/f'zee_all16_{plane}_{suffix}.png'
            draw_all16(plots, nufit, out,
                       plane=plane, view=args.view, juno_mode=args.juno_mode)
            outputs[f'all16_{plane}'] = {'models': list(original.MODEL_ORDER),
                                        'png': str(out), 'pdf': str(out.with_suffix('.pdf'))}
            print(f'Saved all-16 {plane}:', out, flush=True)

    manifest = {
        'format': 't3-saved-zee-sensitivity-solar-atmospheric-figures-v1',
        'source_campaign': str(args.source),
        'view': args.view,
        'juno_mode': args.juno_mode,
        'focus_models': list(focus_models),
        'planes': {name: {'x': xy[0], 'y': xy[1]} for name, xy in PLANES.items()},
        'physics_scope': 'Saved conditional one-scale general T3 Zee T1/T2/T3 sensitivity; not corrected sequential LLSS EFT1 freeze',
        'physics_rerun': False,
        'outputs': outputs,
    }
    (args.output_dir/'provenance.json').write_text(
        json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
