"""Plot saved Zee sensitivity for three focus models and all 16 T3 models.

Thesis typography refresh with nonoverlapping focus-panel y-axis labels.
Reads EXACTLY the saved source points; does not refit, run RGEs, or match.

Outputs: original 2x3 focus and both 4x4 overviews, plus optional readable
2x2 pages (four per plane) for inclusion at thesis text width.

The source is conditional one-scale C5, NOT corrected sequential LLSS.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts import PlotT3GeneralZeeOvals as original

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / 'Reports/output/thesis_general_5_3/zee_solar_atmospheric_split_short_legends_both_v4'
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


def compact_label(raw: str) -> str:
    """Map source plotting labels to short display names, explicitly."""
    label = str(raw).strip()
    if label.startswith('NuFIT'):
        if '3σ' in label or '3\\sigma' in label:
            return 'NuFIT 3σ'
        if '1σ' in label or '1\\sigma' in label:
            return 'NuFIT 1σ'
    if label.startswith('JUNO'):
        return 'JUNO 6-y 1σ' if ('6-y' in label or 'projection' in label) else 'JUNO 1σ'
    if label.startswith('Model response'):
        return 'Model response'
    return label


def unique_legend_handles(groups):
    """Collect legend artists from all panels, notably solar-only JUNO."""
    by_label = {}
    for handles in groups:
        for handle in handles or []:
            label = compact_label(handle.get_label())
            if not label or label.startswith('_'):
                continue
            by_label.setdefault(label, handle)
    # Stable, compact ordering.
    priority = ('NuFIT 1σ', 'NuFIT 3σ', 'JUNO 1σ', 'JUNO 6-y 1σ',
                'Model response', 'T1', 'T2', 'T3')
    ordered = [name for name in priority if name in by_label]
    ordered.extend(name for name in by_label if name not in ordered)
    return [(by_label[name], name) for name in ordered]


def font_scale(figure_width: float, intended_width: float, body_pt: float) -> float:
    """Source point size that prints at body_pt after figure downscaling.

    Use for axis *text*, not for physics coordinates, limits, or markers.
    """
    if intended_width <= 0 or body_pt <= 0:
        raise ValueError('Thesis width/body point size must be positive')
    return figure_width / intended_width * body_pt


def improve_panel(ax, *, source_width: float, thesis_width: float,
                  body_pt: float, grid: str) -> None:
    """Restyle existing axes without changing any plotted point or region."""
    target = font_scale(source_width, thesis_width, body_pt)
    # The 4x4 overview cannot hold full body-sized titles/ticks at text width:
    # retain a compact overview, while 2x2 sheets provide body-sized labels.
    factor = 0.68 if grid == 'overview' else 0.93
    labelsize = target * factor
    ticksize = target * (0.62 if grid == 'overview' else 0.82)
    titlesize = target * (0.83 if grid == 'overview' else 1.02)
    ax.xaxis.label.set_fontsize(labelsize)
    ax.yaxis.label.set_fontsize(labelsize)
    ax.title.set_fontsize(titlesize)
    ax.tick_params(axis='both', which='major', labelsize=ticksize, pad=3)
    ax.tick_params(axis='both', which='minor', labelsize=ticksize)
    # Prevent scientific offset text from staying at old 8-point size.
    ax.xaxis.get_offset_text().set_fontsize(ticksize)
    ax.yaxis.get_offset_text().set_fontsize(ticksize)


def add_legend(fig, handles, *, source_width: float, thesis_width: float,
               body_pt: float, columns: int, overview: bool = False) -> None:
    entries = unique_legend_handles([handles])
    if not entries:
        return
    # CRITICAL: pass labels explicitly. Changing handle.get_label() does not
    # guarantee that Matplotlib's composed legend text has been replaced.
    artists, labels = zip(*entries)
    point = font_scale(source_width, thesis_width, body_pt)
    legend = fig.legend(
        handles=list(artists), labels=list(labels),
        loc='upper center', bbox_to_anchor=(0.5, 0.992),
        ncol=min(len(labels), 4 if columns <= 2 else 7),
        fontsize=point * (0.66 if overview else 0.77),
        frameon=False, columnspacing=1.1, handlelength=1.5,
    )
    # Fail rather than silently generate old long labels again.
    actual = [entry.get_text() for entry in legend.get_texts()]
    if actual != list(labels) or any('measured' in x.lower() or
                                     'joint gaussian' in x.lower() for x in actual):
        raise RuntimeError(f'Legend labels were not shortened: {actual}')
    print('Legend text:', ' | '.join(actual), flush=True)


def draw_panel(ax, model, data, nufit, *, plane, view, juno_mode):
    x, y = PLANES[plane]
    return original._plot_panel(
        ax, model, data, nufit, x, y,
        show_model_oval=True, juno_mode=juno_mode, zoom=view,
    )


def save_figure(fig, destination: Path) -> None:
    original._save_figure(fig, destination, dpi=240)
    import matplotlib.pyplot as plt
    plt.close(fig)


def draw_focus(plots, nufit, models, out, *, view, juno_mode,
               thesis_width, body_pt):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    n = len(models)
    width = 4.15 * n + 0.9
    fig, axes = plt.subplots(2, n, figsize=(width, 9.8), squeeze=False)
    fig.set_facecolor('white')
    legend_groups = []
    for row, plane in enumerate(('solar', 'atmospheric')):
        for col, model in enumerate(models):
            ax = axes[row, col]
            handles = draw_panel(ax, model, plots[model], nufit,
                                 plane=plane, view=view, juno_mode=juno_mode)
            legend_groups.append(handles)
            improve_panel(ax, source_width=width, thesis_width=thesis_width,
                          body_pt=body_pt, grid='focus')
    # The old vertical row headings ("Solar" / "Atmospheric") overlap
    # the leftmost axis label at thesis width. The axis observables already
    # identify each row, so remove the redundant headings entirely.
    # Reserve a wider left margin for the atmospheric mass-squared label.
    fig.subplots_adjust(left=0.135, right=0.985, bottom=0.085,
                        top=0.858, wspace=0.46, hspace=0.47)
    add_legend(fig, [h for group in legend_groups for h in (group or [])],
               source_width=width, thesis_width=thesis_width,
               body_pt=body_pt, columns=n)
    save_figure(fig, out)


def draw_all16(plots, nufit, out, *, plane, view, juno_mode,
               thesis_width, body_pt):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    width = 15.8
    fig, axes = plt.subplots(4, 4, figsize=(width, 14.9), squeeze=False)
    fig.set_facecolor('white')
    handles = None
    for ax, model in zip(axes.flat, original.MODEL_ORDER):
        handles = draw_panel(ax, model, plots[model], nufit,
                             plane=plane, view=view, juno_mode=juno_mode)
        improve_panel(ax, source_width=width, thesis_width=thesis_width,
                      body_pt=body_pt, grid='overview')
    fig.subplots_adjust(left=0.090, right=0.978, bottom=0.080,
                        top=0.875, wspace=0.58, hspace=0.65)
    add_legend(fig, handles, source_width=width, thesis_width=thesis_width,
               body_pt=body_pt, columns=4, overview=True)
    save_figure(fig, out)


def draw_readable_pages(plots, nufit, output_dir, *, plane, view, juno_mode,
                        thesis_width, body_pt) -> list[Path]:
    """Four 2x2 pages: the 16-model data, legible at the thesis text width."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    files = []
    names = original.MODEL_ORDER
    for page, offset in enumerate(range(0, 16, 4), start=1):
        width = 7.2
        fig, axes = plt.subplots(2, 2, figsize=(width, 7.6), squeeze=False)
        fig.set_facecolor('white')
        handles = None
        for ax, model in zip(axes.flat, names[offset:offset + 4]):
            handles = draw_panel(ax, model, plots[model], nufit,
                                 plane=plane, view=view, juno_mode=juno_mode)
            improve_panel(ax, source_width=width, thesis_width=thesis_width,
                          body_pt=body_pt, grid='page')
        fig.subplots_adjust(left=.15, right=.98, bottom=.095, top=.845,
                            wspace=.57, hspace=.64)
        add_legend(fig, handles, source_width=width, thesis_width=thesis_width,
                   body_pt=body_pt, columns=2)
        file = output_dir / f'zee_all16_{plane}_{view}_page{page}.png'
        save_figure(fig, file)
        files.append(file)
    return files


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=original.DEFAULT_SOURCE,
                        help='Directory with saved points.csv and summary.json')
    parser.add_argument('--nufit', type=Path, default=original.DEFAULT_TARGET)
    parser.add_argument('--output-dir', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--focus-models', default=','.join(FOCUS_DEFAULT),
                        help='Default: B_-1,C_-1,A_-4')
    parser.add_argument('--view', choices=('experiment', 'trajectory', 'both'), default='both',
                        help='Plot zoomed-out experiment view, zoomed-in trajectory view, or both')
    parser.add_argument('--juno-mode', choices=('measured', 'future', 'both'), default='measured')
    parser.add_argument('--figures', choices=('all', 'focus', 'all16'), default='all')
    parser.add_argument('--thesis-width-in', type=float, default=6.3)
    parser.add_argument('--thesis-body-pt', type=float, default=11.0)
    parser.add_argument('--no-readable-pages', action='store_true',
                        help='Skip additional 2x2 readable all-16 pages')
    args = parser.parse_args(argv)

    try:
        focus_models = parse_models(args.focus_models)
        font_scale(6.3, args.thesis_width_in, args.thesis_body_pt)
    except ValueError as exc:
        parser.error(str(exc))
    for path in (args.source/'points.csv', args.source/'summary.json', args.nufit):
        if not path.is_file():
            parser.error(f'Missing input: {path}')
    needed = (list(original.MODEL_ORDER) if args.figures != 'focus'
              else list(focus_models))
    plots, _metadata = original.load_points(
        args.source/'points.csv', args.source/'summary.json', models=needed)
    nufit = original.load_nufit(args.nufit)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {}
    views = ('experiment', 'trajectory') if args.view == 'both' else (args.view,)
    for view in views:
        suffix = view
        # The legend text is composed explicitly by add_legend for both views.
        if args.figures in ('all', 'focus'):
            tag = '_'.join(x.replace('-', 'm').replace('+', 'p') for x in focus_models)
            out = args.output_dir/f'zee_focus_solar_atmospheric_{tag}_{suffix}.png'
            draw_focus(plots, nufit, focus_models, out,
                       view=view, juno_mode=args.juno_mode,
                       thesis_width=args.thesis_width_in, body_pt=args.thesis_body_pt)
            outputs[f'{view}_focus_2row'] = {'models': list(focus_models),
                'png': str(out), 'pdf': str(out.with_suffix('.pdf'))}
            print(f'Saved {view} focus figure:', out, flush=True)

        if args.figures in ('all', 'all16'):
            for plane in ('solar', 'atmospheric'):
                out = args.output_dir/f'zee_all16_{plane}_{suffix}.png'
                draw_all16(plots, nufit, out,
                           plane=plane, view=view, juno_mode=args.juno_mode,
                           thesis_width=args.thesis_width_in, body_pt=args.thesis_body_pt)
                outputs[f'{view}_all16_{plane}'] = {'models': list(original.MODEL_ORDER),
                    'png': str(out), 'pdf': str(out.with_suffix('.pdf'))}
                print(f'Saved {view} all-16 {plane} overview:', out, flush=True)
                if not args.no_readable_pages:
                    pages = draw_readable_pages(plots, nufit, args.output_dir,
                        plane=plane, view=view, juno_mode=args.juno_mode,
                        thesis_width=args.thesis_width_in, body_pt=args.thesis_body_pt)
                    outputs[f'{view}_all16_{plane}_readable_pages'] = [
                        {'png': str(q), 'pdf': str(q.with_suffix('.pdf'))} for q in pages]
                    print(f'Saved {len(pages)} {view} readable {plane} pages', flush=True)

    manifest = {
        'format': 't3-saved-zee-sensitivity-solar-atmospheric-figures-v2-thesis-typography',
        'source_campaign': str(args.source),
        'view': args.view, 'views_generated': list(views), 'juno_mode': args.juno_mode,
        'focus_models': list(focus_models),
        'planes': {name: {'x': xy[0], 'y': xy[1]} for name, xy in PLANES.items()},
        'thesis_width_in': args.thesis_width_in,
        'thesis_body_pt': args.thesis_body_pt,
        'physics_scope': 'Saved conditional one-scale general T3 Zee T1/T2/T3 sensitivity; not corrected sequential LLSS EFT1 freeze',
        'physics_rerun': False,
        'outputs': outputs,
    }
    (args.output_dir/'provenance.json').write_text(
        json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
