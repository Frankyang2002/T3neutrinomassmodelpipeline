"""Index selected T3 thesis figures without creating duplicate PNG files.

This is a presentation-only file locator, not a plot generator. The existing
representation-comparison plots must be generated separately from saved study
results. By default no figures are copied; use --copy for a thesis export.

From the repository root:
    python -m Numerical.orchestration.CollectThesisFigures
    python -m Numerical.orchestration.CollectThesisFigures --all
    python -m Numerical.orchestration.CollectThesisFigures --copy
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
COMPARISON = Path('Reports/output/full/representation_comparison/smallY_smallL')
DRAFT = Path('Reports/thesis_draft/figures')

# Filename, priority, description, explanation/limitation.
FIGURES = (
    ('baseline_solar_stage_contributions.png', 'primary',
     'Stage contributions at fitted benchmarks',
     'Individual-model UV, EFT1, Weinberg and total contributions to solar mass splitting at r=1.'),
    ('baseline_eft1_by_model.png', 'primary',
     'Intermediate EFT dependence by model',
     'Highlights the substantial spread of EFT1 effects; fitted couplings differ by model.'),
    ('eft1_by_representation_vs_mass_ratio.png', 'primary',
     'Intermediate EFT dependence on heavy-fermion threshold',
     'Shows variation with M_F and mu_F scaled together; only r=1 is a fitted point.'),
    ('uv_by_representation_vs_mass_ratio.png', 'supporting',
     'UV running versus fermion threshold',
     'Supporting comparison of the first running interval.'),
    ('weinberg_by_representation_vs_mass_ratio.png', 'supporting',
     'Weinberg running versus fermion threshold',
     'Supporting comparison of the final EFT running interval.'),
)


def build_index(*, include_all: bool = False, copy: bool = False,
                destination: Path = Path('Reports/figure_review'),
                root: Path = PROJECT_ROOT) -> int:
    """Index canonical figure paths; copy only with explicit permission.

    Never deletes source files or existing destination figures. An unavailable
    source is reported as missing rather than silently substituted.
    """
    root = Path(root).resolve()
    destination = Path(destination)
    if not destination.is_absolute():
        destination = root / destination
    destination = destination.resolve()
    destination.mkdir(parents=True, exist_ok=True)

    selected = [item for item in FIGURES if include_all or item[1] == 'primary']
    lines = [
        '# T3 numerical figure review', '',
        'This index references existing plots. It does not recalculate RGEs, fits, or plots.',
        'Use `--copy` only when a separate export of PNGs is required.', '',
        '| Priority | Figure | Purpose | Source |',
        '|---|---|---|---|',
    ]
    missing = []
    available = 0
    for filename, priority, title, note in selected:
        candidates = (root / COMPARISON / filename, root / DRAFT / filename)
        source = next((p for p in candidates if p.is_file()), None)
        if source is None:
            missing.append((filename, candidates))
            print(f'MISSING: {filename}')
            continue
        available += 1
        if copy:
            target = destination / filename
            if source.resolve() != target.resolve():
                shutil.copy2(source, target)
            print(f'COPIED: {source.relative_to(root)} -> {target}')
        else:
            print(f'FOUND: {source.relative_to(root)}')
        lines.append(f'| {priority} | `{filename}` | {title} | `{source.relative_to(root).as_posix()}` |')

    lines.extend(['', '## Interpretation notes', ''])
    for filename, _priority, _title, note in selected:
        lines.append(f'- **{filename}:** {note}')
    lines.extend([
        '',
        'All staged percentages are conditional, order-dependent comparisons, '
        'not scheme-independent diagrammatic contributions.',
        'Across independently fitted models, differences cannot be assigned '
        'to representation alone.', '',
    ])
    if missing:
        lines.extend(['## Missing requested figures', ''])
        for filename, candidates in missing:
            locations = ' and '.join(f'`{p.relative_to(root).as_posix()}`' for p in candidates)
            lines.append(f'- `{filename}`: searched {locations}')
        lines.append('')
    index = destination / 'FIGURE_INDEX.md'
    index.write_text('\n'.join(lines), encoding='utf-8')
    print(f'\nIndexed {available}/{len(selected)} figures: {index}')
    return 1 if missing else 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--all', action='store_true',
                        help='Include the two supporting RG figures')
    parser.add_argument('--copy', action='store_true',
                        help='Copy selected PNGs to the destination (off by default)')
    parser.add_argument('--destination', type=Path,
                        default=Path('Reports/figure_review'),
                        help='Index/export location relative to the project root')
    args = parser.parse_args(argv)
    return build_index(include_all=args.all, copy=args.copy, destination=args.destination)


if __name__ == '__main__':
    raise SystemExit(main())
