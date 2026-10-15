"""Copy selected existing T3 numerical figures into a review folder.

Run from anywhere with:
    python -m Numerical.orchestration.CollectThesisFigures
    python -m Numerical.orchestration.CollectThesisFigures --all

No plotting, fitting, RG running, or matching is performed.
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--all', action='store_true',
                        help='Also include the two supporting RG figures.')
    parser.add_argument('--destination', type=Path,
                        default=Path('Reports/figure_review'),
                        help='Destination (relative paths are resolved from the project root).')
    args = parser.parse_args()

    destination = args.destination
    if not destination.is_absolute():
        destination = PROJECT_ROOT / destination
    destination.mkdir(parents=True, exist_ok=True)

    selected = [item for item in FIGURES if args.all or item[1] == 'primary']
    copied = []
    missing = []
    for filename, priority, title, note in selected:
        candidates = [PROJECT_ROOT / COMPARISON / filename,
                      PROJECT_ROOT / DRAFT / filename]
        source = next((path for path in candidates if path.is_file()), None)
        if source is None:
            missing.append((filename, candidates))
            print(f'MISSING: {filename}')
            continue
        target = destination / filename
        # Avoid accidentally copying a file onto itself.
        if source.resolve() != target.resolve():
            shutil.copy2(source, target)
        copied.append((filename, priority, title, note, source.relative_to(PROJECT_ROOT)))
        print(f'COPIED: {source.relative_to(PROJECT_ROOT)} -> {target}')

    lines = [
        '# T3 numerical figure review', '',
        'These are copies of existing figures, not regenerated calculations.', '',
        '| Priority | Figure | Purpose | Source |',
        '|---|---|---|---|',
    ]
    for filename, priority, title, note, source in copied:
        lines.append(f'| {priority} | `{filename}` | {title} | `{source.as_posix()}` |')
    lines.extend(['', '## Interpretation notes', ''])
    for filename, _priority, _title, note, _source in copied:
        lines.append(f'- **{filename}:** {note}')
    lines.extend(['',
                  'All staged percentages are conditional, order-dependent comparisons, '
                  'not scheme-independent diagrammatic contributions.',
                  'Across independently fitted models, differences cannot be assigned '
                  'to representation alone.', ''])
    if missing:
        lines.extend(['## Missing requested figures', ''])
        for filename, candidates in missing:
            lines.append(f'- `{filename}`: searched ' + ' and '.join(
                f'`{p.relative_to(PROJECT_ROOT).as_posix()}`' for p in candidates))
        lines.append('')
    (destination / 'FIGURE_INDEX.md').write_text('\n'.join(lines), encoding='utf-8')
    print(f'\nCopied {len(copied)}/{len(selected)} figures to: {destination}')
    print(f'Index: {destination / "FIGURE_INDEX.md"}')
    return 1 if missing else 0


if __name__ == '__main__':
    raise SystemExit(main())
