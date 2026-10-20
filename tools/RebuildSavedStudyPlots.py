"""Rebuild retained fermion-threshold or Weinberg plots from saved scan JSON.

Project path: tools/RebuildSavedStudyPlots.py

PowerShell, repository root:
  python -m tools.RebuildSavedStudyPlots --study fermion-threshold
  python -m tools.RebuildSavedStudyPlots --study weinberg

This replaces the previous 'staged-rg' / 'rebuilt_studies' output convention.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from Numerical.plotting.FermionThresholdVariationFigures import (
    load_rows, stage_switch_off_plots, weinberg_running_plots,
)

ROOT = Path(__file__).resolve().parents[1]
CONFIG = {
    'fermion-threshold': ('staged_rg_comparison', 'fermion_threshold_variation', stage_switch_off_plots),
    'weinberg': ('weinberg_running_comparison', 'weinberg_running', weinberg_running_plots),
}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--study', required=True, choices=tuple(CONFIG))
    ap.add_argument('--scenario', default='smallY_smallL')
    ap.add_argument('--source', type=Path, help='Override saved scan.json (read only)')
    ap.add_argument('--output-dir', type=Path, help='Override figure destination')
    args = ap.parse_args(argv)
    src_folder, display_folder, draw = CONFIG[args.study]
    source = args.source or ROOT / 'Reports' / 'output' / 'full' / src_folder / args.scenario / 'scan.json'
    output = args.output_dir or ROOT / 'Reports' / 'output' / 'display' / display_folder / args.scenario
    source, output = Path(source).resolve(), Path(output).resolve()
    if output == source.parent or output == ROOT.resolve():
        ap.error('Destination cannot be a numerical source directory or project root')
    if not source.is_file():
        ap.error(f'Saved scan missing: {source}')
    rows = load_rows(source)
    if not any(r.get('status') == 'Success' for r in rows):
        ap.error(f'No successful saved scan rows: {source}')
    draw(rows, output)
    print(f'Read-only numerical source: {source}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
