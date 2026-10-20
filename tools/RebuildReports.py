"""Rebuild only supported threshold/Weinberg plots from saved numerical data.

Project path: tools/RebuildReports.py

Examples (repository root):
  python -m tools.RebuildReports
  python -m tools.RebuildReports --study fermion-threshold
  python -m tools.RebuildReports --study weinberg

Old 'staged-rg', four-scenario comparison, Study A/B display, and standalone
interactive-correlation generation are retired from supported rebuilds.
Does not recompute physics or modify saved scan.json / scan.csv files.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from Numerical.plotting.FermionThresholdVariationFigures import rebuild

ROOT = Path(__file__).resolve().parents[1]
STUDIES = ('fermion-threshold', 'weinberg', 'fitted-mass')


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', choices=STUDIES, action='append',
                        help='Default: fermion-threshold and Weinberg comparisons')
    parser.add_argument('--scenario', default='smallY_smallL',
                        help='Saved-data folder name; not the obsolete four-scenario sweep')
    args = parser.parse_args(argv)
    if not (ROOT / 'pipeline.py').is_file():
        parser.error(f'Run from valid project layout: {ROOT}')
    studies = dict.fromkeys(args.study or ('fermion-threshold', 'weinberg'))
    failed = False
    for study in studies:
        try:
            rebuild(study, args.scenario)
        except (OSError, ValueError, RuntimeError, ImportError, KeyError, TypeError) as exc:
            print(f'[FAIL] {study}: {type(exc).__name__}: {exc}')
            failed = True
    return int(failed)


if __name__ == '__main__':
    raise SystemExit(main())
