"""Safely rename T3 model *presentation* folders to T3-A-m4 etc.

Default dry-run. Only existing model-named directories are renamed, preserving
all their contents. No JSON payloads or physics artifacts are rewritten.

Usage:
    python scripts/RenameT3ModelFolders.py
    python scripts/RenameT3ModelFolders.py --apply

Raw output/ and input configs intentionally retain internal names by default.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from Numerical.orchestration.ModelNames import LONG, shorten

SCENARIOS = ('smallY_smallL', 'smallY_largeL', 'largeY_smallL', 'largeY_largeL')


def locations(root: Path) -> tuple[Path, ...]:
    reports = root / 'Reports' / 'output'
    return (
        reports / 'models',
        reports / 'full' / 'within_model_comparison',
        *(reports / 'full' / 'comparison' / s for s in SCENARIOS),
    )


def plan(root: Path) -> list[tuple[Path, Path]]:
    renames = []
    for parent in locations(root):
        if not parent.is_dir():
            continue
        for old in sorted(parent.iterdir()):
            if old.is_dir() and LONG.fullmatch(old.name):
                renames.append((old, old.with_name(shorten(old.name))))
    destinations = set()
    for old, new in renames:
        if new in destinations or new.exists():
            raise FileExistsError(f'Conflicting destination for {old}: {new}')
        destinations.add(new)
    return renames


def migrate(root: Path, execute: bool = False) -> dict:
    changes = plan(root)
    if execute:
        for old, new in changes:
            old.rename(new)
    return {
        'mode': 'apply' if execute else 'dry-run',
        'count': len(changes),
        'changes': [{'from': str(old), 'to': str(new)} for old, new in changes],
        'note': 'Existing report folders only; source diagnostics and input config paths unchanged.'
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args()
    report = migrate(ROOT, execute=args.apply)
    print(json.dumps(report, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
