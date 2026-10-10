"""Run the existing full study, then consolidate its repeated symbolic PDFs.

Run from repository root:
    python scripts/FullStudyOrganised.py
    python scripts/FullStudyOrganised.py --prune-duplicates

All extra options following -- are forwarded to pipeline.py --full.
The migration is skipped if the full pipeline fails. No physics implementation
is changed; this wrapper only changes the postprocessing workflow.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))
from Numerical.orchestration.ModelReportLayout import migrate


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prune-duplicates', action='store_true',
                        help='Remove scenario-local TEX/PDF after verifying consolidated copies')
    parser.add_argument('pipeline_args', nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    extra = list(args.pipeline_args)
    if extra and extra[0] == '--':
        extra.pop(0)
    command = [sys.executable, str(PROJECT_ROOT / 'pipeline.py'), '--full', *extra]
    result = subprocess.run(command, cwd=PROJECT_ROOT, check=False)
    if result.returncode != 0:
        print(f'Full study exited {result.returncode}; report consolidation skipped.')
        return result.returncode
    report = migrate(PROJECT_ROOT / 'Reports' / 'output', execute=True,
                     prune=args.prune_duplicates)
    print(json.dumps(report, indent=2))
    return 2 if report['conflicts'] else 0

if __name__ == '__main__':
    raise SystemExit(main())
