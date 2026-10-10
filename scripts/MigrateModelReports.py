"""Migrate repeated symbolic T3 PDFs/TEX into Reports/output/models/.

Examples (project root):
 python scripts/MigrateModelReports.py
 python scripts/MigrateModelReports.py --apply
 python scripts/MigrateModelReports.py --apply --prune-duplicates
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))
from Numerical.orchestration.ModelReportLayout import migrate


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report-root', type=Path, default=PROJECT_ROOT / 'Reports' / 'output')
    parser.add_argument('--apply', action='store_true', help='Copy verified reports to shared model directory')
    parser.add_argument('--prune-duplicates', action='store_true',
                        help='After copying and hash verification, delete only original TEX/PDF duplicates')
    args = parser.parse_args(argv)
    if args.prune_duplicates and not args.apply:
        parser.error('--prune-duplicates requires --apply')
    report = migrate(args.report_root, execute=args.apply, prune=args.prune_duplicates)
    print(json.dumps(report, indent=2))
    return 2 if report['conflicts'] else 0

if __name__ == '__main__':
    raise SystemExit(main())
