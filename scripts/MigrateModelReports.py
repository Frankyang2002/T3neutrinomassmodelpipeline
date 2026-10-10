"""Move repeated symbolic report PDFs/TeX to Reports/output/models.

  python scripts/MigrateModelReports.py                # dry-run
  python scripts/MigrateModelReports.py --apply        # copy+verify+remove originals
  python scripts/MigrateModelReports.py --apply --keep-originals
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
    parser.add_argument("--report-root", type=Path,
                        default=PROJECT_ROOT / "Reports" / "output")
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--keep-originals", action="store_true",
                        help="Copy only; do not remove verified source files")
    args = parser.parse_args(argv)
    report = migrate(args.report_root, execute=args.apply, prune=not args.keep_originals)
    print(json.dumps(report, indent=2))
    return 2 if report["conflicts"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
