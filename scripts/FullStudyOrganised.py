"""Run existing full study and relocate symbolic reports to models/ afterwards.

This does not alter or accelerate the physics calculation.
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
    parser.add_argument("--keep-originals", action="store_true")
    parser.add_argument("pipeline_args", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    extra = list(args.pipeline_args)
    if extra and extra[0] == "--":
        extra.pop(0)
    status = subprocess.run([sys.executable, str(PROJECT_ROOT / "pipeline.py"),
                             "--full", *extra], cwd=PROJECT_ROOT, check=False).returncode
    if status:
        print(f"Full study failed ({status}); migration skipped")
        return status
    report = migrate(PROJECT_ROOT / "Reports" / "output", execute=True,
                     prune=not args.keep_originals)
    print(json.dumps(report, indent=2))
    return 2 if report["conflicts"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
