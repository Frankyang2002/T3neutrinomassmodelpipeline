from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
W = ROOT / "wolfram"
OUT = W / "output"
PIPELINE = ROOT / "pipeline.py"
WOLFRAM_REGRESSION = W / "RegressionC5.wl"


def run(cmd: list[str], cwd: Path) -> int:
    proc = subprocess.run(cmd, cwd=cwd, text=True)
    return proc.returncode


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Final regression gate for the T3 UV->matching->C5 stage."
    )
    parser.add_argument(
        "--no-rematch",
        action="store_true",
        help="reuse existing wolfram/output smoke results instead of rerunning matching",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="show full Matchete output during the smoke rematch",
    )
    args = parser.parse_args()

    if not args.no_rematch:
        print("[1/2] Running five-model smoke matching...")
        pipeline_cmd = [sys.executable, str(PIPELINE), "--smoke"]
        if args.verbose:
            pipeline_cmd.append("--verbose")
        rc = run(pipeline_cmd, ROOT)
        if rc != 0:
            print("REGRESSION FAIL: smoke matching failed.")
            return rc
    else:
        print("[1/2] Reusing existing smoke outputs (--no-rematch).")

    print("[2/2] Validating C5 outputs and T3-B scotogenic limit...")
    rc = run(
        ["wolframscript", "-file", str(WOLFRAM_REGRESSION), str(OUT)],
        W,
    )

    report_path = OUT / "t3_regression_report.json"
    if report_path.exists():
        report = json.loads(report_path.read_text(encoding="utf-8"))
        print("\n" + "=" * 72)
        print("FINAL T3 REGRESSION")
        print("=" * 72)
        print(
            f"{report.get('Status', 'UNKNOWN')}: "
            f"{report.get('Passed', 0)} passed, {report.get('Failed', 0)} failed"
        )
        if report.get("Failed", 0):
            for test in report.get("Tests", []):
                if not test.get("Pass", False):
                    detail = test.get("Detail", "")
                    suffix = f" — {detail}" if detail and detail != '""' else ""
                    print(f"  FAIL: {test.get('Test')}{suffix}")
        print(f"Report: {report_path}")
    else:
        print("REGRESSION FAIL: Wolfram regression report was not produced.")
        return rc or 2

    return rc


if __name__ == "__main__":
    raise SystemExit(main())
