from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent

OUTPUT_DIR = ROOT / "output"
SMOKE_OUTPUT_DIR = OUTPUT_DIR / "smoke"
PIPELINE_SCRIPT = ROOT / "pipeline.py"
C5_REGRESSION_SCRIPT = ROOT / "tests" / "wolfram" / "RegressionC5.wl"
REPORT_FILE = SMOKE_OUTPUT_DIR / "t3_regression_report.json"


def run_command(command: list[str], cwd: Path) -> int:
    """Run one regression stage and return its process status."""
    return subprocess.run(
        command,
        cwd=cwd,
        text=True,
        check=False,
    ).returncode


def run_smoke_matching(verbose: bool) -> int:
    """Regenerate the five historical smoke-model outputs."""
    command = [
        sys.executable,
        str(PIPELINE_SCRIPT),
        "--smoke",
    ]

    if verbose:
        command.append("--verbose")

    return run_command(command, ROOT)


def run_c5_regression() -> int:
    """Validate extracted C5 data."""

    # Remove an old report so a failed Wolfram run cannot accidentally
    # reuse results from an earlier successful regression.
    if REPORT_FILE.exists():
        REPORT_FILE.unlink()

    return run_command(
        [
            "wolframscript",
            "-file",
            str(C5_REGRESSION_SCRIPT),
            str(SMOKE_OUTPUT_DIR),
        ],
        ROOT,
    )


def print_report(report_path: Path) -> None:
    """Print the machine-readable regression report in a compact form."""
    report = json.loads(
        report_path.read_text(encoding="utf-8")
    )

    print("\n" + "=" * 72)
    print("FINAL T3 REGRESSION")
    print("=" * 72)

    print(
        f"{report.get('Status', 'UNKNOWN')}: "
        f"{report.get('Passed', 0)} passed, "
        f"{report.get('Failed', 0)} failed"
    )

    for test in report.get("Tests", []):
        if test.get("Pass", False):
            continue

        detail = test.get("Detail", "")
        suffix = (
            f" — {detail}"
            if detail and detail != '""'
            else ""
        )

        print(
            f"  FAIL: {test.get('Test')}{suffix}"
        )

    print(f"Report: {report_path}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Regression gate for T3 UV matching and C5 extraction."
    )

    parser.add_argument(
        "--no-rematch",
        action="store_true",
        help="reuse existing smoke outputs instead of rerunning matching",
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="show full Matchete output during the smoke rematch",
    )

    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Stage 1: UV construction and matching
    # ------------------------------------------------------------------

    if args.no_rematch:
        print("[1/2] Reusing existing smoke outputs (--no-rematch).")

    else:
        print("[1/2] Running five-model smoke matching...")

        rc = run_smoke_matching(args.verbose)

        if rc != 0:
            print("REGRESSION FAIL: smoke matching failed.")
            return rc

    # ------------------------------------------------------------------
    # Stage 2: C5 regression validation
    # ------------------------------------------------------------------

    print("[2/2] Validating C5 outputs and T3-B scotogenic limit...")

    rc = run_c5_regression()

    if rc != 0:
        print("REGRESSION FAIL: Wolfram C5 validation failed.")
        return rc

    if not REPORT_FILE.exists():
        print("REGRESSION FAIL: Wolfram regression report was not produced.")
        return 2

    print_report(REPORT_FILE)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
