"""Run a read-only saved-expression D(-2) restricted-Yukawa diagnostic.

No Matchete Match[] is run; this never edits the production calculation.
"""
from __future__ import annotations
import argparse
import json
from datetime import datetime
from pathlib import Path
import subprocess
import sys
from uuid import uuid4

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / "output/general_sequential_pilot/D_m2_stage1"
SCRIPT = ROOT / "Lagrangian/RunWorkerCRestrictedDMinus2Finite.wl"
REQUIRED = (
    "profile_c5_flavour_candidate_hard_symmetric.wxf",
    "profile_c5_flavour_candidate_inherited_symmetric.wxf",
    "profile_c5_flavour_candidate_running_as_supplied.wxf",
    "profile_c5_flavour_candidate_total_provisional.wxf",
    "profile_c5_flavour_candidate_result.json",
    "profile_c5_pole_cancellation_result.json",
)


def preflight(stage_root: Path, output_dir: Path, wolfram_script: Path) -> None:
    stage_root = stage_root.resolve()
    output_dir = output_dir.resolve()
    if not stage_root.is_dir():
        raise ValueError(f"Missing saved first-threshold directory: {stage_root}")
    if not wolfram_script.is_file():
        raise ValueError(f"Missing Wolfram implementation: {wolfram_script}")
    profile = stage_root / "data/general_d_m2_scalar_threshold/profile_tree"
    missing = [name for name in REQUIRED if not (profile / name).is_file()]
    if missing:
        raise ValueError(f"Saved matching files missing: {missing}")
    if output_dir.exists():
        raise ValueError(f"Refusing to overwrite existing output: {output_dir}")
    if stage_root == output_dir or stage_root in output_dir.parents:
        raise ValueError("Output cannot be inside preserved first-threshold directory")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage-root", type=Path, default=DEFAULT)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--preflight", action="store_true", help="Validate inputs without Wolfram")
    args = parser.parse_args()
    stage_root = args.stage_root.resolve()
    if args.output_dir is None:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S") + "_" + uuid4().hex[:8]
        out = ROOT / "output/audits" / f"worker_c_restricted_D_m2_{stamp}"
    else:
        out = args.output_dir.resolve()
    try:
        preflight(stage_root, out, SCRIPT)
    except ValueError as exc:
        print(f"BLOCKED: {exc}", file=sys.stderr)
        return 2
    if args.preflight:
        print("PASS: saved D(-2) candidate files present; fresh output path")
        return 0
    command = ["wolframscript", "-file", str(SCRIPT), str(stage_root), str(out)]
    print("Running restricted saved-expression diagnostic (not new Match[])")
    try:
        status = subprocess.run(command, check=False, cwd=ROOT).returncode
    except FileNotFoundError:
        print("BLOCKED: wolframscript not found", file=sys.stderr)
        return 3
    if status != 0:
        print(f"BLOCKED: Wolfram returned exit code {status}", file=sys.stderr)
        result = out / "result.json"
        if result.is_file():
            print(f"Diagnostic report: {result}", file=sys.stderr)
        return status
    result = out / "result.json"
    if not result.is_file():
        print(f"BLOCKED: Wolfram exited 0 but did not create {result}", file=sys.stderr)
        return 4
    try:
        report = json.loads(result.read_text(encoding="utf-8-sig"))
    except (OSError, ValueError) as exc:
        print(f"BLOCKED: unreadable Wolfram report: {exc}", file=sys.stderr)
        return 5
    if (report.get("status") != "RestrictedNF1PoleFreeCandidateExported"
            or report.get("restricted_total_pole_zero_exact_NF1") is not True
            or report.get("safe_to_add_to_worker_c_joint") is not False
            or not (out / "restricted_total_finite_candidate.wxf").is_file()):
        print("BLOCKED: Wolfram report or finite WXF failed output validation", file=sys.stderr)
        return 6
    print(f"PASS: finite restricted candidate diagnostic exported: {result}")
    print("Independent correction safe to add: False")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
