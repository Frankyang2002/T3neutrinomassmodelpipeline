"""Run one full comparison-mode model/scenario for fast debugging.

Example:
    python comparison_test.py --dims 3 3 2 --alpha 0 \
        --scenario largeY_largeL --threshold F --threshold S1 S2
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from studies.FullT3Study import (
    COMPARISON_SCENARIOS,
    DEFAULT_TEMPLATE,
    PROJECT_ROOT,
    _fixed_comparison_config,
    _load_json,
    _model_key,
    _write_json,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run one ordinary T3 model at one common comparison point."
    )
    parser.add_argument(
        "--dims",
        nargs=3,
        type=int,
        required=True,
        metavar=("DS1", "DS2", "DF"),
        help="ordinary T3 dimensions DS1 DS2 DF",
    )
    parser.add_argument(
        "--alpha",
        type=int,
        required=True,
        help="T3 hypercharge parameter",
    )
    parser.add_argument(
        "--scenario",
        choices=tuple(COMPARISON_SCENARIOS),
        required=True,
        help="comparison point to test",
    )
    parser.add_argument(
        "--numerical",
        type=Path,
        default=DEFAULT_TEMPLATE,
        help="numerical template (default: standard T3 template)",
    )
    parser.add_argument(
        "--threshold",
        action="append",
        nargs="+",
        metavar="FIELD",
        default=None,
        help="ordered threshold group; repeat for successive thresholds",
    )
    parser.add_argument(
        "--threshold-scale",
        action="append",
        metavar="SCALE",
        default=None,
        help="matching scale corresponding to each --threshold group",
    )
    parser.add_argument("--debug-reports", action="store_true")
    parser.add_argument("--force", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    ds1, ds2, df = args.dims
    alpha = args.alpha
    yukawa, scalar = COMPARISON_SCENARIOS[args.scenario]

    template_path = args.numerical
    if not template_path.is_absolute():
        template_path = PROJECT_ROOT / template_path
    template = _load_json(template_path)

    key = _model_key(ds1, ds2, df, alpha)
    config_dir = (
        PROJECT_ROOT
        / "configs"
        / "generated_model_sets"
        / "comparison_test"
        / args.scenario
    )
    config_path = config_dir / f"{key}.json"

    payload = _fixed_comparison_config(
        template,
        ds1,
        ds2,
        df,
        alpha,
        yukawa,
        scalar,
    )
    _write_json(config_path, payload)

    study = f"comparison_test/{args.scenario}/{key}"
    command = [
        sys.executable,
        str(PROJECT_ROOT / "pipeline.py"),
        "--dims",
        str(ds1),
        str(ds2),
        str(df),
        "--alpha",
        str(alpha),
        "--study",
        study,
        "--numerical",
        str(config_path),
        "--reset-numerical-configs",
    ]

    if args.threshold:
        for group in args.threshold:
            command.append("--threshold")
            command.extend(group)
    else:
        # The numerical project currently cares about F -> (S1,S2).
        command.extend(["--threshold", "F", "--threshold", "S1", "S2"])

    if args.threshold_scale:
        for scale in args.threshold_scale:
            command.extend(["--threshold-scale", str(scale)])
    if args.debug_reports:
        command.append("--debug-reports")
    if args.force:
        command.append("--force")

    print("=" * 72)
    print("SINGLE COMPARISON TEST")
    print("=" * 72)
    print(f"Model: dims=({ds1},{ds2},{df}), alpha={alpha}")
    print(f"Scenario: {args.scenario}")
    print(f"Yukawa scale Y: {yukawa}")
    print(f"Scalar coupling lambdaT3: {scalar}")
    print(f"Config: {config_path}")
    print(f"Study: {study}")
    print("\n> " + " ".join(command), flush=True)

    return int(subprocess.run(command, cwd=PROJECT_ROOT).returncode)


if __name__ == "__main__":
    raise SystemExit(main())
