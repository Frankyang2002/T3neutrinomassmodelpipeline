"""Run the complete collection of standard T3 comparison studies.

``--full`` is intentionally an orchestration layer around ordinary
``pipeline.py`` invocations.  Each child study therefore uses exactly the same
model construction, matching, RGE, neutrino-physics, and reporting code as a
standalone run.

This module owns only subprocess orchestration and the master summary.  Keeping
those details outside ``pipeline.py`` leaves the central pipeline focused on
the physics/data-flow order for one study.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
PIPELINE_SCRIPT = PROJECT_ROOT / "pipeline.py"
OUTPUT_DIR = PROJECT_ROOT / "output"
REPORT_OUTPUT_DIR = PROJECT_ROOT / "Reports" / "output"


FULL_STUDIES = (
    ("smoke", "--smoke"),
    ("hypercharge", "--hypercharge-comparison"),
    ("dimensions", "--dimension-comparison"),
)


def _forward_common_cli_args(args: argparse.Namespace) -> list[str]:
    """Return options that every child study of ``--full`` should inherit."""

    forwarded: list[str] = []

    if args.debug_reports:
        forwarded.append("--debug-reports")

    if args.force:
        forwarded.append("--force")

    if args.allow_truncated_scalar_first:
        forwarded.append("--allow-truncated-scalar-first")

    eft_max_dimension = int(getattr(args, "eft_max_dimension", 5))
    if eft_max_dimension != 5:
        forwarded.extend(["--eft-max-dimension", str(eft_max_dimension)])

    if args.numerical is not None:
        forwarded.extend(["--numerical", str(args.numerical)])

    if args.threshold:
        for group in args.threshold:
            forwarded.append("--threshold")
            forwarded.extend(group)

    if args.threshold_scale:
        for scale in args.threshold_scale:
            forwarded.extend(["--threshold-scale", str(scale)])

    return forwarded


def _study_model_counts(aggregate: Path) -> tuple[int | None, int | None]:
    """Return ``(successful, total)`` from one child aggregate when available."""

    if not aggregate.is_file():
        return None, None

    try:
        payload = json.loads(aggregate.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None, None

    if not isinstance(payload, list):
        return None, None

    total_models = len(payload)
    successful_models = sum(
        1
        for item in payload
        if isinstance(item, dict)
        and item.get("BuildStatus") == "Success"
        and item.get("MatchingStatus") == "Success"
    )
    return successful_models, total_models


def run_full_study(args: argparse.Namespace) -> int:
    """Run smoke, hypercharge, and dimension studies sequentially.

    Each child is a fresh invocation of ``pipeline.py``.  This preserves the
    historical output layout and avoids a second implementation of the physics
    pipeline inside the ``--full`` mode.
    """

    full_output_dir = OUTPUT_DIR / "full"
    full_report_dir = REPORT_OUTPUT_DIR / "full"
    full_output_dir.mkdir(parents=True, exist_ok=True)
    full_report_dir.mkdir(parents=True, exist_ok=True)

    forwarded = _forward_common_cli_args(args)
    started = time.time()
    study_results: list[dict[str, object]] = []
    overall_status = 0

    print("=" * 72)
    print("FULL T3 STUDY")
    print("=" * 72)
    print(f"Raw output root: {full_output_dir}")
    print(f"Report root: {full_report_dir}")

    for study_name, mode_flag in FULL_STUDIES:
        print("\n" + "=" * 72)
        print(f"FULL STUDY: {study_name}")
        print("=" * 72)

        command = [
            sys.executable,
            str(PIPELINE_SCRIPT),
            mode_flag,
            "--study",
            f"full/{study_name}",
            *forwarded,
        ]

        study_started = time.time()
        completed = subprocess.run(command, cwd=PROJECT_ROOT)
        elapsed = time.time() - study_started
        status = int(completed.returncode)
        overall_status = max(overall_status, status)

        child_output_dir = full_output_dir / study_name
        child_report_dir = full_report_dir / study_name
        aggregate = child_output_dir / "t3_model_comparison.json"
        successful_models, total_models = _study_model_counts(aggregate)

        study_results.append(
            {
                "Study": study_name,
                "ModeFlag": mode_flag,
                "Status": "Success" if status == 0 else "Failed",
                "ReturnCode": status,
                "RuntimeSeconds": elapsed,
                "SuccessfulModels": successful_models,
                "TotalModels": total_models,
                "OutputDirectory": str(child_output_dir),
                "ReportDirectory": str(child_report_dir),
                "AggregateFile": str(aggregate) if aggregate.is_file() else None,
            }
        )

    total_elapsed = time.time() - started
    master_summary = {
        "Status": "Success" if overall_status == 0 else "Failed",
        "RuntimeSeconds": total_elapsed,
        "Studies": study_results,
        "OutputRoot": str(full_output_dir),
        "ReportRoot": str(full_report_dir),
    }

    summary_path = full_output_dir / "full_run_summary.json"
    summary_path.write_text(
        json.dumps(master_summary, indent=2),
        encoding="utf-8",
    )

    print("\n" + "=" * 72)
    print("FULL STUDY SUMMARY")
    print("=" * 72)
    for item in study_results:
        model_text = ""
        if item["TotalModels"] is not None:
            model_text = (
                f"; models={item['SuccessfulModels']}/{item['TotalModels']}"
            )
        print(
            f"{item['Study']}: {item['Status']} "
            f"({item['RuntimeSeconds']:.1f}s{model_text})"
        )

    print(f"Total runtime: {total_elapsed:.1f}s")
    print(f"Master summary: {summary_path}")
    print(f"Reports: {full_report_dir}")

    return overall_status
