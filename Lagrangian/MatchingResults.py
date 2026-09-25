"""Interpret Wolfram matching outputs and build Python run records.

This module owns the Python-facing result boundary after ``RunModel.wl`` has
finished. It loads ``comparison_summary.json``, translates shared-scalar
metadata back to physical field language, verifies sequential stage export, and
constructs the stable ``RunRecord``.

It does not launch Wolfram and does not validate the requested T3 model.
"""

from __future__ import annotations

import json
from pathlib import Path

from common.RunRecords import RunRecord
from common.Thresholds import ThresholdPlan
from Lagrangian.WolframRunner import WolframProcessResult


def physicalize_shared_scalar_summary(
    summary: dict,
    *,
    d_s: int,
    d_f: int,
) -> dict:
    """Return a shared-scalar summary expressed in physical field language.

    Wolfram matching always works with the formal T3 roles S1 and S2. For the
    shared-scalar branch those two roles represent one physical scalar S, so
    Python-facing metadata must collapse simultaneous S1/S2 occurrences back
    to S.

    Only summary metadata is rewritten here. Matching expressions and files
    produced by Wolfram are untouched.
    """

    physical_summary = dict(summary)
    physical_summary["SharedScalar"] = True
    physical_summary["PhysicalScalarField"] = "S"
    physical_summary["FormalScalarIdentification"] = "S1=C*S*, S2=S"
    physical_summary["PhysicalDimensions"] = {"dS": d_s, "dF": d_f}

    physical_stages: list[dict] = []

    for raw_stage in summary.get("EFTStages", []):
        stage = dict(raw_stage)
        integrated = list(stage.get("IntegratedFields", []))
        active = list(stage.get("ActiveHeavyFields", []))

        if "S1" in integrated and "S2" in integrated:
            integrated = [
                field
                for field in integrated
                if field not in {"S1", "S2"}
            ] + ["S"]

        if "S1" in active and "S2" in active:
            active = [
                field
                for field in active
                if field not in {"S1", "S2"}
            ] + ["S"]

        stage["IntegratedFields"] = integrated
        stage["ActiveHeavyFields"] = active

        for key in ("Label", "label"):
            if key in stage:
                stage[key] = str(stage[key]).replace("S1_S2", "S")

        physical_stages.append(stage)

    if "EFTStages" in summary:
        physical_summary["EFTStages"] = physical_stages

    return physical_summary


def load_matching_summary(
    output_dir: Path,
    *,
    shared_scalar: bool,
    d_s: int,
    d_f: int,
) -> dict:
    """Load the Wolfram comparison summary or the established failure fallback."""

    summary_path = output_dir / "comparison_summary.json"

    summary = (
        json.loads(summary_path.read_text(encoding="utf-8"))
        if summary_path.exists()
        else {
            "BuildStatus": "ProcessFailed",
            "MatchingStatus": "NotRun",
        }
    )

    if shared_scalar:
        summary = physicalize_shared_scalar_summary(
            summary,
            d_s=d_s,
            d_f=d_f,
        )

    return summary


def validate_sequential_stage_export(
    name: str,
    output_dir: Path,
    summary: dict,
    process: WolframProcessResult,
    threshold_plan: ThresholdPlan,
) -> dict:
    """Check that Wolfram exported every requested sequential threshold stage."""

    requested_stage_count = len(threshold_plan)
    wolfram_stages = summary.get("EFTStages", [])
    sequential_status = summary.get("SequentialMatchingStatus", "Missing")

    print(
        f"  {name}: sequential matching={sequential_status}; "
        f"EFT stages={len(wolfram_stages)}/{requested_stage_count}",
        flush=True,
    )

    if requested_stage_count > 1 and len(wolfram_stages) != requested_stage_count:
        debug_dir = output_dir / "debug"
        debug_dir.mkdir(parents=True, exist_ok=True)

        (debug_dir / "wolfram_stdout.log").write_text(
            process.stdout,
            encoding="utf-8",
        )
        (debug_dir / "wolfram_stderr.log").write_text(
            process.stderr,
            encoding="utf-8",
        )

        print(
            f"ERROR: {name} requested {requested_stage_count} threshold stages "
            f"but Wolfram exported {len(wolfram_stages)}. "
            f"See {debug_dir / 'wolfram_stdout.log'}",
            flush=True,
        )

        summary["MatchingStatus"] = "SequentialStageExportFailed"

    return summary


def build_run_record(
    *,
    name: str,
    alpha: int,
    d_s1: int,
    d_s2: int,
    d_f: int,
    output_dir: Path,
    shared_scalar: bool,
    process: WolframProcessResult,
    threshold_plan: ThresholdPlan,
) -> RunRecord:
    """Load, validate, and package one completed Wolfram model run."""

    summary = load_matching_summary(
        output_dir,
        shared_scalar=shared_scalar,
        d_s=d_s1,
        d_f=d_f,
    )
    summary = validate_sequential_stage_export(
        name,
        output_dir,
        summary,
        process,
        threshold_plan,
    )

    return RunRecord(
        name=name,
        alpha=alpha,
        d_s1=d_s1,
        d_s2=d_s2,
        d_f=d_f,
        return_code=process.returncode,
        summary=summary,
        output_dir=output_dir,
        shared_scalar=shared_scalar,
    )
