# Main script for the T3 model pipeline.
#
# This file:
# 1. Reads the requested T3 representation assignment / benchmark model.
# 2. Builds the UV model and performs EFT matching.
# 3. Organises the matched Weinberg coefficient C5.
# 4. Runs the UV and EFT RGE stages.
# 5. Constructs the neutrino-mass matrix.
# 6. Optionally performs numerical running and calculates neutrino observables.
# 7. Generates the final reports and run summaries.
#
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

from RGE.matching.MatchedEFTRGE import run_matched_eft_rge
from RGE.stages.FlavorMatchedRGEStage import (
    run_flavor_matched_rge,
    run_neutrino_mass_stage,
)
from RGE.stages.NumericalPipelineStage import run_numerical_pipeline_stage
from RGE.phenomenology.NeutrinoObservables import run_neutrino_observables_stage
from RGE.running.RGBetaT3Running import (
    run_rgbeta_t3,
    run_rgbeta_t3_eft1,
)
from RGE.running.EFT1WilsonRGE import run_eft1_wilson_rge
from RGE.running.EFT1WilsonRunning import run_eft1_wilson_transport
from RGE.running.EFT1WilsonFlavorSeed import run_flavor_seed_export
from RGE.running.EFT1DirectWeinberg import (
    build_direct_weinberg_flavor_transport,
    export_direct_weinberg_matchete,
)
from RGE.running.EFT1ThresholdResume import rerun_threshold2_with_running
from RGE.running.FinalWeinbergCoefficient import (
    build_final_weinberg_coefficient,
    normalize_pole_rge_consistency,
)

from common.Records import RunRecord
from common.Thresholds import (
    build_eft_stage_records,
    threshold_plan_label,
    threshold_plan_to_json,
    validate_threshold_plan,
)
from common.T3Model import EXTENDED, INTERESTING, SMOKE
from Lagrangian.Runner import validate_dimensions, validate_shared_dimensions, obtain_class_dimensions
from Reports.ReportGeneration import (
    compile_latex_document,
    write_bsm_uv_field_table,
    write_bsm_matched_field_table,
    write_c5_coefficient_report,
)
from Reports.RGEComparison import (
    write_and_compile_eft1_rge_comparison,
    write_and_compile_final_eft_rge_comparison,
    write_and_compile_rge_comparison,
)
from Reports.GroupFactorReports import (
    write_and_compile_stage_group_factor_reports,
)


PROJECT_ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = PROJECT_ROOT / "output"
REPORT_OUTPUT_DIR = PROJECT_ROOT / "Reports" / "output"
RUN_THRESHOLD_STAGE_SCRIPT = PROJECT_ROOT / "Lagrangian" / "RunThresholdStage.wl"


# ---------------------------------------------------------------------------
# Comparison-study definitions
# ---------------------------------------------------------------------------
# smoke: regression check on the five historical benchmark points
# dimensions: compare the five SU(2) assignments at fixed alpha = 0
# hypercharge: scan alpha while keeping each SU(2) assignment fixed
DIMENSION_COMPARISON = tuple(
    (model_class, 0)
    for model_class in ("A", "B", "C", "D", "E")
)

HYPERCHARGE_ALPHAS = (-2, -1, 0, 1, 2)
HYPERCHARGE_COMPARISON = tuple(
    (model_class, alpha)
    for model_class in ("A", "B", "C", "D", "E")
    for alpha in HYPERCHARGE_ALPHAS
)


def _forward_common_cli_args(args: argparse.Namespace) -> list[str]:
    """Return options that each child study of --full should inherit."""

    forwarded: list[str] = []

    if args.debug_reports:
        forwarded.append("--debug-reports")

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


def run_full_study(args: argparse.Namespace) -> int:
    """Run every comparison study sequentially and write a master summary.

    --full is intentionally only an orchestrator: each child invocation goes
    through the ordinary pipeline, so there is one implementation of the
    physics stages and report generation.
    """

    full_output_dir = OUTPUT_DIR / "full"
    full_report_dir = REPORT_OUTPUT_DIR / "full"
    full_output_dir.mkdir(parents=True, exist_ok=True)
    full_report_dir.mkdir(parents=True, exist_ok=True)

    studies = (
        ("smoke", "--smoke"),
        ("hypercharge", "--hypercharge-comparison"),
        ("dimensions", "--dimension-comparison"),
    )

    forwarded = _forward_common_cli_args(args)
    started = time.time()
    study_results: list[dict[str, object]] = []
    overall_status = 0

    print("=" * 72)
    print("FULL T3 STUDY")
    print("=" * 72)
    print(f"Raw output root: {full_output_dir}")
    print(f"Report root: {full_report_dir}")

    for study_name, mode_flag in studies:
        print("\n" + "=" * 72)
        print(f"FULL STUDY: {study_name}")
        print("=" * 72)

        command = [
            sys.executable,
            str(Path(__file__).resolve()),
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

        successful_models = None
        total_models = None
        if aggregate.is_file():
            try:
                payload = json.loads(aggregate.read_text(encoding="utf-8"))
                if isinstance(payload, list):
                    total_models = len(payload)
                    successful_models = sum(
                        1
                        for item in payload
                        if isinstance(item, dict)
                        and item.get("BuildStatus") == "Success"
                        and item.get("MatchingStatus") == "Success"
                    )
            except (OSError, json.JSONDecodeError):
                pass

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


def run_uv_rgbeta_stage(record: RunRecord) -> bool:
    """ Generate the one-loop beta functions for the renormalisable UV T3 model
    using RGBeta.

    This stage describes running above the heavy-particle matching threshold,
    before the T3 fields are integrated out."""
    summary = record.summary

    # Fails if we dont have an actual build for our model
    if summary.get("BuildStatus") != "Success":
        summary["UVRGEStatus"] = "NotRun"
        return False

    data_dir = record.output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    output_path = data_dir / "uv_rgbeta_rge.json"

    print(f"  {record.name}: starting RGBeta UV-RGE stage...", flush=True)

    # Use run_rgbeta_t3
    try:
        result = run_rgbeta_t3(
            record.d_s1,
            record.d_s2,
            record.d_f,
            record.alpha,
            shared_scalar=record.shared_scalar,
        )
    except Exception as exc:
        summary["UVRGEStatus"] = "Failed"
        summary["UVRGEError"] = str(exc)
        print(f"  {record.name}: RGBeta UV-RGE failed: {exc}")
        return False

    # Pass the T3 SU(2) representation dimensions and hypercharge parameter
    # to the RGBeta Wolfram runner, which constructs the UV model beta functions.
    output_path.write_text(
        json.dumps(result.raw, indent=2),
        encoding="utf-8",
    )

    summary["UVRGEStatus"] = result.status
    summary["UVRGEFile"] = output_path.relative_to(record.output_dir).as_posix()
    summary["UVRGEBetaCount"] = len(result.betas)

    print(
        f"  {record.name}: RGBeta UV-RGE={result.status} "
        f"({len(result.betas)} beta functions) -> {output_path}"
    )

    return result.status == "Success"


def _physical_scalar_fields(record: RunRecord) -> set[str]:
    """Return the physical scalar fields active after integrating out F."""
    return {"S"} if record.shared_scalar else {"S1", "S2"}


def _has_f_first_eft_stage(record: RunRecord) -> bool:
    """Return whether the first matched EFT is obtained by integrating out F."""
    stages = record.summary.get("EFTStages", [])
    return bool(
        stages
        and stages[0].get("IntegratedFields") == ["F"]
        and set(stages[0].get("ActiveHeavyFields", []))
        == _physical_scalar_fields(record)
    )


def run_eft1_rgbeta_stage(record: RunRecord) -> bool:
    """Generate the renormalisable one-loop RGEs in EFT1 = SM + S1 + S2.

    This is the dimension-four part of the intermediate EFT after integrating
    out F.  The higher-dimensional Wilson operators generated at the F
    threshold are deliberately not claimed here; their running is a separate
    stage that must be combined with these beta functions before the EFT1 PDF
    is considered complete.
    """

    summary = record.summary

    if not _has_f_first_eft_stage(record):
        summary["EFT1RenormalisableRGEStatus"] = "NotApplicable"
        return True

    data_dir = record.output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    output_path = data_dir / "eft1_rgbeta_rge.json"

    print(
        f"  {record.name}: starting RGBeta EFT1 renormalisable RGE stage...",
        flush=True,
    )

    try:
        result = run_rgbeta_t3_eft1(
            record.d_s1,
            record.d_s2,
            record.d_f,
            record.alpha,
            shared_scalar=record.shared_scalar,
        )
    except Exception as exc:
        summary["EFT1RenormalisableRGEStatus"] = "Failed"
        summary["EFT1RenormalisableRGEError"] = str(exc)
        print(
            f"  {record.name}: RGBeta EFT1 renormalisable RGE failed: {exc}"
        )
        return False

    output_path.write_text(
        json.dumps(result.raw, indent=2),
        encoding="utf-8",
    )

    summary["EFT1RenormalisableRGEStatus"] = result.status
    summary["EFT1RenormalisableRGEFile"] = (
        output_path.relative_to(record.output_dir).as_posix()
    )
    summary["EFT1RenormalisableRGEBetaCount"] = len(result.betas)

    print(
        f"  {record.name}: RGBeta EFT1 renormalisable RGE={result.status} "
        f"({len(result.betas)} beta functions) -> {output_path}"
    )

    return result.status == "Success"


def run_eft1_wilson_rge_stage(record: RunRecord) -> bool:
    """Generate the one-loop dimension-five Wilson RGE in EFT1 = SM + S1 + S2.

    This stage is applicable when the first threshold integrates out F while
    leaving S1 and S2 active.  It combines the exact tree-level Wilson seed,
    the exact scalar-quartic seed, and the EFT1 renormalisable model metadata
    with the general psi^2 phi^2 master RGE.

    Fixed-order bookkeeping:
        only C^(0), the tree-generated stage-1 Wilson tensor, is evolved.
        The one-loop matching contribution at the F threshold remains a
        boundary term and is not inserted into this one-loop beta function.
    """

    summary = record.summary
    if not _has_f_first_eft_stage(record):
        summary["EFT1WilsonRGEStatus"] = "NotApplicable"
        return True

    data_dir = record.output_dir / "data"
    wilson_seed_path = data_dir / "eft1_after_F_wilson_seed.json"
    quartic_seed_path = data_dir / "eft1_after_F_scalar_quartic_seed.json"
    rgbeta_path = data_dir / "eft1_rgbeta_rge.json"
    output_path = data_dir / "eft1_wilson_rge.json"

    required_inputs = (
        wilson_seed_path,
        quartic_seed_path,
        rgbeta_path,
    )
    missing = [path for path in required_inputs if not path.is_file()]

    if missing:
        summary["EFT1WilsonRGEStatus"] = "Failed"
        summary["EFT1WilsonRGEError"] = (
            "Missing EFT1 Wilson-RGE input(s): "
            + ", ".join(
                path.relative_to(record.output_dir).as_posix()
                for path in missing
            )
        )
        print(
            f"  {record.name}: EFT1 Wilson RGE failed: "
            f"{summary['EFT1WilsonRGEError']}"
        )
        return False

    print(
        f"  {record.name}: starting EFT1 dimension-five Wilson RGE stage...",
        flush=True,
    )

    try:
        result = run_eft1_wilson_rge(
            wilson_seed_path=wilson_seed_path,
            quartic_seed_path=quartic_seed_path,
            rgbeta_path=rgbeta_path,
            output_path=output_path,
            seed_support_only=False,
        )
    except Exception as exc:
        summary["EFT1WilsonRGEStatus"] = "Failed"
        summary["EFT1WilsonRGEError"] = str(exc)
        print(
            f"  {record.name}: EFT1 dimension-five Wilson RGE failed: {exc}"
        )
        return False

    summary["EFT1WilsonRGEStatus"] = result.get("status", "Unknown")
    summary["EFT1WilsonRGEFile"] = (
        output_path.relative_to(record.output_dir).as_posix()
    )
    summary["EFT1WilsonRGEInitialComponentCount"] = result.get(
        "initial_independent_component_count",
        0,
    )
    summary["EFT1WilsonRGEBetaComponentCount"] = result.get(
        "nonzero_beta_component_count",
        0,
    )
    summary["EFT1WilsonRGEGeneratedComponentCount"] = result.get(
        "generated_component_count",
        0,
    )

    weinberg_validation = result.get("weinberg_subspace_validation") or {}
    summary["EFT1WilsonRGEWeinbergSubspaceValid"] = (
        weinberg_validation.get("matches_weinberg_subspace")
    )
    summary["EFT1WilsonRGEWeinbergBeta"] = (
        weinberg_validation.get("beta_kappa_16pi2", "")
    )

    # Mirror the calculation metadata onto the actual EFT1 stage record so
    # stage-aware reporting can consume it without rediscovering files.
    if record.eft_stages:
        eft1_stage = record.first_eft_stage
        eft1_stage.summary.update(
            {
                "RenormalisableRGEStatus": summary.get(
                    "EFT1RenormalisableRGEStatus",
                    "Unknown",
                ),
                "RenormalisableRGEFile": summary.get(
                    "EFT1RenormalisableRGEFile",
                    "",
                ),
                "WilsonRGEStatus": summary["EFT1WilsonRGEStatus"],
                "WilsonRGEFile": summary["EFT1WilsonRGEFile"],
                "WilsonRGEInitialComponentCount": summary[
                    "EFT1WilsonRGEInitialComponentCount"
                ],
                "WilsonRGEBetaComponentCount": summary[
                    "EFT1WilsonRGEBetaComponentCount"
                ],
                "WilsonRGEGeneratedComponentCount": summary[
                    "EFT1WilsonRGEGeneratedComponentCount"
                ],
                "WeinbergSubspaceValid": summary[
                    "EFT1WilsonRGEWeinbergSubspaceValid"
                ],
                "WeinbergBeta16Pi2": summary[
                    "EFT1WilsonRGEWeinbergBeta"
                ],
            }
        )

    print(
        f"  {record.name}: EFT1 Wilson RGE="
        f"{summary['EFT1WilsonRGEStatus']} "
        f"({summary['EFT1WilsonRGEBetaComponentCount']} nonzero beta "
        f"components; "
        f"{summary['EFT1WilsonRGEGeneratedComponentCount']} generated)"
        f" -> {output_path}"
    )

    if weinberg_validation:
        print(
            f"  {record.name}: EFT1 Weinberg-subspace validation="
            f"{weinberg_validation.get('matches_weinberg_subspace')} "
            f"(16*pi^2 beta_kappa="
            f"{weinberg_validation.get('beta_kappa_16pi2', '')})"
        )

    return summary["EFT1WilsonRGEStatus"] == "Success"


def _default_threshold_scale(group: list[str]) -> str:
    """Return a symbolic scale name for one threshold group."""

    fields = set(group)
    if fields == {"F"}:
        return "MF"
    if fields == {"S"}:
        return "MS"
    if fields == {"S1", "S2"}:
        return "MS"
    if fields == {"S1"}:
        return "MS1"
    if fields == {"S2"}:
        return "MS2"

    ordered = "_".join(group)
    return f"M_{ordered}"


def run_eft1_wilson_transport_stage(
    record: RunRecord,
    *,
    mu_high: str,
    mu_low: str,
) -> bool:
    """Transport tree-generated EFT1 Wilson coefficients to threshold 2.

    This is fixed-order one-loop transport only:

        C(mu_low) = C^(0)
                    + hbar log(mu_low/mu_high) beta^(1)[C^(0)]
                    + O(hbar^2).

    The stage-1 one-loop matching boundary contribution is intentionally not
    evolved, because that would first affect the calculation at O(hbar^2).
    """

    summary = record.summary

    if summary.get("EFT1WilsonRGEStatus") != "Success":
        summary["EFT1WilsonTransportStatus"] = "NotRun"
        return False

    data_dir = record.output_dir / "data"
    wilson_seed_path = data_dir / "eft1_after_F_wilson_seed.json"
    wilson_rge_path = data_dir / "eft1_wilson_rge.json"
    rgbeta_path = data_dir / "eft1_rgbeta_rge.json"
    output_path = data_dir / "eft1_wilson_at_S_threshold.json"

    required_inputs = (
        wilson_seed_path,
        wilson_rge_path,
        rgbeta_path,
    )
    missing = [path for path in required_inputs if not path.is_file()]
    if missing:
        summary["EFT1WilsonTransportStatus"] = "Failed"
        summary["EFT1WilsonTransportError"] = (
            "Missing EFT1 Wilson-transport input(s): "
            + ", ".join(
                path.relative_to(record.output_dir).as_posix()
                for path in missing
            )
        )
        print(
            f"  {record.name}: EFT1 Wilson transport failed: "
            f"{summary['EFT1WilsonTransportError']}"
        )
        return False

    print(
        f"  {record.name}: transporting EFT1 Wilson coefficients "
        f"{mu_high} -> {mu_low}...",
        flush=True,
    )

    try:
        result = run_eft1_wilson_transport(
            wilson_seed_path=wilson_seed_path,
            wilson_rge_path=wilson_rge_path,
            rgbeta_path=rgbeta_path,
            mu_high=mu_high,
            mu_low=mu_low,
            output_path=output_path,
        )
    except Exception as exc:
        summary["EFT1WilsonTransportStatus"] = "Failed"
        summary["EFT1WilsonTransportError"] = str(exc)
        print(
            f"  {record.name}: EFT1 Wilson transport failed: {exc}"
        )
        return False

    summary["EFT1WilsonTransportStatus"] = result.get("status", "Unknown")
    summary["EFT1WilsonTransportFile"] = (
        output_path.relative_to(record.output_dir).as_posix()
    )
    summary["EFT1WilsonTransportMuHigh"] = mu_high
    summary["EFT1WilsonTransportMuLow"] = mu_low
    summary["EFT1WilsonTransportLogRatio"] = (
        result.get("scales", {}).get("log_ratio", "")
    )
    summary["EFT1WilsonTransportCorrectionCount"] = result.get(
        "running_correction_component_count",
        0,
    )
    summary["EFT1WilsonTransportGeneratedCount"] = result.get(
        "generated_by_running_component_count",
        0,
    )
    summary["EFT1WilsonTransportEqualScaleCheck"] = result.get(
        "equal_scale_running_vanishes"
    )

    if record.eft_stages:
        record.first_eft_stage.summary.update(
            {
                "WilsonTransportStatus": summary[
                    "EFT1WilsonTransportStatus"
                ],
                "WilsonTransportFile": summary[
                    "EFT1WilsonTransportFile"
                ],
                "WilsonTransportMuHigh": mu_high,
                "WilsonTransportMuLow": mu_low,
                "WilsonTransportLogRatio": summary[
                    "EFT1WilsonTransportLogRatio"
                ],
                "WilsonTransportCorrectionCount": summary[
                    "EFT1WilsonTransportCorrectionCount"
                ],
                "WilsonTransportGeneratedCount": summary[
                    "EFT1WilsonTransportGeneratedCount"
                ],
                "WilsonTransportEqualScaleCheck": summary[
                    "EFT1WilsonTransportEqualScaleCheck"
                ],
            }
        )

    print(
        f"  {record.name}: EFT1 Wilson transport="
        f"{summary['EFT1WilsonTransportStatus']} "
        f"({summary['EFT1WilsonTransportCorrectionCount']} running "
        f"components; "
        f"{summary['EFT1WilsonTransportGeneratedCount']} generated; "
        f"log={summary['EFT1WilsonTransportLogRatio']}) "
        f"-> {output_path}"
    )

    return summary["EFT1WilsonTransportStatus"] == "Success"



def _is_f_then_s1_s2_plan(record: RunRecord) -> bool:
    """Return whether recorded EFT stages realize the F -> scalar hierarchy."""
    stages = record.summary.get("EFTStages", [])
    scalar_fields = _physical_scalar_fields(record)
    return bool(
        _has_f_first_eft_stage(record)
        and len(stages) >= 2
        and set(stages[1].get("IntegratedFields", [])) == scalar_fields
    )


def _is_f_then_scalar_threshold_plan(
    threshold_plan,
    record: RunRecord,
) -> bool:
    """Return whether the requested physical threshold plan is F -> scalar(s)."""
    scalar_fields = _physical_scalar_fields(record)
    return bool(
        len(threshold_plan) >= 2
        and tuple(threshold_plan[0]) == ("F",)
        and set(threshold_plan[1]) == scalar_fields
    )


def run_eft1_full_flavor_threshold_bridge(
    record: RunRecord,
    *,
    mu_high: str,
    mu_low: str,
    debug_reports: bool = False,
) -> bool:
    """Build the authoritative fixed-one-loop hierarchical C5.

    The tree LLSS Wilson boundary is O(hbar^0). Its direct mixing into O5 is
    O(hbar^1). LLSS self-running is also O(hbar^1), so matching that corrected
    heavy operator through the scalar loop would be O(hbar^2). Therefore only
    the direct full-flavor C12 -> Weinberg running enters the authoritative
    one-loop neutrino-mass calculation.
    """

    summary = record.summary

    if not _is_f_then_s1_s2_plan(record):
        summary["EFT1FullFlavorBridgeStatus"] = "NotApplicable"
        return True

    data_dir = record.output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    wilson_seed_path = data_dir / "eft1_after_F_wilson_seed.json"
    component_rge_path = data_dir / "eft1_wilson_rge.json"
    flavor_seed_path = data_dir / "eft1_wilson_flavor_seed.json"
    flavor_transport_path = data_dir / "eft1_wilson_flavor_at_S_threshold.json"
    insertion_path = data_dir / "eft1_wilson_flavor_running_insertion.wl"

    required = (wilson_seed_path, component_rge_path)
    missing = [path for path in required if not path.is_file()]
    if missing:
        summary["EFT1FullFlavorBridgeStatus"] = "Failed"
        summary["EFT1FullFlavorBridgeError"] = (
            "Missing direct-Weinberg bridge input(s): "
            + ", ".join(
                path.relative_to(record.output_dir).as_posix()
                for path in missing
            )
        )
        print(
            f"  {record.name}: full-flavor EFT1 bridge failed: "
            f"{summary['EFT1FullFlavorBridgeError']}"
        )
        return False

    print(
        f"  {record.name}: starting fixed-one-loop direct Weinberg bridge "
        f"{mu_high} -> {mu_low}...",
        flush=True,
    )

    try:
        flavor_seed = run_flavor_seed_export(
            wilson_seed_path=wilson_seed_path,
            output_path=flavor_seed_path,
        )
        if flavor_seed.get("status") != "Success":
            raise RuntimeError("Full-flavor Wilson seed regression failed.")

        flavor_transport = build_direct_weinberg_flavor_transport(
            flavor_seed_path=flavor_seed_path,
            component_rge_path=component_rge_path,
            mu_high=mu_high,
            mu_low=mu_low,
            output_path=flavor_transport_path,
        )
        if flavor_transport.get("status") != "Success":
            raise RuntimeError(
                "Direct full-flavor Weinberg transport regression failed."
            )

        insertion = export_direct_weinberg_matchete(
            transport_path=flavor_transport_path,
            output_path=insertion_path,
        )
        if insertion.get("status") != "Success":
            raise RuntimeError(
                "Direct Weinberg Matchete insertion export failed."
            )

        resume = rerun_threshold2_with_running(
            output_dir=record.output_dir,
            running_insertion=insertion_path,
            run_threshold_script=RUN_THRESHOLD_STAGE_SCRIPT,
            d_s1=record.d_s1,
            d_s2=record.d_s2,
            d_f=record.d_f,
            alpha=record.alpha,
            shared_scalar=record.shared_scalar,
            validation_mode=debug_reports,
        )
        if resume.get("status") != "Success":
            raise RuntimeError(
                "Threshold-2 rerun with direct Weinberg running failed."
            )

        pole_rge_path = data_dir / "c5_pole_rge_consistency.json"
        if pole_rge_path.is_file():
            pole_rge = normalize_pole_rge_consistency(pole_rge_path)
            if not pole_rge.get("ConsistencyValidated", False):
                raise RuntimeError(
                    "C5 pole/RGE consistency failed strict regression: "
                    "direct log coefficient must equal 2 x hard pole residue; "
                    f"reported ratio="
                    f"{pole_rge.get('DirectLogToPoleRatioOneGenerationInputForm')}."
                )

    except Exception as exc:
        summary["EFT1FullFlavorBridgeStatus"] = "Failed"
        summary["EFT1FullFlavorBridgeError"] = str(exc)
        print(f"  {record.name}: full-flavor EFT1 bridge failed: {exc}")
        return False

    authoritative_threshold_c5 = data_dir / "c5_threshold_with_eft1_running.txt"
    authoritative_direct_c5 = data_dir / "c5_direct_eft1_running.txt"
    final_c5_bookkeeping_path = data_dir / "final_weinberg_coefficient.json"

    if not authoritative_threshold_c5.is_file():
        summary["EFT1FullFlavorBridgeStatus"] = "Failed"
        summary["EFT1FullFlavorBridgeError"] = (
            "Authoritative resumed threshold C5 was not exported: "
            f"{authoritative_threshold_c5}"
        )
        print(
            f"  {record.name}: full-flavor EFT1 bridge failed: "
            f"{summary['EFT1FullFlavorBridgeError']}"
        )
        return False

    try:
        final_c5 = build_final_weinberg_coefficient(
            threshold_c5_path=authoritative_threshold_c5,
            flavor_transport_path=flavor_transport_path,
            output_path=final_c5_bookkeeping_path,
        )
    except Exception as exc:
        summary["EFT1FullFlavorBridgeStatus"] = "Failed"
        summary["EFT1FullFlavorBridgeError"] = (
            f"Final Weinberg coefficient bookkeeping failed: {exc}"
        )
        print(
            f"  {record.name}: full-flavor EFT1 bridge failed: "
            f"{summary['EFT1FullFlavorBridgeError']}"
        )
        return False

    if final_c5.get("status") != "Success":
        summary["EFT1FullFlavorBridgeStatus"] = "Failed"
        summary["EFT1FullFlavorBridgeError"] = (
            "Final Weinberg coefficient bookkeeping returned non-success."
        )
        return False

    combined = final_c5.get("combined", {})
    if not combined.get("ready_for_physical_majorana_numerics", False):
        summary["EFT1FullFlavorBridgeStatus"] = "Failed"
        summary["EFT1FullFlavorBridgeError"] = (
            "Final Weinberg coefficient is not ready for physical Majorana "
            f"downstream use: {combined.get('physical_majorana_reason')!r}"
        )
        return False

    # From this point onward every EFT-side stage must consume the authoritative
    # hierarchical full-flavor coefficient, never the preliminary
    # c5_coefficient.txt extracted before the resumed threshold calculation.
    summary["WeinbergCoefficientFile"] = (
        final_c5_bookkeeping_path.relative_to(record.output_dir).as_posix()
    )

    summary.update(
        {
            "EFT1FlavorSeedStatus": flavor_seed.get("status"),
            "EFT1FlavorSeedFile": flavor_seed_path.relative_to(
                record.output_dir
            ).as_posix(),
            "EFT1FlavorSeedOneGenerationCheck": flavor_seed.get(
                "one_generation_reduction_matches"
            ),
            "EFT1FlavorRGEStatus": "DiagnosticNotRequiredForOneLoopC5",
            "EFT1FlavorRGEOneGenerationCheck": None,
            "EFT1FlavorTransportStatus": flavor_transport.get("status"),
            "EFT1FlavorTransportFile": flavor_transport_path.relative_to(
                record.output_dir
            ).as_posix(),
            "EFT1FlavorTransportEqualScaleCheck": flavor_transport.get(
                "equal_scale_running_vanishes"
            ),
            "EFT1FlavorRunningInsertionStatus": insertion.get("status"),
            "EFT1FlavorRunningInsertionFile": insertion_path.relative_to(
                record.output_dir
            ).as_posix(),
            "EFT1HeavySelfRunningIncludedInAuthoritativeC5": False,
            "EFT1DirectWeinbergOnlyAtOneLoop": True,
            "EFT1ThresholdResumeStatus": resume.get("status"),
            "EFT1ThresholdResumeValidationMode": resume.get(
                "validation_mode",
                False,
            ),
            "EFT1ThresholdResumeResultFile": str(
                Path(resume["result_path"]).resolve()
            ),
            "EFT1ThresholdResumeInsertionLoaded": resume.get(
                "running_insertion_loaded"
            ),
            "EFT1ThresholdResumeInsertedInCOnly": resume.get(
                "running_inserted_in_C_only"
            ),
            "EFT1ThresholdResumeDirectWeinbergCarried": resume.get(
                "direct_weinberg_carried_separately"
            ),
            "EFT1ThresholdResumeDirectWeinbergEqualScaleCheck": resume.get(
                "direct_weinberg_equal_scale_vanishes"
            ),
            "EFT1FullFlavorBridgeStatus": "Success",
            "AuthoritativeThresholdC5File":
                authoritative_threshold_c5.relative_to(
                    record.output_dir
                ).as_posix(),
            "AuthoritativeDirectRunningC5File": (
                authoritative_direct_c5.relative_to(
                    record.output_dir
                ).as_posix()
                if authoritative_direct_c5.is_file()
                else ""
            ),
            "FinalWeinbergCoefficientFile": (
                final_c5_bookkeeping_path.relative_to(
                    record.output_dir
                ).as_posix()
            ),
            "AuthoritativeFinalC5Status": "ReadyForDownstream",
            "AuthoritativeFinalC5FullFlavorReady": combined.get(
                "ready_for_full_flavor_numerics"
            ),
            "AuthoritativeFinalC5PhysicalMajoranaReady": combined.get(
                "ready_for_physical_majorana_numerics"
            ),
            "AuthoritativeDownstreamC5File": (
                final_c5_bookkeeping_path.relative_to(
                    record.output_dir
                ).as_posix()
            ),
        }
    )

    print(
        f"  {record.name}: direct one-loop Weinberg bridge=Success "
        f"(heavy self-running excluded at O(hbar); "
        f"equal-scale check="
        f"{summary['EFT1FlavorTransportEqualScaleCheck']})",
        flush=True,
    )

    return True


def print_summary(records: list[RunRecord], aggregate_dir: Path | None = None) -> int:
    """Print the results of all completed T3 runs."""

    print(
        "\n"
        + "=" * 72
        + "\nT3 MODEL SUMMARY\n"
        + "=" * 72
    )

    successful = 0

    for record in records:
        summary = record.summary

        build_ok = summary.get("BuildStatus") == "Success"
        match_ok = summary.get("MatchingStatus") == "Success"
        uv_rge_ok = summary.get("UVRGEStatus") == "Success"

        successful += int(build_ok and match_ok and uv_rge_ok)

        print(
            f"{record.name} alpha={record.alpha}: "
            f"build={summary.get('BuildStatus')}, "
            f"UV-RGE={summary.get('UVRGEStatus')}, "
            f"match={summary.get('MatchingStatus')}, "
            f"T3={summary.get('T3IngredientsPresent')}, "
            f"Weinberg={summary.get('WeinbergOperatorPresent')}"
        )

        # If matching did not generate the Weinberg operator, there is
        # no C5 coefficient to extract.
        if not summary.get("WeinbergOperatorPresent"):
            continue

        extraction = summary.get(
            "WeinbergExtractionStatus",
            "Unknown",
        )

        if extraction == "Success":
            n_holo = summary.get(
                "WeinbergHolomorphicTermCount",
                0,
            )
            n_hc = summary.get(
                "WeinbergConjugateTermCount",
                0,
            )

            print(
                f"  C5 extraction: Success "
                f"({n_holo} holomorphic + {n_hc} HC terms)"
            )

            if coefficient_file := summary.get("WeinbergCoefficientFile"):
                print(
                    f"  C5: "
                    f"{record.output_dir / coefficient_file}"
                )

        else:
            print(
                f"  Weinberg terms: "
                f"{summary.get('WeinbergTermCount', 0)}; "
                f"C5: {extraction}"
            )

    # Save all model summaries together so we can easily compare runs
    # or use them for regression testing.
    aggregate_dir = aggregate_dir or OUTPUT_DIR
    aggregate = aggregate_dir / "t3_model_comparison.json"

    aggregate.write_text(
        json.dumps(
            [record.summary for record in records],
            indent=2,
        ),
        encoding="utf-8",
    )

    print(
        f"\n{successful}/{len(records)} completed build+matching."
        f"\nAggregate: {aggregate}"
    )

    # Standard command-line convention:
    #   0 = everything succeeded
    #   1 = at least one model failed
    return 0 if successful == len(records) else 1


def organise_c5_input(record: RunRecord) -> Path | None:
    """
    Put the matched Weinberg coefficient in its standard machine-readable
    location before any EFT RGE stage starts.
    From this point onward, c5_coefficient.txt represents C5 at the
    matching scale M and is the main input to the EFT-side calculations.
    """
    
    summary = record.summary
    coefficient_file = summary.get("WeinbergCoefficientFile")

    if not coefficient_file:
        return None

    c5_path = record.output_dir / coefficient_file

    if not c5_path.is_file():
        return None

    data_dir = record.output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    organised_path = data_dir / c5_path.name

    if c5_path != organised_path:
        if organised_path.exists():
            raise FileExistsError(
                f"Refusing to overwrite existing RGE input: {organised_path}"
            )

        c5_path.replace(organised_path)
        c5_path = organised_path

    summary["WeinbergCoefficientFile"] = (
        c5_path.relative_to(record.output_dir).as_posix()
    )
    return c5_path


def matched_c5_path(record: RunRecord) -> Path | None:
    """Return the organised matched C5 coefficient for the RGE stages."""

    coefficient_file = record.summary.get("WeinbergCoefficientFile")

    if not coefficient_file:
        return None

    c5_path = record.output_dir / coefficient_file

    if not c5_path.is_file():
        return None

    return c5_path


def run_matched_eft_rge_stage(
    record: RunRecord,
    debug_reports: bool = False,
) -> bool:
    """
    Evaluate the matched Weinberg coefficient using the general psi^2 phi^2
    one-loop RGE machinery and verify that it reduces to the known
    one-generation SMEFT Weinberg RGE. 
    Just for testing as one generation is simple to test
    """

    summary = record.summary
    c5_path = matched_c5_path(record)

    if c5_path is None:
        summary["RGEStatus"] = "NotRun"
        return False

    print(f"  {record.name}: starting matched-EFT RGE stage...", flush=True)

    try:
        rge_summary = run_matched_eft_rge(
            c5_path=c5_path,
            output_dir=record.output_dir,
            debug_outputs=debug_reports,
        )
    except Exception as exc:
        summary["RGEStatus"] = "Failed"
        summary["RGEError"] = str(exc)
        print(
            f"  {record.name}: matched-EFT RGE failed: {exc}"
        )
        return False

    summary.update(rge_summary)

    print(
        f"  {record.name}: matched-EFT RGE=Success"
        f" -> {record.output_dir / rge_summary['C5BetaFile']}"
    )

    return summary.get("RGEStatus") == "Success"


def run_flavor_rge_stage(
    record: RunRecord,
    debug_reports: bool = False,
) -> bool:
    """Run the symbolic full-flavor Weinberg RGE for one model."""

    summary = record.summary
    c5_path = matched_c5_path(record)

    if c5_path is None:
        summary["FlavorRGEStatus"] = "NotRun"
        return False

    print(
        f"  {record.name}: starting symbolic full-flavor RGE stage...",
        flush=True,
    )

    try:
        flavor_summary = run_flavor_matched_rge(
            c5_path=c5_path,
            output_dir=record.output_dir,
            debug_outputs=debug_reports,
        )
    except Exception as exc:
        summary["FlavorRGEStatus"] = "Failed"
        summary["FlavorRGEError"] = str(exc)
        print(
            f"  {record.name}: full-flavor RGE failed: {exc}"
        )
        return False

    summary.update(flavor_summary)

    print(
        f"  {record.name}: full-flavor RGE=Success"
        f" -> "
        f"{record.output_dir / flavor_summary['C5FlavorBetaMatrixFile']}"
    )

    return summary.get("FlavorRGEStatus") == "Success"


def run_symbolic_neutrino_mass_stage(record: RunRecord) -> bool:
    """Construct the symbolic neutrino-mass matrix from the matched C5."""

    summary = record.summary
    c5_path = matched_c5_path(record)

    if c5_path is None:
        summary["NeutrinoMassStatus"] = "NotRun"
        return False

    print(
        f"  {record.name}: starting symbolic neutrino mass stage...",
        flush=True,
    )

    try:
        mass_summary = run_neutrino_mass_stage(
            c5_path=c5_path,
            output_dir=record.output_dir,
        )
    except Exception as exc:
        summary["NeutrinoMassStatus"] = "Failed"
        summary["NeutrinoMassError"] = str(exc)
        print(
            f"  {record.name}: neutrino mass stage failed: {exc}"
        )
        return False

    summary.update(mass_summary)

    print(
        f"  {record.name}: neutrino mass=Success"
        f" -> "
        f"{record.output_dir / mass_summary['NeutrinoMassMatrixFile']}"
    )

    return summary.get("NeutrinoMassStatus") == "Success"


def run_numerical_rge_stage(
    record: RunRecord,
    numerical_config: Path,
) -> bool:
    """Run numerical EFT evolution and obtain the resulting neutrino observables."""

    summary = record.summary
    c5_path = matched_c5_path(record)

    if c5_path is None:
        summary["NumericalRGEStatus"] = "NotRun"
        return False

    print(f"  {record.name}: starting numerical RGE stage...", flush=True)

    try:
        numerical_summary = run_numerical_pipeline_stage(
            c5_path=c5_path,
            output_dir=record.output_dir,
            config_path=numerical_config,
        )
    except Exception as exc:
        summary["NumericalRGEStatus"] = "Failed"
        summary["NumericalRGEError"] = str(exc)
        print(
            f"  {record.name}: numerical RGE failed: {exc}"
        )
        return False

    summary.update(numerical_summary)

    print(
        f"  {record.name}: numerical RGE calculation finished.",
        flush=True,
    )

    mass_matrix_path = (
        record.output_dir
        / numerical_summary["NeutrinoMassMatrixLowScaleFile"]
    )

    numerical_payload = json.loads(
        numerical_config.read_text(encoding="utf-8")
    )
    ordering = numerical_payload.get("ordering", "NO")

    print(f"  {record.name}: starting neutrino observables...", flush=True)

    try:
        observable_summary = run_neutrino_observables_stage(
            mass_matrix_path=mass_matrix_path,
            output_dir=record.output_dir,
            ordering=ordering,
        )
    except Exception as exc:
        summary["NeutrinoObservableStatus"] = "Failed"
        summary["NeutrinoObservableError"] = str(exc)
        print(
            f"  {record.name}: neutrino observables failed: {exc}"
        )
        return False

    summary.update(observable_summary)

    print(
        f"  {record.name}: neutrino observables=Success"
        f" -> "
        f"{record.output_dir / observable_summary['NeutrinoObservablesFile']}"
    )

    print(
        f"  {record.name}: numerical RGE=Success"
        f" -> "
        f"{record.output_dir / numerical_summary['NeutrinoMassMatrixLowScaleFile']}"
    )

    return (
        summary.get("NumericalRGEStatus") == "Success"
        and summary.get("NeutrinoObservableStatus") == "Success"
    )


def finish_runs(
    records: list[RunRecord],
    debug_reports: bool = False,
    physics_failed: bool = False,
    *,
    study_output_dir: Path | None = None,
    study_report_dir: Path | None = None,
) -> int:
    """Print final summaries and generate reports after the physics pipeline."""

    study_output_dir = study_output_dir or OUTPUT_DIR
    study_report_dir = study_report_dir or REPORT_OUTPUT_DIR

    status = print_summary(records, study_output_dir)

    if physics_failed:
        status = 1

    # ------------------------------------------------------------------
    # Human-readable comparison reports
    # ------------------------------------------------------------------
    #
    # Keep the report surface deliberately small.  Calculation data remains
    # under output/, while Reports/output/ contains only the comparison
    # Lagrangian and RGE reports used to inspect the physics across models.

    # Lagrangian reports:
    #   rows    = model configurations
    #   columns = field configurations
    uv_lagrangian_tex = write_bsm_uv_field_table(records, report_root=study_report_dir)
    compile_latex_document(uv_lagrangian_tex)

    # This writer emits and compiles every real sequential EFT stage.
    write_bsm_matched_field_table(records, report_root=study_report_dir)

    # Standalone study-level Weinberg-coefficient report, kept alongside the
    # UV/EFT Lagrangian reports.
    c5_tex = write_c5_coefficient_report(
        records,
        report_root=study_report_dir,
    )
    compile_latex_document(c5_tex)

    # RGE report:
    #   one table per running coupling
    #   rows = model configurations
    #
    # RGEComparison will be refined next so that its columns are the
    # individual full RGE terms rather than one combined beta-function cell.
    write_and_compile_rge_comparison(records, report_root=study_report_dir)
    write_and_compile_eft1_rge_comparison(records, report_root=study_report_dir)
    write_and_compile_final_eft_rge_comparison(records, report_root=study_report_dir)

    # Stage-aware analytic group-factor reports.
    write_and_compile_stage_group_factor_reports(
        records,
        report_root=study_report_dir,
    )

    final_stage_label = (
        records[0].eft_stages[-1].label
        if records and records[0].eft_stages
        else "EFT"
    )
    print(
        "\nStage-aware reports:"
        "\n  Lagrangian/UV"
        "\n  Lagrangian/EFT_1_after_F"
        f"\n  Lagrangian/{final_stage_label}"
        "\n  Lagrangian/C5"
        "\n  RGE/UV"
        "\n  RGE/EFT_1_after_F"
        f"\n  RGE/{final_stage_label}"
        "\n  GroupFactors/GF_UV"
        "\n  GroupFactors/GF_EFT_1_after_F"
    )

    return status


# Function:
# 1. Get all arguments
# 2. Start process
# 3. Get reports
def main() -> int:
    # Get Arguments
    parser = argparse.ArgumentParser(
        description=(
            "Run T3 matching for known benchmark models "
            "or arbitrary valid SU(2) irreps."
        )
    )

    # Mutually exclusive command arguments
    mode = parser.add_mutually_exclusive_group()

    # If we use --smoke we use the 5 T3 models we know
    mode.add_argument(
        "--smoke",
        action="store_true",
        help="five T3 regression models",
    )

    # More running (not used)
    mode.add_argument(
        "--extended",
        action="store_true",
        help="seven historical benchmark points",
    )

    mode.add_argument(
        "--hypercharge-comparison",
        action="store_true",
        help=(
            "scan alpha=-2,-1,0,1,2 for every T3-A...E SU(2) assignment "
            "to isolate hypercharge dependence"
        ),
    )

    mode.add_argument(
        "--dimension-comparison",
        action="store_true",
        help=(
            "compare T3-A...E at fixed alpha=0 to isolate SU(2) "
            "representation dependence"
        ),
    )

    mode.add_argument(
        "--full",
        action="store_true",
        help=(
            "run smoke, hypercharge comparison, and dimension comparison "
            "sequentially under output/full and Reports/output/full"
        ),
    )

    # Dimension input to use a specific diagram
    mode.add_argument(
        "--dims",
        nargs="+",
        type=int,
        metavar="D",
        help=(
            "three numbers DS1 DS2 DF give the ordinary T3 model; "
            "two numbers DS DF give one physical shared scalar, e.g. "
            "--dims 2 1 for the scotogenic singlet-fermion model"
        ),
    )


    parser.add_argument(
        "--study",
        type=str,
        default=None,
        help=(
            "output/report study folder. Defaults to smoke, extended, "
            "interesting, or single depending on the selected run mode; "
            "examples: hypercharge, dimensions"
        ),
    )

    # We can attach a file for UV scale initial values for our couplings etc to RGE down to EFT
    parser.add_argument(
        "--numerical",
        type=Path,
        default=None,
        help=(
            "optional JSON parameter point for numerical "
            "matched-EFT running E.g python pipeline.py --numerical examples/t3_numerical_example.json"
        ),
    )

    # Adds the alpha hypercharge which is added onto the dims
    parser.add_argument(
        "--alpha",
        type=int,
        default=None,
        help=(
            "hypercharge parameter for three-number --dims mode (default 0). "
            "Two-number shared-scalar mode fixes alpha=-1."
        ),
    )

    # Ordered heavy-particle thresholds. Repeat --threshold to create
    # successive EFT levels; fields given in the same occurrence are
    # integrated out together.
    #
    # Examples:
    #   --threshold F --threshold S1 S2
    #   --threshold S1 --threshold F --threshold S2
    #   --threshold F S1 S2
    parser.add_argument(
        "--threshold",
        action="append",
        nargs="+",
        metavar="FIELD",
        default=None,
        help=(
            "ordered threshold group; ordinary mode uses F,S1,S2 and shared-scalar mode uses F,S. "
            "Repeat the option for successive thresholds. "
            "Fields in one group are integrated out together. "
            "If omitted, F S1 S2 are integrated out together."
        ),
    )

    parser.add_argument(
        "--threshold-scale",
        action="append",
        metavar="SCALE",
        default=None,
        help=(
            "matching scale for each --threshold occurrence, in the same "
            "order. Values may be symbolic (MF, MS1, ...) or numeric. "
            "For F -> (S1,S2), defaults are MF and MS."
        ),
    )

    # For logging just in case
    parser.add_argument(
        "--debug-reports",
        action="store_true",
        help=(
            "also keep raw Wolfram logs and generate the full "
            "UV/EFT expression reports"
        ),
    )

    args = parser.parse_args()

    # --full is a pure orchestration mode. It launches the ordinary study
    # modes as child processes, each with its own nested output/report root.
    if args.full:
        return run_full_study(args)

    shared_scalar_mode = False
    if args.dims is not None:
        if len(args.dims) not in (2, 3):
            parser.error("--dims requires either DS DF or DS1 DS2 DF.")
        shared_scalar_mode = len(args.dims) == 2

    if shared_scalar_mode:
        if args.alpha is not None and args.alpha != -1:
            parser.error("Two-number shared-scalar mode requires alpha=-1.")
        args.alpha = -1
    elif args.alpha is None:
        args.alpha = 0

    try:
        threshold_plan = validate_threshold_plan(
            args.threshold, shared_scalar=shared_scalar_mode
        )
    except ValueError as exc:
        parser.error(str(exc))

    if args.threshold_scale is not None:
        if len(args.threshold_scale) != len(threshold_plan):
            parser.error(
                "--threshold-scale must be supplied once for each "
                "--threshold group."
            )
        threshold_scales = list(args.threshold_scale)
    else:
        threshold_scales = [
            _default_threshold_scale(group)
            for group in threshold_plan
        ]

    if args.study:
        # Keep '/' so --full can intentionally create nested study roots such
        # as full/hypercharge. Spaces are still normalised for CLI convenience.
        study_name = args.study.strip().replace(" ", "_")
    elif args.smoke:
        study_name = "smoke"
    elif args.hypercharge_comparison:
        study_name = "hypercharge"
    elif args.dimension_comparison:
        study_name = "dimensions"
    elif args.extended:
        study_name = "extended"
    elif args.dims:
        study_name = "single"
    else:
        study_name = "interesting"

    study_output_dir = OUTPUT_DIR / study_name
    study_report_dir = REPORT_OUTPUT_DIR / study_name

    print(
        "Threshold plan: "
        + threshold_plan_label(threshold_plan)
    )

    study_output_dir.mkdir(parents=True, exist_ok=True)
    study_report_dir.mkdir(parents=True, exist_ok=True)

    print(f"Study: {study_name}")
    print(f"Raw output: {study_output_dir}")
    print(f"Reports: {study_report_dir}")
    # ----------------------------------------------------------------------
    # Lagrangian and Weinberg Coefficient pipeline
    # ----------------------------------------------------------------------

    # If the dimensions happen to match T3-A ... T3-E, validate_dimensions()
    # automatically recognises and labels the model appropriately.
    if args.dims:
        try:
            if shared_scalar_mode:
                d_s, d_f = args.dims
                record = validate_shared_dimensions(
                    d_s, d_f,
                    debug_reports=args.debug_reports,
                    export_rge_tensors=False,
                    threshold_plan=threshold_plan,
                    output_root=study_output_dir,
                )
            else:
                d_s1, d_s2, d_f = args.dims
                record = validate_dimensions(
                    d_s1, d_s2, d_f, args.alpha, args.debug_reports, False,
                    threshold_plan, study_output_dir,
                )
        except ValueError as exc:
            parser.error(str(exc))

        records = [record]

    else:
        # The benchmark lists still use the familiar A-E notation because
        # it is convenient for regression testing and comparison with the paper.
        if args.smoke:
            mode_name = "smoke"
            points = SMOKE

        elif args.hypercharge_comparison:
            mode_name = "hypercharge comparison"
            points = HYPERCHARGE_COMPARISON

        elif args.dimension_comparison:
            mode_name = "dimension comparison"
            points = DIMENSION_COMPARISON

        elif args.extended:
            mode_name = "extended"
            points = EXTENDED

        else:
            mode_name = "interesting"
            points = INTERESTING

        print(
            f"T3 scan mode: {mode_name}; "
            f"{len(points)} model(s)."
        )

        # This line only exist if we use the special interesting, extended and smoke options where we dont input any dimensions
        # Known A-E models are converted to dimensions first and then sent
        # through exactly the same validate_dimensions() path as generalised models.
        records = [
            obtain_class_dimensions(
                model_class,
                alpha,
                args.debug_reports,
                False,
                threshold_plan,
                study_output_dir,
            )
            for model_class, alpha in points
        ]

    # Attach Python-side stage records using the same threshold plan.  The
    # Wolfram summary now contains the actual sequentially matched EFT-stage
    # Lagrangians; do not overwrite those physics results here.
    for record in records:
        record.eft_stages = build_eft_stage_records(
            threshold_plan,
            record.output_dir,
            shared_scalar=record.shared_scalar,
        )

        record.summary.setdefault(
            "ThresholdPlan",
            threshold_plan_to_json(threshold_plan),
        )
        record.summary.setdefault(
            "ThresholdPlanLabel",
            threshold_plan_label(threshold_plan),
        )
        record.summary.setdefault(
            "ThresholdScales",
            threshold_scales,
        )

    # Organise the matched coefficient before starting the RGE pipeline, so all
    # later stages read C5 from the same final machine-readable location.
    for record in records:
        organise_c5_input(record)

    # ----------------------------------------------------------------------
    # RGE pipeline
    # ----------------------------------------------------------------------

    physics_failed = False

    # For every model, only run the UV RGE stage if the model Lagrangian was built successfully. Also remember whether any model failed.
    # We want UV RGE for comparison with EFT RGE
    for record in records:
        if record.summary.get("BuildStatus") == "Success":
            physics_failed |= not run_uv_rgbeta_stage(record)
        else:
            record.summary["UVRGEStatus"] = "NotRun"

    # If F is the first threshold, generate the dimension-four running in the
    # genuine intermediate EFT SM + S1 + S2.  This is intentionally separate
    # from the higher-dimensional Wilson-coefficient running.
    for record in records:
        if record.summary.get("BuildStatus") == "Success":
            eft1_renormalisable_ok = run_eft1_rgbeta_stage(record)
            physics_failed |= not eft1_renormalisable_ok

            # The higher-dimensional EFT1 RGE needs the RGBeta metadata written
            # by the preceding stage, so only run it when that stage succeeded.
            if eft1_renormalisable_ok:
                eft1_wilson_ok = run_eft1_wilson_rge_stage(record)
                physics_failed |= not eft1_wilson_ok

                if (
                    eft1_wilson_ok
                    and _is_f_then_scalar_threshold_plan(
                        threshold_plan,
                        record,
                    )
                ):
                    # Keep the one-generation component transport as a compact
                    # regression/diagnostic, then run the authoritative
                    # full-flavor stop/run/resume threshold bridge.
                    transport_ok = run_eft1_wilson_transport_stage(
                        record,
                        mu_high=threshold_scales[0],
                        mu_low=threshold_scales[1],
                    )
                    physics_failed |= not transport_ok

                    if transport_ok:
                        physics_failed |= not run_eft1_full_flavor_threshold_bridge(
                            record,
                            mu_high=threshold_scales[0],
                            mu_low=threshold_scales[1],
                            debug_reports=args.debug_reports,
                        )

    # After matching has generated the Weinberg coefficient C5, the EFT RGE pipeline starts here. 
    # Each later stage runs only if the previous stage for that model succeeded.
    for record in records:
        summary = record.summary

        # Check if any part works for a model, if a model fails we skip
        if summary.get("BuildStatus") != "Success":
            continue

        if summary.get("MatchingStatus") != "Success":
            continue

        if summary.get("WeinbergExtractionStatus") != "Success":
            continue

        # Hierarchical F -> (S1,S2) runs promote
        # data/final_weinberg_coefficient.json to WeinbergCoefficientFile in
        # run_eft1_full_flavor_threshold_bridge.  Legacy/common-threshold runs
        # continue to use their organised scalar C5 text file.

        # Evaluate the matched Weinberg coefficient with the general
        # one-generation SMEFT RGE machinery and verify the known SMEFT result.
        if not run_matched_eft_rge_stage(
            record,
            args.debug_reports,
        ):
            physics_failed = True
            continue

        # Get the matched one-generation C5 into the full 3x3 lepton-flavor
        # matrix and calculate its symbolic SMEFT beta matrix.
        if not run_flavor_rge_stage(
            record,
            args.debug_reports,
        ):
            physics_failed = True
            continue

        # Convert the matched symbolic flavor coefficient into the symbolic
        # Majorana neutrino-mass matrix using m_nu = -v^2 C5.
        if not run_symbolic_neutrino_mass_stage(record):
            physics_failed = True
            continue

        # If a numerical parameter point and argument is given, evaluate C5 at the
        # matching scale and numerically integrate the coupled one-loop SMEFT RGEs
        # down to the requested low scale. Then calculate neutrino observables.
        if args.numerical is not None:
            if not run_numerical_rge_stage(
                record,
                args.numerical,
            ):
                physics_failed = True
                continue

    # Reports and final summaries happen only after the full physics pipeline.
    return finish_runs(
        records,
        args.debug_reports,
        physics_failed,
        study_output_dir=study_output_dir,
        study_report_dir=study_report_dir,
    )


if __name__ == "__main__":
    raise SystemExit(main())
