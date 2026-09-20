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
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

from RGE.matching.MatchedWeinbergRGE import run_matched_weinberg_rge
from RGE.running.weinberg.FlavorMatchedWeinbergStage import (
    run_flavor_matched_weinberg_rge,
    run_symbolic_neutrino_mass_stage as _run_symbolic_neutrino_mass_stage,
)
from RGE.running.weinberg.NumericalWeinbergStage import (
    run_numerical_weinberg_stage as _run_numerical_weinberg_stage,
)
from RGE.phenomenology.NeutrinoObservables import run_neutrino_observables_stage
from RGE.running.rgbeta.RGBetaT3Running import (
    run_rgbeta_t3,
    run_rgbeta_t3_eft1,
)
from RGE.running.eft1.EFT1WilsonRGE import run_eft1_wilson_rge
from RGE.running.eft1.EFT1WilsonFlow import run_eft1_wilson_transport
from RGE.running.eft1.EFT1WilsonFlow import run_flavor_seed_export
from RGE.running.eft1.EFT1DirectWeinberg import (
    build_direct_weinberg_flavor_transport,
    export_direct_weinberg_matchete,
)
from RGE.running.eft1.EFT1ThresholdResume import rerun_threshold2_with_running
from RGE.running.weinberg.FinalWeinbergCoefficient import (
    build_final_weinberg_coefficient,
    normalize_pole_rge_consistency,
)

from common.RunRecords import RunRecord
from common.Thresholds import (
    build_eft_stage_records,
    threshold_plan_label,
    threshold_plan_to_json,
    validate_threshold_plan,
)
from common.T3Model import SMOKE
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


def _eft1_full_flavor_bridge_paths(record: RunRecord) -> dict[str, Path]:
    """Return all files used by the fixed-one-loop EFT1 Weinberg bridge."""
    data_dir = record.output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    return {
        "data_dir": data_dir,
        "wilson_seed": data_dir / "eft1_after_F_wilson_seed.json",
        "component_rge": data_dir / "eft1_wilson_rge.json",
        "flavor_seed": data_dir / "eft1_wilson_flavor_seed.json",
        "flavor_transport": data_dir / "eft1_wilson_flavor_at_S_threshold.json",
        "insertion": data_dir / "eft1_wilson_flavor_running_insertion.wl",
        "pole_rge": data_dir / "c5_pole_rge_consistency.json",
        "threshold_c5": data_dir / "c5_threshold_with_eft1_running.txt",
        "direct_c5": data_dir / "c5_direct_eft1_running.txt",
        "final_c5": data_dir / "final_weinberg_coefficient.json",
    }


def _fail_eft1_full_flavor_bridge(record: RunRecord, message: str) -> bool:
    """Record and print a failed fixed-one-loop EFT1 Weinberg bridge."""
    record.summary["EFT1FullFlavorBridgeStatus"] = "Failed"
    record.summary["EFT1FullFlavorBridgeError"] = message
    print(f"  {record.name}: full-flavor EFT1 bridge failed: {message}")
    return False


def _validate_eft1_full_flavor_bridge_inputs(
    record: RunRecord,
    paths: dict[str, Path],
) -> bool:
    """Check that the direct-Weinberg bridge inputs have been produced."""
    required = (paths["wilson_seed"], paths["component_rge"])
    missing = [path for path in required if not path.is_file()]
    if not missing:
        return True

    message = (
        "Missing direct-Weinberg bridge input(s): "
        + ", ".join(
            path.relative_to(record.output_dir).as_posix()
            for path in missing
        )
    )
    return _fail_eft1_full_flavor_bridge(record, message)


def _run_direct_weinberg_bridge_stages(
    record: RunRecord,
    paths: dict[str, Path],
    *,
    mu_high: str,
    mu_low: str,
    debug_reports: bool,
) -> tuple[dict, dict, dict, dict]:
    """Run flavor seeding, direct Weinberg transport, and threshold resume."""
    flavor_seed = run_flavor_seed_export(
        wilson_seed_path=paths["wilson_seed"],
        output_path=paths["flavor_seed"],
    )
    if flavor_seed.get("status") != "Success":
        raise RuntimeError("Full-flavor Wilson seed regression failed.")

    flavor_transport = build_direct_weinberg_flavor_transport(
        flavor_seed_path=paths["flavor_seed"],
        component_rge_path=paths["component_rge"],
        mu_high=mu_high,
        mu_low=mu_low,
        output_path=paths["flavor_transport"],
    )
    if flavor_transport.get("status") != "Success":
        raise RuntimeError(
            "Direct full-flavor Weinberg transport regression failed."
        )

    insertion = export_direct_weinberg_matchete(
        transport_path=paths["flavor_transport"],
        output_path=paths["insertion"],
    )
    if insertion.get("status") != "Success":
        raise RuntimeError("Direct Weinberg Matchete insertion export failed.")

    resume = rerun_threshold2_with_running(
        output_dir=record.output_dir,
        running_insertion=paths["insertion"],
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

    if paths["pole_rge"].is_file():
        pole_rge = normalize_pole_rge_consistency(paths["pole_rge"])
        if not pole_rge.get("ConsistencyValidated", False):
            raise RuntimeError(
                "C5 pole/RGE consistency failed strict regression: "
                "direct log coefficient must equal 2 x hard pole residue; "
                f"reported ratio="
                f"{pole_rge.get('DirectLogToPoleRatioOneGenerationInputForm')}."
            )

    return flavor_seed, flavor_transport, insertion, resume


def _build_authoritative_final_weinberg(
    record: RunRecord,
    paths: dict[str, Path],
    *,
    flavor_transport: dict,
) -> dict:
    """Build and validate the authoritative hierarchical Weinberg coefficient."""
    if not paths["threshold_c5"].is_file():
        raise RuntimeError(
            "Authoritative resumed threshold C5 was not exported: "
            f"{paths['threshold_c5']}"
        )

    try:
        final_c5 = build_final_weinberg_coefficient(
            threshold_c5_path=paths["threshold_c5"],
            flavor_transport_path=paths["flavor_transport"],
            output_path=paths["final_c5"],
        )
    except Exception as exc:
        raise RuntimeError(
            f"Final Weinberg coefficient bookkeeping failed: {exc}"
        ) from exc

    if final_c5.get("status") != "Success":
        raise RuntimeError(
            "Final Weinberg coefficient bookkeeping returned non-success."
        )

    combined = final_c5.get("combined", {})
    if not combined.get("ready_for_physical_majorana_numerics", False):
        raise RuntimeError(
            "Final Weinberg coefficient is not ready for physical Majorana "
            f"downstream use: {combined.get('physical_majorana_reason')!r}"
        )

    return final_c5


def _record_eft1_full_flavor_bridge_success(
    record: RunRecord,
    paths: dict[str, Path],
    *,
    flavor_seed: dict,
    flavor_transport: dict,
    insertion: dict,
    resume: dict,
    final_c5: dict,
) -> None:
    """Store successful EFT1 bridge provenance in the run summary."""
    summary = record.summary
    combined = final_c5.get("combined", {})

    # Every EFT-side stage after this point consumes the authoritative
    # hierarchical full-flavor coefficient, never the preliminary C5 exported
    # before the resumed threshold calculation.
    summary["WeinbergCoefficientFile"] = (
        paths["final_c5"].relative_to(record.output_dir).as_posix()
    )

    summary.update(
        {
            "EFT1FlavorSeedStatus": flavor_seed.get("status"),
            "EFT1FlavorSeedFile": paths["flavor_seed"].relative_to(
                record.output_dir
            ).as_posix(),
            "EFT1FlavorSeedOneGenerationCheck": flavor_seed.get(
                "one_generation_reduction_matches"
            ),
            "EFT1FlavorRGEStatus": "DiagnosticNotRequiredForOneLoopC5",
            "EFT1FlavorRGEOneGenerationCheck": None,
            "EFT1FlavorTransportStatus": flavor_transport.get("status"),
            "EFT1FlavorTransportFile": paths["flavor_transport"].relative_to(
                record.output_dir
            ).as_posix(),
            "EFT1FlavorTransportEqualScaleCheck": flavor_transport.get(
                "equal_scale_running_vanishes"
            ),
            "EFT1FlavorRunningInsertionStatus": insertion.get("status"),
            "EFT1FlavorRunningInsertionFile": paths["insertion"].relative_to(
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
            "AuthoritativeThresholdC5File": paths["threshold_c5"].relative_to(
                record.output_dir
            ).as_posix(),
            "AuthoritativeDirectRunningC5File": (
                paths["direct_c5"].relative_to(record.output_dir).as_posix()
                if paths["direct_c5"].is_file()
                else ""
            ),
            "FinalWeinbergCoefficientFile": paths["final_c5"].relative_to(
                record.output_dir
            ).as_posix(),
            "AuthoritativeFinalC5Status": "ReadyForDownstream",
            "AuthoritativeFinalC5FullFlavorReady": combined.get(
                "ready_for_full_flavor_numerics"
            ),
            "AuthoritativeFinalC5PhysicalMajoranaReady": combined.get(
                "ready_for_physical_majorana_numerics"
            ),
            "AuthoritativeDownstreamC5File": paths["final_c5"].relative_to(
                record.output_dir
            ).as_posix(),
        }
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
    if not _is_f_then_s1_s2_plan(record):
        record.summary["EFT1FullFlavorBridgeStatus"] = "NotApplicable"
        return True

    paths = _eft1_full_flavor_bridge_paths(record)
    if not _validate_eft1_full_flavor_bridge_inputs(record, paths):
        return False

    print(
        f"  {record.name}: starting fixed-one-loop direct Weinberg bridge "
        f"{mu_high} -> {mu_low}...",
        flush=True,
    )

    try:
        flavor_seed, flavor_transport, insertion, resume = (
            _run_direct_weinberg_bridge_stages(
                record,
                paths,
                mu_high=mu_high,
                mu_low=mu_low,
                debug_reports=debug_reports,
            )
        )
        final_c5 = _build_authoritative_final_weinberg(
            record,
            paths,
            flavor_transport=flavor_transport,
        )
    except Exception as exc:
        return _fail_eft1_full_flavor_bridge(record, str(exc))

    _record_eft1_full_flavor_bridge_success(
        record,
        paths,
        flavor_seed=flavor_seed,
        flavor_transport=flavor_transport,
        insertion=insertion,
        resume=resume,
        final_c5=final_c5,
    )

    print(
        f"  {record.name}: direct one-loop Weinberg bridge=Success "
        f"(heavy self-running excluded at O(hbar); "
        f"equal-scale check="
        f"{record.summary['EFT1FlavorTransportEqualScaleCheck']})",
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


def run_matched_weinberg_rge_stage(
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
        rge_summary = run_matched_weinberg_rge(
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


def run_symbolic_flavor_weinberg_stage(
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
        flavor_summary = run_flavor_matched_weinberg_rge(
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
        mass_summary = _run_symbolic_neutrino_mass_stage(
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


def run_numerical_weinberg_stage(
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
        numerical_summary = _run_numerical_weinberg_stage(
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



def _build_argument_parser() -> argparse.ArgumentParser:
    """Build the command-line interface for the T3 pipeline."""

    parser = argparse.ArgumentParser(
        description=(
            "Run T3 matching for known benchmark models "
            "or arbitrary valid SU(2) irreps."
        )
    )
    mode = parser.add_mutually_exclusive_group()

    mode.add_argument(
        "--smoke",
        action="store_true",
        help="five T3 regression models",
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
    parser.add_argument(
        "--numerical",
        type=Path,
        default=None,
        help=(
            "optional JSON parameter point for numerical matched-EFT running; "
            "e.g. python pipeline.py --numerical examples/t3_numerical_example.json"
        ),
    )
    parser.add_argument(
        "--alpha",
        type=int,
        default=None,
        help=(
            "hypercharge parameter for three-number --dims mode (default 0). "
            "Two-number shared-scalar mode fixes alpha=-1."
        ),
    )
    parser.add_argument(
        "--threshold",
        action="append",
        nargs="+",
        metavar="FIELD",
        default=None,
        help=(
            "ordered threshold group; ordinary mode uses F,S1,S2 and "
            "shared-scalar mode uses F,S. Repeat the option for successive "
            "thresholds. Fields in one group are integrated out together. "
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
    parser.add_argument(
        "--debug-reports",
        action="store_true",
        help=(
            "also keep raw Wolfram logs and generate the full "
            "UV/EFT expression reports"
        ),
    )
    return parser


def _resolve_model_mode(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
) -> bool:
    """Validate dimension/alpha options and return shared-scalar mode."""

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

    return shared_scalar_mode


def _resolve_threshold_configuration(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
    *,
    shared_scalar_mode: bool,
) -> tuple[list[list[str]], list[str]]:
    """Validate threshold ordering and choose one scale per threshold group."""

    try:
        threshold_plan = validate_threshold_plan(
            args.threshold,
            shared_scalar=shared_scalar_mode,
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

    return threshold_plan, threshold_scales


def _study_name(args: argparse.Namespace) -> str:
    """Return the output/report study directory name for one pipeline run."""

    if args.study:
        # Keep '/' so --full can intentionally create nested study roots such
        # as full/hypercharge. Spaces are normalised for CLI convenience.
        return args.study.strip().replace(" ", "_")
    if args.smoke:
        return "smoke"
    if args.hypercharge_comparison:
        return "hypercharge"
    if args.dimension_comparison:
        return "dimensions"
    if args.dims:
        return "single"
    return


def _scan_definition(
    args: argparse.Namespace,
) -> tuple[str, tuple[tuple[str, int], ...]]:
    """Return the label and benchmark points for a non-``--dims`` scan."""

    if args.smoke:
        return "smoke", SMOKE
    if args.hypercharge_comparison:
        return "hypercharge comparison", HYPERCHARGE_COMPARISON
    if args.dimension_comparison:
        return "dimension comparison", DIMENSION_COMPARISON
    return


def _build_run_records(
    parser: argparse.ArgumentParser,
    args: argparse.Namespace,
    *,
    shared_scalar_mode: bool,
    threshold_plan: list[list[str]],
    study_output_dir: Path,
) -> list[RunRecord]:
    """Build the requested model(s) through the common Lagrangian path."""

    if args.dims:
        try:
            if shared_scalar_mode:
                d_s, d_f = args.dims
                record = validate_shared_dimensions(
                    d_s,
                    d_f,
                    debug_reports=args.debug_reports,
                    export_rge_tensors=False,
                    threshold_plan=threshold_plan,
                    output_root=study_output_dir,
                )
            else:
                d_s1, d_s2, d_f = args.dims
                record = validate_dimensions(
                    d_s1,
                    d_s2,
                    d_f,
                    args.alpha,
                    args.debug_reports,
                    False,
                    threshold_plan,
                    study_output_dir,
                )
        except ValueError as exc:
            parser.error(str(exc))
        return [record]

    mode_name, points = _scan_definition(args)
    print(f"T3 scan mode: {mode_name}; {len(points)} model(s).")
    return [
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


def _attach_threshold_metadata(
    records: list[RunRecord],
    *,
    threshold_plan: list[list[str]],
    threshold_scales: list[str],
) -> None:
    """Attach Python-side stage records and threshold metadata to each run."""

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
        record.summary.setdefault("ThresholdScales", threshold_scales)


def _run_intermediate_eft_stages(
    records: list[RunRecord],
    *,
    threshold_plan: list[list[str]],
    threshold_scales: list[str],
    debug_reports: bool,
) -> bool:
    """Run UV and, when present, intermediate-EFT RGBeta/Wilson stages."""

    physics_failed = False

    for record in records:
        if record.summary.get("BuildStatus") == "Success":
            physics_failed |= not run_uv_rgbeta_stage(record)
        else:
            record.summary["UVRGEStatus"] = "NotRun"

    for record in records:
        if record.summary.get("BuildStatus") != "Success":
            continue

        eft1_renormalisable_ok = run_eft1_rgbeta_stage(record)
        physics_failed |= not eft1_renormalisable_ok
        if not eft1_renormalisable_ok:
            continue

        eft1_wilson_ok = run_eft1_wilson_rge_stage(record)
        physics_failed |= not eft1_wilson_ok
        if not eft1_wilson_ok:
            continue

        if not _is_f_then_scalar_threshold_plan(threshold_plan, record):
            continue

        # Keep the one-generation component transport as a compact diagnostic,
        # then run the authoritative full-flavor stop/run/resume bridge.
        transport_ok = run_eft1_wilson_transport_stage(
            record,
            mu_high=threshold_scales[0],
            mu_low=threshold_scales[1],
        )
        physics_failed |= not transport_ok
        if not transport_ok:
            continue

        bridge_ok = run_eft1_full_flavor_threshold_bridge(
            record,
            mu_high=threshold_scales[0],
            mu_low=threshold_scales[1],
            debug_reports=debug_reports,
        )
        physics_failed |= not bridge_ok

    return physics_failed


def _ready_for_final_weinberg_stages(record: RunRecord) -> bool:
    """Return whether matching produced the inputs needed by final C5 stages."""

    summary = record.summary
    return (
        summary.get("BuildStatus") == "Success"
        and summary.get("MatchingStatus") == "Success"
        and summary.get("WeinbergExtractionStatus") == "Success"
    )


def _run_final_weinberg_stages(
    records: list[RunRecord],
    *,
    debug_reports: bool,
    numerical_input: Path | None,
) -> bool:
    """Run matched, flavor, neutrino-mass, and optional numerical stages."""

    physics_failed = False

    for record in records:
        if not _ready_for_final_weinberg_stages(record):
            continue

        if not run_matched_weinberg_rge_stage(record, debug_reports):
            physics_failed = True
            continue

        if not run_symbolic_flavor_weinberg_stage(record, debug_reports):
            physics_failed = True
            continue

        if not run_symbolic_neutrino_mass_stage(record):
            physics_failed = True
            continue

        if numerical_input is not None:
            if not run_numerical_weinberg_stage(record, numerical_input):
                physics_failed = True

    return physics_failed


def main() -> int:
    """Parse one pipeline run, execute its stages, and generate reports."""

    parser = _build_argument_parser()
    args = parser.parse_args()

    if args.full:
        return run_full_study(args)

    shared_scalar_mode = _resolve_model_mode(parser, args)
    threshold_plan, threshold_scales = _resolve_threshold_configuration(
        parser,
        args,
        shared_scalar_mode=shared_scalar_mode,
    )

    study_name = _study_name(args)
    study_output_dir = OUTPUT_DIR / study_name
    study_report_dir = REPORT_OUTPUT_DIR / study_name
    study_output_dir.mkdir(parents=True, exist_ok=True)
    study_report_dir.mkdir(parents=True, exist_ok=True)

    print(f"Threshold plan: {threshold_plan_label(threshold_plan)}")
    print(f"Study: {study_name}")
    print(f"Raw output: {study_output_dir}")
    print(f"Reports: {study_report_dir}")

    records = _build_run_records(
        parser,
        args,
        shared_scalar_mode=shared_scalar_mode,
        threshold_plan=threshold_plan,
        study_output_dir=study_output_dir,
    )
    _attach_threshold_metadata(
        records,
        threshold_plan=threshold_plan,
        threshold_scales=threshold_scales,
    )

    for record in records:
        organise_c5_input(record)

    physics_failed = _run_intermediate_eft_stages(
        records,
        threshold_plan=threshold_plan,
        threshold_scales=threshold_scales,
        debug_reports=args.debug_reports,
    )
    physics_failed |= _run_final_weinberg_stages(
        records,
        debug_reports=args.debug_reports,
        numerical_input=args.numerical,
    )

    return finish_runs(
        records,
        args.debug_reports,
        physics_failed,
        study_output_dir=study_output_dir,
        study_report_dir=study_report_dir,
    )


if __name__ == "__main__":
    raise SystemExit(main())
