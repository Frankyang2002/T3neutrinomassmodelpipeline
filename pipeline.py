# Central backbone for the T3 neutrino-mass calculation.
#
# The ordinary 16-model symbolic comparison now uses the full general
# BSM interaction basis by default. Restricted matching/numerical production
# is retained explicitly; it is not falsely labelled general matching.
from __future__ import annotations

import argparse
from pathlib import Path

from studies.FullT3Study import run_full_study
from RGE.running.IntermediateEFTRunning import (
    mark_intermediate_running_not_applicable,
    run_intermediate_eft_interval,
)
from validation.IntermediateEFTValidation import run_intermediate_eft_validation
from common.EFT import EFTRunningInterval
from common.PipelinePlan import PipelinePlan
from common.PipelineCLI import (
    build_argument_parser as _build_argument_parser,
    resolve_model_mode as _resolve_model_mode,
    resolve_pipeline_plan as _resolve_pipeline_plan,
    study_name as _study_name,
)
from common.RunRecords import EFTStageRecord, RunRecord
from Numerical.orchestration.NumericalConfig import (
    freeze_resolved_benchmark_config as _freeze_resolved_benchmark_config,
    prepare_model_numerical_config as _prepare_model_numerical_config,
)
from Numerical.orchestration.PipelineNumericalResults import (
    is_pipeline_numerical_results_config,
    run_pipeline_numerical_results,
)
from physics.LowEnergyNeutrino import (
    build_symbolic_neutrino_mass,
    organise_matched_weinberg_coefficient,
    run_full_flavor_weinberg_rge,
    run_numerical_neutrino_observables,
    run_smeft_weinberg_rge,
)
from model.T3Study import select_study_models as _select_study_models
from Lagrangian.T3ModelMatching import (
    build_and_match_t3_models as _build_and_match_t3_models,
)
from RGE.running.UVRunning import run_uv_rge as run_uv_rgbeta_stage
from Reports.PipelineReports import finish_pipeline_run, print_pipeline_summary

PROJECT_ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = PROJECT_ROOT / "output"
REPORT_OUTPUT_DIR = PROJECT_ROOT / "Reports" / "output"


def _stage_for_running_interval(record: RunRecord, interval: EFTRunningInterval) -> EFTStageRecord:
    return record.stage_with_active_fields(*sorted(interval.active_heavy_fields))


def print_summary(records: list[RunRecord], aggregate_dir: Path | None = None) -> int:
    return print_pipeline_summary(records, aggregate_dir)


def finish_runs(
    records: list[RunRecord],
    debug_reports: bool = False,
    physics_failed: bool = False,
    *,
    study_output_dir: Path | None = None,
    study_report_dir: Path | None = None,
) -> int:
    return finish_pipeline_run(
        records,
        debug_reports=debug_reports,
        physics_failed=physics_failed,
        study_output_dir=study_output_dir,
        study_report_dir=study_report_dir,
    )


def _attach_threshold_metadata(records: list[RunRecord], *, pipeline_plan: PipelinePlan) -> None:
    metadata = pipeline_plan.summary_metadata()
    for record in records:
        if record.shared_scalar != pipeline_plan.shared_scalar:
            raise ValueError("RunRecord shared-scalar mode does not match PipelinePlan.")
        record.eft_stages = pipeline_plan.build_stage_records(record.output_dir)
        for key, value in metadata.items():
            record.summary.setdefault(key, value)


def _run_intermediate_eft_stages(
    records: list[RunRecord],
    *,
    pipeline_plan: PipelinePlan,
    debug_reports: bool,
) -> bool:
    physics_failed = False
    for record in records:
        if record.summary.get("BuildStatus") == "Success":
            physics_failed |= not run_uv_rgbeta_stage(record)
        else:
            record.summary["UVRGEStatus"] = "NotRun"

    for record in records:
        if record.summary.get("BuildStatus") != "Success":
            continue
        attempted_any = False
        production_complete = True
        unsupported_intervals: list[dict[str, object]] = []
        for interval in pipeline_plan.running_intervals:
            stage = _stage_for_running_interval(record, interval)
            outcome = run_intermediate_eft_interval(
                record, stage, interval, debug_reports=debug_reports,
            )
            attempted_any |= outcome.attempted
            if not outcome.attempted:
                production_complete = False
                physics_failed = True
                unsupported_intervals.append({
                    "ActiveHeavyFields": sorted(interval.active_heavy_fields),
                    "HighScale": interval.high_scale,
                    "LowScale": interval.low_scale,
                    "EnteredBy": list(interval.entered_by.fields_to_integrate),
                    "ExitedBy": list(interval.exited_by.fields_to_integrate),
                    "Status": outcome.status,
                })
                continue
            if not outcome.success:
                production_complete = False
                physics_failed = True
                continue
            validation = run_intermediate_eft_validation(record, stage, interval)
            if validation.attempted and not validation.success:
                physics_failed = True

        if not pipeline_plan.running_intervals:
            record.summary["IntermediateEFTProductionStatus"] = "NotRequired"
        elif production_complete:
            record.summary["IntermediateEFTProductionStatus"] = "Complete"
        elif unsupported_intervals:
            record.summary["IntermediateEFTProductionStatus"] = "IncompleteUnsupportedFieldContent"
        else:
            record.summary["IntermediateEFTProductionStatus"] = "Failed"
        record.summary["IntermediateEFTUnsupportedIntervals"] = unsupported_intervals
        if not attempted_any:
            mark_intermediate_running_not_applicable(record)
    return physics_failed


def _ready_for_low_energy_neutrino_stages(record: RunRecord) -> bool:
    summary = record.summary
    return (
        summary.get("BuildStatus") == "Success"
        and summary.get("MatchingStatus") == "Success"
        and summary.get("WeinbergExtractionStatus") == "Success"
        and summary.get("IntermediateEFTProductionStatus")
        not in {"IncompleteUnsupportedFieldContent", "Failed"}
    )


def _run_low_energy_neutrino_stages(
    records: list[RunRecord],
    *,
    debug_reports: bool,
    numerical_input: Path | None,
    reset_numerical_configs: bool,
) -> bool:
    physics_failed = False
    for record in records:
        if not _ready_for_low_energy_neutrino_stages(record):
            continue
        if not run_smeft_weinberg_rge(record, debug_reports):
            physics_failed = True
            continue
        if not run_full_flavor_weinberg_rge(record, debug_reports):
            physics_failed = True
            continue
        if not build_symbolic_neutrino_mass(record):
            physics_failed = True
            continue
        if numerical_input is not None:
            if is_pipeline_numerical_results_config(numerical_input):
                model_config = _prepare_model_numerical_config(
                    record, numerical_input, reset=reset_numerical_configs,
                )
                if not run_pipeline_numerical_results(record, model_config):
                    physics_failed = True
                else:
                    _freeze_resolved_benchmark_config(record, model_config)
            elif not run_numerical_neutrino_observables(record, numerical_input):
                physics_failed = True
    return physics_failed


def _use_general_symbolic(args: argparse.Namespace, parser: argparse.ArgumentParser) -> bool:
    """The 16-model symbolic scan defaults to the full interaction basis.

    No general-theory matching result is inferred from legacy restricted
    matching. The opt-out is explicit and preserves old production commands.
    """
    if args.general_symbolic and args.restricted_production:
        parser.error("--general-symbolic and --restricted-production cannot be combined")
    if args.reports_only and args.restricted_production:
        parser.error("--reports-only is not available with --restricted-production")
    requested = args.general_symbolic or args.reports_only or (
        args.hypercharge_comparison and args.numerical is None
        and not args.full and not args.restricted_production
    )
    if requested:
        if not args.hypercharge_comparison:
            parser.error("General symbolic mode requires --hypercharge-comparison (the 16-model catalogue)")
        if args.numerical is not None or args.full or args.dims:
            parser.error("General symbolic mode does not support numerical/full/dims matching")
        if args.force:
            parser.error("--force is not supported by the fixed 16-model general symbolic catalogue")
        if args.reset_numerical_configs:
            parser.error("--reset-numerical-configs requires restricted numerical production")
    return requested


def _run_general_symbolic(args: argparse.Namespace) -> int:
    """Delegate to the already-working 16-model symbolic/report workflow."""
    from scripts.RunT3GeneralSymbolicStudy import main as run_general

    argv = ["--study", args.study or "full16_general_symbolic"]
    if args.reports_only:
        argv.append("--reports-only")
    print("Pipeline: full general BSM symbolic Lagrangian + UV/EFT1 RGBeta RGEs")
    if args.threshold or args.threshold_scale:
        print("Note: threshold arguments do not perform general-theory matching; "
              "symbolic UV/EFT1 reports only.")
    if args.debug_reports:
        print("Note: --debug-reports is a legacy report option; "
              "general symbolic reports are written automatically.")
    result = run_general(argv)
    if result != 0:
        return result
    from Reports.GeneralSymbolicAdapter import render
    try:
        pdfs = render(study=args.study or "full16_general_symbolic")
    except (OSError, ValueError, RuntimeError, KeyError) as exc:
        print(f"General symbolic PDF generation failed: {exc}")
        return 1
    print(f"General symbolic PDF reports: {len(pdfs)}")
    return 0


def main() -> int:
    parser = _build_argument_parser()
    parser.add_argument(
        "--general-symbolic", action="store_true",
        help="Use the complete general-interaction 16-model UV/EFT1 symbolic workflow",
    )
    parser.add_argument(
        "--restricted-production", action="store_true",
        help="Keep the historical restricted Lagrangian, matching, and numerical pipeline",
    )
    parser.add_argument(
        "--reports-only", action="store_true",
        help="Regenerate general-symbolic reports from saved outputs, without Wolfram",
    )
    args = parser.parse_args()

    if _use_general_symbolic(args, parser):
        return _run_general_symbolic(args)
    if args.full:
        return run_full_study(args)
    if args.reports_only:
        parser.error("--reports-only requires --hypercharge-comparison")

    shared_scalar_mode = _resolve_model_mode(parser, args)
    pipeline_plan = _resolve_pipeline_plan(
        parser, args, shared_scalar_mode=shared_scalar_mode,
    )
    study_name = _study_name(args)
    study_output_dir = OUTPUT_DIR / study_name
    study_report_dir = REPORT_OUTPUT_DIR / study_name
    study_output_dir.mkdir(parents=True, exist_ok=True)
    study_report_dir.mkdir(parents=True, exist_ok=True)
    for line in pipeline_plan.description_lines():
        print(line)
    print(f"Study: {study_name}")
    print(f"Raw output: {study_output_dir}")
    print(f"Reports: {study_report_dir}")
    model_requests = _select_study_models(args, shared_scalar_mode=shared_scalar_mode)
    try:
        records = _build_and_match_t3_models(
            model_requests,
            pipeline_plan=pipeline_plan,
            study_output_dir=study_output_dir,
            debug_reports=args.debug_reports,
            force=args.force,
        )
    except ValueError as exc:
        parser.error(str(exc))
        raise AssertionError("argparse.error() should not return") from exc
    _attach_threshold_metadata(records, pipeline_plan=pipeline_plan)
    for record in records:
        organise_matched_weinberg_coefficient(record)
    physics_failed = _run_intermediate_eft_stages(
        records, pipeline_plan=pipeline_plan, debug_reports=args.debug_reports,
    )
    physics_failed |= _run_low_energy_neutrino_stages(
        records,
        debug_reports=args.debug_reports,
        numerical_input=args.numerical,
        reset_numerical_configs=args.reset_numerical_configs,
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
