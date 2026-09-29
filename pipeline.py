# Central backbone for the T3 neutrino-mass calculation.
#
# Read this file to see the project flow and execution order:
# 1. Parse the requested T3 study and threshold plan.
# 2. Build the UV model(s) and perform matching.
# 3. Organise the matched Weinberg coefficient C5.
# 4. Run the UV theory and each intermediate EFT region.
# 5. Run the low-energy Weinberg theory and construct the neutrino mass.
# 6. Optionally run numerical evolution and calculate observables/figures.
# 7. Generate summaries and reports.
#
# Detailed model selection, beta functions, matching algebra, validation, and
# report formatting live in specialised modules; this file owns their order.
from __future__ import annotations

import argparse
from copy import deepcopy
import json
from pathlib import Path
from typing import Any

from studies.FullT3Study import run_full_study

from RGE.running.IntermediateEFTRunning import (
    mark_intermediate_running_not_applicable,
    run_intermediate_eft_interval,
)
from validation.IntermediateEFTValidation import (
    run_intermediate_eft_validation,
)

from common.EFT import EFTRunningInterval
from common.PipelinePlan import PipelinePlan
from common.PipelineCLI import (
    build_argument_parser as _build_argument_parser,
    resolve_model_mode as _resolve_model_mode,
    resolve_pipeline_plan as _resolve_pipeline_plan,
    study_name as _study_name,
)
from common.RunRecords import EFTStageRecord, RunRecord
from Numerical.PipelineNumericalResults import (
    is_pipeline_numerical_results_config,
    run_pipeline_numerical_results,
)
from Numerical.SobolBenchmarkSearch import (
    _parameters_from_config,
    apply_benchmark_parameters,
    model_benchmark_key,
    retarget_payload_to_record,
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
from Reports.PipelineReports import (
    finish_pipeline_run,
    print_pipeline_summary,
)


PROJECT_ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = PROJECT_ROOT / "output"
REPORT_OUTPUT_DIR = PROJECT_ROOT / "Reports" / "output"
GENERATED_NUMERICAL_CONFIG_DIR = PROJECT_ROOT / "configs" / "generated_models"


def _stage_for_running_interval(
    record: RunRecord,
    interval: EFTRunningInterval,
) -> EFTStageRecord:
    """Return the recorded EFT stage corresponding to one running interval."""
    return record.stage_with_active_fields(
        *sorted(interval.active_heavy_fields)
    )


def print_summary(records: list[RunRecord], aggregate_dir: Path | None = None) -> int:
    """Compatibility wrapper for the pipeline-level run summary."""
    return print_pipeline_summary(records, aggregate_dir)


def finish_runs(
    records: list[RunRecord],
    debug_reports: bool = False,
    physics_failed: bool = False,
    *,
    study_output_dir: Path | None = None,
    study_report_dir: Path | None = None,
) -> int:
    """Compatibility wrapper for summary and report generation."""
    return finish_pipeline_run(
        records,
        debug_reports=debug_reports,
        physics_failed=physics_failed,
        study_output_dir=study_output_dir,
        study_report_dir=study_report_dir,
    )


def _load_json_object(path: Path) -> dict[str, Any]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path} must contain a JSON object.")
    return payload


def _atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def _generated_numerical_config_path(record: RunRecord) -> Path:
    return GENERATED_NUMERICAL_CONFIG_DIR / f"{model_benchmark_key(record)}.json"


def _prepare_model_numerical_config(
    record: RunRecord,
    template_path: Path,
    *,
    reset: bool,
) -> Path:
    """Return the persistent model-specific numerical config for one record.

    Existing generated configs are reused by default.  With reset=True the
    supplied numerical template is retargeted to the current model and replaces
    the generated config before the numerical run.
    """
    generated_path = _generated_numerical_config_path(record)

    if generated_path.is_file() and not reset:
        print(
            f"  {record.name}: using saved numerical config -> {generated_path}",
            flush=True,
        )
        return generated_path

    template = _load_json_object(Path(template_path))
    retargeted = retarget_payload_to_record(record, template)
    retargeted["generated_config"] = {
        "model_key": model_benchmark_key(record),
        "source_template": str(Path(template_path).resolve()),
        "resolved_benchmark": False,
    }
    _atomic_write_json(generated_path, retargeted)
    print(
        f"  {record.name}: saved model numerical config -> {generated_path}",
        flush=True,
    )
    return generated_path


def _freeze_resolved_benchmark_config(
    record: RunRecord,
    generated_path: Path,
) -> None:
    """Freeze the benchmark actually used into the persistent model config.

    The integrated numerical run writes benchmark provenance into
    running_diagnostics.json.  When that provenance points to a Sobol cache,
    apply the cache's best parameters to the generated config and disable
    benchmark_search so future runs reproduce that exact point.
    """
    diagnostics_relative = record.summary.get("RunningDiagnosticsFile")
    if not diagnostics_relative:
        return

    diagnostics_path = record.output_dir / str(diagnostics_relative)
    if not diagnostics_path.is_file():
        return

    diagnostics = _load_json_object(diagnostics_path)
    metadata = diagnostics.get("benchmark_search")
    if not isinstance(metadata, dict):
        return

    cache_value = metadata.get("cache_path")
    if not cache_value:
        return

    cache_path = Path(str(cache_value))
    if not cache_path.is_absolute():
        cache_path = PROJECT_ROOT / cache_path
    if not cache_path.is_file():
        return

    cache = _load_json_object(cache_path)
    best = cache.get("best")
    if not isinstance(best, dict) or not isinstance(best.get("parameters"), dict):
        return

    payload = _load_json_object(generated_path)
    definitions = _parameters_from_config(payload)
    resolved = apply_benchmark_parameters(
        payload,
        best["parameters"],
        definitions,
    )

    search = resolved.get("benchmark_search")
    if isinstance(search, dict):
        search = deepcopy(search)
        search["enabled"] = False
        resolved["benchmark_search"] = search

    resolved["generated_config"] = {
        "model_key": model_benchmark_key(record),
        "source_template": payload.get("generated_config", {}).get(
            "source_template"
        ),
        "resolved_benchmark": True,
        "benchmark_cache": str(cache_path.resolve()),
        "benchmark_status": metadata.get("benchmark_status"),
        "benchmark_chi2": metadata.get("chi2"),
        "benchmark_method": metadata.get("method"),
    }
    _atomic_write_json(generated_path, resolved)
    print(
        f"  {record.name}: froze resolved benchmark config -> {generated_path}",
        flush=True,
    )


def _attach_threshold_metadata(
    records: list[RunRecord],
    *,
    pipeline_plan: PipelinePlan,
) -> None:
    """Attach stage records and stable threshold metadata to each model run.

    Existing summary keys and report-stage labels are preserved.  The plan now
    also records the explicit EFT truncation (currently ``d<=5``).
    """

    metadata = pipeline_plan.summary_metadata()

    for record in records:
        if record.shared_scalar != pipeline_plan.shared_scalar:
            raise ValueError(
                "RunRecord shared-scalar mode does not match PipelinePlan."
            )

        record.eft_stages = pipeline_plan.build_stage_records(record.output_dir)

        for key, value in metadata.items():
            record.summary.setdefault(key, value)


def _run_intermediate_eft_stages(
    records: list[RunRecord],
    *,
    pipeline_plan: PipelinePlan,
    debug_reports: bool,
) -> bool:
    """Run UV RGEs and every configured intermediate heavy-field region.

    ``pipeline.py`` owns the order of operations.  Production RGE/matching is
    delegated to ``run_intermediate_eft_interval`` and independent regression
    checks are then delegated to ``run_intermediate_eft_validation``.  Both
    dispatch from physical EFT content rather than an EFT stage number.
    """

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
                record,
                stage,
                interval,
                debug_reports=debug_reports,
            )

            attempted_any |= outcome.attempted

            if not outcome.attempted:
                production_complete = False
                physics_failed = True
                unsupported_intervals.append(
                    {
                        "ActiveHeavyFields": sorted(interval.active_heavy_fields),
                        "HighScale": interval.high_scale,
                        "LowScale": interval.low_scale,
                        "EnteredBy": list(interval.entered_by.fields_to_integrate),
                        "ExitedBy": list(interval.exited_by.fields_to_integrate),
                        "Status": outcome.status,
                    }
                )
                continue

            if not outcome.success:
                production_complete = False
                physics_failed = True
                continue

            validation = run_intermediate_eft_validation(
                record,
                stage,
                interval,
            )
            if validation.attempted and not validation.success:
                physics_failed = True

        if not pipeline_plan.running_intervals:
            record.summary["IntermediateEFTProductionStatus"] = "NotRequired"
        elif production_complete:
            record.summary["IntermediateEFTProductionStatus"] = "Complete"
        elif unsupported_intervals:
            record.summary["IntermediateEFTProductionStatus"] = (
                "IncompleteUnsupportedFieldContent"
            )
        else:
            record.summary["IntermediateEFTProductionStatus"] = "Failed"

        record.summary["IntermediateEFTUnsupportedIntervals"] = unsupported_intervals

        if not attempted_any:
            mark_intermediate_running_not_applicable(record)

    return physics_failed


def _ready_for_low_energy_neutrino_stages(record: RunRecord) -> bool:
    """Return whether matching produced the inputs needed by final C5 stages."""

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
    """Run the low-energy Weinberg, neutrino-mass, and optional numerical stages."""

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
                    record,
                    numerical_input,
                    reset=reset_numerical_configs,
                )
                if not run_pipeline_numerical_results(
                    record,
                    model_config,
                ):
                    physics_failed = True
                else:
                    _freeze_resolved_benchmark_config(record, model_config)
            elif not run_numerical_neutrino_observables(
                record,
                numerical_input,
            ):
                physics_failed = True

    return physics_failed


def main() -> int:
    """Parse one pipeline run, execute its stages, and generate reports."""

    parser = _build_argument_parser()
    args = parser.parse_args()

    if args.full:
        return run_full_study(args)

    shared_scalar_mode = _resolve_model_mode(parser, args)
    pipeline_plan = _resolve_pipeline_plan(
        parser,
        args,
        shared_scalar_mode=shared_scalar_mode,
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

    model_requests = _select_study_models(
        args,
        shared_scalar_mode=shared_scalar_mode,
    )
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
    _attach_threshold_metadata(
        records,
        pipeline_plan=pipeline_plan,
    )

    for record in records:
        organise_matched_weinberg_coefficient(record)

    physics_failed = _run_intermediate_eft_stages(
        records,
        pipeline_plan=pipeline_plan,
        debug_reports=args.debug_reports,
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
