"""Pipeline-level run summaries and human-readable report generation.

This module owns presentation and report-file orchestration only. Physics
stages remain in ``pipeline.py`` and the dedicated physics/RGE modules.

The public functions preserve the existing command-line output, aggregate JSON
schema, report paths, and report-generation order while keeping those details
out of the central pipeline backbone.
"""

from __future__ import annotations

import json
from pathlib import Path

from common.RunRecords import RunRecord
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


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "output"
DEFAULT_REPORT_OUTPUT_DIR = PROJECT_ROOT / "Reports" / "output"


def print_pipeline_summary(records: list[RunRecord], aggregate_dir: Path | None = None) -> int:
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
    aggregate_dir = aggregate_dir or DEFAULT_OUTPUT_DIR
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



def _unique_eft_stage_labels(records: list[RunRecord]) -> list[str]:
    """Return recorded EFT-stage labels in first-seen order."""
    labels: list[str] = []
    for record in records:
        for stage in record.eft_stages:
            if stage.is_uv or stage.label in labels:
                continue
            labels.append(stage.label)
    return labels


def _report_stem(path: Path, report_root: Path) -> str:
    """Return one generated report path relative to the study report root."""
    try:
        relative = path.relative_to(report_root)
    except ValueError:
        relative = path
    return relative.with_suffix("").as_posix()


def _print_generated_report_summary(
    records: list[RunRecord],
    *,
    report_root: Path,
    uv_lagrangian_tex: Path,
    c5_tex: Path,
    uv_rge_tex: Path,
    intermediate_rge_tex: Path | None,
    final_rge_tex: Path,
    group_factor_tex: list[Path],
) -> None:
    """Print only report surfaces that were actually generated.

    Older pipeline output hard-coded ``EFT_1_after_F`` even for common-threshold
    studies where that stage did not exist.  Report filenames themselves are
    unchanged; this helper makes the terminal summary follow the actual
    threshold plan.
    """
    lines = ["", "Stage-aware reports:"]
    lines.append(f"  {_report_stem(uv_lagrangian_tex, report_root)}")

    for stage_label in _unique_eft_stage_labels(records):
        lines.append(f"  Lagrangian/{stage_label}")

    lines.append(f"  {_report_stem(c5_tex, report_root)}")
    lines.append(f"  {_report_stem(uv_rge_tex, report_root)}")

    if intermediate_rge_tex is not None:
        lines.append(f"  {_report_stem(intermediate_rge_tex, report_root)}")

    lines.append(f"  {_report_stem(final_rge_tex, report_root)}")
    lines.extend(
        f"  {_report_stem(path, report_root)}"
        for path in group_factor_tex
    )

    print("\n".join(lines))

def finish_pipeline_run(
    records: list[RunRecord],
    debug_reports: bool = False,
    physics_failed: bool = False,
    *,
    study_output_dir: Path | None = None,
    study_report_dir: Path | None = None,
) -> int:
    """Print final summaries and generate reports after the physics pipeline."""

    study_output_dir = study_output_dir or DEFAULT_OUTPUT_DIR
    study_report_dir = study_report_dir or DEFAULT_REPORT_OUTPUT_DIR

    status = print_pipeline_summary(records, study_output_dir)

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
    # Historical report-writer names are retained for compatibility, but the
    # physical applicability of the intermediate report is decided from the
    # recorded stage content inside RGEComparison.
    uv_rge_tex = write_and_compile_rge_comparison(
        records, report_root=study_report_dir
    )
    intermediate_rge_tex = write_and_compile_eft1_rge_comparison(
        records, report_root=study_report_dir
    )
    final_rge_tex = write_and_compile_final_eft_rge_comparison(
        records, report_root=study_report_dir
    )

    # Stage-aware analytic group-factor reports.
    group_factor_tex = write_and_compile_stage_group_factor_reports(
        records,
        report_root=study_report_dir,
    )

    _print_generated_report_summary(
        records,
        report_root=study_report_dir,
        uv_lagrangian_tex=uv_lagrangian_tex,
        c5_tex=c5_tex,
        uv_rge_tex=uv_rge_tex,
        intermediate_rge_tex=intermediate_rge_tex,
        final_rge_tex=final_rge_tex,
        group_factor_tex=group_factor_tex,
    )

    return status
