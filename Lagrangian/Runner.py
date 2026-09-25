"""Public Python boundary for T3 model construction and matching.

The historical entry points in this module are retained as compatibility API:

- ``run_model``;
- ``validate_dimensions``;
- ``validate_shared_dimensions``;
- ``obtain_class_dimensions``.

Detailed responsibilities are delegated to:

- ``Lagrangian.ModelValidation`` for representation checks and model naming;
- ``Lagrangian.WolframRunner`` for external ``wolframscript`` execution;
- ``Lagrangian.MatchingResults`` for summary ingestion and ``RunRecord``
  construction.

No matching formula is implemented in Python here.
"""

from __future__ import annotations

from pathlib import Path

from common.RunRecords import RunRecord
from common.Thresholds import ThresholdPlan
from Lagrangian.MatchingResults import (
    build_run_record,
    physicalize_shared_scalar_summary as _physicalize_shared_scalar_summary,
)
from Lagrangian.ModelValidation import (
    dimensions_for_class,
    ordinary_run_specification,
    shared_run_specification,
)
from Lagrangian.WolframRunner import (
    EFT_ORDER,
    LOOP_ORDER,
    OUTPUT_DIR,
    PROJECT_ROOT,
    RUN_MODEL_SCRIPT,
    run_wolfram_model,
)


def run_model(
    name: str,
    alpha: int,
    d_s1: int,
    d_s2: int,
    d_f: int,
    output_dir: Path,
    model_args: list[str],
    debug_reports: bool = False,
    export_rge_tensors: bool = False,
    threshold_plan: ThresholdPlan | None = None,
    shared_scalar: bool = False,
) -> RunRecord:
    """Run one already-specified T3 model through the Wolfram matching boundary."""

    process, resolved_plan = run_wolfram_model(
        output_dir,
        model_args,
        debug_reports=debug_reports,
        export_rge_tensors=export_rge_tensors,
        threshold_plan=threshold_plan,
        shared_scalar=shared_scalar,
    )

    return build_run_record(
        name=name,
        alpha=alpha,
        d_s1=d_s1,
        d_s2=d_s2,
        d_f=d_f,
        output_dir=output_dir,
        shared_scalar=shared_scalar,
        process=process,
        threshold_plan=resolved_plan,
    )


def validate_dimensions(
    d_s1: int,
    d_s2: int,
    d_f: int,
    alpha: int,
    debug_reports: bool = False,
    export_rge_tensors: bool = False,
    threshold_plan: ThresholdPlan | None = None,
    output_root: Path | None = None,
    force: bool = False,
) -> RunRecord:
    """Validate, name, and run one ordinary T3 representation."""

    specification = ordinary_run_specification(
        d_s1,
        d_s2,
        d_f,
        alpha,
        output_root=output_root or OUTPUT_DIR,
        force=force,
    )

    print(
        f"Running {specification.name}, "
        f"dims=({d_s1}, {d_s2}, {d_f}), "
        f"alpha={alpha} ...",
        flush=True,
    )

    return run_model(
        specification.name,
        specification.alpha,
        specification.d_s1,
        specification.d_s2,
        specification.d_f,
        specification.output_dir,
        list(specification.model_args),
        debug_reports,
        export_rge_tensors,
        threshold_plan,
    )


def validate_shared_dimensions(
    d_s: int,
    d_f: int,
    *,
    debug_reports: bool = False,
    export_rge_tensors: bool = False,
    threshold_plan: ThresholdPlan | None = None,
    output_root: Path | None = None,
    force: bool = False,
) -> RunRecord:
    """Validate, name, and run the one-physical-scalar scotogenic branch."""

    specification = shared_run_specification(
        d_s,
        d_f,
        output_root=output_root or OUTPUT_DIR,
        force=force,
    )

    print(
        f"Running {specification.name}, "
        f"shared scalar dims=({d_s}, {d_f}), alpha=-1 ...",
        flush=True,
    )

    return run_model(
        specification.name,
        specification.alpha,
        specification.d_s1,
        specification.d_s2,
        specification.d_f,
        specification.output_dir,
        list(specification.model_args),
        debug_reports,
        export_rge_tensors,
        threshold_plan,
        shared_scalar=True,
    )


def obtain_class_dimensions(
    model_class: str,
    alpha: int,
    debug_reports: bool = False,
    export_rge_tensors: bool = False,
    threshold_plan: ThresholdPlan | None = None,
    output_root: Path | None = None,
    force: bool = False,
) -> RunRecord:
    """Convert a known A-E class into dimensions and then run normally."""

    d_s1, d_s2, d_f = dimensions_for_class(model_class)

    return validate_dimensions(
        d_s1,
        d_s2,
        d_f,
        alpha,
        debug_reports,
        export_rge_tensors,
        threshold_plan,
        output_root,
        force,
    )
