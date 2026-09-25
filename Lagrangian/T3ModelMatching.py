"""Build T3 UV Lagrangians and run the configured EFT matching sequence.

This module is the Python boundary to the existing Wolfram/Matchete model
builder in :mod:`Lagrangian.Runner`.  Study selection is intentionally handled
elsewhere; each function here receives an already selected
:class:`model.T3Study.T3ModelRequest` and returns the corresponding
:class:`common.RunRecords.RunRecord`.

No matching formula is reimplemented here.  The current, verified Runner entry
points and their output/report conventions are retained exactly.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

from common.PipelinePlan import PipelinePlan
from common.RunRecords import RunRecord
from Lagrangian.Runner import (
    obtain_class_dimensions,
    validate_dimensions,
    validate_shared_dimensions,
)
from model.T3Study import T3ModelRequest


def build_and_match_t3_model(
    request: T3ModelRequest,
    *,
    pipeline_plan: PipelinePlan,
    study_output_dir: Path,
    debug_reports: bool = False,
    force: bool = False,
) -> RunRecord:
    """Build and match one selected T3 model using the existing Runner.

    ``PipelinePlan`` remains the Python source of truth.  Its historical tuple
    threshold representation is exposed only here because the Wolfram layer
    still consumes that compatibility format.
    """

    threshold_plan = pipeline_plan.threshold_plan

    if request.kind == "class":
        assert request.model_class is not None
        return obtain_class_dimensions(
            request.model_class,
            request.alpha,
            debug_reports,
            False,
            threshold_plan,
            study_output_dir,
            force,
        )

    if request.kind == "shared_dimensions":
        d_s, d_f = request.dimensions
        return validate_shared_dimensions(
            d_s,
            d_f,
            debug_reports=debug_reports,
            export_rge_tensors=False,
            threshold_plan=threshold_plan,
            output_root=study_output_dir,
            force=force,
        )

    d_s1, d_s2, d_f = request.dimensions
    return validate_dimensions(
        d_s1,
        d_s2,
        d_f,
        request.alpha,
        debug_reports,
        False,
        threshold_plan,
        study_output_dir,
        force,
    )


def build_and_match_t3_models(
    requests: Iterable[T3ModelRequest],
    *,
    pipeline_plan: PipelinePlan,
    study_output_dir: Path,
    debug_reports: bool = False,
    force: bool = False,
) -> list[RunRecord]:
    """Build and match all selected T3 models in their requested order."""

    return [
        build_and_match_t3_model(
            request,
            pipeline_plan=pipeline_plan,
            study_output_dir=study_output_dir,
            debug_reports=debug_reports,
            force=force,
        )
        for request in requests
    ]
