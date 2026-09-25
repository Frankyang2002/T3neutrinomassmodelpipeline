"""Production-scope checks for scalar-first threshold ordering.

Scalar-first matching can require leading intermediate operators above the
project's dimension-five production truncation. The current honours-project
pipeline therefore keeps the generic threshold data model but rejects
scalar-first execution before launching matching/RGE production.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from common.PipelineCLI import (
    build_argument_parser,
    resolve_model_mode,
    resolve_pipeline_plan,
)
from common.PipelinePlan import PipelinePlan


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_fermion_first_d5_remains_supported_production() -> None:
    plan = PipelinePlan.from_threshold_configuration(
        [["F"], ["S1", "S2"]],
        ["MF", "MS"],
    )

    assert plan.is_verified_fermion_first_hierarchy is True
    assert plan.is_supported_production_order is True
    assert plan.truncation.max_operator_dimension == 5


def test_scalar_first_is_rejected_by_the_user_facing_pipeline() -> None:
    parser = build_argument_parser()
    args = parser.parse_args(
        [
            "--dims", "2", "2", "1",
            "--threshold", "S1",
            "--threshold", "F", "S2",
        ]
    )
    shared = resolve_model_mode(parser, args)

    with pytest.raises(SystemExit) as exc_info:
        resolve_pipeline_plan(parser, args, shared_scalar_mode=shared)

    assert exc_info.value.code == 2


def test_generic_plan_can_describe_scalar_first_without_authorising_production() -> None:
    plan = PipelinePlan.from_threshold_configuration(
        [["S1"], ["F", "S2"]],
        ["MS1", "M_F_S2"],
    )

    assert plan.threshold_plan == (("S1",), ("F", "S2"))
    assert plan.is_supported_production_order is False
    assert "Scalar-first" in plan.production_scope_error()


def test_retired_scalar_first_experimental_modules_are_absent() -> None:
    retired = (
        PROJECT_ROOT / "RGE" / "general" / "Psi2Phi3RGE.py",
        PROJECT_ROOT / "RGE" / "running" / "intermediate" / "ScalarFirstDimensionSixSeed.py",
        PROJECT_ROOT / "tests" / "test_psi2phi3_rge.py",
        PROJECT_ROOT / "tests" / "test_scalar_first_dimension_six_foundation.py",
    )

    assert all(not path.exists() for path in retired)
