"""Physics-safety checks for scalar-first threshold ordering.

The current authoritative pipeline is truncated at dimension five.  A scalar
threshold before F can generate a dimension-six intermediate operator that later
feeds the Weinberg coefficient, so scalar-first d<=5 runs must remain explicitly
marked as truncated and cannot silently reach the authoritative low-energy path.
"""

from __future__ import annotations

from common.EFT import EFTTruncation
from common.PipelinePlan import PipelinePlan


def test_fermion_first_d5_remains_authoritative() -> None:
    plan = PipelinePlan.from_threshold_configuration(
        [["F"], ["S1", "S2"]],
        ["MF", "MS"],
    )

    assert plan.has_scalar_first_threshold is False
    assert plan.scalar_first_truncates_leading_path is False
    assert plan.summary_metadata()["ThresholdOrderingPhysics"] == {
        "ScalarFirst": False,
        "ScalarFirstDimensionSixPathRetained": False,
        "ScalarFirstLeadingPathTruncated": False,
        "AuthoritativeAtConfiguredTruncation": True,
    }


def test_scalar_first_d5_is_explicitly_non_authoritative() -> None:
    plan = PipelinePlan.from_threshold_configuration(
        [["S1"], ["F"], ["S2"]],
        ["MS1", "MF", "MS2"],
    )

    assert plan.has_scalar_first_threshold is True
    assert plan.scalar_first_truncates_leading_path is True

    physics = plan.summary_metadata()["ThresholdOrderingPhysics"]
    assert physics["ScalarFirst"] is True
    assert physics["ScalarFirstDimensionSixPathRetained"] is False
    assert physics["ScalarFirstLeadingPathTruncated"] is True
    assert physics["AuthoritativeAtConfiguredTruncation"] is False
    assert any("experimental/non-authoritative" in line for line in plan.physics_warnings())


def test_scalar_first_flag_clears_when_dimension_six_is_retained() -> None:
    plan = PipelinePlan.from_threshold_configuration(
        [["S1"], ["F"], ["S2"]],
        ["MS1", "MF", "MS2"],
        truncation=EFTTruncation(max_operator_dimension=6),
    )

    assert plan.has_scalar_first_threshold is True
    assert plan.scalar_first_truncates_leading_path is False
    physics = plan.summary_metadata()["ThresholdOrderingPhysics"]
    assert physics["ScalarFirstDimensionSixPathRetained"] is True
    assert physics["AuthoritativeAtConfiguredTruncation"] is True
