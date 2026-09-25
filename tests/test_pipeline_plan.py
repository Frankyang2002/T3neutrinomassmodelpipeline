from __future__ import annotations

from pathlib import Path

import pytest

from common.EFT import EFTTruncation
from common.PipelinePlan import PipelinePlan
from common.Thresholds import (
    default_threshold_scale,
    resolve_threshold_scales,
)


def test_default_common_threshold_plan_preserves_existing_scale_convention() -> None:
    plan = PipelinePlan.from_threshold_configuration(None)

    assert plan.threshold_plan == (("F", "S1", "S2"),)
    assert plan.threshold_scales == ("M_F_S1_S2",)
    assert plan.label == "(F,S1,S2)"


def test_fermion_first_plan_retains_existing_default_scales() -> None:
    plan = PipelinePlan.from_threshold_configuration(
        [["F"], ["S1", "S2"]]
    )

    assert plan.threshold_plan == (("F",), ("S1", "S2"))
    assert plan.threshold_scales == ("MF", "MS")
    assert plan.label == "F -> (S1,S2)"


def test_scalar_first_plan_uses_same_generic_plan_object() -> None:
    plan = PipelinePlan.from_threshold_configuration(
        [["S1"], ["F"], ["S2"]]
    )

    assert plan.threshold_plan == (("S1",), ("F",), ("S2",))
    assert plan.threshold_scales == ("MS1", "MF", "MS2")

    lines = plan.description_lines()
    assert lines[0] == "EFT truncation: d<=5"
    assert lines[1] == "Threshold plan: S1 -> F -> S2"
    assert "active heavy fields: F, S2" in lines[3]
    assert lines[-1].endswith("active heavy fields: none")


def test_explicit_threshold_scales_are_preserved() -> None:
    plan = PipelinePlan.from_threshold_configuration(
        [["S2"], ["S1"], ["F"]],
        [3.0e10, "MS1_running", 1.0e8],
    )

    assert plan.threshold_scales == (3.0e10, "MS1_running", 1.0e8)


def test_pipeline_plan_exposes_dimension_five_truncation_in_summary_metadata() -> None:
    plan = PipelinePlan.from_threshold_configuration(
        [["F"], ["S1", "S2"]]
    )

    metadata = plan.summary_metadata()

    # Historical keys are retained for report/output compatibility.
    assert metadata["ThresholdPlan"] == [["F"], ["S1", "S2"]]
    assert metadata["ThresholdPlanLabel"] == "F -> (S1,S2)"
    assert metadata["ThresholdScales"] == ["MF", "MS"]

    # The approximation that used to be implicit is now explicit metadata.
    assert metadata["EFTTruncation"] == {
        "MaxOperatorDimension": 5,
        "Label": "d<=5",
    }


def test_pipeline_plan_can_carry_a_different_truncation_without_changing_order() -> None:
    plan = PipelinePlan.from_threshold_configuration(
        [["S1"], ["F"], ["S2"]],
        truncation=EFTTruncation(max_operator_dimension=6),
    )

    assert plan.truncation.max_operator_dimension == 6
    assert plan.threshold_plan == (("S1",), ("F",), ("S2",))


def test_stage_records_preserve_historical_report_labels() -> None:
    plan = PipelinePlan.from_threshold_configuration(
        [["S1"], ["F"], ["S2"]]
    )
    stages = plan.build_stage_records(Path("output/model"))

    assert [stage.label for stage in stages] == [
        "UV",
        "EFT_1_after_S1",
        "EFT_2_after_F",
        "EFT_3_after_S2",
    ]
    assert stages[1].output_dir == Path("output/model/stages/EFT_1_after_S1")
    assert stages[-1].active_field_set == frozenset()


def test_shared_scalar_plan_uses_physical_scalar_name() -> None:
    plan = PipelinePlan.from_threshold_configuration(
        [["S"], ["F"]],
        shared_scalar=True,
    )

    assert plan.threshold_plan == (("S",), ("F",))
    assert plan.threshold_scales == ("MS", "MF")
    stages = plan.build_stage_records()
    assert stages[0].active_field_set == frozenset({"F", "S"})
    assert stages[1].active_field_set == frozenset({"F"})


def test_pipeline_plan_rejects_missing_heavy_field_even_when_constructed_directly() -> None:
    from common.EFT import ThresholdStep

    with pytest.raises(ValueError, match="missing field"):
        PipelinePlan(
            threshold_steps=(
                ThresholdStep(("S1",), "MS1"),
                ThresholdStep(("F",), "MF"),
            )
        )


def test_scale_count_mismatch_is_rejected_before_execution() -> None:
    with pytest.raises(ValueError, match="once for each threshold group"):
        PipelinePlan.from_threshold_configuration(
            [["F"], ["S1", "S2"]],
            ["MF"],
        )


def test_default_threshold_scale_matches_existing_named_cases() -> None:
    assert default_threshold_scale(("F",)) == "MF"
    assert default_threshold_scale(("S",)) == "MS"
    assert default_threshold_scale(("S1",)) == "MS1"
    assert default_threshold_scale(("S2",)) == "MS2"
    assert default_threshold_scale(("S1", "S2")) == "MS"
    assert default_threshold_scale(("F", "S1")) == "M_F_S1"


def test_resolve_threshold_scales_preserves_explicit_values() -> None:
    resolved = resolve_threshold_scales(
        (("S1",), ("F",), ("S2",)),
        (1.0e9, "MF_running", 5.0e7),
    )
    assert resolved == (1.0e9, "MF_running", 5.0e7)
