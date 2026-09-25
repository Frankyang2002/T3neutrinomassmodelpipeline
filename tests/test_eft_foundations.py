from __future__ import annotations

from itertools import permutations
from pathlib import Path

import pytest

from common.EFT import DEFAULT_EFT_TRUNCATION, EFTTruncation, ThresholdStep
from common.RunRecords import RunRecord
from common.Thresholds import (
    build_eft_stage_records,
    build_threshold_steps,
    threshold_plan_for_wolfram,
    threshold_plan_label,
    validate_threshold_plan,
)


def test_default_eft_truncation_is_explicitly_dimension_five() -> None:
    assert DEFAULT_EFT_TRUNCATION.max_operator_dimension == 5
    assert DEFAULT_EFT_TRUNCATION.label == "d<=5"
    assert DEFAULT_EFT_TRUNCATION.keeps(4)
    assert DEFAULT_EFT_TRUNCATION.keeps(5)
    assert not DEFAULT_EFT_TRUNCATION.keeps(6)


def test_invalid_operator_truncation_is_rejected() -> None:
    with pytest.raises(ValueError, match="positive"):
        EFTTruncation(max_operator_dimension=0)


def test_threshold_step_is_model_independent() -> None:
    step = ThresholdStep(fields_to_integrate=("S1", "F"), scale="M_S1F")
    assert step.fields_to_integrate == ("S1", "F")
    assert step.scale == "M_S1F"
    assert step.label == "(S1,F)"


def test_every_single_field_threshold_order_is_valid() -> None:
    for order in permutations(("F", "S1", "S2")):
        plan = validate_threshold_plan([[field] for field in order])
        assert tuple(group[0] for group in plan) == order


def test_grouped_scalar_first_plan_is_valid() -> None:
    plan = validate_threshold_plan([["S1", "S2"], ["F"]])
    assert plan == (("S1", "S2"), ("F",))
    assert threshold_plan_label(plan) == "(S1,S2) -> F"


def test_threshold_steps_preserve_order_and_scale() -> None:
    plan = validate_threshold_plan([["S1"], ["F"], ["S2"]])
    steps = build_threshold_steps(plan, ["MS1", "MF", "MS2"])

    assert [step.fields_to_integrate for step in steps] == [
        ("S1",),
        ("F",),
        ("S2",),
    ]
    assert [step.scale for step in steps] == ["MS1", "MF", "MS2"]


def test_threshold_step_count_must_match_scale_count() -> None:
    plan = validate_threshold_plan([["F"], ["S1", "S2"]])
    with pytest.raises(ValueError, match="once for each threshold group"):
        build_threshold_steps(plan, ["MF"])


def test_scalar_first_stage_metadata_tracks_physical_content() -> None:
    plan = validate_threshold_plan([["S1"], ["F"], ["S2"]])
    stages = build_eft_stage_records(plan)

    assert stages[0].label == "UV"
    assert stages[0].active_heavy_fields == ("F", "S1", "S2")

    assert stages[1].label == "EFT_1_after_S1"
    assert stages[1].integrated_fields == ("S1",)
    assert stages[1].active_field_set == frozenset({"F", "S2"})

    assert stages[2].label == "EFT_2_after_F"
    assert stages[2].integrated_fields == ("F",)
    assert stages[2].active_field_set == frozenset({"S2"})

    assert stages[3].label == "EFT_3_after_S2"
    assert stages[3].integrated_fields == ("S2",)
    assert stages[3].active_field_set == frozenset()


def test_run_record_can_find_stages_by_physical_content() -> None:
    plan = validate_threshold_plan([["S1"], ["F"], ["S2"]])
    stages = build_eft_stage_records(plan)
    record = RunRecord(
        name="test",
        alpha=0,
        d_s1=1,
        d_s2=3,
        d_f=2,
        return_code=0,
        summary={},
        output_dir=Path("output/test"),
        eft_stages=stages,
    )

    assert record.uv_stage.label == "UV"
    assert record.stage_after_integrating("S1").label == "EFT_1_after_S1"
    assert record.stage_with_active_fields("S2").label == "EFT_2_after_F"
    assert record.final_eft_stage.label == "EFT_3_after_S2"


def test_existing_first_eft_stage_property_is_preserved() -> None:
    plan = validate_threshold_plan([["F"], ["S1", "S2"]])
    record = RunRecord(
        name="test",
        alpha=-1,
        d_s1=2,
        d_s2=2,
        d_f=1,
        return_code=0,
        summary={},
        output_dir=Path("output/test"),
        eft_stages=build_eft_stage_records(plan),
    )

    assert record.first_eft_stage.label == "EFT_1_after_F"


def test_shared_scalar_plan_uses_physical_S_but_wolfram_gets_formal_roles() -> None:
    plan = validate_threshold_plan([["S"], ["F"]], shared_scalar=True)
    assert plan == (("S",), ("F",))
    assert threshold_plan_for_wolfram(plan, shared_scalar=True) == (
        ("S1", "S2"),
        ("F",),
    )

    stages = build_eft_stage_records(plan, shared_scalar=True)
    assert stages[0].active_field_set == frozenset({"F", "S"})
    assert stages[1].active_field_set == frozenset({"F"})
    assert stages[2].active_field_set == frozenset()
