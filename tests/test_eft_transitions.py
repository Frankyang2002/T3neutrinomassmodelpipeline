from __future__ import annotations

import pytest

from common.EFT import EFTContent, EFTTransition, EFTTruncation, ThresholdStep
from common.PipelinePlan import PipelinePlan


def test_eft_content_is_identified_by_physical_active_fields() -> None:
    content = EFTContent.from_fields(["S2", "F"])

    assert content.active_heavy_fields == frozenset({"F", "S2"})
    assert content.has_exactly("F", "S2")
    assert not content.has_exactly("F", "S1")
    assert content.heavy_field_label == "F, S2"


def test_integrating_out_field_returns_new_content_without_stage_number() -> None:
    before = EFTContent.from_fields(["F", "S1", "S2"])
    after = before.after_integrating(["S1"])

    assert after.has_exactly("F", "S2")
    assert before.has_exactly("F", "S1", "S2")


def test_cannot_integrate_out_inactive_field() -> None:
    content = EFTContent.from_fields(["F", "S2"])

    with pytest.raises(ValueError, match="not active"):
        content.after_integrating(["S1"])


def test_scalar_first_pipeline_transitions_follow_physical_content() -> None:
    plan = PipelinePlan.from_threshold_configuration(
        [["S1"], ["F"], ["S2"]]
    )

    transitions = plan.transitions

    assert transitions[0].before.has_exactly("F", "S1", "S2")
    assert transitions[0].fields_to_integrate == ("S1",)
    assert transitions[0].after.has_exactly("F", "S2")

    assert transitions[1].before.has_exactly("F", "S2")
    assert transitions[1].fields_to_integrate == ("F",)
    assert transitions[1].after.has_exactly("S2")

    assert transitions[2].before.has_exactly("S2")
    assert transitions[2].fields_to_integrate == ("S2",)
    assert transitions[2].after.is_fully_decoupled


def test_fermion_first_uses_same_transition_type() -> None:
    plan = PipelinePlan.from_threshold_configuration(
        [["F"], ["S1", "S2"]]
    )

    first, second = plan.transitions

    assert isinstance(first, EFTTransition)
    assert first.before.has_exactly("F", "S1", "S2")
    assert first.after.has_exactly("S1", "S2")
    assert second.after.is_fully_decoupled


def test_grouped_scalar_first_transition_is_supported() -> None:
    plan = PipelinePlan.from_threshold_configuration(
        [["S1", "S2"], ["F"]]
    )

    first = plan.transitions[0]
    assert first.fields_to_integrate == ("S1", "S2")
    assert first.after.has_exactly("F")


def test_shared_scalar_transitions_use_physical_S() -> None:
    plan = PipelinePlan.from_threshold_configuration(
        [["S"], ["F"]],
        shared_scalar=True,
    )

    first, second = plan.transitions
    assert first.before.has_exactly("F", "S")
    assert first.after.has_exactly("F")
    assert second.after.is_fully_decoupled


def test_truncation_is_carried_through_every_transition() -> None:
    truncation = EFTTruncation(max_operator_dimension=5)
    plan = PipelinePlan.from_threshold_configuration(
        [["S2"], ["F"], ["S1"]],
        truncation=truncation,
    )

    assert plan.initial_content.truncation == truncation
    for transition in plan.transitions:
        assert transition.before.truncation == truncation
        assert transition.after.truncation == truncation


def test_transition_rejects_inconsistent_after_content() -> None:
    before = EFTContent.from_fields(["F", "S1", "S2"])
    wrong_after = EFTContent.from_fields(["S1", "S2"])

    with pytest.raises(ValueError, match="does not match"):
        EFTTransition(
            index=1,
            step=ThresholdStep(("S1",), "MS1"),
            before=before,
            after=wrong_after,
        )


def test_final_content_is_fully_decoupled_for_every_valid_plan() -> None:
    plans = (
        [["F", "S1", "S2"]],
        [["F"], ["S1", "S2"]],
        [["S1"], ["F"], ["S2"]],
        [["S2"], ["S1"], ["F"]],
        [["S1", "S2"], ["F"]],
    )

    for threshold_groups in plans:
        plan = PipelinePlan.from_threshold_configuration(threshold_groups)
        assert plan.final_content.is_fully_decoupled
