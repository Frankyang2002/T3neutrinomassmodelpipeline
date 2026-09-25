from __future__ import annotations

import pytest

from common.EFT import EFTRunningInterval
from common.PipelinePlan import PipelinePlan


def test_common_threshold_has_no_intermediate_running_interval() -> None:
    plan = PipelinePlan.from_threshold_configuration([["F", "S1", "S2"]])

    assert plan.running_intervals == ()


def test_fermion_first_interval_is_scalar_only() -> None:
    plan = PipelinePlan.from_threshold_configuration(
        [["F"], ["S1", "S2"]],
        ["MF", "MS"],
    )

    (interval,) = plan.running_intervals

    assert isinstance(interval, EFTRunningInterval)
    assert interval.has_exactly("S1", "S2")
    assert interval.high_scale == "MF"
    assert interval.low_scale == "MS"
    assert interval.entered_by.fields_to_integrate == ("F",)
    assert interval.exited_by.fields_to_integrate == ("S1", "S2")


def test_scalar_first_plan_exposes_each_running_region() -> None:
    plan = PipelinePlan.from_threshold_configuration(
        [["S1"], ["F"], ["S2"]],
        ["MS1", "MF", "MS2"],
    )

    first, second = plan.running_intervals

    assert first.has_exactly("F", "S2")
    assert first.high_scale == "MS1"
    assert first.low_scale == "MF"

    assert second.has_exactly("S2")
    assert second.high_scale == "MF"
    assert second.low_scale == "MS2"


def test_shared_scalar_plan_uses_physical_scalar_content() -> None:
    plan = PipelinePlan.from_threshold_configuration(
        [["F"], ["S"]],
        ["MF", "MS"],
        shared_scalar=True,
    )

    (interval,) = plan.running_intervals
    assert interval.has_exactly("S")


def test_running_interval_rejects_non_adjacent_thresholds() -> None:
    plan = PipelinePlan.from_threshold_configuration(
        [["S1"], ["F"], ["S2"]],
    )
    first, second, third = plan.transitions

    with pytest.raises(ValueError, match="adjacent"):
        EFTRunningInterval(
            index=first.index,
            content=first.after,
            high_scale=first.scale,
            low_scale=third.scale,
            entered_by=first,
            exited_by=third,
        )
