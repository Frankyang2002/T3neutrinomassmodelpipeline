"""Physical field-identity tests for ordinary and shared-scalar T3 models."""

from __future__ import annotations

from pathlib import Path

import pytest

from common.PipelinePlan import PipelinePlan
from common.RunRecords import RunRecord
from common.T3Fields import (
    FORMAL_SCALAR_FIELDS,
    ORDINARY_T3_FIELDS,
    SHARED_SCALAR_T3_FIELDS,
    t3_heavy_field_scheme,
)
from common.Thresholds import normalise_threshold_field, threshold_plan_for_wolfram


def test_ordinary_scheme_has_three_distinct_physical_heavy_fields() -> None:
    scheme = ORDINARY_T3_FIELDS
    assert scheme.heavy_fields == ("F", "S1", "S2")
    assert scheme.scalar_fields == ("S1", "S2")
    assert scheme.scalar_field_set == frozenset({"S1", "S2"})


def test_shared_scheme_has_one_physical_scalar() -> None:
    scheme = SHARED_SCALAR_T3_FIELDS
    assert scheme.heavy_fields == ("F", "S")
    assert scheme.scalar_fields == ("S",)
    assert scheme.scalar_field_set == frozenset({"S"})


def test_shared_physical_scalar_expands_only_at_formal_matching_boundary() -> None:
    scheme = SHARED_SCALAR_T3_FIELDS
    assert FORMAL_SCALAR_FIELDS == ("S1", "S2")
    assert scheme.formal_roles_for("S") == ("S1", "S2")
    assert scheme.formal_roles_for("F") == ("F",)

    plan = (("F",), ("S",))
    assert threshold_plan_for_wolfram(plan, shared_scalar=True) == (
        ("F",),
        ("S1", "S2"),
    )


def test_shared_cli_aliases_resolve_to_physical_s() -> None:
    for alias in ("S", "scalar", "S1", "S2", "scalar_1", "scalar_2"):
        assert normalise_threshold_field(alias, shared_scalar=True) == "S"

    with pytest.raises(ValueError, match="F or S"):
        normalise_threshold_field("X", shared_scalar=True)


def test_ordinary_scheme_does_not_invent_shared_s_alias() -> None:
    with pytest.raises(ValueError, match="F, S1, or S2"):
        normalise_threshold_field("S", shared_scalar=False)


def test_pipeline_plan_and_run_record_share_one_field_scheme() -> None:
    plan = PipelinePlan.from_threshold_configuration(
        [["F"], ["S"]],
        ["MF", "MS"],
        shared_scalar=True,
    )
    record = RunRecord(
        name="shared",
        alpha=-1,
        d_s1=2,
        d_s2=2,
        d_f=1,
        return_code=0,
        summary={},
        output_dir=Path("output/shared"),
        shared_scalar=True,
    )

    assert plan.field_scheme is t3_heavy_field_scheme(shared_scalar=True)
    assert record.field_scheme is plan.field_scheme
    assert record.physical_heavy_fields == ("F", "S")
    assert record.physical_scalar_fields == frozenset({"S"})


def test_ordinary_plan_uses_same_central_scalar_identity() -> None:
    plan = PipelinePlan.from_threshold_configuration(
        [["F"], ["S1", "S2"]],
        ["MF", "MS"],
    )
    assert plan.field_scheme is t3_heavy_field_scheme(shared_scalar=False)
    assert plan.initial_content.active_heavy_fields == frozenset({"F", "S1", "S2"})
