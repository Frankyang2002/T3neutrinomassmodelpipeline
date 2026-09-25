"""Physical threshold-plan handling for the T3 EFT sequence.

This module defines *ordering*, not matching physics.  Any valid ordered
partition of the physical heavy fields is accepted.  Downstream calculation
code is responsible for running to each threshold and matching the fields in
that threshold group.

The current production EFT approximation retains operators through dimension
five.  That approximation is represented explicitly by ``EFTTruncation`` in
``common.EFT`` rather than by restricting the threshold order.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence, TypeAlias

from common.EFT import ThresholdStep
from common.RunRecords import EFTStageRecord
from common.T3Fields import (
    ORDINARY_T3_FIELDS,
    SHARED_SCALAR_T3_FIELDS,
    t3_heavy_field_scheme,
)


ThresholdField: TypeAlias = str
ThresholdGroup: TypeAlias = tuple[ThresholdField, ...]
ThresholdPlan: TypeAlias = tuple[ThresholdGroup, ...]

# Compatibility constants.  The physical field identity itself is centralised
# in ``common.T3Fields`` so all threshold/RGE code uses the same mapping.
T3_HEAVY_FIELDS: tuple[str, ...] = ORDINARY_T3_FIELDS.heavy_fields
SHARED_HEAVY_FIELDS: tuple[str, ...] = SHARED_SCALAR_T3_FIELDS.heavy_fields


def heavy_fields(*, shared_scalar: bool = False) -> tuple[str, ...]:
    """Return the physical heavy fields that must be removed by the plan."""
    return t3_heavy_field_scheme(shared_scalar=shared_scalar).heavy_fields


def default_threshold_plan(*, shared_scalar: bool = False) -> ThresholdPlan:
    """Return the common-threshold plan in which all heavy fields are removed together."""
    return (heavy_fields(shared_scalar=shared_scalar),)


def _threshold_aliases(*, shared_scalar: bool) -> dict[str, str]:
    """Compatibility wrapper for accepted threshold aliases."""
    return dict(
        t3_heavy_field_scheme(shared_scalar=shared_scalar).threshold_aliases
    )


def normalise_threshold_field(
    value: str,
    *,
    shared_scalar: bool = False,
) -> str:
    """Convert one user-facing field token to its physical pipeline name."""
    return t3_heavy_field_scheme(
        shared_scalar=shared_scalar
    ).normalise_threshold_field(value)


def validate_threshold_plan(
    threshold_groups: Sequence[Sequence[str]] | None,
    *,
    shared_scalar: bool = False,
) -> ThresholdPlan:
    """Validate and canonicalise an arbitrary physical threshold ordering.

    A threshold plan is an ordered partition of the physical heavy fields.
    Every heavy field must appear exactly once, but there is deliberately no
    restriction on which field is integrated out first.
    """
    fields = heavy_fields(shared_scalar=shared_scalar)

    if not threshold_groups:
        return default_threshold_plan(shared_scalar=shared_scalar)

    canonical_groups: list[ThresholdGroup] = []
    flattened: list[str] = []

    for raw_group in threshold_groups:
        if not raw_group:
            raise ValueError("A threshold group cannot be empty.")

        group = tuple(
            normalise_threshold_field(
                field,
                shared_scalar=shared_scalar,
            )
            for field in raw_group
        )

        if len(set(group)) != len(group):
            raise ValueError(
                "A field cannot appear twice inside one threshold group: "
                + " ".join(group)
            )

        canonical_groups.append(group)
        flattened.extend(group)

    duplicates = sorted(
        field
        for field in fields
        if flattened.count(field) > 1
    )
    missing = [field for field in fields if field not in flattened]
    unknown = sorted(set(flattened) - set(fields))

    errors: list[str] = []
    if duplicates:
        errors.append("repeated field(s): " + ", ".join(duplicates))
    if missing:
        errors.append("missing field(s): " + ", ".join(missing))
    if unknown:
        errors.append("unknown field(s): " + ", ".join(unknown))

    if errors:
        raise ValueError(
            "Invalid threshold plan; "
            + "; ".join(errors)
            + ". Every one of "
            + ", ".join(fields)
            + " must be integrated out exactly once."
        )

    return tuple(canonical_groups)


def default_threshold_scale(group: Sequence[str]) -> str:
    """Return the established symbolic matching scale for one field group.

    This naming convention already exists in ``pipeline.py``.  Keeping it next
    to threshold-plan handling makes scale resolution reusable without making
    the generic ``ThresholdStep`` object T3-specific.
    """
    fields = frozenset(group)

    if fields == {"F"}:
        return "MF"
    if fields == {"S"}:
        return "MS"
    if fields == {"S1", "S2"}:
        return "MS"
    if fields == {"S1"}:
        return "MS1"
    if fields == {"S2"}:
        return "MS2"

    # This fallback is mainly useful for grouped thresholds such as (F,S1).
    # Preserve the user-specified group order in the readable symbolic name.
    return "M_" + "_".join(group)


def resolve_threshold_scales(
    plan: ThresholdPlan,
    supplied_scales: Sequence[str | float] | None = None,
) -> tuple[str | float, ...]:
    """Return one matching scale for every threshold in ``plan``.

    Explicit scales are preserved exactly.  When they are omitted, use the
    existing project convention implemented by :func:`default_threshold_scale`.
    """
    if supplied_scales is None:
        return tuple(default_threshold_scale(group) for group in plan)

    if len(supplied_scales) != len(plan):
        raise ValueError(
            "Threshold scales must be supplied once for each threshold group."
        )

    return tuple(supplied_scales)


def build_threshold_steps(
    plan: ThresholdPlan,
    scales: Sequence[str | float],
) -> tuple[ThresholdStep, ...]:
    """Combine a validated threshold plan with its matching scales.

    The returned objects are generic threshold transitions and carry no
    fermion-first or scalar-first assumption.
    """
    if len(plan) != len(scales):
        raise ValueError(
            "Threshold scales must be supplied once for each threshold group."
        )

    return tuple(
        ThresholdStep(
            fields_to_integrate=group,
            scale=scale,
        )
        for group, scale in zip(plan, scales, strict=True)
    )


def _formal_fields_for_physical_field(
    field: str,
    *,
    shared_scalar: bool,
) -> tuple[str, ...]:
    """Map a physical field to the formal T3 roles used by Wolfram matching."""
    return t3_heavy_field_scheme(
        shared_scalar=shared_scalar
    ).formal_roles_for(field)


def threshold_plan_for_wolfram(
    plan: ThresholdPlan,
    *,
    shared_scalar: bool = False,
) -> ThresholdPlan:
    """Convert a physical threshold plan to the formal T3 matching roles."""
    if not shared_scalar:
        return plan

    formal_groups: list[ThresholdGroup] = []

    for group in plan:
        expanded: list[str] = []
        for field in group:
            expanded.extend(
                _formal_fields_for_physical_field(
                    field,
                    shared_scalar=True,
                )
            )
        formal_groups.append(tuple(expanded))

    return tuple(formal_groups)


def threshold_plan_to_json(plan: ThresholdPlan) -> list[list[str]]:
    """Return the physical threshold plan in JSON-safe form."""
    return [list(group) for group in plan]


def threshold_plan_label(plan: ThresholdPlan) -> str:
    """Return a compact human-readable threshold ordering."""
    labels: list[str] = []

    for group in plan:
        if len(group) == 1:
            labels.append(group[0])
        else:
            labels.append("(" + ",".join(group) + ")")

    return " -> ".join(labels)


def build_eft_stage_records(
    plan: ThresholdPlan,
    model_output_dir: Path | None = None,
    *,
    shared_scalar: bool = False,
) -> list[EFTStageRecord]:
    """Build stage metadata for every theory level in ``plan``.

    Stage labels intentionally retain the existing ``EFT_<n>_after_<fields>``
    format because report paths and historical output currently depend on it.
    New calculation logic should identify a stage from its field content rather
    than from this compatibility label.
    """
    stages: list[EFTStageRecord] = []
    active = list(heavy_fields(shared_scalar=shared_scalar))

    uv_dir = (
        model_output_dir / "stages" / "uv"
        if model_output_dir is not None
        else None
    )

    stages.append(
        EFTStageRecord(
            level=0,
            integrated_fields=(),
            active_heavy_fields=tuple(active),
            label="UV",
            output_dir=uv_dir,
        )
    )

    for level, group in enumerate(plan, start=1):
        for field in group:
            active.remove(field)

        label = f"EFT_{level}_after_{'_'.join(group)}"
        stage_dir = (
            model_output_dir / "stages" / label
            if model_output_dir is not None
            else None
        )

        stages.append(
            EFTStageRecord(
                level=level,
                integrated_fields=group,
                active_heavy_fields=tuple(active),
                label=label,
                output_dir=stage_dir,
            )
        )

    return stages
