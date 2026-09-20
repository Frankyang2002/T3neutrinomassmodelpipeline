# This file tells us how heavy T3 fields are integrated across
# EFT thresholds

from __future__ import annotations

from pathlib import Path
from typing import TypeAlias

from common.RunRecords import EFTStageRecord


ThresholdField: TypeAlias = str # The heavy field name, like S, S1, F etc
ThresholdGroup: TypeAlias = tuple[ThresholdField, ...] # Set of fields integrated out, like ("S1","S2")
ThresholdPlan: TypeAlias = tuple[ThresholdGroup, ...] # The sequence of integrating out, Eg: ((F),(S1,S2))

# Convention for Field or Shared integration out
T3_HEAVY_FIELDS: tuple[str, ...] = ("F", "S1", "S2") 
SHARED_HEAVY_FIELDS: tuple[str, ...] = ("F", "S")


def heavy_fields(*, shared_scalar: bool = False) -> tuple[str, ...]:
    """Based on settings give the fields either shared or split."""
    return SHARED_HEAVY_FIELDS if shared_scalar else T3_HEAVY_FIELDS


def default_threshold_plan(*, shared_scalar: bool = False) -> ThresholdPlan:
    """We just go back to default without threshold structure, just the fields.
    This helps doing UV -> Final EFT as backup"""
    return (heavy_fields(shared_scalar=shared_scalar),)


def _threshold_aliases(*, shared_scalar: bool) -> dict[str, str]:
    """We map different names to the specific field names,
     helps with shared and split scalar differences ."""
    if shared_scalar:
        return {
            "F": "F",
            "FERMION": "F",
            "S": "S",
            "SCALAR": "S",
            "S1": "S",
            "S2": "S",
            "SCALAR1": "S",
            "SCALAR2": "S",
        }

    return {
        "F": "F",
        "FERMION": "F",
        "S1": "S1",
        "SCALAR1": "S1",
        "S2": "S2",
        "SCALAR2": "S2",
    }


def normalise_threshold_field(
    value: str,
    *,
    shared_scalar: bool = False,
) -> str:
    """We normalise different names to the actual physical fields
    We use _threshold_aliases to do so
    """
    token = value.strip().upper().replace("_", "")
    aliases = _threshold_aliases(shared_scalar=shared_scalar)
    allowed = "F or S" if shared_scalar else "F, S1, or S2"

    try:
        return aliases[token]
    except KeyError as exc:
        raise ValueError(
            f"Unknown heavy field {value!r}. Use {allowed}."
        ) from exc


def validate_threshold_plan(
    threshold_groups: list[list[str]] | None,
    *,
    shared_scalar: bool = False,
) -> ThresholdPlan:
    """Make sure our threshold plans are valid
    
    So we cant integrate out a field more than once
    """
    fields = heavy_fields(shared_scalar=shared_scalar)

    if not threshold_groups:
        return default_threshold_plan(shared_scalar=shared_scalar)

    canonical_groups: list[ThresholdGroup] = []
    flattened: list[str] = []

    # Loop through our groups and check if threshold is empty
    # or has duplicated fields across groups or within a group
    # Or if a field is not what we want, as we want fermion/scalar
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
    missing = [
        field
        for field in fields
        if field not in flattened
    ]
    unknown = sorted(set(flattened) - set(fields))

    # Error message
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


def _formal_fields_for_physical_field(
    field: str,
    *,
    shared_scalar: bool,
) -> tuple[str, ...]:
    """For shared scalar, S -> S1,S2, for split nothing really changes"""
    if shared_scalar and field == "S":
        return "S1", "S2"

    return (field,)


def threshold_plan_for_wolfram(
    plan: ThresholdPlan,
    *,
    shared_scalar: bool = False,
) -> ThresholdPlan:
    """We prepare the plan for thresholds for wolfram and fix shared scalar
    """
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
    """Return a JSON-safe representation of a physical threshold plan.
    It just converts our tuples into lists as JSON doesnt preserve Python Tuples"""
    return [list(group) for group in plan]


def threshold_plan_label(plan: ThresholdPlan) -> str:
    """We just make it readable for us, like if we have ((F),(S1,S2))
    It becomes F -> (S1,S2) for us to read well"""
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
    """
    We convert our threshold plan into EFTStageRecord objects metadata
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
