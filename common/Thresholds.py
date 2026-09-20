from __future__ import annotations

from pathlib import Path
from typing import TypeAlias

from common.RunRecords import EFTStageRecord


ThresholdField: TypeAlias = str
ThresholdGroup: TypeAlias = tuple[ThresholdField, ...]
ThresholdPlan: TypeAlias = tuple[ThresholdGroup, ...]

T3_HEAVY_FIELDS: tuple[str, ...] = ("F", "S1", "S2")
SHARED_HEAVY_FIELDS: tuple[str, ...] = ("F", "S")


def heavy_fields(*, shared_scalar: bool = False) -> tuple[str, ...]:
    """Return the physical heavy fields that may appear in a threshold plan."""
    return SHARED_HEAVY_FIELDS if shared_scalar else T3_HEAVY_FIELDS


def default_threshold_plan(*, shared_scalar: bool = False) -> ThresholdPlan:
    """Return the historical common-threshold plan for one physical model."""
    return (heavy_fields(shared_scalar=shared_scalar),)


def _threshold_aliases(*, shared_scalar: bool) -> dict[str, str]:
    """Return CLI aliases mapped to canonical physical heavy-field names."""
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
    """Normalise one CLI field name to a canonical physical heavy-field name.

    Ordinary T3 uses the physical fields F, S1 and S2.

    Shared-scalar mode uses the physical fields F and S.  In that mode S1/S2
    are accepted only as input aliases and both normalise to the single
    physical scalar S.
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
    """Validate and canonicalise an ordered physical decoupling plan.

    Every physical heavy field must occur exactly once across the plan.
    Fields in the same group are integrated out at the same threshold.

    If no explicit plan is supplied, all physical heavy fields are integrated
    out together, preserving the historical common-threshold behaviour.
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
    missing = [
        field
        for field in fields
        if field not in flattened
    ]
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


def _formal_fields_for_physical_field(
    field: str,
    *,
    shared_scalar: bool,
) -> tuple[str, ...]:
    """Map one physical threshold field to its formal T3 matching role(s)."""
    if shared_scalar and field == "S":
        return "S1", "S2"

    return (field,)


def threshold_plan_for_wolfram(
    plan: ThresholdPlan,
    *,
    shared_scalar: bool = False,
) -> ThresholdPlan:
    """Convert a physical threshold plan to formal T3 matching field names.

    Ordinary T3 is unchanged.

    In shared-scalar mode the one physical scalar S occupies both formal
    matching roles, so a physical S threshold is exported to Wolfram as the
    simultaneous formal threshold (S1, S2).
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
    """Return a JSON-safe representation of a physical threshold plan."""
    return [list(group) for group in plan]


def threshold_plan_label(plan: ThresholdPlan) -> str:
    """Return a compact human-readable label for a physical threshold plan."""
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
    """Construct UV/EFT metadata using physical heavy-field content.

    The resulting ``integrated_fields`` and ``active_heavy_fields`` always
    refer to physical fields.  Therefore shared-scalar records contain S, not
    the formal matching roles S1 and S2.

    Level 0 is the UV theory.  Each subsequent level is the theory after one
    threshold group from ``plan`` has been integrated out.
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
