from __future__ import annotations

from pathlib import Path

from common.Records import EFTStageRecord


T3_HEAVY_FIELDS: tuple[str, ...] = ("F", "S1", "S2")


def normalise_threshold_field(value: str) -> str:
    """Normalise CLI aliases to the canonical T3 heavy-field labels."""

    token = value.strip().upper().replace("_", "")

    aliases = {
        "F": "F",
        "FERMION": "F",
        "S1": "S1",
        "SCALAR1": "S1",
        "S2": "S2",
        "SCALAR2": "S2",
    }

    try:
        return aliases[token]
    except KeyError as exc:
        raise ValueError(
            f"Unknown heavy field {value!r}. "
            "Use F, S1, or S2."
        ) from exc


def validate_threshold_plan(
    threshold_groups: list[list[str]] | None,
) -> tuple[tuple[str, ...], ...]:
    """Validate and canonicalise an ordered decoupling plan.

    Each heavy T3 field must occur exactly once.  Fields in the same inner
    tuple are integrated out together at one threshold.

    No explicit plan keeps the historical common-threshold behaviour:
        (("F", "S1", "S2"),)
    """

    if not threshold_groups:
        return (T3_HEAVY_FIELDS,)

    canonical_groups: list[tuple[str, ...]] = []
    flattened: list[str] = []

    for raw_group in threshold_groups:
        if not raw_group:
            raise ValueError("A threshold group cannot be empty.")

        group = tuple(normalise_threshold_field(field) for field in raw_group)

        if len(set(group)) != len(group):
            raise ValueError(
                "A field cannot appear twice inside one threshold group: "
                + " ".join(group)
            )

        canonical_groups.append(group)
        flattened.extend(group)

    duplicates = sorted(
        field
        for field in T3_HEAVY_FIELDS
        if flattened.count(field) > 1
    )
    missing = [
        field
        for field in T3_HEAVY_FIELDS
        if field not in flattened
    ]
    unknown = sorted(set(flattened) - set(T3_HEAVY_FIELDS))

    errors: list[str] = []

    if duplicates:
        errors.append("repeated field(s): " + ", ".join(duplicates))
    if missing:
        errors.append("missing field(s): " + ", ".join(missing))
    if unknown:
        errors.append("unknown field(s): " + ", ".join(unknown))

    if errors:
        raise ValueError(
            "Invalid threshold plan; " + "; ".join(errors) + ". "
            "Every one of F, S1, S2 must be integrated out exactly once."
        )

    return tuple(canonical_groups)


def threshold_plan_to_json(
    plan: tuple[tuple[str, ...], ...],
) -> list[list[str]]:
    """Convert an immutable plan into a JSON-safe representation."""

    return [list(group) for group in plan]


def threshold_plan_label(
    plan: tuple[tuple[str, ...], ...],
) -> str:
    """Compact human-readable representation such as F -> (S1,S2)."""

    labels: list[str] = []

    for group in plan:
        if len(group) == 1:
            labels.append(group[0])
        else:
            labels.append("(" + ",".join(group) + ")")

    return " -> ".join(labels)


def build_eft_stage_records(
    plan: tuple[tuple[str, ...], ...],
    model_output_dir: Path | None = None,
) -> list[EFTStageRecord]:
    """Create the UV + EFT stage metadata implied by a threshold plan."""

    stages: list[EFTStageRecord] = []

    active = list(T3_HEAVY_FIELDS)

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

        group_slug = "_".join(group)
        label = f"EFT_{level}_after_{group_slug}"

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
