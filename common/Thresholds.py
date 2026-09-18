from __future__ import annotations

from pathlib import Path
from common.Records import EFTStageRecord

T3_HEAVY_FIELDS = ("F", "S1", "S2")
SHARED_HEAVY_FIELDS = ("F", "S")


def heavy_fields(*, shared_scalar: bool = False) -> tuple[str, ...]:
    return SHARED_HEAVY_FIELDS if shared_scalar else T3_HEAVY_FIELDS


def normalise_threshold_field(value: str, *, shared_scalar: bool = False) -> str:
    token = value.strip().upper().replace("_", "")
    if shared_scalar:
        aliases = {
            "F": "F", "FERMION": "F", "S": "S", "SCALAR": "S",
            "S1": "S", "S2": "S", "SCALAR1": "S", "SCALAR2": "S",
        }
        allowed = "F or S"
    else:
        aliases = {
            "F": "F", "FERMION": "F", "S1": "S1", "SCALAR1": "S1",
            "S2": "S2", "SCALAR2": "S2",
        }
        allowed = "F, S1, or S2"
    try:
        return aliases[token]
    except KeyError as exc:
        raise ValueError(f"Unknown heavy field {value!r}. Use {allowed}.") from exc


def validate_threshold_plan(threshold_groups, *, shared_scalar: bool = False):
    fields = heavy_fields(shared_scalar=shared_scalar)
    if not threshold_groups:
        return (fields,)

    canonical_groups = []
    flattened = []
    for raw_group in threshold_groups:
        if not raw_group:
            raise ValueError("A threshold group cannot be empty.")
        group = tuple(
            normalise_threshold_field(x, shared_scalar=shared_scalar)
            for x in raw_group
        )
        if len(set(group)) != len(group):
            raise ValueError("A field cannot appear twice inside one threshold group: " + " ".join(group))
        canonical_groups.append(group)
        flattened.extend(group)

    duplicates = sorted(f for f in fields if flattened.count(f) > 1)
    missing = [f for f in fields if f not in flattened]
    unknown = sorted(set(flattened) - set(fields))
    errors = []
    if duplicates:
        errors.append("repeated field(s): " + ", ".join(duplicates))
    if missing:
        errors.append("missing field(s): " + ", ".join(missing))
    if unknown:
        errors.append("unknown field(s): " + ", ".join(unknown))
    if errors:
        raise ValueError(
            "Invalid threshold plan; " + "; ".join(errors) + ". Every one of "
            + ", ".join(fields) + " must be integrated out exactly once."
        )
    return tuple(canonical_groups)


def threshold_plan_for_wolfram(plan, *, shared_scalar: bool = False):
    if not shared_scalar:
        return plan
    out = []
    for group in plan:
        expanded = []
        for field in group:
            expanded.extend(("S1", "S2") if field == "S" else (field,))
        out.append(tuple(expanded))
    return tuple(out)


def threshold_plan_to_json(plan):
    return [list(group) for group in plan]


def threshold_plan_label(plan):
    labels = []
    for group in plan:
        labels.append(group[0] if len(group) == 1 else "(" + ",".join(group) + ")")
    return " -> ".join(labels)


def build_eft_stage_records(plan, model_output_dir: Path | None = None, *, shared_scalar: bool = False):
    stages = []
    active = list(heavy_fields(shared_scalar=shared_scalar))
    uv_dir = model_output_dir / "stages" / "uv" if model_output_dir else None
    stages.append(EFTStageRecord(0, (), tuple(active), "UV", uv_dir))
    for level, group in enumerate(plan, start=1):
        for field in group:
            active.remove(field)
        label = f"EFT_{level}_after_{'_'.join(group)}"
        stage_dir = model_output_dir / "stages" / label if model_output_dir else None
        stages.append(EFTStageRecord(level, group, tuple(active), label, stage_dir))
    return stages
