"""Extract the exact scalar-first ``psi^2 phi^3`` Matchete boundary terms.

When one T3 scalar is integrated out before the heavy fermion ``F``, tree-level
matching can generate a dimension-six operator with field content schematically
``L F H H S``.  This module does not guess its Clebsch-Gordan normalization or
field conjugations.  Instead it reads the exact Matchete ``InputForm`` tree EFT
already exported by ``Lagrangian/RunMatching.wl`` and records candidate terms
verbatim for the later tensor adapter.

The output is intentionally a boundary/audit artifact, not a production RGE
result.  Scalar-first running remains unsupported until those exact terms are
mapped into the general ``psi^2 phi^3`` tensor basis and all required one-loop
mixing is verified.
"""

from __future__ import annotations

import json
from pathlib import Path

from common.EFT import EFTRunningInterval
from common.RunRecords import EFTStageRecord, RunRecord
from RGE.running.intermediate.IntermediateMatcheteParsing import (
    matching_bracket,
    split_top_level,
)


OPERATOR_CLASS = "psi2phi3"
OPERATOR_DIMENSION = 6


def _split_top_level_sum(expression: str) -> list[str]:
    """Split a Mathematica ``InputForm`` sum without touching nested signs.

    ``InputForm`` normally prints algebraic sums infix, but the functional
    ``Plus[...]`` form is accepted as well so the extractor is robust to the
    exact serializer used by a particular Mathematica version.
    """

    text = expression.strip()
    if not text:
        return []

    if text.startswith("Plus["):
        open_index = text.index("[")
        if matching_bracket(text, open_index) == len(text) - 1:
            return [
                item
                for item in split_top_level(text[open_index + 1 : -1])
                if item
            ]

    terms: list[str] = []
    start = 0
    square = curly = paren = 0

    for pos, char in enumerate(text):
        if char == "[":
            square += 1
        elif char == "]":
            square -= 1
        elif char == "{":
            curly += 1
        elif char == "}":
            curly -= 1
        elif char == "(":
            paren += 1
        elif char == ")":
            paren -= 1
        elif (
            pos > start
            and char in "+-"
            and square == 0
            and curly == 0
            and paren == 0
        ):
            term = text[start:pos].strip()
            if term:
                terms.append(term)
            start = pos

    final = text[start:].strip()
    if final:
        terms.append(final)
    return terms


def _field_count(term: str, field_name: str, field_type: str) -> int:
    """Count serialized Matchete fields by their stable InputForm prefix."""

    return term.count(f"Field[{field_name}, {field_type}")


def _candidate_metadata(
    term: str,
    *,
    integrated_scalar: str,
    remaining_scalar: str,
) -> dict[str, object] | None:
    """Return structural metadata if ``term`` has the scalar-first d=6 fields."""

    formal_scalar_name = {
        "S1": "NewScalar1",
        "S2": "NewScalar2",
    }
    integrated_formal = formal_scalar_name[integrated_scalar]
    remaining_formal = formal_scalar_name[remaining_scalar]

    counts = {
        "LeptonDoublet": _field_count(term, "l", "Fermion"),
        "HeavyFermion": _field_count(term, "NewFermion", "Fermion"),
        "Higgs": _field_count(term, "H", "Scalar"),
        "RemainingScalar": _field_count(term, remaining_formal, "Scalar"),
        "IntegratedScalar": _field_count(term, integrated_formal, "Scalar"),
    }
    total_fermions = term.count(", Fermion,")
    total_scalars = term.count(", Scalar,")

    is_candidate = (
        counts["LeptonDoublet"] == 1
        and counts["HeavyFermion"] == 1
        and counts["Higgs"] == 2
        and counts["RemainingScalar"] == 1
        and counts["IntegratedScalar"] == 0
        and total_fermions == 2
        and total_scalars == 3
    )
    if not is_candidate:
        return None

    return {
        "FieldCounts": counts,
        "TotalFermionFields": total_fermions,
        "TotalScalarFields": total_scalars,
        "ContainsPLProjector": "Proj[-1]" in term,
        "ContainsPRProjector": "Proj[1]" in term,
        "TermInputForm": term,
    }


def supports_scalar_first_dimension_six_seed(
    record: RunRecord,
    stage: EFTStageRecord,
    interval: EFTRunningInterval,
) -> bool:
    """Return whether this interval is the ordinary one-scalar-first T3 EFT."""

    if record.shared_scalar or not interval.content.truncation.keeps(6):
        return False

    integrated_scalars = stage.integrated_field_set & record.physical_scalar_fields
    active_scalars = interval.active_heavy_fields & record.physical_scalar_fields

    return bool(
        len(integrated_scalars) == 1
        and len(active_scalars) == 1
        and "F" in interval.active_heavy_fields
        and stage.integrated_field_set == integrated_scalars
    )


def extract_scalar_first_dimension_six_seed(
    record: RunRecord,
    stage: EFTStageRecord,
    interval: EFTRunningInterval,
) -> dict[str, object]:
    """Extract exact candidate ``L F H H S`` terms from the matched tree EFT."""

    if not supports_scalar_first_dimension_six_seed(record, stage, interval):
        raise ValueError(
            "Scalar-first dimension-six seed extraction was requested for an "
            "incompatible EFT interval."
        )

    integrated_scalar = next(iter(stage.integrated_field_set))
    remaining_scalar = next(
        iter(interval.active_heavy_fields & record.physical_scalar_fields)
    )

    prefix = (
        f"eft{stage.level}_after_"
        + "_".join(interval.entered_by.fields_to_integrate)
    )
    data_dir = record.output_dir / "data"
    tree_path = data_dir / f"{prefix}_tree_inputform.txt"
    output_path = data_dir / f"{prefix}_psi2phi3_seed.json"

    if not tree_path.is_file():
        payload: dict[str, object] = {
            "status": "MissingTreeEFT",
            "operator_class": OPERATOR_CLASS,
            "operator_dimension": OPERATOR_DIMENSION,
            "source_tree_file": tree_path.name,
            "candidate_term_count": 0,
            "candidates": [],
        }
    else:
        expression = tree_path.read_text(encoding="utf-8")
        candidates: list[dict[str, object]] = []

        for index, term in enumerate(_split_top_level_sum(expression), start=1):
            metadata = _candidate_metadata(
                term,
                integrated_scalar=integrated_scalar,
                remaining_scalar=remaining_scalar,
            )
            if metadata is None:
                continue
            candidates.append({"Index": index, **metadata})

        payload = {
            "status": "Success" if candidates else "NoCandidateFound",
            "operator_class": OPERATOR_CLASS,
            "operator_dimension": OPERATOR_DIMENSION,
            "coefficient_symmetry_target": "C6_(ij)(abc)",
            "integrated_scalar": integrated_scalar,
            "remaining_scalar": remaining_scalar,
            "active_heavy_fields": sorted(interval.active_heavy_fields),
            "source_tree_file": tree_path.name,
            "candidate_term_count": len(candidates),
            "candidates": candidates,
            "normalization_status": "RawMatcheteTermsOnly",
            "tensor_adapter_status": "NotImplemented",
            "production_rge_status": "NotImplemented",
            "note": (
                "Terms are copied verbatim from the Matchete tree EFT. No "
                "Clebsch-Gordan, conjugation, or coefficient normalization "
                "has been inferred by this extractor."
            ),
        }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    stage.summary["ScalarFirstDimensionSixSeedStatus"] = payload["status"]
    stage.summary["ScalarFirstDimensionSixSeedFile"] = (
        output_path.relative_to(record.output_dir).as_posix()
    )
    stage.summary["ScalarFirstDimensionSixCandidateCount"] = payload.get(
        "candidate_term_count",
        0,
    )
    return payload
