"""High-level execution plan for one T3 EFT calculation.

The pipeline should be readable as a sequence of physical operations. This
module therefore collects the configuration of that sequence in one place:
which physical fields are removed at each threshold, the scale of each
threshold, and the EFT operator-dimension truncation.

No matching formula, beta function, numerical integration, report generation,
or validation calculation belongs here. ``pipeline.py`` consumes this object
as its central description of what should happen, while detailed physics
remains in specialised modules.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from common.EFT import (
    DEFAULT_EFT_TRUNCATION,
    EFTContent,
    EFTRunningInterval,
    EFTTransition,
    EFTTruncation,
    ThresholdStep,
)
from common.RunRecords import EFTStageRecord
from common.T3Fields import T3HeavyFieldScheme, t3_heavy_field_scheme
from common.Thresholds import (
    ThresholdPlan,
    build_eft_stage_records,
    build_threshold_steps,
    resolve_threshold_scales,
    threshold_plan_label,
    threshold_plan_to_json,
    validate_threshold_plan,
)


@dataclass(frozen=True, slots=True)
class PipelinePlan:
    """Physical threshold sequence and EFT truncation for one pipeline run.

    ``threshold_steps`` uses physical field names. For ordinary T3 these are
    ``F``, ``S1`` and ``S2``. Shared-scalar mode instead uses the physical
    scalar ``S`` together with ``F``.

    The data model remains generic enough to represent an arbitrary ordered
    partition of the heavy fields. Production execution is intentionally
    narrower: the current project supports either one common threshold or the
    verified fermion-first hierarchy in which ``F`` is removed first and the
    complete physical scalar sector is removed together at the next threshold.
    """

    threshold_steps: tuple[ThresholdStep, ...]
    shared_scalar: bool = False
    truncation: EFTTruncation = DEFAULT_EFT_TRUNCATION

    def __post_init__(self) -> None:
        if not self.threshold_steps:
            raise ValueError("A pipeline plan must contain at least one threshold.")

        if not isinstance(self.truncation, EFTTruncation):
            raise TypeError("truncation must be an EFTTruncation instance.")

        if self.truncation != DEFAULT_EFT_TRUNCATION:
            raise ValueError(
                "The production T3 pipeline currently supports only the "
                f"{DEFAULT_EFT_TRUNCATION.label} EFT truncation."
            )

        # Reuse the canonical threshold validator so direct construction of a
        # PipelinePlan cannot bypass the same field-completeness rules used by
        # the CLI configuration path.
        validate_threshold_plan(
            [list(step.fields_to_integrate) for step in self.threshold_steps],
            shared_scalar=self.shared_scalar,
        )

    @classmethod
    def from_threshold_configuration(
        cls,
        threshold_groups: Sequence[Sequence[str]] | None,
        threshold_scales: Sequence[str | float] | None = None,
        *,
        shared_scalar: bool = False,
        truncation: EFTTruncation = DEFAULT_EFT_TRUNCATION,
    ) -> "PipelinePlan":
        """Build a complete plan from user-facing threshold configuration.

        ``threshold_groups`` is canonicalised and validated first. If scales
        are omitted, the established symbolic scale convention is retained:
        ``MF``, ``MS``, ``MS1`` and ``MS2`` where applicable.

        The ``truncation`` argument remains explicit so the approximation is
        visible in architecture/tests, but production currently accepts only
        the default dimension-five truncation.
        """
        plan = validate_threshold_plan(
            threshold_groups,
            shared_scalar=shared_scalar,
        )
        scales = resolve_threshold_scales(plan, threshold_scales)
        steps = build_threshold_steps(plan, scales)

        return cls(
            threshold_steps=steps,
            shared_scalar=shared_scalar,
            truncation=truncation,
        )

    @property
    def field_scheme(self) -> T3HeavyFieldScheme:
        """Return the physical heavy-field naming scheme for this plan."""
        return t3_heavy_field_scheme(shared_scalar=self.shared_scalar)

    @property
    def threshold_plan(self) -> ThresholdPlan:
        """Return the compatibility tuple representation of the threshold plan."""
        return tuple(step.fields_to_integrate for step in self.threshold_steps)

    @property
    def threshold_scales(self) -> tuple[str | float, ...]:
        """Return threshold scales in execution order."""
        return tuple(step.scale for step in self.threshold_steps)

    @property
    def label(self) -> str:
        """Return the existing human-readable threshold-order label."""
        return threshold_plan_label(self.threshold_plan)

    @property
    def initial_content(self) -> EFTContent:
        """Return the UV heavy-field content before the first threshold."""
        return EFTContent.from_fields(
            self.field_scheme.heavy_fields,
            truncation=self.truncation,
        )

    @property
    def transitions(self) -> tuple[EFTTransition, ...]:
        """Return every physical threshold as a content-aware transition.

        Dispatch should depend on physical field content rather than ordinal
        stage labels such as ``EFT1`` or ``EFT2``.
        """
        transitions: list[EFTTransition] = []
        current = self.initial_content

        for index, step in enumerate(self.threshold_steps, start=1):
            next_content = current.after_integrating(step.fields_to_integrate)
            transitions.append(
                EFTTransition(
                    index=index,
                    step=step,
                    before=current,
                    after=next_content,
                )
            )
            current = next_content

        return tuple(transitions)

    @property
    def running_intervals(self) -> tuple[EFTRunningInterval, ...]:
        """Return the heavy-field EFT regions between adjacent thresholds.

        A common-threshold calculation has no intermediate running interval.
        Sequential plans produce one interval for every gap between adjacent
        matching thresholds. Each interval is identified only by its active
        field content and matching scales.
        """
        transitions = self.transitions
        intervals: list[EFTRunningInterval] = []

        for entered_by, exited_by in zip(
            transitions[:-1],
            transitions[1:],
            strict=True,
        ):
            intervals.append(
                EFTRunningInterval(
                    index=entered_by.index,
                    content=entered_by.after,
                    high_scale=entered_by.scale,
                    low_scale=exited_by.scale,
                    entered_by=entered_by,
                    exited_by=exited_by,
                )
            )

        return tuple(intervals)

    @property
    def is_common_threshold(self) -> bool:
        """Return whether all physical heavy fields are removed together."""
        return len(self.threshold_steps) == 1

    @property
    def is_verified_fermion_first_hierarchy(self) -> bool:
        """Return whether this is the implemented hierarchical production path.

        Ordinary T3 requires ``F -> (S1,S2)``. Shared-scalar mode requires
        ``F -> S``. The remaining physical scalar sector must be removed
        together because that is the field content implemented by the current
        intermediate-EFT backend.
        """
        if len(self.threshold_steps) != 2:
            return False

        first, second = self.threshold_steps
        return bool(
            frozenset(first.fields_to_integrate) == frozenset({"F"})
            and frozenset(second.fields_to_integrate)
            == self.field_scheme.scalar_field_set
        )

    @property
    def is_supported_production_order(self) -> bool:
        """Return whether the threshold ordering has a verified production path."""
        return self.is_common_threshold or self.is_verified_fermion_first_hierarchy

    def production_scope_error(self) -> str:
        """Explain the currently supported threshold-order scope."""
        scalar_label = "S" if self.shared_scalar else "(S1,S2)"
        return (
            "Current production threshold support is limited to one common "
            f"threshold or the verified fermion-first hierarchy F -> {scalar_label}. "
            "Scalar-first and partially split scalar hierarchies are outside "
            "the production scope used for this project; scalar-first matching "
            "can require leading intermediate operators above dimension five."
        )

    @property
    def final_content(self) -> EFTContent:
        """Return the EFT content after all requested thresholds."""
        return self.transitions[-1].after

    def build_stage_records(
        self,
        model_output_dir: Path | None = None,
    ) -> list[EFTStageRecord]:
        """Build runtime metadata for the UV theory and every matched EFT.

        Historical report labels are preserved by ``build_eft_stage_records``;
        calculation code should use the recorded active fields instead of the
        ordinal stage number.
        """
        return build_eft_stage_records(
            self.threshold_plan,
            model_output_dir,
            shared_scalar=self.shared_scalar,
        )

    def summary_metadata(self) -> dict[str, object]:
        """Return JSON-safe run metadata without changing historical keys.

        The established threshold keys are retained. ``EFTTruncation`` keeps
        the dimension-five approximation explicit rather than hiding it in the
        matching implementation.
        """
        return {
            "ThresholdPlan": threshold_plan_to_json(self.threshold_plan),
            "ThresholdPlanLabel": self.label,
            "ThresholdScales": list(self.threshold_scales),
            "EFTTruncation": {
                "MaxOperatorDimension": self.truncation.max_operator_dimension,
                "Label": self.truncation.label,
            },
        }

    def description_lines(self) -> list[str]:
        """Return a compact, physics-facing description of the EFT sequence."""
        lines = [
            f"EFT truncation: {self.truncation.label}",
            f"Threshold plan: {self.label}",
            "  UV active heavy fields: " + self.initial_content.heavy_field_label,
        ]

        for transition in self.transitions:
            integrated = ", ".join(transition.fields_to_integrate)
            lines.append(
                f"  [{transition.index}] integrate {integrated} at "
                f"{transition.scale} -> active heavy fields: "
                f"{transition.after.heavy_field_label}"
            )

        return lines
