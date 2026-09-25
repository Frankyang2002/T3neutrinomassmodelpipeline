"""High-level execution plan for one T3 EFT calculation.

The pipeline should be readable as a sequence of physical operations.  This
module therefore collects the *configuration* of that sequence in one place:
which physical fields are removed at each threshold, the scale of each
threshold, and the EFT operator-dimension truncation.

No matching formula, beta function, numerical integration, report generation,
or validation calculation belongs here.  ``pipeline.py`` can consume this
object as its central description of what should happen, while the detailed
physics remains in specialised modules.
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

    ``threshold_steps`` uses physical field names.  For ordinary T3 these are
    ``F``, ``S1`` and ``S2``.  Shared-scalar mode instead uses the physical
    scalar ``S`` together with ``F``.

    The plan deliberately has no concept of ``EFT1`` or a preferred first
    threshold.  Fermion-first, scalar-first, grouped and common-threshold runs
    are represented by the same object.
    """

    threshold_steps: tuple[ThresholdStep, ...]
    shared_scalar: bool = False
    truncation: EFTTruncation = DEFAULT_EFT_TRUNCATION

    def __post_init__(self) -> None:
        if not self.threshold_steps:
            raise ValueError("A pipeline plan must contain at least one threshold.")

        if not isinstance(self.truncation, EFTTruncation):
            raise TypeError("truncation must be an EFTTruncation instance.")

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

        ``threshold_groups`` is canonicalised and validated first.  If scales
        are omitted, the existing symbolic scale convention is retained:
        ``MF``, ``MS``, ``MS1`` and ``MS2`` where applicable.
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

        This is the preferred execution view for new pipeline code.  Dispatch
        should depend on ``transition.before.active_heavy_fields`` and
        ``transition.after.active_heavy_fields`` rather than on labels such as
        ``EFT1`` or on a hard-coded fermion-first route.
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
        matching thresholds.  Each interval is identified only by its active
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
    def has_scalar_first_threshold(self) -> bool:
        """Return whether a scalar is removed before the heavy fermion ``F``.

        This is a physical ordering property, not a statement about whether a
        numerical backend exists for the resulting EFT.  It is kept explicit
        because scalar-first T3 matching can generate an intermediate
        dimension-six operator before the final Weinberg operator appears.
        """
        first_fields = frozenset(self.threshold_steps[0].fields_to_integrate)
        return (
            "F" not in first_fields
            and bool(first_fields & self.field_scheme.scalar_field_set)
        )

    @property
    def scalar_first_truncates_leading_path(self) -> bool:
        """Return whether the configured truncation drops the scalar-first d=6 path.

        Integrating out a scalar before ``F`` can first produce a dimension-six
        operator of schematic form ``L F H H S / M_S^2``.  Later thresholds can
        feed that operator into the dimension-five Weinberg coefficient.  A
        scalar-first run truncated at ``d<=5`` therefore cannot be treated as
        physically equivalent to the verified fermion-first calculation.
        """
        return self.has_scalar_first_threshold and not self.truncation.keeps(6)

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
        """Return JSON-safe run metadata without changing existing key names.

        The three historical threshold keys are retained exactly.  The new
        ``EFTTruncation`` entry makes the dimension-five approximation explicit
        for future reports and downstream tools.
        """
        return {
            "ThresholdPlan": threshold_plan_to_json(self.threshold_plan),
            "ThresholdPlanLabel": self.label,
            "ThresholdScales": list(self.threshold_scales),
            "EFTTruncation": {
                "MaxOperatorDimension": self.truncation.max_operator_dimension,
                "Label": self.truncation.label,
            },
            "ThresholdOrderingPhysics": {
                "ScalarFirst": self.has_scalar_first_threshold,
                "ScalarFirstDimensionSixPathRetained": (
                    self.has_scalar_first_threshold and self.truncation.keeps(6)
                ),
                "ScalarFirstLeadingPathTruncated": (
                    self.scalar_first_truncates_leading_path
                ),
                "AuthoritativeAtConfiguredTruncation": (
                    not self.scalar_first_truncates_leading_path
                ),
            },
        }

    def physics_warnings(self) -> tuple[str, ...]:
        """Return explicit physics limitations implied by this plan.

        Warnings are kept separate from ``description_lines`` so the historical
        threshold-plan display remains stable while the pipeline can still make
        non-authoritative approximations impossible to miss.
        """
        if self.scalar_first_truncates_leading_path:
            return (
                "Scalar-first d<=5 matching omits an intermediate dimension-six "
                "operator that can feed the leading Weinberg coefficient; this "
                "ordering is experimental/non-authoritative.",
            )
        return ()

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
