"""Runtime records passed through the Python pipeline.

The objects in this module contain metadata only.  They describe what happened
at each EFT level and where outputs were written; they do not perform matching,
RGE running, or validation calculations.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from common.T3Fields import T3HeavyFieldScheme, t3_heavy_field_scheme


@dataclass
class EFTStageRecord:
    """Metadata for one physical theory level in a threshold sequence.

    ``level`` is an ordering index used for report compatibility.  Physics code
    should prefer ``active_heavy_fields`` and ``integrated_fields`` when it
    needs to identify a stage, because those fields remain meaningful for any
    threshold ordering.

    Field names are physical names.  Ordinary T3 uses ``F``, ``S1`` and ``S2``;
    shared-scalar mode uses ``F`` and ``S``.
    """

    level: int
    integrated_fields: tuple[str, ...]
    active_heavy_fields: tuple[str, ...]
    label: str
    output_dir: Path | None = None
    summary: dict = field(default_factory=dict)

    @property
    def active_field_set(self) -> frozenset[str]:
        """Return active heavy fields as an order-independent set."""
        return frozenset(self.active_heavy_fields)

    @property
    def integrated_field_set(self) -> frozenset[str]:
        """Return fields removed at this threshold as an order-independent set."""
        return frozenset(self.integrated_fields)

    @property
    def is_uv(self) -> bool:
        """Return whether this record describes the UV theory."""
        return self.level == 0 and not self.integrated_fields

    def has_active_fields(self, *fields: str) -> bool:
        """Return whether the stage has exactly the requested active fields."""
        return self.active_field_set == frozenset(fields)

    def integrates(self, *fields: str) -> bool:
        """Return whether this threshold removes exactly the requested fields."""
        return self.integrated_field_set == frozenset(fields)


@dataclass
class RunRecord:
    """Bookkeeping for one complete T3 pipeline run.

    ``shared_scalar`` records whether the formal topology roles ``S1`` and
    ``S2`` represent one physical scalar ``S``.  The stage records themselves
    always use physical field names.
    """

    name: str
    alpha: int
    d_s1: int
    d_s2: int
    d_f: int
    return_code: int
    summary: dict
    output_dir: Path
    eft_stages: list[EFTStageRecord] = field(default_factory=list)
    shared_scalar: bool = False

    @property
    def field_scheme(self) -> T3HeavyFieldScheme:
        """Return the physical heavy-field naming scheme for this run."""
        return t3_heavy_field_scheme(shared_scalar=self.shared_scalar)

    @property
    def physical_heavy_fields(self) -> tuple[str, ...]:
        """Return all physical heavy-field names for this run."""
        return self.field_scheme.heavy_fields

    @property
    def physical_scalar_fields(self) -> frozenset[str]:
        """Return the physical T3 scalar sector for this run."""
        return self.field_scheme.scalar_field_set

    @property
    def physical_dimensions(self) -> tuple[int, int, int] | tuple[int, int]:
        """Return the physical heavy-field dimensions for this run."""
        # Imported lazily to keep this metadata module independent of the T3
        # representation implementation at import time.
        from common.T3Model import physical_dimensions_from_formal

        return physical_dimensions_from_formal(
            self.d_s1,
            self.d_s2,
            self.d_f,
            shared_scalar=self.shared_scalar,
        )

    @property
    def d_s(self) -> int | None:
        """Return the shared-scalar dimension, or ``None`` for ordinary T3."""
        if not self.shared_scalar:
            return None

        d_s, _ = self.physical_dimensions
        return d_s

    def stage(self, level: int) -> EFTStageRecord:
        """Return the unique stage with the requested compatibility level."""
        matches = [stage for stage in self.eft_stages if stage.level == level]

        if len(matches) == 1:
            return matches[0]
        if not matches:
            raise KeyError(f"{self.name} has no EFT stage with level {level}.")

        raise ValueError(
            f"{self.name} has multiple EFT stages with level {level}."
        )

    def stage_with_active_fields(self, *fields: str) -> EFTStageRecord:
        """Return the unique stage whose active heavy fields match ``fields``."""
        requested = frozenset(fields)
        matches = [
            stage
            for stage in self.eft_stages
            if stage.active_field_set == requested
        ]

        if len(matches) == 1:
            return matches[0]
        if not matches:
            field_label = ", ".join(fields) if fields else "none"
            raise KeyError(
                f"{self.name} has no EFT stage with active fields: {field_label}."
            )

        raise ValueError(
            f"{self.name} has multiple EFT stages with active fields {sorted(requested)}."
        )

    def stage_after_integrating(self, *fields: str) -> EFTStageRecord:
        """Return the unique stage produced by integrating exactly ``fields``."""
        requested = frozenset(fields)
        matches = [
            stage
            for stage in self.eft_stages
            if stage.integrated_field_set == requested
        ]

        if len(matches) == 1:
            return matches[0]
        if not matches:
            field_label = ", ".join(fields) if fields else "none"
            raise KeyError(
                f"{self.name} has no EFT stage after integrating: {field_label}."
            )

        raise ValueError(
            f"{self.name} has multiple EFT stages integrating {sorted(requested)}."
        )

    @property
    def uv_stage(self) -> EFTStageRecord:
        """Return the recorded UV theory."""
        matches = [stage for stage in self.eft_stages if stage.is_uv]
        if len(matches) == 1:
            return matches[0]
        if not matches:
            raise KeyError(f"{self.name} has no UV stage.")
        raise ValueError(f"{self.name} has multiple UV stages.")

    @property
    def final_eft_stage(self) -> EFTStageRecord:
        """Return the fully decoupled EFT stage with no active heavy fields."""
        return self.stage_with_active_fields()

    @property
    def first_eft_stage(self) -> EFTStageRecord:
        """Return level 1 for compatibility with existing fermion-first code.

        New production code should identify stages by field content instead of
        relying on this ordering-specific convenience property.
        """
        return self.stage(1)
