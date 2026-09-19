from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from common.T3Model import (
    FormalT3Dimensions,
    SharedScalarDimensions,
    physical_dimensions_from_formal,
)


@dataclass
class EFTStageRecord:
    """Metadata for one physical theory level in the threshold sequence.

    ``level == 0`` is the UV theory.

    ``integrated_fields`` and ``active_heavy_fields`` use physical field names.
    Therefore shared-scalar mode uses ``S`` rather than the formal matching
    roles ``S1`` and ``S2``.
    """

    level: int
    integrated_fields: tuple[str, ...]
    active_heavy_fields: tuple[str, ...]
    label: str
    output_dir: Path | None = None
    summary: dict = field(default_factory=dict)

    @property
    def is_uv(self) -> bool:
        """Return whether this record describes the UV theory."""
        return self.level == 0

    @property
    def is_fully_decoupled(self) -> bool:
        """Return whether no physical heavy fields remain active."""
        return not self.active_heavy_fields


@dataclass
class RunRecord:
    """Bookkeeping for one complete T3 pipeline run.

    ``d_s1``, ``d_s2`` and ``d_f`` always store the formal T3 topology
    dimensions used by the matching code.

    ``shared_scalar`` determines whether the two formal scalar roles S1 and S2
    correspond to one physical scalar S.
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
    def formal_dimensions(self) -> FormalT3Dimensions:
        """Return formal matching-topology dimensions (dS1, dS2, dF)."""
        return self.d_s1, self.d_s2, self.d_f

    @property
    def physical_dimensions(
        self,
    ) -> FormalT3Dimensions | SharedScalarDimensions:
        """Return physical heavy-field dimensions for this run.

        Ordinary T3 returns ``(dS1, dS2, dF)``.
        Shared-scalar mode returns ``(dS, dF)``.
        """
        return physical_dimensions_from_formal(
            self.d_s1,
            self.d_s2,
            self.d_f,
            shared_scalar=self.shared_scalar,
        )

    @property
    def d_s(self) -> int | None:
        """Return the physical shared-scalar dimension, if applicable."""
        if not self.shared_scalar:
            return None

        d_s, _ = self.physical_dimensions
        return d_s

    def stage(self, level: int) -> EFTStageRecord:
        """Return the stage with the requested physical EFT level.

        This lookup is intentionally based on ``EFTStageRecord.level`` rather
        than list position.  Level 0 is UV, level 1 is the first EFT after one
        threshold, and so on.
        """
        matches = [
            stage
            for stage in self.eft_stages
            if stage.level == level
        ]

        if len(matches) == 1:
            return matches[0]

        if not matches:
            raise KeyError(
                f"{self.name} has no EFT stage with level {level}."
            )

        raise ValueError(
            f"{self.name} has multiple EFT stages with level {level}."
        )

    @property
    def uv_stage(self) -> EFTStageRecord:
        """Return the UV stage (level 0)."""
        return self.stage(0)

    @property
    def first_eft_stage(self) -> EFTStageRecord:
        """Return the first EFT after the first threshold (level 1)."""
        return self.stage(1)

    @property
    def final_stage(self) -> EFTStageRecord:
        """Return the highest-level recorded theory stage."""
        if not self.eft_stages:
            raise KeyError(f"{self.name} has no EFT stage records.")

        return max(self.eft_stages, key=lambda stage: stage.level)
