from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class EFTStageRecord:
    """One theory level in an ordered heavy-particle decoupling plan.

    level = 0 is the UV theory.  Levels >= 1 are EFTs after the corresponding
    threshold group has been integrated out.
    """

    level: int
    integrated_fields: tuple[str, ...]
    active_heavy_fields: tuple[str, ...]
    label: str
    output_dir: Path | None = None
    summary: dict = field(default_factory=dict)


# For each completed model run we have this object.
@dataclass
class RunRecord:
    name: str
    alpha: int
    d_s1: int
    d_s2: int
    d_f: int
    return_code: int
    summary: dict
    output_dir: Path
    eft_stages: list[EFTStageRecord] = field(default_factory=list)
