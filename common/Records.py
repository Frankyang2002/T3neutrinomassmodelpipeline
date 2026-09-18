from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class EFTStageRecord:
    level: int
    integrated_fields: tuple[str, ...]
    active_heavy_fields: tuple[str, ...]
    label: str
    output_dir: Path | None = None
    summary: dict = field(default_factory=dict)


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
    shared_scalar: bool = False

    @property
    def d_s(self) -> int | None:
        return self.d_s1 if self.shared_scalar else None
