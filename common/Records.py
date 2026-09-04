from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


# For each completed run we have this object
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
