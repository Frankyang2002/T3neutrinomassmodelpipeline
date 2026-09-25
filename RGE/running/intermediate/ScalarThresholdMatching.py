"""Lower scalar-threshold matching after intermediate-EFT running."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from RGE.running.intermediate.LegacyEFT1Compatibility import (
    LegacyScalarThresholdContinuation,
    legacy_resume_scalar_threshold_with_running,
)


def resume_scalar_threshold_with_running(
    *,
    output_dir: Path,
    running_insertion: Path,
    run_threshold_script: Path,
    d_s1: int,
    d_s2: int,
    d_f: int,
    alpha: int,
    eft_order: int = 5,
    loop_order: int = 1,
    continuation: LegacyScalarThresholdContinuation | None = None,
    result_path: Path | None = None,
    shared_scalar: bool = False,
    validation_mode: bool = False,
) -> dict[str, Any]:
    """Resume the existing scalar threshold with the calculated running term."""
    return legacy_resume_scalar_threshold_with_running(
        output_dir=output_dir,
        running_insertion=running_insertion,
        run_threshold_script=run_threshold_script,
        d_s1=d_s1,
        d_s2=d_s2,
        d_f=d_f,
        alpha=alpha,
        eft_order=eft_order,
        loop_order=loop_order,
        continuation=continuation,
        result_path=result_path,
        shared_scalar=shared_scalar,
        validation_mode=validation_mode,
    )
