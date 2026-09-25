"""Generic post-production validation dispatch for intermediate EFT regions.

Validation is separate from the production calculation.  This module selects an
independent validation backend from the physical EFT content and records the
result; detailed regression calculations live under ``validation/backends``.
"""

from __future__ import annotations

from dataclasses import dataclass

from common.EFT import EFTRunningInterval
from common.RunRecords import EFTStageRecord, RunRecord
from RGE.running.IntermediateEFTRunning import (
    SCALAR_ONLY_AFTER_FERMION_BACKEND,
    intermediate_running_backend,
)
from validation.backends.ScalarOnlyAfterFermionValidation import (
    BACKEND_NAME as SCALAR_ONLY_AFTER_FERMION_VALIDATION,
    run as run_scalar_only_after_fermion_validation,
)


@dataclass(frozen=True, slots=True)
class IntermediateEFTValidationOutcome:
    """Result of validating one completed production EFT interval."""

    attempted: bool
    success: bool
    backend: str | None
    status: str

    def __post_init__(self) -> None:
        if not self.attempted and not self.success:
            raise ValueError(
                "An unattempted validation interval cannot be marked failed."
            )
        if self.attempted and self.backend is None:
            raise ValueError(
                "An attempted validation interval must identify its backend."
            )


def intermediate_validation_backend(
    record: RunRecord,
    interval: EFTRunningInterval,
) -> str | None:
    """Return the validation backend corresponding to production support."""
    production_backend = intermediate_running_backend(record, interval)
    if production_backend == SCALAR_ONLY_AFTER_FERMION_BACKEND:
        return SCALAR_ONLY_AFTER_FERMION_VALIDATION
    return None


def run_intermediate_eft_validation(
    record: RunRecord,
    stage: EFTStageRecord,
    interval: EFTRunningInterval,
) -> IntermediateEFTValidationOutcome:
    """Validate one completed intermediate EFT interval.

    Validation may read production artifacts and write diagnostic artifacts, but
    it must not construct or replace the production EFT result.
    """
    backend = intermediate_validation_backend(record, interval)

    if backend is None:
        status = "NotImplementedForFieldContent"
        stage.summary.setdefault("IntermediateValidationStatus", status)
        stage.summary.setdefault("IntermediateValidationBackend", None)
        return IntermediateEFTValidationOutcome(
            attempted=False,
            success=True,
            backend=None,
            status=status,
        )

    stage.summary["IntermediateValidationBackend"] = backend
    stage.summary["IntermediateValidationStatus"] = "Running"
    record.summary["IntermediateEFTValidationBackend"] = backend
    record.summary["IntermediateEFTValidationStatus"] = "Running"

    if backend == SCALAR_ONLY_AFTER_FERMION_VALIDATION:
        success = run_scalar_only_after_fermion_validation(
            record,
            stage,
            interval,
        )
    else:  # Defensive guard for future validation registry edits.
        raise RuntimeError(f"No validation runner registered for backend {backend!r}.")

    status = "Success" if success else "Failed"
    stage.summary["IntermediateValidationStatus"] = status
    record.summary["IntermediateEFTValidationStatus"] = status

    return IntermediateEFTValidationOutcome(
        attempted=True,
        success=success,
        backend=backend,
        status=status,
    )
