"""Generic dispatch for intermediate heavy-field EFT running.

``pipeline.py`` owns the order of the calculation and passes each physical
:class:`~common.EFT.EFTRunningInterval` here.  This module decides which
production backend can evaluate that interval; it does not contain the detailed
RGE or matching algebra itself.

This keeps the central architecture independent of historical labels such as
``EFT1``.  Backend implementations are selected from the actual fields active
between two thresholds.
"""

from __future__ import annotations

from dataclasses import dataclass

from common.EFT import EFTRunningInterval
from common.RunRecords import EFTStageRecord, RunRecord
from RGE.running.backends.ScalarOnlyAfterFermion import (
    BACKEND_NAME as SCALAR_ONLY_AFTER_FERMION_BACKEND,
    run as run_scalar_only_after_fermion,
    supports as supports_scalar_only_after_fermion,
)
from RGE.running.intermediate.ScalarFirstDimensionSixSeed import (
    extract_scalar_first_dimension_six_seed,
    supports_scalar_first_dimension_six_seed,
)


@dataclass(frozen=True, slots=True)
class IntermediateEFTRunningOutcome:
    """Result of attempting one intermediate heavy-field running interval.

    ``attempted`` is false when no production backend exists for the interval's
    field content.  That state is explicit rather than silently routing the
    interval through an incorrect backend.  The central pipeline treats such a
    required-but-unimplemented interval as incomplete production and blocks the
    authoritative low-energy path.
    """

    attempted: bool
    success: bool
    backend: str | None
    status: str

    def __post_init__(self) -> None:
        if not self.attempted and not self.success:
            raise ValueError(
                "An unattempted intermediate running interval cannot be marked failed."
            )
        if self.attempted and self.backend is None:
            raise ValueError(
                "An attempted intermediate running interval must identify its backend."
            )


def intermediate_running_backend(
    record: RunRecord,
    interval: EFTRunningInterval,
) -> str | None:
    """Return the production backend for an interval, if one is implemented."""
    if supports_scalar_only_after_fermion(record, interval):
        return SCALAR_ONLY_AFTER_FERMION_BACKEND
    return None


def run_intermediate_eft_interval(
    record: RunRecord,
    stage: EFTStageRecord,
    interval: EFTRunningInterval,
    *,
    debug_reports: bool = False,
) -> IntermediateEFTRunningOutcome:
    """Run one intermediate EFT region using its physical field content.

    Unsupported field contents remain represented by the generic EFT plan but
    return ``NotImplementedForFieldContent``.  ``pipeline.py`` records that as an
    incomplete production route and does not let it silently feed the final
    neutrino result.  This is preferable to pretending that a fermion-first
    calculation also implements scalar-first running.
    """
    backend = intermediate_running_backend(record, interval)

    if backend is None:
        status = "NotImplementedForFieldContent"

        # A d=6 scalar-first run already contains the exact Matchete tree EFT.
        # Preserve the candidate psi^2 phi^3 boundary terms now, without
        # pretending that their tensor adapter/RGE is production-complete.
        if supports_scalar_first_dimension_six_seed(record, stage, interval):
            try:
                extract_scalar_first_dimension_six_seed(record, stage, interval)
            except Exception as exc:
                stage.summary["ScalarFirstDimensionSixSeedStatus"] = "Failed"
                stage.summary["ScalarFirstDimensionSixSeedError"] = str(exc)

        stage.summary.setdefault("IntermediateRunningStatus", status)
        stage.summary.setdefault("IntermediateRunningBackend", None)
        return IntermediateEFTRunningOutcome(
            attempted=False,
            success=True,
            backend=None,
            status=status,
        )

    stage.summary["IntermediateRunningBackend"] = backend
    stage.summary["IntermediateRunningStatus"] = "Running"

    if backend == SCALAR_ONLY_AFTER_FERMION_BACKEND:
        success = run_scalar_only_after_fermion(
            record,
            stage,
            interval,
            debug_reports=debug_reports,
        )
    else:  # Defensive guard for future registry edits.
        raise RuntimeError(f"No runner registered for backend {backend!r}.")

    status = "Success" if success else "Failed"
    stage.summary["IntermediateRunningStatus"] = status

    return IntermediateEFTRunningOutcome(
        attempted=True,
        success=success,
        backend=backend,
        status=status,
    )


def mark_intermediate_running_not_applicable(record: RunRecord) -> None:
    """Preserve legacy report keys when no intermediate backend was attempted.

    These keys are compatibility outputs only.  They must never be used to
    identify or dispatch an EFT in new production code.
    """
    record.summary.setdefault("EFT1RenormalisableRGEStatus", "NotApplicable")
    record.summary.setdefault("EFT1WilsonRGEStatus", "NotApplicable")
