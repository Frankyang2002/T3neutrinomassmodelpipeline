"""Verified production backend for the scalar-only T3 intermediate EFT.

Physics content
---------------
This backend applies after the heavy T3 fermion ``F`` has been integrated out
while the complete physical T3 scalar sector remains active.  The interval ends
when those remaining scalars are integrated out together.

The implementation is now housed directly under ``RGE.running.intermediate``
using physical names.  Historical ``EFT1`` names remain only in serialized
keys, filenames and Wolfram/report contracts while production orchestration
identifies the theory by physical field content.

Fixed-order convention
----------------------
The production neutrino-mass path is the verified one-loop calculation already
used by the project.  The tree-generated LLSS dimension-five coefficient is run
only through the pieces that contribute to the authoritative one-loop Weinberg
coefficient.  Heavy-operator self-running followed by scalar-loop matching would
first enter at two-loop order and is therefore excluded from the one-loop result.

Independent regressions and alternative checks do not live here; see
``validation/backends/ScalarOnlyAfterFermionValidation.py``.
"""

from __future__ import annotations

import json
from pathlib import Path

from common.EFT import EFTRunningInterval
from common.RunRecords import EFTStageRecord, RunRecord
from RGE.running.intermediate.DirectWeinbergRunning import (
    build_direct_weinberg_running,
    export_direct_weinberg_insertion,
)
from RGE.running.intermediate.FlavorWilsonBoundary import (
    export_full_flavor_wilson_boundary,
)
from RGE.running.intermediate.ScalarOnlyRenormalisableRunning import (
    run_scalar_only_renormalisable_rge,
)
from RGE.running.intermediate.ScalarOnlyWilsonRunning import (
    run_scalar_only_dimension_five_wilson_rge,
)
from RGE.running.intermediate.ScalarThresholdMatching import (
    resume_scalar_threshold_with_running,
)
from RGE.running.weinberg.FinalWeinbergCoefficient import (
    build_final_weinberg_coefficient,
)


PROJECT_ROOT = Path(__file__).resolve().parents[3]
RUN_THRESHOLD_STAGE_SCRIPT = PROJECT_ROOT / "Lagrangian" / "RunThresholdStage.wl"

BACKEND_NAME = "scalar_only_after_fermion_d5"


def supports(
    record: RunRecord,
    interval: EFTRunningInterval,
) -> bool:
    """Return whether this backend matches the interval's physical content.

    The condition is entirely transition based: ``F`` is removed on entry, all
    physical T3 scalars remain active, and those scalars are removed together on
    exit.  No EFT stage number is consulted.
    """
    scalar_fields = record.physical_scalar_fields
    return bool(
        interval.active_heavy_fields == scalar_fields
        and interval.entered_by.fields_to_integrate == ("F",)
        and frozenset(interval.exited_by.fields_to_integrate) == scalar_fields
        and interval.exited_by.after.is_fully_decoupled
    )


def _run_renormalisable_rge(
    record: RunRecord,
    stage: EFTStageRecord,
) -> bool:
    """Generate renormalisable one-loop RGEs for the active scalar EFT.

    The historical output files and summary keys still use ``EFT1`` for
    compatibility, but the production interface is defined by physical field
    content: F has been removed while the complete T3 scalar sector is active.
    """

    summary = record.summary
    data_dir = record.output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    output_path = data_dir / "eft1_rgbeta_rge.json"

    print(
        f"  {record.name}: starting scalar-only intermediate-EFT renormalisable RGE...",
        flush=True,
    )

    try:
        result = run_scalar_only_renormalisable_rge(
            record.d_s1,
            record.d_s2,
            record.d_f,
            record.alpha,
            shared_scalar=record.shared_scalar,
        )
    except Exception as exc:
        summary["EFT1RenormalisableRGEStatus"] = "Failed"
        summary["EFT1RenormalisableRGEError"] = str(exc)
        print(
            f"  {record.name}: scalar-only renormalisable RGE failed: {exc}"
        )
        return False

    output_path.write_text(
        json.dumps(result.raw, indent=2),
        encoding="utf-8",
    )

    summary["EFT1RenormalisableRGEStatus"] = result.status
    summary["EFT1RenormalisableRGEFile"] = (
        output_path.relative_to(record.output_dir).as_posix()
    )
    summary["EFT1RenormalisableRGEBetaCount"] = len(result.betas)
    stage.summary.update(
        {
            "RenormalisableRGEStatus": result.status,
            "RenormalisableRGEFile": summary["EFT1RenormalisableRGEFile"],
            "RenormalisableRGEBetaCount": len(result.betas),
        }
    )

    print(
        f"  {record.name}: scalar-only renormalisable RGE={result.status} "
        f"({len(result.betas)} beta functions) -> {output_path}"
    )

    return result.status == "Success"


def _run_dimension_five_wilson_rge(
    record: RunRecord,
    stage: EFTStageRecord,
) -> bool:
    """Generate the one-loop dimension-five Wilson RGE in the scalar EFT.

    The active heavy fields are the physical T3 scalars left after removing F.
    The calculation combines the matched tree-level Wilson seed, scalar-quartic
    seed and renormalisable scalar-EFT metadata with the general
    ``psi^2 phi^2`` master RGE.

    Fixed-order bookkeeping:
        only C^(0), the tree-generated Wilson tensor, is evolved.  The one-loop
        matching contribution at the F threshold remains a boundary term and is
        not inserted into this one-loop beta function.
    """

    summary = record.summary
    data_dir = record.output_dir / "data"
    wilson_seed_path = data_dir / "eft1_after_F_wilson_seed.json"
    quartic_seed_path = data_dir / "eft1_after_F_scalar_quartic_seed.json"
    rgbeta_path = data_dir / "eft1_rgbeta_rge.json"
    output_path = data_dir / "eft1_wilson_rge.json"

    required_inputs = (
        wilson_seed_path,
        quartic_seed_path,
        rgbeta_path,
    )
    missing = [path for path in required_inputs if not path.is_file()]

    if missing:
        summary["EFT1WilsonRGEStatus"] = "Failed"
        summary["EFT1WilsonRGEError"] = (
            "Missing scalar-only Wilson-RGE input(s): "
            + ", ".join(
                path.relative_to(record.output_dir).as_posix()
                for path in missing
            )
        )
        print(
            f"  {record.name}: scalar-only Wilson RGE failed: "
            f"{summary['EFT1WilsonRGEError']}"
        )
        return False

    print(
        f"  {record.name}: starting scalar-only dimension-five Wilson RGE stage...",
        flush=True,
    )

    try:
        result = run_scalar_only_dimension_five_wilson_rge(
            wilson_seed_path=wilson_seed_path,
            quartic_seed_path=quartic_seed_path,
            rgbeta_path=rgbeta_path,
            output_path=output_path,
            seed_support_only=False,
        )
    except Exception as exc:
        summary["EFT1WilsonRGEStatus"] = "Failed"
        summary["EFT1WilsonRGEError"] = str(exc)
        print(
            f"  {record.name}: scalar-only dimension-five Wilson RGE failed: {exc}"
        )
        return False

    summary["EFT1WilsonRGEStatus"] = result.get("status", "Unknown")
    summary["EFT1WilsonRGEFile"] = (
        output_path.relative_to(record.output_dir).as_posix()
    )
    summary["EFT1WilsonRGEInitialComponentCount"] = result.get(
        "initial_independent_component_count",
        0,
    )
    summary["EFT1WilsonRGEBetaComponentCount"] = result.get(
        "nonzero_beta_component_count",
        0,
    )
    summary["EFT1WilsonRGEGeneratedComponentCount"] = result.get(
        "generated_component_count",
        0,
    )

    # Mirror the calculation metadata onto the physical intermediate-stage record so
    # stage-aware reporting can consume it without rediscovering files.
    if record.eft_stages:
        stage.summary.update(
            {
                "RenormalisableRGEStatus": summary.get(
                    "EFT1RenormalisableRGEStatus",
                    "Unknown",
                ),
                "RenormalisableRGEFile": summary.get(
                    "EFT1RenormalisableRGEFile",
                    "",
                ),
                "WilsonRGEStatus": summary["EFT1WilsonRGEStatus"],
                "WilsonRGEFile": summary["EFT1WilsonRGEFile"],
                "WilsonRGEInitialComponentCount": summary[
                    "EFT1WilsonRGEInitialComponentCount"
                ],
                "WilsonRGEBetaComponentCount": summary[
                    "EFT1WilsonRGEBetaComponentCount"
                ],
                "WilsonRGEGeneratedComponentCount": summary[
                    "EFT1WilsonRGEGeneratedComponentCount"
                ],
            }
        )

    print(
        f"  {record.name}: scalar-only Wilson RGE="
        f"{summary['EFT1WilsonRGEStatus']} "
        f"({summary['EFT1WilsonRGEBetaComponentCount']} nonzero beta "
        f"components; "
        f"{summary['EFT1WilsonRGEGeneratedComponentCount']} generated)"
        f" -> {output_path}"
    )

    return summary["EFT1WilsonRGEStatus"] == "Success"


def _production_paths(record: RunRecord) -> dict[str, Path]:
    """Return all files used by the fixed-one-loop scalar-EFT Weinberg bridge."""
    data_dir = record.output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    return {
        "data_dir": data_dir,
        "wilson_seed": data_dir / "eft1_after_F_wilson_seed.json",
        "component_rge": data_dir / "eft1_wilson_rge.json",
        "flavor_seed": data_dir / "eft1_wilson_flavor_seed.json",
        "flavor_transport": data_dir / "eft1_wilson_flavor_at_S_threshold.json",
        "insertion": data_dir / "eft1_wilson_flavor_running_insertion.wl",
        "pole_rge": data_dir / "c5_pole_rge_consistency.json",
        "threshold_c5": data_dir / "c5_threshold_with_eft1_running.txt",
        "direct_c5": data_dir / "c5_direct_eft1_running.txt",
        "final_c5": data_dir / "final_weinberg_coefficient.json",
    }


def _fail_full_flavor_matching(record: RunRecord, message: str) -> bool:
    """Record and print a failed fixed-one-loop scalar-EFT Weinberg bridge."""
    record.summary["EFT1FullFlavorBridgeStatus"] = "Failed"
    record.summary["EFT1FullFlavorBridgeError"] = message
    print(f"  {record.name}: full-flavor scalar-EFT bridge failed: {message}")
    return False


def _require_full_flavor_inputs(
    record: RunRecord,
    paths: dict[str, Path],
) -> bool:
    """Check that the direct-Weinberg bridge inputs have been produced."""
    required = (paths["wilson_seed"], paths["component_rge"])
    missing = [path for path in required if not path.is_file()]
    if not missing:
        return True

    message = (
        "Missing direct-Weinberg bridge input(s): "
        + ", ".join(
            path.relative_to(record.output_dir).as_posix()
            for path in missing
        )
    )
    return _fail_full_flavor_matching(record, message)


def _run_direct_weinberg_transport_and_resume(
    record: RunRecord,
    paths: dict[str, Path],
    *,
    mu_high: str,
    mu_low: str,
    debug_reports: bool,
) -> tuple[dict, dict, dict, dict]:
    """Run flavor seeding, direct Weinberg transport, and threshold resume."""
    flavor_seed = export_full_flavor_wilson_boundary(
        wilson_seed_path=paths["wilson_seed"],
        output_path=paths["flavor_seed"],
    )
    if flavor_seed.get("status") != "Success":
        raise RuntimeError("Full-flavor Wilson boundary construction failed.")

    flavor_transport = build_direct_weinberg_running(
        flavor_boundary_path=paths["flavor_seed"],
        component_rge_path=paths["component_rge"],
        mu_high=mu_high,
        mu_low=mu_low,
        output_path=paths["flavor_transport"],
    )
    if flavor_transport.get("status") != "Success":
        raise RuntimeError(
            "Direct full-flavor Weinberg running construction failed."
        )

    insertion = export_direct_weinberg_insertion(
        transport_path=paths["flavor_transport"],
        output_path=paths["insertion"],
    )
    if insertion.get("status") != "Success":
        raise RuntimeError("Direct Weinberg Matchete insertion export failed.")

    resume = resume_scalar_threshold_with_running(
        output_dir=record.output_dir,
        running_insertion=paths["insertion"],
        run_threshold_script=RUN_THRESHOLD_STAGE_SCRIPT,
        d_s1=record.d_s1,
        d_s2=record.d_s2,
        d_f=record.d_f,
        alpha=record.alpha,
        shared_scalar=record.shared_scalar,
        validation_mode=debug_reports,
    )
    if resume.get("status") != "Success":
        raise RuntimeError(
            "Threshold-2 rerun with direct Weinberg running failed."
        )

    return flavor_seed, flavor_transport, insertion, resume


def _build_authoritative_threshold_c5(
    record: RunRecord,
    paths: dict[str, Path],
    *,
    flavor_transport: dict,
) -> dict:
    """Build the authoritative hierarchical Weinberg coefficient."""
    if not paths["threshold_c5"].is_file():
        raise RuntimeError(
            "Authoritative resumed threshold C5 was not exported: "
            f"{paths['threshold_c5']}"
        )

    try:
        final_c5 = build_final_weinberg_coefficient(
            threshold_c5_path=paths["threshold_c5"],
            flavor_transport_path=paths["flavor_transport"],
            output_path=paths["final_c5"],
        )
    except Exception as exc:
        raise RuntimeError(
            f"Final Weinberg coefficient bookkeeping failed: {exc}"
        ) from exc

    if final_c5.get("status") != "Success":
        raise RuntimeError(
            "Final Weinberg coefficient bookkeeping returned non-success."
        )

    combined = final_c5.get("combined", {})
    if not combined.get("ready_for_physical_majorana_numerics", False):
        raise RuntimeError(
            "Final Weinberg coefficient is not ready for physical Majorana "
            f"downstream use: {combined.get('physical_majorana_reason')!r}"
        )

    return final_c5


def _record_full_flavor_success(
    record: RunRecord,
    paths: dict[str, Path],
    *,
    flavor_seed: dict,
    flavor_transport: dict,
    insertion: dict,
    resume: dict,
    final_c5: dict,
) -> None:
    """Store successful scalar-EFT bridge provenance in the run summary."""
    summary = record.summary
    combined = final_c5.get("combined", {})

    # Every EFT-side stage after this point consumes the authoritative
    # hierarchical full-flavor coefficient, never the preliminary C5 exported
    # before the resumed threshold calculation.
    summary["WeinbergCoefficientFile"] = (
        paths["final_c5"].relative_to(record.output_dir).as_posix()
    )

    summary.update(
        {
            "EFT1FlavorSeedStatus": flavor_seed.get("status"),
            "EFT1FlavorSeedFile": paths["flavor_seed"].relative_to(
                record.output_dir
            ).as_posix(),
            "EFT1FlavorRGEStatus": "DiagnosticNotRequiredForOneLoopC5",
            "EFT1FlavorRGEOneGenerationCheck": None,
            "EFT1FlavorTransportStatus": flavor_transport.get("status"),
            "EFT1FlavorTransportFile": paths["flavor_transport"].relative_to(
                record.output_dir
            ).as_posix(),
            "EFT1FlavorRunningInsertionStatus": insertion.get("status"),
            "EFT1FlavorRunningInsertionFile": paths["insertion"].relative_to(
                record.output_dir
            ).as_posix(),
            "EFT1HeavySelfRunningIncludedInAuthoritativeC5": False,
            "EFT1DirectWeinbergOnlyAtOneLoop": True,
            "EFT1ThresholdResumeStatus": resume.get("status"),
            "EFT1ThresholdResumeValidationMode": resume.get(
                "validation_mode",
                False,
            ),
            "EFT1ThresholdResumeResultFile": str(
                Path(resume["result_path"]).resolve()
            ),
            "EFT1ThresholdResumeInsertionLoaded": resume.get(
                "running_insertion_loaded"
            ),
            "EFT1ThresholdResumeInsertedInCOnly": resume.get(
                "running_inserted_in_C_only"
            ),
            "EFT1ThresholdResumeDirectWeinbergCarried": resume.get(
                "direct_weinberg_carried_separately"
            ),
            "EFT1ThresholdResumeDirectWeinbergEqualScaleCheck": resume.get(
                "direct_weinberg_equal_scale_vanishes"
            ),
            "EFT1FullFlavorBridgeStatus": "Success",
            "AuthoritativeThresholdC5File": paths["threshold_c5"].relative_to(
                record.output_dir
            ).as_posix(),
            "AuthoritativeDirectRunningC5File": (
                paths["direct_c5"].relative_to(record.output_dir).as_posix()
                if paths["direct_c5"].is_file()
                else ""
            ),
            "FinalWeinbergCoefficientFile": paths["final_c5"].relative_to(
                record.output_dir
            ).as_posix(),
            "AuthoritativeFinalC5Status": "ReadyForDownstream",
            "AuthoritativeFinalC5FullFlavorReady": combined.get(
                "ready_for_full_flavor_numerics"
            ),
            "AuthoritativeFinalC5PhysicalMajoranaReady": combined.get(
                "ready_for_physical_majorana_numerics"
            ),
            "AuthoritativeDownstreamC5File": paths["final_c5"].relative_to(
                record.output_dir
            ).as_posix(),
        }
    )


def _run_full_flavor_weinberg_matching(
    record: RunRecord,
    stage: EFTStageRecord,
    *,
    mu_high: str | float,
    mu_low: str | float,
    debug_reports: bool = False,
) -> bool:
    """Build the authoritative fixed-one-loop hierarchical C5.

    The tree LLSS Wilson boundary is O(hbar^0). Its direct mixing into O5 is
    O(hbar^1). LLSS self-running is also O(hbar^1), so matching that corrected
    heavy operator through the scalar loop would be O(hbar^2). Therefore only
    the direct full-flavor C12 -> Weinberg running enters the authoritative
    one-loop neutrino-mass calculation.
    """
    paths = _production_paths(record)
    if not _require_full_flavor_inputs(record, paths):
        stage.summary["FullFlavorThresholdBridgeStatus"] = "Failed"
        return False

    print(
        f"  {record.name}: starting fixed-one-loop direct Weinberg bridge "
        f"{mu_high} -> {mu_low}...",
        flush=True,
    )

    try:
        flavor_seed, flavor_transport, insertion, resume = (
            _run_direct_weinberg_transport_and_resume(
                record,
                paths,
                mu_high=mu_high,
                mu_low=mu_low,
                debug_reports=debug_reports,
            )
        )
        final_c5 = _build_authoritative_threshold_c5(
            record,
            paths,
            flavor_transport=flavor_transport,
        )
    except Exception as exc:
        stage.summary["FullFlavorThresholdBridgeStatus"] = "Failed"
        return _fail_full_flavor_matching(record, str(exc))

    _record_full_flavor_success(
        record,
        paths,
        flavor_seed=flavor_seed,
        flavor_transport=flavor_transport,
        insertion=insertion,
        resume=resume,
        final_c5=final_c5,
    )

    stage.summary["FullFlavorThresholdBridgeStatus"] = "Success"
    print(
        f"  {record.name}: direct one-loop Weinberg bridge=Success "
        f"(heavy self-running excluded at O(hbar))",
        flush=True,
    )
    return True

def run(
    record: RunRecord,
    stage: EFTStageRecord,
    interval: EFTRunningInterval,
    *,
    debug_reports: bool = False,
) -> bool:
    """Run the verified scalar-only intermediate EFT calculation.

    The order here is part of the retained production behaviour:

    1. renormalisable RGBeta running for the active scalar theory;
    2. dimension-five LLSS Wilson RGE;
    3. direct full-flavor LLSS -> Weinberg running and the lower-threshold
       matching needed to construct the authoritative one-loop ``C5``.
    """
    if not supports(record, interval):
        raise ValueError(
            "ScalarOnlyAfterFermion backend was called for incompatible EFT content."
        )

    if not _run_renormalisable_rge(record, stage):
        return False

    if not _run_dimension_five_wilson_rge(record, stage):
        return False

    return _run_full_flavor_weinberg_matching(
        record,
        stage,
        mu_high=interval.high_scale,
        mu_low=interval.low_scale,
        debug_reports=debug_reports,
    )
