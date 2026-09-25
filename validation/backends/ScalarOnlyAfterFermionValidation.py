"""Validation backend for the scalar-only-after-fermion EFT interval.

This module is intentionally outside the production RGE package.  It consumes
artifacts already produced by the scalar-only production backend and performs
independent regression/consistency checks.  It must not construct or modify the
authoritative EFT result.

The historical ``EFT1...`` summary keys and diagnostic filename are preserved
for report compatibility only.  Their validation values are populated here,
not by the production backend.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from common.EFT import EFTRunningInterval
from common.RunRecords import EFTStageRecord, RunRecord
from RGE.running.intermediate.LegacyEFT1Compatibility import (
    legacy_component_wilson_transport,
)
from RGE.running.weinberg.FinalWeinbergCoefficient import (
    normalize_pole_rge_consistency,
)


BACKEND_NAME = "scalar_only_after_fermion_d5_validation"


def _run_component_transport_regression(
    record: RunRecord,
    stage: EFTStageRecord,
    *,
    mu_high: str | float,
    mu_low: str | float,
) -> bool:
    """Run the historical one-generation Wilson-transport regression.

    This calculation is not used to construct the authoritative full-flavor
    one-loop C5.  It is retained as an independent compact check of the
    component RGE and fixed-order leading-log transport.
    """
    summary = record.summary

    if summary.get("EFT1WilsonRGEStatus") != "Success":
        summary["EFT1WilsonTransportStatus"] = "NotRun"
        return False

    data_dir = record.output_dir / "data"
    wilson_seed_path = data_dir / "eft1_after_F_wilson_seed.json"
    wilson_rge_path = data_dir / "eft1_wilson_rge.json"
    rgbeta_path = data_dir / "eft1_rgbeta_rge.json"
    output_path = data_dir / "eft1_wilson_at_S_threshold.json"

    required_inputs = (
        wilson_seed_path,
        wilson_rge_path,
        rgbeta_path,
    )
    missing = [path for path in required_inputs if not path.is_file()]
    if missing:
        summary["EFT1WilsonTransportStatus"] = "Failed"
        summary["EFT1WilsonTransportError"] = (
            "Missing EFT1 Wilson-transport validation input(s): "
            + ", ".join(
                path.relative_to(record.output_dir).as_posix()
                for path in missing
            )
        )
        print(
            f"  {record.name}: Wilson-transport validation failed: "
            f"{summary['EFT1WilsonTransportError']}"
        )
        return False

    print(
        f"  {record.name}: validating one-generation Wilson transport "
        f"{mu_high} -> {mu_low}...",
        flush=True,
    )

    try:
        result = legacy_component_wilson_transport(
            wilson_seed_path=wilson_seed_path,
            wilson_rge_path=wilson_rge_path,
            rgbeta_path=rgbeta_path,
            mu_high=mu_high,
            mu_low=mu_low,
            output_path=output_path,
        )
    except Exception as exc:
        summary["EFT1WilsonTransportStatus"] = "Failed"
        summary["EFT1WilsonTransportError"] = str(exc)
        print(
            f"  {record.name}: Wilson-transport validation failed: {exc}"
        )
        return False

    summary["EFT1WilsonTransportStatus"] = result.get("status", "Unknown")
    summary["EFT1WilsonTransportFile"] = (
        output_path.relative_to(record.output_dir).as_posix()
    )
    summary["EFT1WilsonTransportMuHigh"] = mu_high
    summary["EFT1WilsonTransportMuLow"] = mu_low
    summary["EFT1WilsonTransportLogRatio"] = (
        result.get("scales", {}).get("log_ratio", "")
    )
    summary["EFT1WilsonTransportCorrectionCount"] = result.get(
        "running_correction_component_count",
        0,
    )
    summary["EFT1WilsonTransportGeneratedCount"] = result.get(
        "generated_by_running_component_count",
        0,
    )
    summary["EFT1WilsonTransportEqualScaleCheck"] = result.get(
        "equal_scale_running_vanishes"
    )

    stage.summary.update(
        {
            "WilsonTransportStatus": summary["EFT1WilsonTransportStatus"],
            "WilsonTransportFile": summary["EFT1WilsonTransportFile"],
            "WilsonTransportMuHigh": mu_high,
            "WilsonTransportMuLow": mu_low,
            "WilsonTransportLogRatio": summary[
                "EFT1WilsonTransportLogRatio"
            ],
            "WilsonTransportCorrectionCount": summary[
                "EFT1WilsonTransportCorrectionCount"
            ],
            "WilsonTransportGeneratedCount": summary[
                "EFT1WilsonTransportGeneratedCount"
            ],
            "WilsonTransportEqualScaleCheck": summary[
                "EFT1WilsonTransportEqualScaleCheck"
            ],
        }
    )

    return bool(
        summary["EFT1WilsonTransportStatus"] == "Success"
        and summary["EFT1WilsonTransportEqualScaleCheck"] is True
    )


def _load_json_diagnostic(
    path: Path,
    *,
    label: str,
    failures: list[str],
) -> dict[str, Any] | None:
    """Load one production JSON artifact for independent validation."""
    if not path.is_file():
        failures.append(f"{label} payload is missing")
        return None

    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        failures.append(f"{label} payload is invalid: {exc}")
        return None

    if not isinstance(payload, dict):
        failures.append(f"{label} payload is not a JSON object")
        return None

    return payload


def _validate_production_diagnostics(
    record: RunRecord,
    stage: EFTStageRecord,
    failures: list[str],
) -> None:
    """Interpret diagnostic metadata emitted by production calculations.

    The production backend writes the raw RGE/flavor artifacts because they are
    required for the calculation.  Whether their diagnostic reductions pass is
    a scientific cross-check and is recorded here.  Historical summary keys are
    filled here so existing reports remain unchanged.
    """
    summary = record.summary
    data_dir = record.output_dir / "data"

    component_rge = _load_json_diagnostic(
        data_dir / "eft1_wilson_rge.json",
        label="component Wilson RGE",
        failures=failures,
    )
    if component_rge is not None:
        projection = component_rge.get("weinberg_subspace_validation") or {}
        subspace_ok = projection.get("matches_weinberg_subspace") is True
        beta_kappa = projection.get("beta_kappa_16pi2", "")
        summary["EFT1WilsonRGEWeinbergSubspaceValid"] = subspace_ok
        summary["EFT1WilsonRGEWeinbergBeta"] = beta_kappa
        stage.summary["WeinbergSubspaceValid"] = subspace_ok
        stage.summary["WeinbergBeta16Pi2"] = beta_kappa
        if not subspace_ok:
            failures.append("component RGE does not lie in the Weinberg subspace")

    flavor_seed = _load_json_diagnostic(
        data_dir / "eft1_wilson_flavor_seed.json",
        label="full-flavor Wilson boundary",
        failures=failures,
    )
    if flavor_seed is not None:
        one_generation_ok = (
            flavor_seed.get("one_generation_reduction_matches") is True
        )
        summary["EFT1FlavorSeedOneGenerationCheck"] = one_generation_ok
        if not one_generation_ok:
            failures.append("full-flavor seed fails the one-generation reduction")

    flavor_transport = _load_json_diagnostic(
        data_dir / "eft1_wilson_flavor_at_S_threshold.json",
        label="direct Weinberg running",
        failures=failures,
    )
    if flavor_transport is not None:
        equal_scale_ok = (
            flavor_transport.get("equal_scale_running_vanishes") is True
        )
        summary["EFT1FlavorTransportEqualScaleCheck"] = equal_scale_ok
        if not equal_scale_ok:
            failures.append(
                "full-flavor direct running does not vanish at equal scales"
            )


def _validate_consistency_metadata(
    record: RunRecord,
    stage: EFTStageRecord,
) -> bool:
    """Validate independent consistency metadata from production outputs."""
    summary = record.summary
    failures: list[str] = []

    _validate_production_diagnostics(record, stage, failures)

    pole_path = record.output_dir / "data" / "c5_pole_rge_consistency.json"
    if not pole_path.is_file():
        failures.append("pole/RGE consistency payload is missing")
        pole_payload = None
    else:
        try:
            pole_payload = normalize_pole_rge_consistency(pole_path)
        except Exception as exc:
            failures.append(f"pole/RGE consistency payload is invalid: {exc}")
            pole_payload = None

    if pole_payload is not None:
        summary["PoleRGEConsistencyValidated"] = pole_payload.get(
            "ConsistencyValidated"
        )
        summary["PoleRGEConsistencyReason"] = pole_payload.get(
            "ConsistencyReason"
        )
        summary["PoleRGEConsistencyRatio"] = pole_payload.get(
            "DirectLogToPoleRatioOneGenerationInputForm"
        )
        if pole_payload.get("ConsistencyValidated") is not True:
            failures.append(
                "direct-running log coefficient is not twice the hard pole residue"
            )

    final_c5_path = record.output_dir / "data" / "final_weinberg_coefficient.json"
    if not final_c5_path.is_file():
        failures.append("authoritative final C5 payload is missing")

    success = not failures
    status = "Success" if success else "Failed"
    summary["IntermediateEFTValidationStatus"] = status
    summary["IntermediateEFTValidationFailures"] = failures
    stage.summary["IntermediateValidationStatus"] = status
    stage.summary["IntermediateValidationFailures"] = failures

    if failures:
        print(
            f"  {record.name}: intermediate EFT validation failed: "
            + "; ".join(failures)
        )

    return success


def run(
    record: RunRecord,
    stage: EFTStageRecord,
    interval: EFTRunningInterval,
) -> bool:
    """Run all independent checks for the verified scalar-only interval."""
    transport_ok = _run_component_transport_regression(
        record,
        stage,
        mu_high=interval.high_scale,
        mu_low=interval.low_scale,
    )
    metadata_ok = _validate_consistency_metadata(record, stage)
    return transport_ok and metadata_ok
