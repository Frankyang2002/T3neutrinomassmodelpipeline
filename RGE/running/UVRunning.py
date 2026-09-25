"""One-loop renormalisable running in the ultraviolet T3 theory."""

from __future__ import annotations

import json

from common.RunRecords import RunRecord
from RGE.running.rgbeta.RGBetaT3Running import run_rgbeta_t3


def run_uv_rge(record: RunRecord) -> bool:
    """Generate the one-loop beta functions for the renormalisable UV T3 model.

    This stage describes running above the heavy-particle matching thresholds,
    before any T3 field has been integrated out.
    """
    summary = record.summary

    # Fails if we dont have an actual build for our model
    if summary.get("BuildStatus") != "Success":
        summary["UVRGEStatus"] = "NotRun"
        return False

    data_dir = record.output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    output_path = data_dir / "uv_rgbeta_rge.json"

    print(f"  {record.name}: starting RGBeta UV-RGE stage...", flush=True)

    # Use run_rgbeta_t3
    try:
        result = run_rgbeta_t3(
            record.d_s1,
            record.d_s2,
            record.d_f,
            record.alpha,
            shared_scalar=record.shared_scalar,
        )
    except Exception as exc:
        summary["UVRGEStatus"] = "Failed"
        summary["UVRGEError"] = str(exc)
        print(f"  {record.name}: RGBeta UV-RGE failed: {exc}")
        return False

    # Pass the T3 SU(2) representation dimensions and hypercharge parameter
    # to the RGBeta Wolfram runner, which constructs the UV model beta functions.
    output_path.write_text(
        json.dumps(result.raw, indent=2),
        encoding="utf-8",
    )

    summary["UVRGEStatus"] = result.status
    summary["UVRGEFile"] = output_path.relative_to(record.output_dir).as_posix()
    summary["UVRGEBetaCount"] = len(result.betas)

    print(
        f"  {record.name}: RGBeta UV-RGE={result.status} "
        f"({len(result.betas)} beta functions) -> {output_path}"
    )

    return result.status == "Success"
