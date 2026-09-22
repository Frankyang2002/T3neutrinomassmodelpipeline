
# Python wrapper for Wolfram RGBeta runners for python pipeline
# It just calls our wolfram scripts in python
# We run RGBeta for UV and RGBeta running

from __future__ import annotations

import json
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class RGBetaT3Result:
    """Result returned by the Wolfram RGBeta UV T3 runner."""

    status: str
    metadata: dict[str, Any]
    betas: dict[str, str]
    report_betas: dict[str, str]
    report_latex_betas: dict[str, str]
    raw: dict[str, Any]


@dataclass(frozen=True)
class RGBetaT3IntermediateResult:
    """Result returned by the Wolfram RGBeta intermediate-EFT runner."""

    status: str
    metadata: dict[str, Any]
    betas: dict[str, str]
    report_betas: dict[str, str]
    report_latex_betas: dict[str, str]
    raw: dict[str, Any]


def _wolfram_integer_token(value: int) -> str:
    """Give -3 -> m3 for wolfram to read"""
    return f"m{abs(value)}" if value < 0 else str(value)


def _run_rgbeta(
    d_s1: int,
    d_s2: int,
    d_f: int,
    alpha: int,
    *,
    shared_scalar: bool,
    runner_path: Path,
    output_name: str,
    failure_message: str,
    result_type: type[RGBetaT3Result] | type[RGBetaT3IntermediateResult],
    wolframscript: str,
):
    """Run one of the UV/EFT1 RGBeta Wolfram front ends and parse its JSON output."""

    if any(d not in {1, 2, 3} for d in (d_s1, d_s2, d_f)):
        raise ValueError(
            "RGBeta T3 running currently supports only SU(2) dimensions 1, 2 and 3."
        )

    runner = Path(runner_path)
    if not runner.exists():
        raise FileNotFoundError(f"RGBeta Wolfram runner not found: {runner}")

    with tempfile.TemporaryDirectory(prefix="t3_rgbeta_") as tmpdir:
        output_path = Path(tmpdir) / output_name

        if shared_scalar:
            if alpha != -1 or d_s1 != d_s2:
                raise ValueError(
                    "Shared-scalar RGBeta mode requires alpha=-1 and dS1=dS2."
                )
            command = [
                wolframscript,
                "-file",
                str(runner),
                "SHARED",
                _wolfram_integer_token(d_s1),
                _wolfram_integer_token(d_f),
                str(output_path),
            ]
        else:
            command = [
                wolframscript,
                "-file",
                str(runner),
                _wolfram_integer_token(d_s1),
                _wolfram_integer_token(d_s2),
                _wolfram_integer_token(d_f),
                _wolfram_integer_token(alpha),
                str(output_path),
            ]

        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False,
        )

        if not output_path.exists():
            raise RuntimeError(
                "RGBeta runner did not produce its JSON output.\n"
                f"return code: {completed.returncode}\n"
                f"stdout:\n{completed.stdout}\n"
                f"stderr:\n{completed.stderr}"
            )

        raw_text = output_path.read_text(encoding="utf-8-sig").strip()
        if not raw_text:
            raise RuntimeError(
                "RGBeta runner created an empty JSON output file.\n"
                f"return code: {completed.returncode}\n"
                f"stdout:\n{completed.stdout}\n"
                f"stderr:\n{completed.stderr}"
            )

        try:
            payload = json.loads(raw_text)
        except json.JSONDecodeError as exc:
            raise RuntimeError(
                "RGBeta runner produced invalid JSON.\n"
                f"JSON error: {exc}\n"
                f"output contents:\n{raw_text}\n"
                f"return code: {completed.returncode}\n"
                f"stdout:\n{completed.stdout}\n"
                f"stderr:\n{completed.stderr}"
            ) from exc

        if payload.get("status") != "Success":
            raise RuntimeError(
                failure_message + "\n"
                f"payload: {payload}\n"
                f"stdout:\n{completed.stdout}\n"
                f"stderr:\n{completed.stderr}"
            )

        return result_type(
            status=payload["status"],
            metadata=dict(payload.get("metadata", {})),
            betas=dict(payload.get("betas", {})),
            report_betas=dict(payload.get("report_betas", {})),
            report_latex_betas=dict(payload.get("report_beta_latex", {})),
            raw=payload,
        )


def run_rgbeta_t3(
    d_s1: int,
    d_s2: int,
    d_f: int,
    alpha: int,
    *,
    shared_scalar: bool = False,
    runner_path: Path | None = None,
    wolframscript: str = "wolframscript",
) -> RGBetaT3Result:
    """Run the one-loop renormalisable RGEs for the UV T3 model."""

    runner = (
        Path(runner_path)
        if runner_path is not None
        else Path(__file__).resolve().parent / "RunT3RGBeta.wl"
    )

    return _run_rgbeta(
        d_s1,
        d_s2,
        d_f,
        alpha,
        shared_scalar=shared_scalar,
        runner_path=runner,
        output_name="rgbeta_t3_uv_rge.json",
        failure_message="RGBeta T3 UV-RGE generation failed.",
        result_type=RGBetaT3Result,
        wolframscript=wolframscript,
    )


def run_rgbeta_t3_eft1(
    d_s1: int,
    d_s2: int,
    d_f: int,
    alpha: int,
    *,
    shared_scalar: bool = False,
    runner_path: Path | None = None,
    wolframscript: str = "wolframscript",
) -> RGBetaT3IntermediateResult:
    """Run the renormalisable RGEs in EFT1 after the heavy fermion is removed."""

    runner = (
        Path(runner_path)
        if runner_path is not None
        else Path(__file__).resolve().parent / "RunT3EFT1RGBeta.wl"
    )

    return _run_rgbeta(
        d_s1,
        d_s2,
        d_f,
        alpha,
        shared_scalar=shared_scalar,
        runner_path=runner,
        output_name="rgbeta_t3_eft1_rge.json",
        failure_message="RGBeta T3 EFT1 renormalisable RGE generation failed.",
        result_type=RGBetaT3IntermediateResult,
        wolframscript=wolframscript,
    )
