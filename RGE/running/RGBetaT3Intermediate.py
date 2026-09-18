from __future__ import annotations

import json
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

# Note that frozen means results cannot be modified after creation
@dataclass(frozen=True)
class RGBetaT3IntermediateResult:
    """Result returned by the Wolfram RGBeta intermediate-EFT runner."""

    status: str
    metadata: dict[str, Any]
    betas: dict[str, str]
    report_betas: dict[str, str]
    report_latex_betas: dict[str, str]
    raw: dict[str, Any]


def _default_runner_path() -> Path:
    """This finds the Wolfram Runner automatically"""
    return Path(__file__).resolve().parent / "wolfram" / "RunT3EFT1RGBeta.wl"


def _wolfram_integer_token(value: int) -> str:
    """Make negative values be -1 -> m1 instead."""
    return f"m{abs(value)}" if value < 0 else str(value)


def run_rgbeta_t3_eft1(
    d_s1: int,
    d_s2: int,
    d_f: int,
    alpha: int,
    *,
    shared_scalar: bool = False,
    runner_path: Path | None = None,
    wolframscript: str = "wolframscript",
) -> RGBetaT3Result:
    """Run the renormalisable RGEs in EFT1 = SM + S1 + S2 after F is removed."""

    # RGBeta seems to not support higher dimensions so we restrict dimensions for RGBeta
    if any(d not in {1, 2, 3} for d in (d_s1, d_s2, d_f)):
        raise ValueError(
            "RGBeta T3 EFT1 running currently supports only SU(2) dimensions 1, 2 and 3."
        )

    # Note that RunT3RGBeta.wl is the real runner
    runner = Path(runner_path) if runner_path is not None else _default_runner_path()

    if not runner.exists():
        raise FileNotFoundError(f"RGBeta Wolfram runner not found: {runner}")

    
    with tempfile.TemporaryDirectory(prefix="t3_rgbeta_") as tmpdir:
        output_path = Path(tmpdir) / "rgbeta_t3_eft1_rge.json"

        # We use the runner to run 
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
                "RGBeta T3 EFT1 renormalisable RGE generation failed.\n"
                f"payload: {payload}\n"
                f"stdout:\n{completed.stdout}\n"
                f"stderr:\n{completed.stderr}"
            )

        return RGBetaT3IntermediateResult(
            status=payload["status"],
            metadata=dict(payload.get("metadata", {})),
            betas=dict(payload.get("betas", {})),
            report_betas=dict(payload.get("report_betas", {})),
            report_latex_betas=dict(payload.get("report_beta_latex", {})),
            raw=payload,
        )
