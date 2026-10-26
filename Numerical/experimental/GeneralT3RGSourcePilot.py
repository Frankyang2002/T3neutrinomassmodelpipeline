"""One-step leading-log pilot for special RGBeta quartics sourced from zero.

This does not replace a coupled, scale-dependent RGE integration and it does
not identify RGBeta coupling normalisations with Matchete couplings.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any
import numpy as np

from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression

RGBETA_ACTIVE = {"lambdaHHdagS2S2", "lambdaHHdagS1barS1bar"}


def _complex(raw: Any) -> complex:
    if not isinstance(raw, dict) or set(raw) != {"real", "imag"}:
        raise ValueError("RGBeta lambdaT3 requires real/imag")
    z = complex(float(raw["real"]), float(raw["imag"]))
    if not np.isfinite([z.real, z.imag]).all():
        raise ValueError("Nonfinite lambdaT3")
    return z


def one_step_source(report: Path, model: str, parameters: dict[str, Any]) -> dict[str, Any]:
    if model not in {"B_m1", "C_m1"}:
        raise ValueError("Only B/C alpha=-1 have checked RGBeta special quartic sources")
    obj = json.loads(Path(report).read_text(encoding="utf-8-sig"))
    if obj.get("status") != "OneLoopSourcesVerified":
        raise ValueError("Formal one-loop RGBeta sources not verified")
    row = next((r for r in obj.get("results", []) if r.get("model") == model), None)
    if row is None or row.get("status") != "Analyzed":
        raise ValueError("Model RGBeta source proof not present")
    lam = _complex(parameters["rgbeta_lambdaT3"])
    yr = {}
    for species in ("y1", "y2"):
        real = np.asarray(parameters[f"{species}_real"], dtype=float)
        imag = np.asarray(parameters[f"{species}_imag"], dtype=float)
        if real.shape != (3, 3) or imag.shape != real.shape:
            raise ValueError(f"{species} must be a 3x3 matrix for this UV RGBeta model")
        yr[species] = real + 1j * imag
        if not np.isfinite(yr[species]).all():
            raise ValueError(f"{species} contains nonfinite entries")
    uv = float(parameters["scale_uv_gev"])
    ir = float(parameters["scale_ir_gev"])
    if uv <= 0 or ir <= 0 or not np.isfinite([uv, ir]).all():
        raise ValueError("RGE scales must be finite and positive")
    hbar = 1 / (16 * np.pi**2)
    env = {"lambdaT3": lam, **yr}
    results = {}
    seen = set()
    for entry in row.get("couplings", []):
        coupling = entry.get("coupling")
        if coupling not in RGBETA_ACTIVE or coupling in seen:
            raise ValueError("Unexpected or duplicate special RGBeta quartic")
        if entry.get("status") != "Analyzed" or entry.get("nonzero_formal_source") is not True or entry.get("source_independent_of_special_quartics") is not True:
            raise ValueError(f"Unverified zero-coupling source: {coupling}")
        source = entry.get("source_at_zero_special_quartics_inputform")
        if not isinstance(source, str) or not source:
            raise ValueError("Empty beta source expression")
        val = evaluate_inputform_expression(source, env)
        if isinstance(val, np.ndarray):
            raise ValueError("Expected scalar RGBeta beta source")
        beta = complex(val)
        delta = hbar * np.log(ir / uv) * beta
        if not np.isfinite([delta.real, delta.imag]).all():
            raise ValueError("Nonfinite induced coupling")
        results[coupling] = {
            "beta_one_loop": {"real": float(beta.real), "imag": float(beta.imag)},
            "delta_lambda_leading_log": {"real": float(delta.real), "imag": float(delta.imag)},
            "initial_lambda": {"real": 0.0, "imag": 0.0},
        }
        seen.add(coupling)
    if seen != RGBETA_ACTIVE:
        raise ValueError("Missing RGBeta source expression")
    return {
        "status": "ExperimentalFrozenBackgroundLeadingLog",
        "model": model,
        "scale_uv_gev": uv,
        "scale_ir_gev": ir,
        "beta_convention": "16*pi^2*d(lambda)/dln(mu)=beta^(1)",
        "results": results,
        "warning": "First-order log with frozen UV backgrounds. No coupled running or Matchete coupling conversion is claimed.",
    }
