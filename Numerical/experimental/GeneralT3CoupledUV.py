"""Coupled UV one-loop RGBeta evolution for ordinary B/C(alpha=-1).

Runs the *existing* Numerical.running.UVRunner on the model's full RGBeta
report_betas, including all six already-supported complex special quartics.
It DOES NOT apply an unverified RGBeta-to-Matchete coupling conversion and
DOES NOT continue through the fermion or scalar thresholds.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping
import json

import numpy as np

MODELS = {"B_m1": 1, "C_m1": 3}
SPECIAL = ("lambdaHHdagS2S2", "lambdaHHdagS1barS1bar")
ALL_COMPLEX = (
    "lambdaHHdagS2S2", "lambdaHHdagS1barS1bar", "lambdaS1bar2S2bar2",
    "lambdaS1barS2S2bar2", "lambdaS1S1bar2S2bar", "lambdaHHdagS1barS2barCross",
)


def _complex(value: Any, label: str) -> complex:
    if not isinstance(value, Mapping) or set(value) != {"real", "imag"}:
        raise ValueError(f"{label} must be a real/imag object")
    z = complex(float(value["real"]), float(value["imag"]))
    if not np.isfinite([z.real, z.imag]).all():
        raise ValueError(f"{label} must be finite")
    return z


def _matrix(pars: Mapping[str, Any], key: str) -> np.ndarray:
    a = np.asarray(pars[key+"_real"], dtype=float)
    b = np.asarray(pars[key+"_imag"], dtype=float)
    if a.shape != (3,3) or b.shape != a.shape:
        raise ValueError(f"{key} must be a 3x3 complex matrix")
    z = a+1j*b
    if not np.isfinite(z).all():
        raise ValueError(f"{key} must have finite entries")
    return z


def construct_uv_state(model: str, parameters: Mapping[str, Any]):
    """Explicitly construct toy UV state using existing numerical conventions."""
    from Numerical.core.State import SMNumericalState, T3Representation, T3UVState
    if model not in MODELS:
        raise ValueError("Only ordinary B/C alpha=-1 are supported")
    heavy = np.asarray(parameters["MF"], dtype=float)
    if heavy.shape != (3,) or np.any(heavy <= 0) or not np.isfinite(heavy).all():
        raise ValueError("Three positive finite heavy masses are required")
    uv = float(parameters["scale_uv_gev"])
    if not np.isfinite(uv) or uv <= 0:
        raise ValueError("Positive finite UV scale required")
    ms1 = float(parameters["MS1"])
    ms2 = float(parameters["MS2"])
    if min(ms1, ms2) <= 0 or not np.isfinite([ms1, ms2]).all():
        raise ValueError("Positive scalar masses required")
    # Explicit illustrative SM inputs at the UV scale, not a fitted SM run.
    sm_config = parameters.get("uv_sm", {})
    if not isinstance(sm_config, Mapping):
        raise ValueError("uv_sm must be an object")
    def diag(name: str, default: tuple[float, float, float]):
        values = np.asarray(sm_config.get(name, default), dtype=float)
        if values.shape != (3,) or not np.isfinite(values).all():
            raise ValueError(f"uv_sm.{name} must have three finite diagonal entries")
        return np.diag(values.astype(complex))
    sm = SMNumericalState(
        gY=float(sm_config.get("gY", 0.36)),
        g2=float(sm_config.get("g2", 0.63)),
        g3=float(sm_config.get("g3", 1.01)),
        lambdaH=float(sm_config.get("lambdaH", 0.13)),
        yu=diag("yu", (1e-5, 0.004, 0.83)),
        yd=diag("yd", (2e-5, 0.0004, 0.018)),
        ye=diag("ye", (3e-6, 0.0006, 0.010)),
    )
    background = parameters.get("uv_scalar_quartics", {})
    if not isinstance(background, Mapping):
        raise ValueError("uv_scalar_quartics must be an object")
    def quartic(name: str, default: float):
        return float(background.get(name, default))
    optional = {name: 0.0j for name in ALL_COMPLEX}
    supplied = parameters.get("uv_rgbeta_special_quartics", {})
    if not isinstance(supplied, Mapping) or set(supplied)-set(ALL_COMPLEX):
        raise ValueError("Unknown optional uv_rgbeta_special_quartics")
    for name, value in supplied.items():
        optional[name] = _complex(value, name)
    state = T3UVState(
        mu_gev=uv,
        representation=T3Representation(2, 2, MODELS[model], -1),
        sm=sm,
        y1=_matrix(parameters,"y1"), y2=_matrix(parameters,"y2"),
        MF=np.diag(heavy.astype(complex)),
        mS1Sq=ms1**2, mS2Sq=ms2**2,
        lambdaS1=quartic("lambdaS1", 0.08),
        lambdaS2=quartic("lambdaS2", 0.08),
        lambdaH1=quartic("lambdaH1", 0.06),
        lambdaH2=quartic("lambdaH2", 0.06),
        lambda12=quartic("lambda12", 0.03),
        lambdaT3=_complex(parameters["rgbeta_lambdaT3"], "rgbeta_lambdaT3"),
        lambdaH1Adj=quartic("lambdaH1Adj", 0.02),
        lambdaH2Adj=quartic("lambdaH2Adj", 0.02),
        lambda12Adj=quartic("lambda12Adj", 0.01),
        **optional,
    )
    return state.validated()


def _encode(z: complex) -> dict[str, float]:
    v = complex(z)
    return {"real": float(v.real), "imag": float(v.imag)}


def run_coupled_uv(
    model: str, parameters: Mapping[str, Any], rgbeta_path: Path,
    source_proof_path: Path, *, n_save: int = 11,
) -> dict[str, Any]:
    """Evaluate all UV RGBeta beta functions via the production ODE runner."""
    from Numerical.running.UVRunner import run_uv_segment
    from Numerical.experimental.GeneralT3RGSourcePilot import one_step_source
    from Numerical.core.RGBetaEvaluator import load_rgbeta_payload

    state = construct_uv_state(model, parameters)
    endpoint = float(parameters.get("uv_coupled_endpoint_gev", 3500.0))
    if not np.isfinite(endpoint) or endpoint <= 0:
        raise ValueError("Invalid UV-only endpoint")
    if not (state.mu_gev > endpoint > 1.001*np.max(np.real(np.diag(state.MF)))):
        raise ValueError("UV-only endpoint must be above every heavy initial fermion mass and below UV scale")
    if n_save < 2 or n_save > 2001:
        raise ValueError("n_save must be between 2 and 2001")
    payload = load_rgbeta_payload(Path(rgbeta_path))
    # Need exact exported model/representation and beta convention. The existing
    # production evaluator validates metadata and the entire state block layout.
    scales = np.geomspace(state.mu_gev, endpoint, n_save)
    trajectory = run_uv_segment(
        state, payload, endpoint, save_scales_gev=scales,
        rtol=3e-8, atol=3e-12,
    )
    final = trajectory.final_state
    if min(final.mS1Sq, final.mS2Sq) <= 0:
        raise ValueError("Negative scalar mass-squared at UV endpoint")

    # Re-evaluate the original frozen-background calculation over the EXACT
    # same interval. No Matchete values enter this RG-only comparison.
    leading_params = dict(parameters)
    leading_params["scale_ir_gev"] = endpoint
    leading = one_step_source(Path(source_proof_path), model, leading_params)

    rows = []
    for name in SPECIAL:
        z0 = complex(getattr(state, name))
        z1 = complex(getattr(final, name))
        row = {
            "coupling": name,
            "uv_initial_rgbeta": _encode(z0),
            "uv_endpoint_rgbeta": _encode(z1),
            "change_coupled": _encode(z1-z0),
            "leading_log_zero_special_background": leading["results"][name]["delta_lambda_leading_log"],
        }
        rows.append(row)
    others = {name: _encode(getattr(final, name)) for name in ALL_COMPLEX if name not in SPECIAL}
    sampled = []
    for i in range(trajectory.n_points):
        point = trajectory.state_at_index(i)
        sampled.append({
            "mu_gev": float(trajectory.mu_gev[i]),
            **{name: _encode(getattr(point, name)) for name in SPECIAL},
        })
    return {
        "status": "CoupledUVOneLoopRGBetaOnly",
        "model": model,
        "source": str(rgbeta_path),
        "uv_scale_gev": float(state.mu_gev),
        "endpoint_gev": endpoint,
        "solver_success": trajectory.solver_success,
        "solver_nfev": trajectory.nfev,
        "sampled": sampled,
        "active_couplings": rows,
        "other_special_quartics_at_endpoint": others,
        "scalar_masses_sq_at_endpoint_gev2": {
            "mS1Sq": float(final.mS1Sq), "mS2Sq": float(final.mS2Sq)
        },
        "full_matchete_rgbeta_coupling_mapping_verified": False,
        "includes_mass_threshold_matching": False,
        "includes_neutrino_mass_matching": False,
        "note": "This runs the EXISTING complete B/C(-1) UV RGBeta beta system with coupled backgrounds down to a scale above all heavy fermions. The leading-log row freezes UV backgrounds. Neither output is inserted into the Matchete C5 without a vertex-level coupling-convention check.",
    }
