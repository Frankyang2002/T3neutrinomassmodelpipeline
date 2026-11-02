"""Numerical integrator for an explicitly supplied EFT1 Wilson gamma matrix.

This is a *solver*, not a derivation of gamma. The source must supply the
complete ordered-coordinate basis (including flavour axes and all components)
and its anomalous-dimension samples with compatible normalization. Without
that source no Wilson evolution is claimed. All reported runs are conditional
until the matched Matchete-RGBeta operator/coupling dictionary is certified.
"""
from __future__ import annotations

import math
from typing import Any, Mapping

import numpy as np
from scipy.integrate import solve_ivp

LOOP = 1 / (16 * math.pi ** 2)


def _matrix(raw: Mapping[str, Any], n: int) -> np.ndarray:
    r = np.asarray(raw["real"], dtype=float)
    im = np.asarray(raw["imag"], dtype=float)
    if r.shape != (n, n) or im.shape != (n, n):
        raise ValueError("Expected separate real and imaginary NxN gamma matrices")
    a = r + 1j * im
    if not np.isfinite(a).all():
        raise ValueError("Expected finite NxN anomalous-dimension matrix")
    return a


def validate_gamma(payload: Mapping[str, Any], coordinates: list[str],
                   mu_f: float, mu_s: float) -> tuple[np.ndarray, np.ndarray]:
    """Validate explicit source and sample order; no guessed missing entries."""
    if payload.get("format") != "t3-explicit-ordered-wilson-gamma-v1":
        raise ValueError("Unrecognized Wilson anomalous-dimension source format")
    if payload.get("coordinates") != coordinates or len(coordinates) != len(set(coordinates)):
        raise ValueError("Gamma basis does not exactly match ordered tensor coordinates")
    if payload.get("equation") != "16pi2*dC/dlnmu=Gamma(mu)@C":
        raise ValueError("Unverified gamma normalization/sign convention")
    if payload.get("basis") != "ordered_Matchete_real_components_with_flavour":
        raise ValueError("Gamma not in the same ordered scalar/lepton/flavour basis")
    raw = payload.get("samples")
    if not isinstance(raw, list) or len(raw) < 2:
        raise ValueError("At least two gamma samples required")
    mus = np.array([float(x["mu_gev"]) for x in raw])
    if not np.isfinite(mus).all() or np.any(mus <= 0):
        raise ValueError("Invalid gamma sampling scales")
    if not np.all(np.diff(mus) > 0):
        raise ValueError("Gamma sample scales must be strictly increasing")
    if not (mus[0] <= mu_s < mu_f <= mus[-1]):
        raise ValueError("Gamma samples do not cover both EFT1 thresholds")
    arr = np.stack([_matrix(x["gamma"], len(coordinates)) for x in raw])
    return np.log(mus), arr


def evolve_explicit_gamma(boundary: np.ndarray, coordinates: list[str],
                          gamma: Mapping[str, Any], mu_f: float, mu_s: float,
                          *, points: int = 41, rtol: float = 1e-8,
                          atol: float = 1e-11) -> dict:
    """Solve dC/dln(mu)=Gamma(mu) C/(16pi²) on an explicit complete basis."""
    if not (math.isfinite(mu_f) and math.isfinite(mu_s) and mu_f > mu_s > 0):
        raise ValueError("Require finite mu_F > mu_S > 0")
    if points < 2 or points > 10001:
        raise ValueError("2 <= points <= 10001 required")
    v = np.asarray(boundary, dtype=complex)
    if v.shape != (len(coordinates),) or not np.isfinite(v).all():
        raise ValueError("Boundary vector does not match gamma coordinates")
    times, matrices = validate_gamma(gamma, coordinates, mu_f, mu_s)
    initial, final = math.log(mu_f), math.log(mu_s)

    def rhs(t: float, c: np.ndarray) -> np.ndarray:
        hi = int(np.searchsorted(times, t, side="right"))
        hi = min(max(hi, 1), len(times) - 1)
        lo = hi - 1
        alpha = (t - times[lo]) / (times[hi] - times[lo])
        g = matrices[lo] * (1 - alpha) + matrices[hi] * alpha
        return LOOP * (g @ c)

    t_eval = np.linspace(initial, final, points)
    sol = solve_ivp(rhs, (initial, final), v, t_eval=t_eval,
                    rtol=rtol, atol=atol, method="DOP853")
    if not sol.success or sol.y.shape != (len(coordinates), points) or not np.isfinite(sol.y).all():
        raise RuntimeError(f"Wilson ODE solver failed: {sol.message}")
    return {
        "status": "ConditionalExplicitGammaNumericalTransport",
        "physical_matching_certified": False,
        "gamma_generated_by_this_module": False,
        "complete_gamma_verified_by_this_module": False,
        "coordinates": coordinates,
        "scales_gev": np.exp(sol.t).tolist(),
        "start_real": sol.y[:, 0].real.tolist(),
        "start_imag": sol.y[:, 0].imag.tolist(),
        "end_real": sol.y[:, -1].real.tolist(),
        "end_imag": sol.y[:, -1].imag.tolist(),
        "nfev": int(sol.nfev),
        "relative_change": float(np.linalg.norm(sol.y[:, -1] - v) / max(np.linalg.norm(v), 1e-300)),
    }
