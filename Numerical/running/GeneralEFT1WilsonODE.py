"""Numerically integrate the existing EFT1 psi^2 phi^2 master Wilson RGE.

This is a numerical engine for the *specified* master-basis operator class,
NOT a physical T3 boundary adapter and NOT a proof of the complete EFT1 operator
basis. It reuses RGE.running.intermediate.GeneralEFT1C5Mixing.calculate_master_beta
without altering any of its seven terms.

Inputs are canonical master-Weyl boundary JSON objects accepted by
load_boundary. A sequence of canonical coupling snapshots, each with the same
field/component dictionary but its own mu_gev, supplies running EFT1
renormalisable couplings. Couplings are linearly interpolated in ln(mu); supply
sufficient snapshots for an accurate interpolation. The actual coefficient
state is integrated jointly, including components induced radiatively.

The output C5 is the Weinberg component generated during EFT1, *not* an
additional contribution to append to an existing fixed-order LL result.
Finite scalar matching is deliberately NOT supplied by this module.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from math import log, pi
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.integrate import solve_ivp

from RGE.running.intermediate.GeneralEFT1C5Mixing import (
    BoundaryError, Context, TERMS, calculate_master_beta, encode_matrix,
    load_boundary, project_c5,
)

LOOP = 1.0 / (16.0 * pi * pi)
FORMAT = "t3-master-psi2phi2-numerical-running-v1"


def _check_context_geometry(reference: Context, other: Context) -> None:
    if (reference.model != other.model or reference.nf != other.nf or
        reference.ns != other.ns or reference.leptons != other.leptons or
        reference.higgs != other.higgs or len(reference.gauge) != len(other.gauge)):
        raise BoundaryError("Changing EFT1 operator or field dictionary across snapshots")
    for (g0, a), (g1, b) in zip(reference.gauge, other.gauge):
        if len(a) != len(b):
            raise BoundaryError("Changing gauge generator count across snapshots")
        for (f0, s0), (f1, s1) in zip(a, b):
            if not (np.array_equal(f0, f1) and np.array_equal(s0, s1)):
                raise BoundaryError("Gauge generator basis changes across snapshots")
    if [name for name, _ in reference.operators] != [name for name, _ in other.operators]:
        raise BoundaryError("Changing Wilson operator source dictionary across snapshots")


def _complex_norm(x: np.ndarray) -> float:
    return float(np.linalg.norm(x.ravel()))


def solve_coupled_wilson(boundary: Mapping[str, Any],
                         snapshots: Sequence[Mapping[str, Any]], *,
                         allow_synthetic: bool = False,
                         rtol: float = 1e-7,
                         atol: float = 1e-16,
                         max_step_log: float = 0.2,
                         save_points: int = 21) -> dict[str, Any]:
    """Integrate coupled C_ijab in log(mu), including the generated HH entries.

    `boundary` and each `snapshot` must be a full canonical input accepted by
    load_boundary; each snapshot additionally contains `mu_gev`. Snapshot
    operator coefficients are checked for identity only, never used as new
    Wilson initial conditions. Ordered snapshots must span [muS, muF].

    At one-loop in the Wilson beta itself the generated Weinberg component
    automatically includes gamma_55*C5; adding a separate leading log is wrong.
    """
    ctx = load_boundary(boundary, allow_synthetic=allow_synthetic)
    if not (len(snapshots) >= 2 and 2 <= save_points <= 1001 and
            np.isfinite(rtol) and rtol > 0 and np.isfinite(atol) and atol > 0 and
            np.isfinite(max_step_log) and max_step_log > 0):
        raise ValueError("Require >=2 coupling snapshots and positive solver controls")

    points = []
    for k, raw in enumerate(snapshots):
        if not isinstance(raw, Mapping):
            raise BoundaryError("Each trajectory snapshot must be a JSON object")
        mu = float(raw.get("mu_gev", float("nan")))
        if not np.isfinite(mu) or mu <= 0:
            raise BoundaryError(f"Invalid coupling snapshot scale at {k}")
        other = load_boundary(raw, allow_synthetic=allow_synthetic)
        _check_context_geometry(ctx, other)
        points.append((log(mu), other))
    points.sort(key=lambda pair: pair[0])
    ts = np.array([p[0] for p in points], float)
    if np.any(np.diff(ts) <= 0):
        raise BoundaryError("Duplicate coupling scales")
    t0, t1 = log(ctx.mu_f), log(ctx.mu_s)
    if ts[0] > t1 + 1e-12 or ts[-1] < t0 - 1e-12:
        raise BoundaryError("Coupling trajectory does not span mu_F to mu_S")
    if abs(ts[0]-t1)>1e-10 or abs(ts[-1]-t0)>1e-10:
        raise BoundaryError("Require endpoint coupling snapshots exactly at mu_S and mu_F")

    # The interpolation is on the actual running tensor data, not endpoint
    # freezing. The fixed representation generators must not be interpolated.
    yukawas = np.stack([item.yukawa for _, item in points])
    quartics = np.stack([item.quartic for _, item in points])
    gs = np.asarray([[g for g, _ in item.gauge] for _, item in points], float)
    source_tensors = [c for _, c in ctx.operators]
    c0 = sum(source_tensors, np.zeros((ctx.nf, ctx.nf, ctx.ns, ctx.ns), complex))
    shape = c0.shape
    n = c0.size
    if not np.isfinite(c0).all():
        raise BoundaryError("Nonfinite initial Wilson coefficients")
    if not np.allclose(c0, c0.swapaxes(0,1), rtol=1e-11, atol=1e-25) or not np.allclose(
            c0, c0.swapaxes(2,3), rtol=1e-11, atol=1e-25):
        raise BoundaryError("Master Wilson boundary must be pairwise symmetric")

    def evaluate(t: float) -> Context:
        i = min(max(int(np.searchsorted(ts, t, side="right") - 1), 0), len(ts)-2)
        fraction = (t-ts[i])/(ts[i+1]-ts[i])
        y = yukawas[i] + fraction*(yukawas[i+1]-yukawas[i])
        q = quartics[i] + fraction*(quartics[i+1]-quartics[i])
        couplings = gs[i] + fraction*(gs[i+1]-gs[i])
        sectors = tuple((float(g), ctx.gauge[j][1]) for j, g in enumerate(couplings))
        return replace(ctx, yukawa=y, quartic=q, gauge=sectors)

    def unpack(v: np.ndarray) -> np.ndarray:
        return (v[:n] + 1j*v[n:]).reshape(shape)

    def rhs(t: float, v: np.ndarray) -> np.ndarray:
        c = unpack(v)
        beta = sum(calculate_master_beta(evaluate(t), c).values(), np.zeros(shape, complex))
        if not np.isfinite(beta).all():
            raise FloatingPointError("Nonfinite Wilson beta during EFT1 evolution")
        deriv = LOOP*beta.ravel()
        return np.concatenate((deriv.real, deriv.imag))

    y0 = np.concatenate((c0.real.ravel(), c0.imag.ravel()))
    teval = np.linspace(t0, t1, save_points)
    solution = solve_ivp(rhs, (t0, t1), y0, method="DOP853", t_eval=teval,
                         rtol=rtol, atol=atol, max_step=max_step_log)
    if not solution.success or not np.isfinite(solution.y).all():
        raise RuntimeError(f"EFT1 Wilson ODE failed: {solution.message}")

    final = unpack(solution.y[:,-1])
    c5_initial, _, _ = project_c5(c0, ctx)
    c5_final, ward, symmetric = project_c5(final, ctx)
    beta0 = sum(calculate_master_beta(evaluate(t0), c0).values(), np.zeros(shape, complex))
    c5_beta0, _, _ = project_c5(beta0, ctx)
    leading_log = LOOP * (t1-t0) * c5_beta0
    trajectory = []
    for k, t in enumerate(solution.t):
        matrix, projection_residual, symmetry_residual = project_c5(
            unpack(solution.y[:,k]), ctx)
        trajectory.append({"mu_gev": float(np.exp(t)),
                           "C5_component_gev_inverse": encode_matrix(matrix),
                           "weinberg_projection_residual": projection_residual,
                           "flavour_symmetry_residual": symmetry_residual})
    return {
        "format": FORMAT, "model": ctx.model,
        "authority": ctx.authority,
        "status": "MasterPsi2Phi2NumericalSolved_NotCompleteEFT1OperatorBasis",
        "scales_gev": {"mu_F":ctx.mu_f, "mu_S":ctx.mu_s},
        "solver": {"success":bool(solution.success), "message":str(solution.message),
                   "method":"DOP853", "nfev":int(solution.nfev),
                   "rtol":rtol, "atol":atol, "max_step_log":max_step_log},
        "operator_shape":list(shape), "operator_sources":[name for name,_ in ctx.operators],
        "running_couplings_from_snapshots":len(points),
        "wilson_boundary_norm":_complex_norm(c0),
        "wilson_endpoint_norm":_complex_norm(final),
        "wilson_pair_symmetry_residual":max(
            float(np.max(np.abs(final-final.swapaxes(0,1)))),
            float(np.max(np.abs(final-final.swapaxes(2,3))))),
        "C5_initial_component_gev_inverse":encode_matrix(c5_initial),
        "C5_at_muS_component_gev_inverse":encode_matrix(c5_final),
        "C5_change_numerical_gev_inverse":encode_matrix(c5_final-c5_initial),
        "C5_fixed_couplings_leading_log_gev_inverse":encode_matrix(leading_log),
        "C5_numerical_minus_LL_norm":_complex_norm((c5_final-c5_initial)-leading_log),
        "weinberg_projection_residual":ward,
        "flavour_symmetry_residual":symmetric,
        "trajectory":trajectory,
        "scope": {
            "rge_equation":"Aebischer et al. 2512.16890 eq (4.85), existing seven-term code",
            "integrates_entire_psi2phi2_master_tensor":True,
            "renormalisable_couplings_are_running":True,
            "other_dimension_five_operator_classes_closed":False,
            "physical_matchete_to_master_normalization_certified":False,
            "scalar_matching_performed":False,
            "one_loop_hard_matching_performed":False,
            "add_existing_fixed_order_leading_log":False,
        },
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--boundary", type=Path, required=True)
    ap.add_argument("--coupling-trajectory", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--allow-synthetic", action="store_true")
    ap.add_argument("--rtol", type=float, default=1e-7)
    ap.add_argument("--atol", type=float, default=1e-16)
    ap.add_argument("--max-step-log", type=float, default=0.2)
    a = ap.parse_args(argv)
    if a.output.exists():
        raise FileExistsError(f"Refusing to overwrite historical output: {a.output}")
    boundary = json.loads(a.boundary.read_text(encoding="utf-8-sig"))
    payload = json.loads(a.coupling_trajectory.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, list):
        raise BoundaryError("Coupling trajectory JSON must be a list of canonical snapshots")
    result = solve_coupled_wilson(boundary, payload, allow_synthetic=a.allow_synthetic,
                                  rtol=a.rtol, atol=a.atol, max_step_log=a.max_step_log)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(f"{result['model']}: master Wilson ODE succeeded; nfev={result['solver']['nfev']}; {a.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
