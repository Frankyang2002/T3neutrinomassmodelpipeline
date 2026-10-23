"""Order-by-order expansion of the *existing* exploratory T3-B coupled solver.

This is bookkeeping of the one-loop Wilson RGEs on a fixed one-loop
renormalisable background, not a complete O(hbar**2) matching calculation.

Write X = X[0] + X[1] + ... where X[n] contains n powers of the Wilson
loop factor hbar=1/(16*pi**2). The one-loop linear RGE implies
   dX[0]/dt = 0,   dX[n]/dt = hbar * M(t) @ X[n-1], n >= 1.
The C5(MF) heavy matching-log boundary is assigned to X[1].
With C5[0]=0, reverse C5->C12 first enters C12[2] and C5[3].

Renormalisable background dependence remains unexpanded; consequently
comparison against frozen-coupling fixed-order matching can differ at an
order formally beyond the first Wilson-loop order. No physical-NLO claim.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from Numerical.diagnostics.T3BHeavyLogAlignment import (
    flavor_kernel_by_heavy, matrix_from_json, matrix_to_json,
)
from Numerical.diagnostics.T3BExploratoryCoupledBenchmark import (
    FWD, LOOP, _background, _rhs, _sym, fixed_order_log_from_common_boundary,
)


def _ratio(a: np.ndarray, b: np.ndarray, *, floor: float = 1.0e-35) -> float:
    return float(np.linalg.norm(a) / max(float(np.linalg.norm(b)), floor))


def _solve_expansion_pair(samples, y1, y2, masses, *, max_order: int = 4,
                          diagonal: bool = True) -> tuple[tuple[np.ndarray, ...],
                                                           tuple[np.ndarray, ...]]:
    """Evolve reverse-on and reverse-off orders on ONE adaptive time mesh.

    Independent solve_ivp runs may select different meshes because the O(hbar^3)
    reverse term changes the adaptive error estimate. Subtracting their tiny
    O(hbar) and O(hbar^2) outputs can then falsely appear to violate the exact
    selection rule. This combined ODE compares both systems at identical stages.
    """
    if not isinstance(max_order, int) or not 1 <= max_order <= 8:
        raise ValueError("max_order must be an integer from 1 to 8")
    t_low, t_high, bg = _background(samples)
    mf = float(samples[0]["mu_gev"])
    ms = float(samples[-1]["mu_gev"])
    masses = np.asarray(masses, dtype=float)
    if (masses.shape != (3,) or not np.all(np.isfinite(masses))
            or np.any(masses <= ms)):
        raise ValueError("Require three finite heavy masses above MS")
    source = flavor_kernel_by_heavy(y1, y2, masses)
    if not np.all(np.isfinite(source)):
        raise ValueError("Non-finite heavy-state kernel")

    lam = complex(bg(t_low)["lambdaT3"])
    boundary = np.zeros((max_order + 1, 18), dtype=complex)
    boundary[0, :9] = _sym(source.sum(axis=0)).reshape(9)
    # Matching-log term starts at Wilson order one, in ordered C5 convention.
    boundary[1, 9:] = (
        _sym(FWD * lam / LOOP *
             np.einsum("r,rij->ij", np.log(mf / masses), source))
    ).reshape(9)
    # Two coupled Taylor towers, sharing one ODE integration grid.
    initial = np.stack((boundary, boundary))

    def rhs(t: float, flat: np.ndarray) -> np.ndarray:
        state = flat.reshape(2, max_order + 1, 18)
        deriv = np.zeros_like(state)
        for n in range(1, max_order + 1):
            deriv[0, n] = _rhs(t, state[0, n - 1], bg,
                                diagonal=diagonal, reverse=True)
            deriv[1, n] = _rhs(t, state[1, n - 1], bg,
                                diagonal=diagonal, reverse=False)
        return deriv.reshape(-1)

    solution = solve_ivp(
        rhs, (t_high, t_low), initial.reshape(-1),
        method="DOP853", rtol=2e-11, atol=1e-25,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    endpoint = solution.y[:, -1].reshape(2, max_order + 1, 18)
    return tuple(
        tuple(2.0 * _sym(item[9:].reshape(3, 3)) for item in tower)
        for tower in endpoint
    )


def _solve_expansion(samples, y1, y2, masses, *, max_order: int = 4,
                     reverse: bool = True, diagonal: bool = True) -> tuple[np.ndarray, ...]:
    """Backward-compatible single-branch interface for existing tests/callers."""
    both = _solve_expansion_pair(
        samples, y1, y2, masses, max_order=max_order, diagonal=diagonal,
    )
    return both[0 if reverse else 1]

def perturbative_accuracy(control: dict, coupled: dict, y1: np.ndarray,
                          y2: np.ndarray, masses: np.ndarray, *,
                          max_order: int = 4) -> dict:
    """Check hierarchy of Wilson-loop orders against the saved full solve.

    The only asserted physics equality is the O(hbar) *frozen coupling*
    direct matching/log result. The resummed full solve is compared solely
    as a numerical convergence test, not as a physical correctness proof.
    """
    if control.get("status") != "FixedOrderControlComplete":
        raise ValueError("Expected successful fixed-order control")
    if coupled.get("status") != "ExploratoryCoupledMatrixComplete":
        raise ValueError("Expected exploratory coupled benchmark output")
    if coupled.get("threshold_scheme") != "common_scale_log_boundary":
        raise ValueError("Only common_scale_log_boundary is supported")
    samples = control["intermediate_running"]["samples"]
    if not np.isclose(float(samples[0]["mu_gev"]), float(control["mu_f_gev"]), rtol=1e-12):
        raise ValueError("Saved MF scale disagrees with the control")
    if not np.isclose(float(samples[-1]["mu_gev"]), float(control["mu_s_gev"]), rtol=1e-12):
        raise ValueError("Saved MS scale disagrees with the control")
    direct = matrix_from_json(coupled["production_direct_c5_at_ms"])
    predicted = fixed_order_log_from_common_boundary(samples, y1, y2, masses)
    direct_anchor_error = _ratio(predicted - direct, direct)
    if direct_anchor_error > 1e-8:
        raise ValueError("Frozen O(hbar) direct term does not agree with production")

    orders, no_reverse = _solve_expansion_pair(
        samples, y1, y2, masses, max_order=max_order,
    )
    full = matrix_from_json(coupled["coupled_eft1_c5_at_ms"])
    partial_sums = np.cumsum(np.asarray(orders), axis=0)
    baseline = matrix_from_json(control["fixed_order_c5_at_ms"])
    if np.linalg.norm(full) == 0 or np.linalg.norm(baseline) == 0:
        raise ValueError("Zero Weinberg norm")

    # Mathematically C5[1] and C5[2] cannot depend on reverse mixing
    # when C5[0]=0; reverse influence starts in C5[3].
    discrepancy_first_two = [
        _ratio(orders[n] - no_reverse[n], orders[1])
        for n in (1, 2) if n <= max_order
    ]
    if any(x > 1e-8 for x in discrepancy_first_two):
        raise RuntimeError("Reverse mixing contaminated C5 at order < 3")
    reconstruction = _ratio(partial_sums[-1] - full, full)
    return {
        "status": "PerturbativeLLBookkeepingComplete",
        "truncated_series_convergence_diagnostic_passed": bool(reconstruction < 1e-5),
        "truncated_series_convergence_tolerance": 1e-5,
        "scheme": "common_scale_simultaneous_fermion_matching",
        "powers_counted": "Wilson loop factor hbar=1/(16*pi^2), with one-loop renormalisable backgrounds held fixed",
        "max_wilson_loop_order": max_order,
        "one_loop_frozen_direct_production_residual": direct_anchor_error,
        "order1_vs_production_direct_relative_to_direct": _ratio(orders[1] - direct, direct),
        "order1_background_effect": "Difference arises from evaluating running lambdaT3(t) rather than freezing it at MS; higher order in renormalisable running",
        "terms": {
            f"order_{n}": {
                "c5_physical_majorana": matrix_to_json(orders[n]),
                "norm_relative_to_fixed_c5": _ratio(orders[n], baseline),
                "norm_relative_to_order1": _ratio(orders[n], orders[1]),
                "partial_sum_relative_to_full_coupled": _ratio(partial_sums[n] - full, full),
                "reverse_contribution_relative_to_order1": _ratio(orders[n] - no_reverse[n], orders[1]),
            } for n in range(1, max_order + 1)
        },
        "resummed_c5_reconstruction_residual": reconstruction,
        "order1_plus_order2_vs_resummed_relative_to_full": _ratio(partial_sums[2] - full, full) if max_order >= 2 else None,
        "reverse_can_affect_c5_first_at_wilson_order": 3,
        "reverse_numerical_difference_at_orders_1_and_2": discrepancy_first_two,
        "order2_complete_nlo_matching": False,
        "validated_physical_prediction_available": False,
        "production_replacement_allowed": False,
        "missing_at_same_formal_second_order": [
            "One-loop corrections to C12 threshold matching (then one-loop C12->C5)",
            "Running of the one-loop hard C5 matching term below MF",
            "Two-loop C12->C5 anomalous dimension and finite scalar-threshold terms",
            "Consistent matching-scale dependence and boundary coupling convention at higher order",
            "Complete three-generation closure beyond the retained two-operator system",
        ],
        "interpretation": (
            "One-loop Wilson anomalous dimensions with leading logarithms resummed; "
            "not a complete O(hbar**2) fixed-order or NLO neutrino-mass prediction"
        ),
    }
