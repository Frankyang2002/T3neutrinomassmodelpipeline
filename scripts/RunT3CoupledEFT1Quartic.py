"""Exploratory coupled EFT1 LLSS <-> LLHH scalar-quartic RGE experiment.

Intended repository location: scripts/RunT3CoupledEFT1Quartic.py

Uses the SAME fermion-tree LLSS source, SAME endpoint HHSS tensor, SAME
finite simultaneous-scalar bubble, SAME renormalisable trajectories, and
SAME final SMEFT evolution as the corrected LLSS ablation.

New: integrate the *coupled scalar-quartic sub-block* of the dimension-5
operator system using scipy.integrate.solve_ivp(DOP853):

   16 pi^2 d C_LLHH[fgijAB] / d ln(mu)
       = sum_ab lambda_HHSS[ABab] C_LLSS[fgijab]

   16 pi^2 d C_LLSS[fgijab] / d ln(mu)
       = sum_AB lambda_HHSS[ABab] C_LLHH[fgijAB].

No complete all-seven-term Wilson RGE, gauge/Yukawa contributions to the
Wilson anomalous dimension, HHHH/SSSS quartics, or independently matched
hard Wilson coefficients are claimed. These omitted terms may matter for
a physical result. All scalar tensors remain fixed at their existing
EFT1 endpoint values, exactly as in the previous leading-log contraction.
This runner is diagnostic only and writes to a NEW output directory.

Source-resolved Matchete GammaCC/Weyl factor and B/C scalar-basis
conventions are inherited, without additional certification. There is
no modification of the authoritative runner or existing output files.

From repository root in PowerShell:
  python -X utf8 -m scripts.RunT3CoupledEFT1Quartic --self-test
  python -X utf8 -m scripts.RunT3CoupledEFT1Quartic --source-run thesis16_merged_v1 --run-id trial_A_m4_v1 --models A_-4 --check-only
  python -X utf8 -m scripts.RunT3CoupledEFT1Quartic --source-run thesis16_merged_v1 --run-id trial_A_m4_v1 --models A_-4
"""
from __future__ import annotations

import math
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import solve_ivp

LOOP = 1.0 / (16.0 * math.pi**2)
SCHEMA = "t3-coupled-quartic-only-llss-llhh-ablation-v1"
RTOL = 1.0e-9
ATOL = 1.0e-12  # dimensionless: the state is scaled by ||C_LLSS(muF)||


def integrate_coupled_quartic(
    c_llss: np.ndarray,
    hhss: np.ndarray,
    mu_f: float,
    mu_s: float,
) -> dict:
    """Couple ALL flavour, SU(2) and scalar-component indices in the HHSS block.

    Does not project onto Weinberg at intermediate steps: the entire LLHH
    tensor is retained to avoid throwing away any generated components.
    """
    source = np.asarray(c_llss, dtype=np.complex128)
    q = np.asarray(hhss, dtype=np.float64)
    if source.ndim != 6 or source.shape[:4] != (3, 3, 2, 2):
        raise ValueError(f"Expected source (3,3,2,2,n,n), got {source.shape}")
    n = source.shape[-1]
    if source.shape[-2] != n or q.shape != (4, 4, n, n):
        raise ValueError("LLSS and HHSS scalar dimensions disagree")
    if not (np.isfinite(source).all() and np.isfinite(q).all()):
        raise ValueError("Nonfinite source or scalar quartic")
    if not (math.isfinite(mu_f) and math.isfinite(mu_s) and mu_f > mu_s > 0):
        raise ValueError("Must have muF > muS > 0")

    norm = float(np.linalg.norm(source))
    if norm <= 0.0:
        raise ValueError("Vanishing LLSS boundary: no nontrivial experiment")

    # Dimensionless state avoids a ~1e-13 GeV^-1 absolute solver tolerance issue.
    ss_shape = source.shape
    hh_shape = (3, 3, 2, 2, 4, 4)
    n_ss, n_hh = source.size, int(np.prod(hh_shape))
    initial = np.concatenate(
        ((source / norm).ravel(), np.zeros(n_hh, dtype=np.complex128))
    )

    def rhs(t: float, v: np.ndarray) -> np.ndarray:
        ss = v[:n_ss].reshape(ss_shape)
        hh = v[n_ss:].reshape(hh_shape)
        beta_ss = np.einsum("fgijAB,ABab->fgijab", hh, q, optimize=True)
        beta_hh = np.einsum("fgijab,ABab->fgijAB", ss, q, optimize=True)
        return LOOP * np.concatenate((beta_ss.ravel(), beta_hh.ravel()))

    # Complex-valued solve_ivp is supported by explicit Runge--Kutta methods.
    a, b = math.log(mu_f), math.log(mu_s)
    sol = solve_ivp(rhs, (a, b), initial, method="DOP853",
                    rtol=RTOL, atol=ATOL)
    if not sol.success:
        raise RuntimeError(f"EFT1 DOP853 failed: {sol.message}")
    ss_final = norm * sol.y[:n_ss, -1].reshape(ss_shape)
    hh_final = norm * sol.y[n_ss:, -1].reshape(hh_shape)
    hh_first_order = LOOP * (b - a) * np.einsum(
        "fgijab,ABab->fgijAB", source, q, optimize=True
    )
    return {
        "c_llss_end": ss_final,
        "c_llhh_end": hh_final,
        "c_llhh_first_order": hh_first_order,
        "c_llss_boundary_norm": norm,
        "c_llss_end_norm": float(np.linalg.norm(ss_final)),
        "c_llss_drift_rel": float(np.linalg.norm(ss_final - source) / norm),
        "c_llhh_resummation_rel": float(
            np.linalg.norm(hh_final - hh_first_order)
            / max(np.linalg.norm(hh_first_order), 1e-300)
        ),
        "nfev": int(sol.nfev),
        "steps": int(len(sol.t)),
        "solver_success": bool(sol.success),
        "operator_terms": ["LLSS <-> LLHH from HHSS quartics"],
        "running_quartics": False,
    }


def self_test() -> None:
    """Closed-form two-component subsystem validates ODE, signs and scaling."""
    n = 1
    ss = np.zeros((3, 3, 2, 2, n, n), dtype=complex)
    ss[0, 0, 0, 0, 0, 0] = 1.0e-13
    q = np.zeros((4, 4, n, n), dtype=float)
    q[2, 2, 0, 0] = 0.5
    result = integrate_coupled_quartic(ss, q, 1.0e5, 1.0e3)
    x = LOOP * q[2, 2, 0, 0] * math.log(1.0e3 / 1.0e5)
    got_h = result["c_llhh_end"][0, 0, 0, 0, 2, 2] / 1.0e-13
    got_s = result["c_llss_end"][0, 0, 0, 0, 0, 0] / 1.0e-13
    if not (np.isclose(got_h, math.sinh(x), rtol=2e-8, atol=1e-12)
            and np.isclose(got_s, math.cosh(x), rtol=2e-8, atol=1e-12)):
        raise AssertionError(f"Coupled ODE differs from analytic solution: "
                             f"hh={got_h} vs {math.sinh(x)}, "
                             f"ss={got_s} vs {math.cosh(x)}")
    print("PASS: coupled quartic ODE matches sinh/cosh analytic solution")
    print(f"  leading log={x:.12e}, coupled LLHH={got_h.real:.12e}")
    print(f"  relative LLHH higher-order correction="
          f"{result['c_llhh_resummation_rel']:.6e}")
    print(f"  nfev={result['nfev']}")


def _independent_coupled_cases(
    model, c5_one, source, scalar, target, *,
    uv_on, eft_on, uv_nfev, eft_nfev, spinor_factor, vev,
):
    """Drop-in replacement for RunT3GeneralSequentialFitAblation._final_cases."""
    from Numerical.joint_scalar_matching.core import project_weinberg
    from Numerical.running.GeneralT3SequentialLeadingLog import (
        response, decode, encode,
    )
    from scripts import RunT3GeneralSequentialFitAblation as original

    mu_f, mu_s = original.SCALES[1:3]
    tree = np.asarray(source.raw, complex) * complex(spinor_factor)
    q = np.asarray(scalar.quartic, float)
    full = integrate_coupled_quartic(tree, q, mu_f, mu_s)

    integrated, ward = project_weinberg(full["c_llhh_end"], strict=False)
    fixed, fixed_ward = project_weinberg(full["c_llhh_first_order"], strict=False)
    independent = response(
        source, scalar, mu_f=mu_f, mu_s=mu_s,
        spinor_factor=complex(spinor_factor),
        check_scale_cancellation=True,
    )
    reference_log = decode(independent["c5_eft1_leading_log"])
    finite = decode(independent["c5_scalar_finite_bubble"])
    one_scale = original.c5_checked(c5_one)

    # Normalisation/sign should reproduce the old scalar-quartic beta exactly.
    calibration_rel = float(
        np.linalg.norm(fixed - reference_log)
        / max(np.linalg.norm(reference_log), 1e-300)
    )
    if calibration_rel > 1.0e-9:
        raise ValueError(
            "Coupled EFT1 first-order limit disagrees with existing LLSS log: "
            f"relative discrepancy {calibration_rel:.3e}"
        )

    # Keep the SAME finite bubble in both calculations. This isolates what
    # changed when the Wilson coefficients were evolved self-consistently.
    cases = {}
    for w_on, c5 in ((False, finite), (True, finite + integrated)):
        c5 = original.c5_checked(c5)
        for sm_on in (False, True):
            key = f"{int(uv_on)}{int(eft_on)}{int(w_on)}{int(sm_on)}"
            data = original.low_energy(
                original.scalar_snapshot_state(scalar), c5, on=sm_on,
                mu_s=mu_s, mu_low=original.SCALES[3],
                vev=vev, target=target,
            )
            data.update({
                "case": key, "uv_on": uv_on, "eft1_on": eft_on,
                "wilson_log_on": w_on, "smeft_on": sm_on,
                "uv_nfev": uv_nfev, "eft1_nfev": eft_nfev,
                "c5_at_scalar_threshold": encode(c5),
                "c5_joint_scalar_bubble": encode(finite),
                "c5_induced_wilson_log": encode(integrated),
                "c5_fixed_order_llss_log": encode(reference_log),
                "c5_coupled_minus_fixed_order_log":
                    encode(integrated - reference_log),
                "c5_matching_anchored_reference": encode(one_scale),
                "c5_bookkeeping_residual_matching":
                    encode(np.zeros((3, 3), complex)),
                "c5_log_frobenius": float(np.linalg.norm(integrated)),
                "c5_one_scale_frobenius": float(np.linalg.norm(one_scale)),
                # Historical CSV field: denominator is independent LLSS full
                # coefficient, NOT the old one-scale Matchete anchor.
                "c5_log_over_one_scale_norm": float(
                    np.linalg.norm(integrated)
                    / max(np.linalg.norm(finite + integrated), 1e-300)
                ),
                "coupled_log_over_fixed_order_log_rel": float(
                    np.linalg.norm(integrated - reference_log)
                    / max(np.linalg.norm(reference_log), 1e-300)
                ),
                "first_order_reference_agreement_rel": calibration_rel,
                "coupled_weinberg_projection_residual": ward,
                "fixed_weinberg_projection_residual": fixed_ward,
                "c_llss_boundary_norm": full["c_llss_boundary_norm"],
                "c_llss_endpoint_norm": full["c_llss_end_norm"],
                "c_llss_endpoint_relative_drift": full["c_llss_drift_rel"],
                "solver_coupled_nfev": full["nfev"],
                "solver_coupled_steps": full["steps"],
                "per_general_interaction": independent["per_channel"],
                "quartic_provenance": independent["quartic_provenance"],
                "matching_scale_cancellation":
                    independent["matching_scale_cancellation_residual"],
                "matching_approximation":
                    "exploratory_coupled_quartic_LLSS_LLHH_plus_fixed_bubble",
                "complete_eft1_dimension5_rge": False,
                "hard_fermion_matching_included": False,
                "other_eft1_dim5_operator_classes_included": False,
                "physical_one_loop_c5_certified": False,
            })
            cases[key] = data
    print(
        f"{model} U={int(uv_on)} E={int(eft_on)} : "
        f"coupled/fixed relative change="
        f"{cases[f'{int(uv_on)}{int(eft_on)}11']['coupled_log_over_fixed_order_log_rel']:.6e}, "
        f"LLSS drift={full['c_llss_drift_rel']:.6e}, "
        f"DOP853 nfev={full['nfev']}",
        flush=True,
    )
    return cases


def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    if "--self-test" in args:
        if len(args) != 1:
            raise SystemExit("--self-test cannot be combined with other options")
        self_test()
        return 0

    # Avoid changing any legacy output. Only one new experimental module.
    if "--resume" in args:
        raise ValueError("Use a fresh --run-id (resume disabled)")
    if "--run-id" not in args:
        raise ValueError("Explicit fresh --run-id required")
    from scripts import RunT3GeneralSequentialFitAblation as original

    original._final_cases = _independent_coupled_cases
    original.OUT = (
        original.ROOT
        / "output/studies/general_scalar_extended/"
          "experimental_coupled_eft1_quartic"
    )
    original_save = original.save

    def save_with_provenance(path: Path, data: dict):
        if isinstance(data, dict):
            data = dict(data)
            if data.get("format") == "t3-full-general-16_U_E_W_S_ablation-v1":
                data["format"] = SCHEMA
                data["matching_scope"] = (
                    "coupled_HHSS_quartic_only_in_LLSS_LLHH_subspace"
                )
                data["complete_eft1_dim5_operator_basis"] = False
                data["hard_fermion_matching_included"] = False
            elif data.get("status") == "Computed" and "cases" in data:
                data["matching_approximation"] = (
                    "experimental_coupled_LLSS_LLHH_plus_finite_bubble"
                )
                data["full_sequential_general_C5_certified"] = False
        original_save(path, data)

    original.save = save_with_provenance
    return original.main(args)


if __name__ == "__main__":
    raise SystemExit(main())
