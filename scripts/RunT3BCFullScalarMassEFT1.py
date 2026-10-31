"""Opt-in B(-1)/C(-1) coupled EFT1 running with the full scalar mass sector.

Reuse Stage33 fermion-threshold candidate boundaries and the COMPLETE
all-quartic EFT1 RGBeta beta export, replacing all older scalar mass betas
with separately verified native Higgs+BSM mass betas. All dimensionless
couplings and five real mass coordinates evolve in the SAME solve_ivp state.

This is an illustrative native-RGBeta numerical study, NOT certified
Matchete/RGBeta matching, scalar eigenstate decoupling, or neutrino physics.
No source beta exports are created, regenerated, or modified.

Run from repository root:
  python -X utf8 -m scripts.RunT3BCFullScalarMassEFT1 --all \
      --higgs-mass-squared -7800 --mu2-fraction 0.01 --mu2-phase 0.6
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np
from scipy.integrate import solve_ivp

from scripts import RunGeneralT3QuadraticEFT1 as old
from scripts import FermionThresholdHandoff as s33
from scripts import IntermediateEFTRunning as s34
from scripts import VerifyT3BCFullScalarMassClosure as closure
from Numerical.core.BetaVector import beta_values_to_derivative
from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression
from Numerical.experimental.GeneralT3ScalarRGState import (
    GeneralT3ScalarRGState, load_registry, pack_general_scalar_state,
    registry_model, unpack_general_scalar_state,
)
from Numerical.experimental.GeneralT3ScalarQuadraticMasses import (
    analytic_mass_squared, load_manifest, quadratic_kind, quadratic_mass_hessian,
)
from Numerical.running.IntermediateScalarRGBetaEvaluator import (
    intermediate_scalar_state_environment,
)

ROOT = old.ROOT
MODELS = ("B_-1", "C_-1")
DIMENSIONLESS_DIR = old.OUTPUT
MASS_DIR = ROOT / "output/audits/workstream_b_bc_full_mass"
OUTPUT_DIR = ROOT / "output/workstream_b_bc_coupled_full_mass_eft1"
FORMAT = "workstream-b-bc-coupled-five-real-mass-eft1-v1"
MASS_NAMES = closure.NAMES
MIXED = "muT3GeneralS1S2"
LOOP = 16.0 * math.pi**2


def require(condition: bool, reason: str) -> None:
    if not condition:
        raise ValueError("B/C coupled five-real-mass EFT1: " + reason)


def merge_beta_reports(old_betas: Mapping[str, str],
                       new_mass_betas: Mapping[str, str],
                       layout_names: tuple[str, ...] | list[str]) -> dict[str, str]:
    """Retain ALL old dimensionless betas, replace ALL four mass betas.

    The baseline layout includes the two diagonal BSM masses but not the
    independent complex mixed bilinear or the Higgs quadratic mass.
    """
    expected_old = set(layout_names) | {MIXED}
    require(set(old_betas) == expected_old,
            f"old beta inventory mismatch: {sorted(set(old_betas) ^ expected_old)}")
    require(set(new_mass_betas) == set(MASS_NAMES),
            "new live mass beta inventory must be exactly four symbols")
    require({"mS1Sq", "mS2Sq"}.issubset(set(layout_names)),
            "baseline numerical state is missing BSM diagonal masses")
    require("mHSq" not in layout_names and MIXED not in layout_names,
            "baseline layout unexpectedly contains an extra mass coordinate")
    require(all(isinstance(v, str) and v.strip() for v in
                list(old_betas.values()) + list(new_mass_betas.values())),
            "non-string or empty beta expression")
    merged = {key: expr for key, expr in old_betas.items()
              if key not in MASS_NAMES}
    merged.update(new_mass_betas)
    require(set(merged) == set(layout_names) | {"mHSq", MIXED},
            "merged beta inventory does not cover the coupled ODE state")
    for key in MASS_NAMES:
        require(merged[key] == new_mass_betas[key],
                f"outdated mass beta survived for {key}")
    return merged


def _real(name: str, value: Any) -> float:
    z = complex(value)
    require(math.isfinite(z.real) and math.isfinite(z.imag),
            f"nonfinite {name} beta")
    # Exact result is real by Hermiticity; tolerate only roundoff generated
    # while evaluating terms of the form z*conjugate(w)+conjugate(z)*w.
    require(abs(z.imag) <= 1.0e-9 * max(1.0, abs(z.real)),
            f"unexpected imaginary component in {name} beta: {z.imag}")
    return float(z.real)


def _complex(name: str, value: Any) -> complex:
    z = complex(value)
    require(math.isfinite(z.real) and math.isfinite(z.imag),
            f"nonfinite complex beta for {name}")
    return z


def beta_derivative(vector: np.ndarray, *, mu_gev: float, layout: Any,
                    registry: dict, betas: Mapping[str, str]) -> np.ndarray:
    """Pack [existing baseline, Re(mu), Im(mu), mHSq] into one derivative."""
    y = np.asarray(vector, dtype=float)
    n = layout.size
    require(y.shape == (n + 3,) and np.isfinite(y).all(),
            "wrong size or nonfinite coupled ODE vector")
    require(math.isfinite(mu_gev) and mu_gev > 0, "invalid RG scale")
    state = unpack_general_scalar_state(y[:n], layout, registry,
                                        uv=False, mu_gev=mu_gev)
    env = intermediate_scalar_state_environment(state.baseline)
    env.update(state.extra)
    env[MIXED] = complex(float(y[n]), float(y[n+1]))
    env["mHSq"] = float(y[n+2])
    # One common environment, including the live Higgs/mixed masses, is used
    # for every dimensionless and every dimensionful beta at this ODE step.
    values = {name: evaluate_inputform_expression(expr, env)
              for name, expr in betas.items()}
    for name in ("mHSq", "mS1Sq", "mS2Sq"):
        values[name] = _real(name, values[name])
    mixed_beta = _complex(MIXED, values.pop(MIXED))
    higgs_beta = values.pop("mHSq")
    base = beta_values_to_derivative(values, layout,
                                     divide_by_loop_factor=True)
    dy = np.concatenate((base, np.array([
        mixed_beta.real/LOOP, mixed_beta.imag/LOOP,
        higgs_beta/LOOP], dtype=float)))
    require(dy.shape == y.shape and np.isfinite(dy).all(),
            "nonfinite or mispacked derivative")
    return dy


def _spectrum(manifest: dict, model: str, base: Any,
              mixed: complex) -> dict[str, Any]:
    m1, m2 = float(base.mS1Sq), float(base.mS2Sq)
    eig = analytic_mass_squared(manifest, model, m1, m2, mixed)
    middle = 0.5 * (m1 + m2)
    split = math.hypot(0.5 * (m1 - m2), abs(mixed))
    lo, hi = middle - split, middle + split
    require(np.allclose(eig, [lo]*4 + [hi]*4,
                        rtol=2e-12, atol=2e-9),
            "analytic eigenvalue degeneracy check failed")
    return {"m_minus_squared_gev2": lo, "m_plus_squared_gev2": hi,
            "real_degeneracy_per_eigenvalue": 4,
            "tachyonic": bool(lo < 0),
            "mass_squared_eigenvalues_gev2": eig.tolist()}


def _snapshot(vector: np.ndarray, *, scale: float, layout: Any,
              registry: dict, manifest: dict, model: str) -> dict:
    n = layout.size
    state = unpack_general_scalar_state(vector[:n], layout, registry,
                                        uv=False, mu_gev=scale)
    mixed = complex(float(vector[n]), float(vector[n+1]))
    return {
        "mu_gev": scale,
        "mHSq_gev2": float(vector[n+2]),
        "mS1Sq_gev2": float(state.baseline.mS1Sq),
        "mS2Sq_gev2": float(state.baseline.mS2Sq),
        "muT3GeneralS1S2_gev2": s33.json_state(mixed),
        "scalar_spectrum": _spectrum(manifest, model, state.baseline, mixed),
        # Complete gauge, Yukawa and quartic trajectory, not only mass plots.
        "dimensionless_and_baseline_state": s33.json_state(state.baseline),
        "additional_quartics": {k:s33.json_state(v) for k,v in state.extra.items()},
    }


def run_model(model: str, *, stage33_dir: Path = s34.DATA,
              dimensionless_dir: Path = DIMENSIONLESS_DIR,
              mass_dir: Path = MASS_DIR, outdir: Path = OUTPUT_DIR,
              higgs_mass_squared: float = -7800.0,
              mu2_fraction: float = 0.01, mu2_phase: float = 0.6,
              points: int = 25, rtol: float = 1e-8,
              atol: float = 1e-11) -> dict[str, Any]:
    require(model in MODELS, "only B_-1 and C_-1 are supported")
    require(math.isfinite(higgs_mass_squared), "invalid initial Higgs mass squared")
    require(math.isfinite(mu2_fraction) and 0 <= mu2_fraction <= 0.45,
            "mu2_fraction must lie in [0, 0.45]")
    require(math.isfinite(mu2_phase), "invalid mixed mass phase")
    require(points >= 2, "require at least two saved trajectory points")
    require(0 < rtol < 1 and 0 < atol < 1, "invalid ODE tolerances")

    registry = load_registry(old.REGISTRY)
    manifest = load_manifest(old.MANIFEST)
    _, baseline, sidecar, reused, hashes = s34.verify_boundary(model, stage33_dir)
    info = registry_model(registry, baseline.representation)
    require(quadratic_kind(manifest, model) == "mixed", "unexpected mass structure")
    require(set(sidecar) == {x["name"] for x in info["new_rgbeta_symbols"]}
            and set(reused) == {x["existing_name"]
                                for x in info["preexisting_rgbeta_symbols"]},
            "Stage33 sidecars and general quartic registry disagree")
    initial = GeneralT3ScalarRGState(
        baseline, {k: old.decode_complex(z, k) for k,z in sidecar.items()}
    ).validated(registry)
    y_base, layout = pack_general_scalar_state(initial, registry)
    require(baseline.mu_gev > s34.STOP > 0, "unexpected EFT1 scale interval")
    scale_mass = min(float(baseline.mS1Sq), float(baseline.mS2Sq))
    require(scale_mass > 0, "initial extra scalar masses must be positive")
    mixed0 = mu2_fraction * scale_mass * complex(
        math.cos(mu2_phase), math.sin(mu2_phase))
    y0 = np.concatenate((y_base, [mixed0.real, mixed0.imag,
                                   higgs_mass_squared]))

    tag = s34.safe_model(model)
    destination = outdir / f"{tag}_eft1_coupled_five_mass_v1.json"
    require(not destination.exists(),
            f"refusing to overwrite an existing candidate result: {destination}")
    dim_path = dimensionless_dir / f"{tag}_eft1_quadratic_live_rgbeta.json"
    mass_path = mass_dir / f"{tag}_eft1_full_mass_v1.json"
    require(dim_path.is_file() and mass_path.is_file(),
            f"missing input beta exports: {dim_path} / {mass_path}")
    old_betas = old.checked_payload(
        model, baseline.representation, "mixed", layout, registry, dim_path)
    mass_payload = closure.read_payload(mass_path)
    mass_betas = closure.check_payload(mass_payload, model, "EFT1")
    betas = merge_beta_reports(old_betas, mass_betas, layout.names)

    mu0, mu1 = float(baseline.mu_gev), float(s34.STOP)
    test_dy = beta_derivative(y0, mu_gev=mu0, layout=layout,
                              registry=registry, betas=betas)
    require(test_dy.shape == y0.shape, "initial beta shape mismatch")
    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        return beta_derivative(y, mu_gev=math.exp(t),
                               layout=layout, registry=registry, betas=betas)
    t_points = np.linspace(math.log(mu0), math.log(mu1), points)
    sol = solve_ivp(rhs, (float(t_points[0]), float(t_points[-1])), y0,
                    t_eval=t_points, method="DOP853", rtol=rtol, atol=atol)
    require(sol.success and sol.y.shape == (len(y0), points) and
            np.isfinite(sol.y).all(),
            "coupled integration failed: " + str(sol.message))
    snapshots = [
        _snapshot(np.asarray(sol.y[:,i]),
                  scale=(mu0 if i == 0 else mu1 if i == points-1
                         else float(math.exp(t))),
                  layout=layout, registry=registry, manifest=manifest,
                  model=model)
        for i,t in enumerate(sol.t)
    ]
    final_y = np.asarray(sol.y[:,-1], dtype=float)
    final_state = unpack_general_scalar_state(final_y[:layout.size],
                        layout, registry, uv=False, mu_gev=mu1)
    mixed_final = complex(float(final_y[layout.size]),
                          float(final_y[layout.size+1]))
    # Independent real-Hessian diagnostic at both endpoints.
    for row, state, mu in ((snapshots[0], baseline, mixed0),
                            (snapshots[-1], final_state.baseline, mixed_final)):
        hess = quadratic_mass_hessian(manifest, model,
                 state.mS1Sq, state.mS2Sq, mu)
        require(np.allclose(hess["mass_squared_eigenvalues_gev2"],
                            row["scalar_spectrum"]["mass_squared_eigenvalues_gev2"],
                            rtol=2e-12, atol=2e-9),
                "independent Cartesian Hessian spectrum mismatch")

    data = {
        "format": FORMAT,
        "status": "CandidateFullDimensionlessAndFiveRealMassEFT1Integrated",
        "model": model, "stage": "EFT1", "solver_success": True,
        "mu_initial_gev": mu0, "mu_final_gev": mu1,
        "solver_method": "DOP853", "solver_nfev": int(sol.nfev),
        "rtol": rtol, "atol": atol,
        "one_loop_beta_convention": closure.CONVENTION,
        "loop_factor_applied_exactly_once": True,
        "five_real_mass_coordinates": ["mHSq", "mS1Sq", "mS2Sq",
                                       "Re(muT3GeneralS1S2)",
                                       "Im(muT3GeneralS1S2)"],
        "ode_vector_order": list(layout.names) +
            ["Re(muT3GeneralS1S2)", "Im(muT3GeneralS1S2)", "mHSq"],
        "mass_beta_source": "LiveNativeRGBetaFullMassEFT1",
        "old_mass_betas_all_replaced": True,
        "beta_keys": sorted(betas),
        "old_beta_keys": sorted(old_betas),
        "initial_higgs_mass_squared_gev2": higgs_mass_squared,
        "initial_mu2_fraction": mu2_fraction,
        "initial_mu2_phase_radians": mu2_phase,
        "initial": snapshots[0], "final": snapshots[-1],
        "trajectory": snapshots,
        "input_sha256": {**hashes,
            str(old.REGISTRY): old.fingerprint(old.REGISTRY),
            str(old.MANIFEST): old.fingerprint(old.MANIFEST),
            str(dim_path): old.fingerprint(dim_path),
            str(mass_path): old.fingerprint(mass_path)},
        "physical_matchete_rgbeta_mapping_certified": False,
        "mass_eigenstate_threshold_matching_performed": False,
        "wilson_coefficient_matching_performed": False,
        "physical_neutrino_mass_prediction": False,
    }
    outdir.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n",
                           encoding="utf-8")
    return data


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    selected = ap.add_mutually_exclusive_group(required=True)
    selected.add_argument("--model", choices=MODELS)
    selected.add_argument("--all", action="store_true")
    ap.add_argument("--stage33-dir", type=Path, default=s34.DATA)
    ap.add_argument("--dimensionless-dir", type=Path,
                    default=DIMENSIONLESS_DIR)
    ap.add_argument("--mass-dir", type=Path, default=MASS_DIR)
    ap.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    ap.add_argument("--higgs-mass-squared", type=float, default=-7800.)
    ap.add_argument("--mu2-fraction", type=float, default=0.01)
    ap.add_argument("--mu2-phase", type=float, default=0.6)
    ap.add_argument("--points", type=int, default=25)
    opts = ap.parse_args()
    for model in (MODELS if opts.all else (opts.model,)):
        row = run_model(model, stage33_dir=opts.stage33_dir,
                        dimensionless_dir=opts.dimensionless_dir,
                        mass_dir=opts.mass_dir, outdir=opts.output_dir,
                        higgs_mass_squared=opts.higgs_mass_squared,
                        mu2_fraction=opts.mu2_fraction,
                        mu2_phase=opts.mu2_phase, points=opts.points)
        initial, final = row["initial"], row["final"]
        print(f"{model}: EFT1 {row['mu_initial_gev']:.6g} -> "
              f"{row['mu_final_gev']:.6g} GeV; "
              f"mH^2 {initial['mHSq_gev2']:.6g} -> "
              f"{final['mHSq_gev2']:.6g} GeV^2; "
              f"m_-^2 {initial['scalar_spectrum']['m_minus_squared_gev2']:.6g} "
              f"-> {final['scalar_spectrum']['m_minus_squared_gev2']:.6g} "
              f"GeV^2; nfev={row['solver_nfev']}; "
              f"physical matching UNCERTIFIED", flush=True)


if __name__ == "__main__":
    main()
