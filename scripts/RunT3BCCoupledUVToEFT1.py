"""Workstream B: opt-in coupled native-RGBeta UV -> EFT1 B(-1)/C(-1) benchmark.

Run from repository root:
  python -X utf8 -m scripts.RunT3BCCoupledUVToEFT1 --all \
      --higgs-mass-squared -7800 --mu2-fraction 0.01 --mu2-phase 0.6

The model's existing Stage29/31 candidate UV starting state is evolved from
1e7 to 1e5 GeV including gauge/Yukawa/quartic parameters, the Majorana MF,
and five real scalar masses. A tree-level renormalisable F removal then hands
the running scalar/quartic state to EFT1, evolved down to 1e3 GeV.

All four scalar-mass beta expressions in each stage come from the previously
extracted native full-mass RGBeta reports, NOT the older incomplete exports.
No Matchete/RGBeta physical map, finite threshold correction, Wilson coefficient
matching, or physical neutrino mass is claimed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np
from scipy.integrate import solve_ivp

from scripts import CoupledUVRunning as s31
from scripts import CoupledA0UVRunning as s30
from scripts import FermionThresholdHandoff as s33
from scripts import RunGeneralT3QuadraticEFT1 as legacy
from scripts import RunT3BCFullScalarMassEFT1 as eft_runner
from scripts import VerifyT3BCFullScalarMassClosure as closure
from Numerical.core.BetaVector import beta_values_to_derivative
from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression, state_environment
from Numerical.core.StateVector import unpack_uv_state
from Numerical.experimental.GeneralT3ScalarRGState import (
    GeneralT3ScalarRGState, load_registry, pack_general_scalar_state,
    unpack_general_scalar_state, project_general_after_fermion,
)
from Numerical.experimental.GeneralT3ScalarQuadraticMasses import (
    load_manifest, quadratic_kind, quadratic_mass_hessian,
)
from Numerical.running.FermionThresholdBoundary import fermion_threshold_masses

ROOT = legacy.ROOT
MODELS = eft_runner.MODELS
STAGE29_DIR = ROOT / "output/general_t3_stage22"
DIMENSIONLESS_DIR = legacy.OUTPUT
MASS_DIR = eft_runner.MASS_DIR
OUTPUT_DIR = ROOT / "output/workstream_b_bc_coupled_uv_eft1"
FORMAT = "workstream-b-bc-coupled-uv-through-eft1-five-mass-v1"
MIXED = eft_runner.MIXED
LOOP = eft_runner.LOOP


def require(condition: bool, reason: str) -> None:
    if not condition:
        raise ValueError("B/C coupled UV->EFT1: " + reason)


def _sha(path: Path) -> str:
    return legacy.fingerprint(path)


def load_uv_start(model: str, *, stage29_dir: Path, registry: dict):
    """Recreate the Stage31 UV input; do not reverse-evolve a Stage33 endpoint."""
    require(model in MODELS, "unsupported UV model")
    entry, new, reused = s31.inventories(registry, model)
    require(len(new) == 1 and len(reused) == 6,
            "B/C extra-quartic inventory is not 1 new and 6 reused")
    stage29_path = stage29_dir / "stage29_convention_uv_handoff.json"
    stage29 = s30.read_json(stage29_path)
    s30.verify_stage29_registry_provenance(stage29, legacy.REGISTRY)
    values = s31.stage29_values(stage29, entry, new, reused)
    config_path = s31.config_path(entry)
    config = s30.read_json(config_path)
    old_vector, old_layout = s31.prepare_state(config, entry, reused, values)
    uv = unpack_uv_state(old_vector, old_layout, mu_gev=s31.START)
    state = GeneralT3ScalarRGState(uv, {name: values[name] for name in new}).validated(registry)
    y, layout = pack_general_scalar_state(state, registry)
    require(tuple(layout.names[:len(old_layout.names)]) == tuple(old_layout.names),
            "packed Stage31 UV baseline layout changed")
    require(np.array_equal(y[:old_layout.size], old_vector),
            "packed UV input differs from the Stage31 fitted benchmark")
    require(tuple(layout.names[ len(old_layout.names):]) == tuple(sorted(new)),
            "UV new-quartic packing order changed")
    require(uv.representation.majorana_fermion and uv.mu_gev == s31.START,
            "expected Majorana B/C UV state at Stage31 start")
    return state, y, layout, {
        stage29_path: _sha(stage29_path),
        config_path: _sha(config_path),
    }


def checked_uv_dimensionless(model: str, *, registry: dict, layout: Any,
                             path: Path) -> dict[str, str]:
    """Fail closed on complete native UV beta inventory and metadata."""
    payload = s30.read_json(path)
    rep = layout.representation
    entry = registry["models"][model]
    wanted_new = {r["name"] for r in entry["new_rgbeta_symbols"]}
    expected = {
        "SharedScalar": False, "dS1": 2, "dS2": 2,
        "dF": 1 if model == "B_-1" else 3, "alpha": -1,
        "GeneralScalarModel": model, "GeneralScalarStage": "UV",
        "GeneralScalarBasisStatus": "CandidateNativeRGBetaUnmappedToMatchete",
        "GeneralScalarQuarticBetaBasis": "RGBetaFullQuarticProjection",
        "GeneralScalarQuadraticMassMixing": "IncludedCandidate",
        "GeneralScalarAdditionalQuartics": 7,
        "GeneralScalarExistingStructuralCount": 6,
        "QuadraticKind": "mixed", "QuadraticMassSymbol": MIXED,
        "QuadraticPrefactorConvention": "V2 = mu2 * invariant + h.c.; no 1/2",
        "PhysicalMatcheteRGBetaMassBasisCertified": False,
        "ReportBetaConvention": closure.CONVENTION,
        "MajoranaFermion": True,
    }
    meta = payload.get("metadata")
    require(payload.get("format") == "general-t3-quadratic-rgbeta-experimental-v1"
            and payload.get("status") == "SuccessCandidateMassBeta"
            and payload.get("physical_mapping_certified") is False
            and isinstance(meta, dict), "unknown UV RGBeta report format/status")
    for key, val in expected.items():
        require(meta.get(key) == val,
                f"UV RGBeta metadata {key}: {meta.get(key)!r} != {val!r}")
    require(set(meta.get("GeneralScalarNewBetaKeys", [])) == wanted_new,
            "UV new-quartic metadata differs from registry")
    expected_keys = set(layout.names) | {MIXED}
    report = payload.get("report_betas")
    require(isinstance(report, dict) and set(report) == expected_keys
            and set(payload.get("beta_symbols", [])) == expected_keys,
            "UV beta dictionary incomplete or inconsistent with packed state")
    require(all(isinstance(v, str) and v.strip() and not any(
        err in v for err in ("BetaTerm[", "$Failed", "$Aborted", "Missing[", "Indeterminate"))
        for v in report.values()), "unresolved UV beta expression")
    return dict(report)


def _symmetrize_majorana(y: np.ndarray, layout: Any) -> np.ndarray:
    """Remove numerical antisymmetry in the Majorana MF coordinate only."""
    arr = np.array(y, dtype=float, copy=True)
    block = layout.block("MF")
    require(block.kind == "complex_matrix" and block.shape == (3, 3),
            "unexpected MF state-vector layout")
    n = 9
    for start in (block.start, block.start + n):
        m = arr[start:start+n].reshape(3, 3)
        mismatch = np.max(np.abs(m-m.T))
        require(mismatch <= 1e-7 * max(1., float(np.max(np.abs(m)))),
                "MF evolved a substantial antisymmetric component")
        m[:] = 0.5*(m+m.T)
    return arr


def uv_derivative(y: np.ndarray, *, mu_gev: float, layout: Any,
                  registry: dict, betas: Mapping[str, str]) -> np.ndarray:
    """UV: same complete environment for all dimensionless + five mass betas."""
    v = np.asarray(y, dtype=float)
    n = layout.size
    require(v.shape == (n+3,) and np.isfinite(v).all(),
            "incorrect or nonfinite UV ODE vector")
    require(math.isfinite(mu_gev) and mu_gev > 0, "invalid UV RG scale")
    v = _symmetrize_majorana(v, layout)
    state = unpack_general_scalar_state(v[:n], layout, registry,
                                        uv=True, mu_gev=mu_gev)
    env = state_environment(state.baseline)
    env.update(state.extra)
    env[MIXED] = complex(v[n], v[n+1])
    env["mHSq"] = float(v[n+2])
    values = {name: evaluate_inputform_expression(expr, env)
              for name, expr in betas.items()}
    require(set(values) == set(layout.names) | {MIXED, "mHSq"},
            "incomplete evaluated UV beta set")
    for key in ("mS1Sq", "mS2Sq", "mHSq"):
        values[key] = eft_runner._real(key, values[key])
    mu_beta = eft_runner._complex(MIXED, values.pop(MIXED))
    higgs_beta = values.pop("mHSq")
    base = beta_values_to_derivative(values, layout, divide_by_loop_factor=True)
    result = np.concatenate([base, [mu_beta.real/LOOP, mu_beta.imag/LOOP,
                                    higgs_beta/LOOP]])
    require(result.shape == v.shape and np.isfinite(result).all(),
            "nonfinite UV derivative or incorrect packing")
    return result


def _snapshot(y: np.ndarray, *, mu_gev: float, layout: Any, registry: dict,
              manifest: dict, model: str, uv: bool) -> dict[str, Any]:
    v = _symmetrize_majorana(y, layout) if uv else np.asarray(y, dtype=float)
    n = layout.size
    state = unpack_general_scalar_state(v[:n], layout, registry,
                                        uv=uv, mu_gev=mu_gev)
    mu_mix = complex(v[n], v[n+1])
    masses = eft_runner._spectrum(manifest, model, state.baseline, mu_mix)
    out = {
        "mu_gev": float(mu_gev),
        "mHSq_gev2": float(v[n+2]),
        "mS1Sq_gev2": float(state.baseline.mS1Sq),
        "mS2Sq_gev2": float(state.baseline.mS2Sq),
        "muT3GeneralS1S2_gev2": s33.json_state(mu_mix),
        "scalar_spectrum": masses,
        "state": s33.json_state(state.baseline),
        "additional_quartics": {k: s33.json_state(z) for k, z in state.extra.items()},
    }
    if uv:
        out["heavy_fermion_singular_masses_gev"] = fermion_threshold_masses(
            state.baseline).masses_gev.tolist()
    return out


def project_uv_endpoint(uv_y: np.ndarray, *, uv_layout: Any,
                       registry: dict, mu_gev: float):
    """Project F/y1/y2/MF away and transport all remaining mass coordinates."""
    n = uv_layout.size
    v = _symmetrize_majorana(np.asarray(uv_y), uv_layout)
    full_uv = unpack_general_scalar_state(v[:n], uv_layout, registry,
                                          uv=True, mu_gev=mu_gev)
    after_f = project_general_after_fermion(full_uv, registry)
    assert after_f.baseline.mu_gev == mu_gev
    eft_base, eft_layout = pack_general_scalar_state(after_f, registry)
    eft_y = np.concatenate((eft_base, v[n:n+3]))
    require(not any(x in eft_layout.names for x in ("MF", "y1", "y2")),
            "fermionic couplings leaked into EFT1")
    require("mS1Sq" in eft_layout.names and "mS2Sq" in eft_layout.names,
            "projected EFT1 diagonal masses absent")
    for name in ("mS1Sq", "mS2Sq"):
        np.testing.assert_allclose(eft_y[eft_layout.block(name).start],
            v[uv_layout.block(name).start], rtol=0, atol=0,
            err_msg=f"{name} discontinuity at candidate F threshold")
    require(np.array_equal(eft_y[-3:], v[-3:]),
            "Higgs or mixed mass discontinuity at candidate F threshold")
    return full_uv, after_f, eft_y, eft_layout


def _integrate(y0: np.ndarray, rhs, start: float, stop: float, *,
               points: int, rtol: float, atol: float, max_step: float | None = None):
    require(start > stop > 0 and points >= 2 and points <= 200,
            "invalid running interval or sample count")
    kw = {} if max_step is None else {"max_step": max_step}
    times = np.linspace(math.log(start), math.log(stop), points)
    sol = solve_ivp(rhs, (float(times[0]), float(times[-1])),
                    np.asarray(y0, dtype=float), t_eval=times,
                    method="DOP853", rtol=rtol, atol=atol, **kw)
    require(sol.success and sol.y.shape == (len(y0), points) and
            np.isfinite(sol.y).all(),
            "coupled integration failed: " + str(sol.message))
    return sol


def run_model(model: str, *, stage29_dir: Path = STAGE29_DIR,
              dimensionless_dir: Path = DIMENSIONLESS_DIR,
              mass_dir: Path = MASS_DIR, output_dir: Path = OUTPUT_DIR,
              higgs_mass_squared: float = -7800., mu2_fraction: float = 0.01,
              mu2_phase: float = 0.6, points: int = 25,
              rtol: float = 1e-8, atol: float = 1e-11) -> dict[str, Any]:
    require(model in MODELS, "only B_-1 and C_-1 are supported")
    require(math.isfinite(higgs_mass_squared), "nonfinite initial Higgs mass")
    require(math.isfinite(mu2_fraction) and 0 <= mu2_fraction <= 0.45,
            "mu2_fraction must lie in [0, 0.45]")
    require(math.isfinite(mu2_phase), "nonfinite complex mass phase")
    require(0 < rtol < 1 and 0 < atol < 1, "invalid ODE tolerances")
    require(quadratic_kind(load_manifest(legacy.MANIFEST), model) == "mixed",
            "unexpected B/C mixed mass invariant")
    registry = load_registry(legacy.REGISTRY)
    manifest = load_manifest(legacy.MANIFEST)
    uv_state, base_y, uv_layout, seed_hashes = load_uv_start(
        model, stage29_dir=stage29_dir, registry=registry)
    mass_scale = min(uv_state.baseline.mS1Sq, uv_state.baseline.mS2Sq)
    require(mass_scale > 0, "UV BSM diagonal masses must be positive")
    mixed0 = (mu2_fraction*mass_scale)*complex(math.cos(mu2_phase),
                                               math.sin(mu2_phase))
    uv_y0 = np.concatenate((base_y, [mixed0.real, mixed0.imag,
                                     higgs_mass_squared]))
    tag = model.replace("-", "m").replace("+", "p")
    uv_dim_path = dimensionless_dir / f"{tag}_uv_quadratic_live_rgbeta.json"
    eft_dim_path = dimensionless_dir / f"{tag}_eft1_quadratic_live_rgbeta.json"
    uv_mass_path = mass_dir / f"{tag}_uv_full_mass_v1.json"
    eft_mass_path = mass_dir / f"{tag}_eft1_full_mass_v1.json"
    for p in (uv_dim_path, eft_dim_path, uv_mass_path, eft_mass_path):
        require(p.is_file(), "missing native RGBeta export: " + str(p))
    uv_old = checked_uv_dimensionless(model, registry=registry,
                                     layout=uv_layout, path=uv_dim_path)
    uv_mass = closure.check_payload(closure.read_payload(uv_mass_path), model, "UV")
    uv_betas = eft_runner.merge_beta_reports(uv_old, uv_mass, uv_layout.names)
    require({"y1", "y2", "MF"}.issubset(uv_betas),
            "UV heavy-fermion running absent")
    uv_start, uv_stop, eft_stop = float(s31.START), float(s31.STOP), float(legacy.s34.STOP)
    require((uv_start, uv_stop, eft_stop) == (1e7, 1e5, 1e3),
            "unexpected Stage31/Stage34 scale conventions")
    uv_rhs = lambda t, y: uv_derivative(y, mu_gev=math.exp(t),
                              layout=uv_layout, registry=registry, betas=uv_betas)
    uv_rhs(math.log(uv_start), uv_y0)
    uv_sol = _integrate(uv_y0, uv_rhs, uv_start, uv_stop,
                        points=points, rtol=rtol, atol=atol, max_step=0.35)
    uv_final_y = _symmetrize_majorana(uv_sol.y[:, -1], uv_layout)
    full_uv, eft_state, eft_y0, eft_layout = project_uv_endpoint(
        uv_final_y, uv_layout=uv_layout, registry=registry, mu_gev=uv_stop)
    eft_old = legacy.checked_payload(model, eft_state.baseline.representation,
                           "mixed", eft_layout, registry, eft_dim_path)
    eft_mass = closure.check_payload(closure.read_payload(eft_mass_path), model, "EFT1")
    eft_betas = eft_runner.merge_beta_reports(eft_old, eft_mass, eft_layout.names)
    eft_rhs = lambda t, y: eft_runner.beta_derivative(
        y, mu_gev=math.exp(t), layout=eft_layout,
        registry=registry, betas=eft_betas)
    eft_rhs(math.log(uv_stop), eft_y0)
    eft_sol = _integrate(eft_y0, eft_rhs, uv_stop, eft_stop,
                         points=points, rtol=rtol, atol=atol)
    uv_rows = [_snapshot(uv_sol.y[:, i],
                       mu_gev=uv_start if i == 0 else uv_stop if i == points-1
                       else math.exp(uv_sol.t[i]),
                       layout=uv_layout, registry=registry,
                       manifest=manifest, model=model, uv=True)
               for i in range(points)]
    eft_rows = [_snapshot(eft_sol.y[:, i],
                        mu_gev=uv_stop if i == 0 else eft_stop if i == points-1
                        else math.exp(eft_sol.t[i]),
                        layout=eft_layout, registry=registry,
                        manifest=manifest, model=model, uv=False)
                for i in range(points)]
    # The projection is continuous in all five scalar-mass coordinates.
    uv_end = uv_rows[-1]
    eft_begin = eft_rows[0]
    for field in ("mHSq_gev2", "mS1Sq_gev2", "mS2Sq_gev2"):
        require(uv_end[field] == eft_begin[field],
                f"unexpected tree-level threshold discontinuity in {field}")
    require(uv_end["muT3GeneralS1S2_gev2"] ==
            eft_begin["muT3GeneralS1S2_gev2"],
            "mixed bilinear discontinuity at fermion projection")
    for row, is_uv in ((uv_rows[0], True), (uv_rows[-1], True),
                       (eft_rows[0], False), (eft_rows[-1], False)):
        base = row["state"]
        mu = row["muT3GeneralS1S2_gev2"]
        hessian = quadratic_mass_hessian(manifest, model, base["mS1Sq"],
                           base["mS2Sq"], complex(mu["real"], mu["imag"]))
        require(np.allclose(hessian["mass_squared_eigenvalues_gev2"],
                row["scalar_spectrum"]["mass_squared_eigenvalues_gev2"],
                rtol=2e-12, atol=2e-9),
                "Cartesian Hessian does not agree with analytic scalar spectrum")
    mass_fermions = uv_rows[-1]["heavy_fermion_singular_masses_gev"]
    # Treat a grouped 100 TeV removal as an illustrative choice; it is not an
    # exact on-shell or sequential mass-eigenstate threshold matching scheme.
    provenance = {str(p): h for p, h in seed_hashes.items()}
    for p in (legacy.REGISTRY, legacy.MANIFEST, uv_dim_path,
              eft_dim_path, uv_mass_path, eft_mass_path):
        provenance[str(p)] = _sha(p)
    data = {
        "format": FORMAT,
        "status": "CandidateCoupledUVAndEFT1FiveRealMassesIntegrated",
        "model": model, "solver_success": True,
        "one_loop_convention": closure.CONVENTION,
        "loop_factor_applied_exactly_once": True,
        "uv_range_gev": [uv_start, uv_stop],
        "eft1_range_gev": [uv_stop, eft_stop],
        "initial_higgs_mass_squared_gev2_at_uv": higgs_mass_squared,
        "initial_mu2_fraction_at_uv": mu2_fraction,
        "initial_mu2_phase_at_uv": mu2_phase,
        "uv_solver_nfev": int(uv_sol.nfev),
        "eft1_solver_nfev": int(eft_sol.nfev),
        "solver_method": "DOP853", "rtol": rtol, "atol": atol,
        "uv_beta_keys": sorted(uv_betas),
        "eft1_beta_keys": sorted(eft_betas),
        "mass_beta_source": "LiveNativeRGBetaFullMassSeparateUVAndEFT1",
        "previous_mass_betas_replaced_in_both_stages": True,
        "five_real_mass_coordinates": ["mHSq", "mS1Sq", "mS2Sq",
            "Re(muT3GeneralS1S2)", "Im(muT3GeneralS1S2)"],
        "uv_vector_order": list(uv_layout.names) +
            ["Re(muT3GeneralS1S2)", "Im(muT3GeneralS1S2)", "mHSq"],
        "eft1_vector_order": list(eft_layout.names) +
            ["Re(muT3GeneralS1S2)", "Im(muT3GeneralS1S2)", "mHSq"],
        "threshold": {
            "scale_gev": uv_stop,
            "projection": "F/y1/y2/MF removed; other renormalisable couplings, all new quartics, and all five scalar masses continuous",
            "tree_level_renormalisable_only": True,
            "finite_mass_matching_performed": False,
            "wilson_coefficient_matching_performed": False,
            "heavy_fermion_singular_masses_gev": mass_fermions,
        },
        "uv_initial": uv_rows[0],
        "uv_threshold": uv_rows[-1],
        "eft1_threshold": eft_rows[0],
        "eft1_final": eft_rows[-1],
        "uv_trajectory": uv_rows,
        "eft1_trajectory": eft_rows,
        "input_sha256": provenance,
        "physical_matchete_rgbeta_mapping_certified": False,
        "scalar_mass_eigenstate_threshold_matching_performed": False,
        "physical_neutrino_mass_prediction": False,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"{tag}_uv_to_eft1_coupled_five_mass_v1.json"
    require(not path.exists(), "refusing to overwrite an existing report: " + str(path))
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    return data


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument("--model", choices=MODELS)
    choice.add_argument("--all", action="store_true")
    parser.add_argument("--stage29-dir", type=Path, default=STAGE29_DIR)
    parser.add_argument("--dimensionless-dir", type=Path, default=DIMENSIONLESS_DIR)
    parser.add_argument("--mass-dir", type=Path, default=MASS_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--higgs-mass-squared", type=float, default=-7800.)
    parser.add_argument("--mu2-fraction", type=float, default=0.01)
    parser.add_argument("--mu2-phase", type=float, default=0.6)
    parser.add_argument("--points", type=int, default=25)
    args = parser.parse_args()
    for model in MODELS if args.all else (args.model,):
        report = run_model(model, stage29_dir=args.stage29_dir,
            dimensionless_dir=args.dimensionless_dir,
            mass_dir=args.mass_dir, output_dir=args.output_dir,
            higgs_mass_squared=args.higgs_mass_squared,
            mu2_fraction=args.mu2_fraction, mu2_phase=args.mu2_phase,
            points=args.points)
        initial, uv_end, end = (report[k] for k in
                             ("uv_initial", "uv_threshold", "eft1_final"))
        print(f"{model}: UV 1e7 -> 1e5 -> EFT1 1e3 GeV; "
              f"mH^2 {initial['mHSq_gev2']:.6g} -> "
              f"{uv_end['mHSq_gev2']:.6g} -> {end['mHSq_gev2']:.6g} GeV^2; "
              f"m_-^2 {initial['scalar_spectrum']['m_minus_squared_gev2']:.6g} "
              f"-> {uv_end['scalar_spectrum']['m_minus_squared_gev2']:.6g} "
              f"-> {end['scalar_spectrum']['m_minus_squared_gev2']:.6g} GeV^2; "
              f"nfev UV/EFT1={report['uv_solver_nfev']}/{report['eft1_solver_nfev']}; "
              "physical matching UNCERTIFIED", flush=True)


if __name__ == "__main__":
    main()
