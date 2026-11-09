"""Rough UV-Yukawa flavour fit of the D(-2), A(-4) seven-term EFT1 pipeline.

Requires the earlier two-model run with --save-tensors:
  python -X utf8 scripts/RunEFT1SevenTermTwoBenchmarks.py --steps 2 --save-tensors

Fit at the F threshold, with the scalar/gauge/SM-Yukawa running background held
fixed. Start from the source-backed diagonal 3-family UV Yukawas and apply the
same lepton-flavour map to all heavy-fermion Yukawa matrices. This produces a
new LLSS coefficient at MF, which is THEN evolved by the full seven-term ODE.
Both heavy scalar multiplets are integrated out JOINTLY at MS, exactly once.
No old c5_eft1_leading_log is added.

Analytic fit seed: For an initially diagonal C5 (and mnu), choose T with
T mnu_old T^T = U_PMNS^* diag(m_target) U_PMNS^dagger. For barred Matchete
UV Yukawas, Y_new = T^* Y_old on the external lepton-flavour index. This
constitutes a real Yukawa change, NOT a post-hoc neutrino mass rotation.
The full ODE is rerun at the transformed F-threshold boundary to measure the
accuracy of the frozen-background analytic seed. There is no repeated UV RGE,
and a saved UV-scale Yukawa set is only a first approximation.

Run from repository root:
  python -X utf8 scripts/FitEFT1SevenTermTwoBenchmarks.py --steps 2

Scope: rough, conditional numerical fitting (not independently certified).
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys
from typing import Any, Mapping

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

MODELS = ("D_-2", "A_-4")
FORMAT = "t3-two-model-flavour-yukawa-fit-seven-term-eft1-v1"
DEFAULT_BASELINE = Path("output/general_observable/eft1_seven_term_Dm2_Am4_rough_v1")
DEFAULT_OUTPUT = Path("output/general_observable/eft1_seven_term_Dm2_Am4_fit_v1")


def require(condition: Any, reason: str) -> None:
    if not condition:
        raise ValueError(reason)


def locate(root: Path, path: Path) -> Path:
    return path if path.is_absolute() else root / path


def read_json(path: Path) -> dict:
    require(path.is_file(), "Missing JSON input: " + str(path))
    result = json.loads(path.read_text(encoding="utf-8-sig"))
    require(isinstance(result, dict), "Expected JSON object: " + str(path))
    return result


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def from_complex_matrix(obj: Any, name: str) -> np.ndarray:
    """Read either an encoded matrix (real/imag arrays) or native cellwise JSON."""
    if isinstance(obj, Mapping) and set(obj) >= {"real", "imag"}:
        mat = np.asarray(obj["real"], dtype=float) + 1j * np.asarray(obj["imag"], dtype=float)
    else:
        try:
            mat = np.asarray([[complex(float(x["real"]), float(x["imag"]))
                               if isinstance(x, Mapping) else complex(x)
                               for x in row] for row in obj], dtype=complex)
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"Invalid {name} complex matrix") from exc
    require(mat.shape == (3, 3) and np.isfinite(mat).all(), f"Expected finite 3x3 {name}")
    return mat


def encode(z: np.ndarray) -> dict:
    z = np.asarray(z, complex)
    require(np.isfinite(z).all(), "Non-finite complex output")
    return {"real": z.real.tolist(), "imag": z.imag.tolist()}


def cellwise(z: np.ndarray) -> list[list[dict]]:
    z = np.asarray(z, complex)
    require(z.shape == (3, 3), "Expected a 3x3 Yukawa matrix")
    return [[{"real": float(v.real), "imag": float(v.imag)} for v in row] for row in z]


def pmns(s12sq: float, s13sq: float, s23sq: float, delta_deg: float) -> np.ndarray:
    for x, key in ((s12sq, "sin2_theta12"), (s13sq, "sin2_theta13"),
                   (s23sq, "sin2_theta23")):
        require(math.isfinite(x) and 0 < x < 1, f"{key} must be strictly in (0,1)")
    require(math.isfinite(delta_deg), "Dirac phase must be finite")
    s12, s13, s23 = np.sqrt([s12sq, s13sq, s23sq])
    c12, c13, c23 = np.sqrt([1-s12sq, 1-s13sq, 1-s23sq])
    ed = np.exp(1j * math.radians(delta_deg))
    r12 = np.array([[c12, s12, 0], [-s12, c12, 0], [0, 0, 1]], complex)
    r13 = np.array([[c13, 0, s13/ed], [0, 1, 0], [-s13*ed, 0, c13]], complex)
    r23 = np.array([[1, 0, 0], [0, c23, s23], [0, -s23, c23]], complex)
    u = r23 @ r13 @ r12
    require(np.linalg.norm(u.conj().T @ u - np.eye(3)) < 1e-12,
            "Constructed PMNS matrix is not unitary")
    return u


def target_masses(m1: float, dm21: float, dm31: float) -> np.ndarray:
    require(np.isfinite([m1, dm21, dm31]).all(), "Mass targets must be finite")
    require(m1 > 0 and dm21 > 0 and dm31 > dm21,
            "This rough fit requires m1>0 and normal ordering: dm31>dm21>0")
    return np.array([m1, math.sqrt(m1*m1 + dm21), math.sqrt(m1*m1 + dm31)])


def flavour_map(old_c5: np.ndarray, masses: np.ndarray, u: np.ndarray,
                vev: float) -> tuple[np.ndarray, dict]:
    """Explicit Yukawa mixing that yields target masses in linear response."""
    old_c5 = np.asarray(old_c5, complex)
    require(old_c5.shape == (3,3) and np.isfinite(old_c5).all(), "Invalid baseline C5")
    m0 = -0.5 * vev**2 * 1e9 * old_c5  # eV, not GeV
    diagonal = np.diag(m0)
    floor = max(float(np.max(abs(m0))), 1e-40)
    require(np.max(abs(m0 - np.diag(diagonal))) < 1e-8 * floor,
            "Baseline C5 is not diagonal: cannot use the explicit diagonal Yukawa fit")
    require(min(abs(diagonal)) > floor*1e-10, "Baseline has a zero mass eigenvalue")
    # Principal complex square root accounts also for a baseline Majorana sign.
    ratios = np.sqrt(np.asarray(masses, float).astype(complex) / diagonal)
    transfer = u.conj() @ np.diag(ratios)
    target = u.conj() @ np.diag(masses) @ u.conj().T
    algebraic = transfer @ m0 @ transfer.T
    mismatch = float(np.linalg.norm(algebraic-target)/max(np.linalg.norm(target), 1e-60))
    require(mismatch < 1e-10, "Yukawa-congruence fit fails algebraic closure")
    return transfer, {
        "baseline_diagonal_masses_complex_eV": encode(np.diag(diagonal)),
        "baseline_singular_masses_eV": np.sort(np.linalg.svd(m0, compute_uv=False)).tolist(),
        "mass_ratio_square_roots": encode(np.diag(ratios)),
        "target_mnu_matrix_eV": encode(target),
        "algebraic_relative_residual": mismatch,
    }


def observables(mass_matrix_ev: np.ndarray) -> dict:
    """Extract masses + PDG mixing angles in the diagonal charged-lepton basis."""
    mass_matrix_ev = np.asarray(mass_matrix_ev, complex)
    require(mass_matrix_ev.shape == (3,3) and np.isfinite(mass_matrix_ev).all(),
            "Expected finite 3x3 neutrino mass matrix in eV")
    denom = max(float(np.linalg.norm(mass_matrix_ev)), 1e-60)
    require(np.linalg.norm(mass_matrix_ev-mass_matrix_ev.T)/denom < 1e-8,
            "Majorana mass matrix not symmetric")
    # mnu = U* D U^dagger -> mnu^dagger mnu = U D^2 U^dagger.
    h = mass_matrix_ev.conj().T @ mass_matrix_ev
    eigen, u = np.linalg.eigh(h)
    masses = np.sqrt(np.maximum(eigen, 0.0))
    s13 = float(abs(u[0,2])**2)
    c13 = max(1-s13, 1e-15)
    s12 = float(abs(u[0,1])**2 / c13)
    s23 = float(abs(u[1,2])**2 / c13)
    v = {"m1_eV": float(masses[0]), "m2_eV": float(masses[1]), "m3_eV": float(masses[2]),
         "delta_m21_sq_eV2": float(masses[1]**2-masses[0]**2),
         "delta_m31_sq_eV2": float(masses[2]**2-masses[0]**2),
         "sin2_theta12": s12, "sin2_theta13": s13, "sin2_theta23": s23,
         "sum_masses_eV": float(sum(masses))}
    # Dirac CP is determined only where all 3 mixing angles are nonzero.
    eps = float(np.imag(u[0,0]*u[1,1]*np.conj(u[0,1])*np.conj(u[1,0])))
    jden = (np.sqrt(max(s12*(1-s12)*s23*(1-s23),0.0)) * (1-s13) * np.sqrt(max(s13,0.0)))
    if jden > 1e-12:
        sin_delta = float(np.clip(eps/jden,-1,1))
        # |U_mu1|^2 identity in the PDG convention.
        cos_delta_den = 2*np.sqrt(max(s12*(1-s12)*s23*(1-s23)*s13,0.0))
        cos_delta = ((abs(u[1,0])**2 - s12*(1-s23)
                      - (1-s12)*s23*s13)/cos_delta_den)
        v["delta_cp_deg"] = float(np.degrees(np.arctan2(sin_delta, cos_delta)) % 360)
        v["jarlskog_J"] = eps
    else:
        v["delta_cp_deg"] = None
        v["jarlskog_J"] = eps
    return v


def compare_observables(actual: Mapping[str, float], target: Mapping[str, float]) -> dict:
    keys = ("delta_m21_sq_eV2", "delta_m31_sq_eV2", "sin2_theta12",
            "sin2_theta13", "sin2_theta23")
    errors = {k: float((actual[k]-target[k])/target[k]) for k in keys}
    return {"relative_error_per_target": errors,
            "max_absolute_fractional_error": max(abs(x) for x in errors.values())}


def transform_uv_yukawas(state: Mapping[str, Any], transfer: np.ndarray) -> dict:
    """Barred-source convention: C_LLSS ~ (Y1*) (Y2*) / MF.

    Therefore multiplying each UV Yukawa externally by T.conj() multiplies
    its Wilson coefficient by T on each external flavour index. Includes
    D(-2)'s independent y2Conj Yukawa as its own third coupling.
    """
    out = {}
    for name in ("y1", "y2", "y2Conj"):
        if name in state:
            y0 = from_complex_matrix(state[name], name)
            out[name] = cellwise(transfer.conj() @ y0)
    require("y1" in out and "y2" in out, "Missing baseline UV Yukawas")
    return out


def uv_fit_guess(running: Mapping[str, Any], transfer: np.ndarray) -> tuple[dict, dict, dict]:
    src_uv = running["uv_initial"]
    src_mf = running["uv_at_fermion_threshold"]
    initial_y = transform_uv_yukawas(src_uv, transfer)
    threshold_y = transform_uv_yukawas(src_mf, transfer)
    init = {**src_uv, **initial_y}
    return init, threshold_y, initial_y


def run_one(root: Path, model: str, args: argparse.Namespace, u: np.ndarray,
            masses: np.ndarray, target: Mapping[str, float]) -> dict:
    from Numerical.running import All16ConditionalFullWilsonTrial as legacy
    from scripts import RunEFT1SevenTermTwoBenchmarks as seven

    require(model in MODELS, "Model is not one of the requested two benchmarks")
    base_dir = locate(root, args.baseline_root)
    tag = legacy.slug(model)
    baseline_path = base_dir / "per_model" / (tag + ".json")
    tensors_path = base_dir / "tensors" / (tag + ".npz")
    baseline = read_json(baseline_path)
    require(baseline.get("model") == model and baseline.get("format") == seven.FORMAT,
            "Baseline file is not the previous seven-term two-benchmark output")
    require(tensors_path.is_file(),
            "Missing saved full Wilson tensors: rerun prior script with --save-tensors")

    release = legacy.load_release_manifest(locate(root, args.release_dir))
    source = legacy.load_source(release[model], expected_model=model)
    bg_path = locate(root, args.v6_dir) / "numeric_backgrounds" / (tag + ".json")
    background = legacy.load_background(bg_path, source)
    run_path = locate(root, args.running_dir) / tag / "running.json"
    running = legacy.read_json(run_path)
    require(running.get("model") == model, "Running input is for a different model")
    require(source.hash == baseline["source_sha256"], "F-threshold source changed since baseline")
    require(legacy.digest(run_path) == baseline["running_json_sha256"],
            "Renormalisable coupling inputs changed since baseline")
    require(legacy.digest(bg_path) == baseline["scalar_background_sha256"],
            "Scalar background changed since baseline")
    require(background.mass2.shape == (8,8), "Need all eight heavy real scalar components")
    state_f = running["eft1_at_fermion_threshold"]
    state_s = running["eft1_at_scalar_threshold"]
    mf, ms = float(source.mu_f_gev), float(background.mu_s)
    require(0 < ms < mf, "Wrong F-first threshold ordering")
    require(abs(float(running["scales_gev"]["F"])-mf) < 1e-6 and
            abs(float(running["scales_gev"]["S"])-ms) < 1e-6,
            "Running points do not agree with threshold scales")

    original = from_complex_matrix(baseline["matching"]["c5_smeft_conditional_gev_inv"],
                                   "baseline C5")
    with np.load(tensors_path, allow_pickle=False) as data:
        require({"full_C_at_MF_gev_inv", "C5_SMEFT_gev_inv"}.issubset(set(data.files)),
                "Wrong baseline saved NPZ schema")
        saved_c0 = np.asarray(data["full_C_at_MF_gev_inv"], complex)
        saved_c5 = np.asarray(data["C5_SMEFT_gev_inv"], complex)
    require(np.allclose(original,saved_c5,rtol=1e-11,atol=1e-30),
            "Saved baseline NPZ and JSON disagree")
    source_c0 = legacy.embedded_boundary(source.raw, background.U)
    relative_c0 = np.linalg.norm(source_c0-saved_c0)/max(np.linalg.norm(source_c0),1e-60)
    require(relative_c0 < 1e-9, "Source boundary disagrees with saved baseline full C")
    transfer, fit_seed = flavour_map(original,masses,u,args.vev)

    # Transform the actual F-matched LLSS tensor, NOT the final neutrino matrix.
    c_source_fit = np.einsum("pr,qs,rsijab->pqijab", transfer, transfer,
                             source.raw, optimize=True)
    c0_fit = legacy.embedded_boundary(c_source_fit, background.U)
    require(np.linalg.norm(c0_fit)>0 and np.isfinite(c0_fit).all(), "New F-threshold C is invalid")

    if model in legacy.SPECIAL:
        snap = legacy.load_running_snapshot(run_path,model)
        native_to_source = legacy.build_special_basis(model,snap,source)["R_original_from_boundary"]
    else:
        native_to_source = np.eye(8)
    scalar_generators = legacy.rotated_scalar_generators(model,native_to_source @ background.U)
    fermionic = {g: legacy.eft.gauge_generators()[g]["fermion"] for g in ("SU3","SU2","U1")}
    source_fit_matching = np.einsum("pqijab,ak,bl->pqijkl",c_source_fit,
                                    background.U,background.U,optimize=True)
    ward = legacy.source_ward_residual(source_fit_matching,scalar_generators)
    require(max(ward.values()) < 2e-8, "Fitted LLSS source fails gauge Ward invariance")

    yf = legacy.yukawa_at_endpoint(state_f,12)
    ys = legacy.yukawa_at_endpoint(state_s,12)
    c_ms, evolution = seven.seven_term_midpoint(
        c0_fit,yf,ys,state_f,state_s,background.quartic,
        fermionic,scalar_generators,mf,ms,args.steps)
    matching = seven.match_once(c_ms,background,args.vev)
    c5 = from_complex_matrix(matching["c5_smeft_conditional_gev_inv"], "fitted SMEFT C5")
    mnu = -0.5*args.vev**2*1e9*c5
    predicted_m = observables(mnu)
    errors = compare_observables(predicted_m,target)
    require(max(matching["direct_weinberg_projection_residual"],
                matching["bubble_weinberg_projection_residual"]) < 1e-6,
            "Full ODE or scalar bubble fails Weinberg tensor projection")

    # The analytic fit is an initial guess, not an exact identity for the
    # running beta because Ye is flavour diagonal and the master RGE has C*.
    c5_analytic = transfer @ original @ transfer.T
    deviation = float(np.linalg.norm(c5-c5_analytic)/max(np.linalg.norm(c5_analytic),1e-60))
    init_state, mf_yuk, uv_yuk = uv_fit_guess(running,transfer)
    all_yuk = [from_complex_matrix(mf_yuk[k],k) for k in mf_yuk]
    largest_y = float(max(np.max(abs(v)) for v in all_yuk))
    require(largest_y < math.sqrt(4*math.pi), "Trial Yukawa coupling is nonperturbatively large")
    return {
        "format":FORMAT,"model":model,
        "status":"ROUGH_FIT_WITH_EXPLICIT_FULL_SEVEN_TERM_RERUN" if
            errors["max_absolute_fractional_error"] <= args.tolerance else "ROUGH_FIT_OUTSIDE_TOLERANCE",
        "mu_F_gev":mf,"mu_S_gev":ms,
        "threshold_order":"UV -> all F removed -> EFT1(H,S1,S2) -> both S removed jointly -> SMEFT",
        "baseline_seven_term_report":str(baseline_path.relative_to(root)) if baseline_path.is_relative_to(root) else str(baseline_path),
        "source_boundary_sha256":source.hash,
        "target_oscillation_parameters":dict(target),
        "analytic_fit_seed":fit_seed,
        "lepton_flavour_transform_T":encode(transfer),
        "Yukawa_rule":"Y1,Y2 and (D_-2) independent Y2Conj -> conj(T) @ originalY at fixed F threshold, as both source factors are barred in Matchete",
        "uv_yukawas_at_MF":mf_yuk,
        "uv_yukawas_at_initial_scale_guess":uv_yuk,
        "uv_initial_fit_guess":init_state,
        "max_abs_uv_yukawa_at_MF":largest_y,
        "C5_exact_from_rerun_gev_inv":encode(c5),
        "C5_analytic_flavour_seed_gev_inv":encode(c5_analytic),
        "C5_exact_vs_analytic_relative":deviation,
        "mnu_exact_from_rerun_eV":encode(mnu),
        "observables_from_rerun":predicted_m,
        "fit_errors":errors,
        "F_source_Ward_residual":ward,
        "full_C_at_MF_norm":float(np.linalg.norm(c0_fit)),
        "seven_term_ODE_evolution":evolution,
        "joint_scalar_matching":matching,
        "fixed_renormalisable_background":True,
        "UV_RGE_reintegrated_with_fitted_yukawas":False,
        "no_legacy_leading_log_added":True,
        "number_joint_heavy_scalar_eigenpairs":matching["joint_scalar_eigenpairs"],
        "signed_Matchete_GammaCC_to_Weyl_factor_certified":False,
        "complete_quartic_tensor_included":False,
        "independently_fitted_physics":False,
    }


def synthetic_self_test() -> None:
    # Does not require any private repository results or RGE modules.
    mref = np.array([0.018378643793548994,0.020420702125149508,0.024504842492942585])
    v = 246.22
    c0 = -2*np.diag(mref)/(1e9*v*v)
    m = target_masses(0.005,7.5e-5,2.5e-3)
    for delta in (0., 195., 310.):
        u = pmns(0.31,0.022,0.55,delta)
        t,info = flavour_map(c0,m,u,v)
        mfit = -v*v*1e9/2*(t@c0@t.T)
        obs = observables(mfit)
        trg = {"delta_m21_sq_eV2":7.5e-5,"delta_m31_sq_eV2":2.5e-3,
               "sin2_theta12":0.31,"sin2_theta13":0.022,"sin2_theta23":0.55}
        err = compare_observables(obs,trg)
        require(err["max_absolute_fractional_error"]<1e-10,
                f"Synthetic fit algebra invalid for delta={delta}: {err}")
        if delta:
            dd = abs((obs["delta_cp_deg"]-delta+180)%360-180)
            require(dd<1e-6, f"CP extraction disagrees at delta={delta}: {obs['delta_cp_deg']}")
    print("PASS: synthetic Majorana/PMNS Yukawa congruence and observable extraction (3 CP phases)")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-root",type=Path,default=DEFAULT_BASELINE)
    parser.add_argument("--output-root",type=Path,default=DEFAULT_OUTPUT)
    parser.add_argument("--release-dir",type=Path,default=None)
    parser.add_argument("--v6-dir",type=Path,default=None)
    parser.add_argument("--running-dir",type=Path,default=None)
    parser.add_argument("--steps",type=int,default=2)
    parser.add_argument("--vev",type=float,default=246.22)
    parser.add_argument("--m1-ev",type=float,default=0.005)
    parser.add_argument("--dm21",type=float,default=7.5e-5)
    parser.add_argument("--dm31",type=float,default=2.5e-3)
    parser.add_argument("--s12",type=float,default=0.31,help="sin^2 theta12")
    parser.add_argument("--s13",type=float,default=0.022,help="sin^2 theta13")
    parser.add_argument("--s23",type=float,default=0.55,help="sin^2 theta23")
    parser.add_argument("--delta-deg",type=float,default=0.,help="optional CP benchmark angle")
    parser.add_argument("--tolerance",type=float,default=0.005,
                        help="maximum relative error among the five oscillation targets")
    parser.add_argument("--self-test",action="store_true",help="test flavour algebra only")
    args=parser.parse_args(argv)
    if args.self_test:
        synthetic_self_test()
        return 0

    from Numerical.running import All16ConditionalFullWilsonTrial as legacy
    args.release_dir = args.release_dir or legacy.DEFAULT_RELEASE
    args.v6_dir = args.v6_dir or legacy.DEFAULT_V6
    args.running_dir = args.running_dir or legacy.DEFAULT_RUNNING
    require(1 <= args.steps <= 64 and 0 < args.vev < 1000 and
            0 < args.tolerance < 1 and math.isfinite(args.tolerance), "Invalid numerical controls")
    target_mass = target_masses(args.m1_ev,args.dm21,args.dm31)
    target_u = pmns(args.s12,args.s13,args.s23,args.delta_deg)
    target = {"delta_m21_sq_eV2":args.dm21,
              "delta_m31_sq_eV2":args.dm31,
              "sin2_theta12":args.s12,"sin2_theta13":args.s13,
              "sin2_theta23":args.s23,
              "m1_eV":args.m1_ev,"delta_cp_deg":args.delta_deg%360}
    out = locate(ROOT,args.output_root)
    require(not out.exists(),f"Refusing to overwrite existing fit: {out}")
    out.mkdir(parents=True)
    rows=[]
    for model in MODELS:
        print("[rough Yukawa fit] "+model,flush=True)
        try:
            result=run_one(ROOT,model,args,target_u,target_mass,target)
            tag=legacy.slug(model)
            path=out/"per_model"/(tag+".json")
            write_json(path,result)
            # Separate immediately usable Yukawa inputs (not a new UV RGE).
            write_json(out/"fitted_yukawas_at_MF"/(tag+".json"),result["uv_yukawas_at_MF"])
            write_json(out/"uv_initial_guesses"/(tag+".json"),result["uv_initial_fit_guess"])
            obs=result["observables_from_rerun"]
            rows.append({"model":model,"status":result["status"],
                         "m1_eV":obs["m1_eV"],"m2_eV":obs["m2_eV"],"m3_eV":obs["m3_eV"],
                         "dm21_eV2":obs["delta_m21_sq_eV2"],"dm31_eV2":obs["delta_m31_sq_eV2"],
                         "s12":obs["sin2_theta12"],"s13":obs["sin2_theta13"],
                         "s23":obs["sin2_theta23"],"delta_deg":obs["delta_cp_deg"],
                         "max_fractional_error":result["fit_errors"]["max_absolute_fractional_error"],
                         "error":""})
            print("  fitted masses (eV):",[obs[k] for k in ("m1_eV","m2_eV","m3_eV")],flush=True)
            print("  maximum fractional oscillation error:",
                  result["fit_errors"]["max_absolute_fractional_error"],flush=True)
        except Exception as exc:
            print("  FAILED:",type(exc).__name__,str(exc),file=sys.stderr,flush=True)
            rows.append({"model":model,"status":"FAILED", "m1_eV":None,"m2_eV":None,
                         "m3_eV":None,"dm21_eV2":None,"dm31_eV2":None,
                         "s12":None,"s13":None,"s23":None,"delta_deg":None,
                         "max_fractional_error":None,"error":f"{type(exc).__name__}: {exc}"})
    with (out/"comparison.csv").open("w",newline="",encoding="utf-8") as fh:
        writer=csv.DictWriter(fh,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    good=all(x["status"]=="ROUGH_FIT_WITH_EXPLICIT_FULL_SEVEN_TERM_RERUN" for x in rows)
    write_json(out/"summary.json",{
        "format":FORMAT,"status":"COMPLETE" if good else "PARTIAL",
        "models":rows,"targets":target,"masses_target_eV":target_mass.tolist(),
        "fit_uses_UV_yukawa_flavour_rotation":True,
        "source_is_transformed_before_eft1_ODE":True,
        "full_seven_term_ODE_rerun_once_per_model":True,
        "scalar_threshold_is_joint_and_single":True,
        "fixed_renormalisable_background":True,
        "uv_reintegrated":False,"independent_physics_certification":False,
        "note":"Conditional rough F-threshold fit; UV-scale Yukawa output is a guess until UV running is repeated.",
    })
    print("Summary:",out/"summary.json",flush=True)
    return 0 if good else 2


if __name__=="__main__":
    raise SystemExit(main())
