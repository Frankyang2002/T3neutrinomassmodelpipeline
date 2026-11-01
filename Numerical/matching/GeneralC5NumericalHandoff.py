"""Provisional A(-2) direct-matching C5 -> low-energy neutrino benchmark.

Uses original matched symbolic C5, evaluates per-heavy-flavour loop kernels via
Wolfram, then reuses the repository's SM+Weinberg runner and neutrino observables.
NOT certified: Matchete/RGBeta quartic-basis mapping, real-triplet CG normalization,
or EFT1 higher-dimensional Wilson evolution. Intended for rough numerical tests.
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from Numerical.running.SMWeinbergEvolution import (
    SMWeinbergInitialConditions, evolve_sm_weinberg,
)
from physics.NeutrinoObservables import calculate_neutrino_observables
from physics.NeutrinoMass import numerical_neutrino_mass_matrix
from scipy.integrate import solve_ivp
from Numerical.running import SMWeinbergEvolution as sm_rge


def _c(x):
    if isinstance(x, dict):
        return complex(float(x["real"]), float(x["imag"]))
    return complex(x)


def _mat(x):
    a = np.asarray([[_c(y) for y in row] for row in x], dtype=complex)
    if a.shape != (3, 3):
        raise ValueError(f"Expected 3x3 matrix, got {a.shape}")
    return a


def _jsonc(z):
    return {"real": float(complex(z).real), "imag": float(complex(z).imag)}


def _dump_matrix(a):
    return [[_jsonc(x) for x in row] for row in a]


def run(args):
    data = json.loads(args.running_json.read_text(encoding="utf-8-sig"))
    if data.get("model") != "A_-2":
        raise ValueError("This matching coefficient was calculated for A_-2 only")
    uv = data["uv_at_fermion_threshold"]
    eft = data["eft1_at_scalar_threshold"]
    mf = np.asarray(_mat(uv["MF"]))
    if not np.allclose(mf, np.diag(np.diag(mf)), rtol=1e-10, atol=1e-7):
        raise ValueError("Non-diagonal heavy-fermion MF requires a flavour rotation")
    masses_f = np.diag(mf)
    if np.max(abs(masses_f.imag)) > 1e-7 or np.any(masses_f.real <= 0):
        raise ValueError("Require positive, real, diagonal heavy fermion masses")
    y1 = _mat(uv["y1"])
    y2 = _mat(uv["y2"])
    y2c = _mat(uv["y2Conj"])
    m1sq = float(eft["mS1Sq"])
    m2sq = float(eft["mS2Sq"])
    kappa = _c(args.kappa_real) + 1j * args.kappa_imag
    # Convention of the transformed Wolfram matching model; CG normalization unverified.
    mplus2 = m2sq + 2 * abs(kappa)
    mminus2 = m2sq - 2 * abs(kappa)
    if min(m1sq, mplus2, mminus2) <= 0:
        raise ValueError("Nonpositive scalar mass squared")
    if abs(kappa) < 1e-15:
        raise ValueError("Use nonzero kappa for the split basis; degeneracy requires a limit")
    phase = np.exp(-0.5j * np.angle(kappa))
    scalar_quartic = _c(args.quartic_real) + 1j * args.quartic_imag
    lam = _c(eft.get("lambdaT3", uv["lambdaT3"]))
    if args.lambda_t3_real is not None:
        lam = complex(args.lambda_t3_real, args.lambda_t3_imag)
    point_data = []
    for index, fermion_mass in enumerate(masses_f.real):
        couplings = {
            "MF": float(fermion_mass),
            "MS1": math.sqrt(m1sq),
            "T32PlusMass": math.sqrt(mplus2),
            "T32MinusMass": math.sqrt(mminus2),
            "T3NeutralPhase": _jsonc(phase),
            "lambdaT3": _jsonc(lam),
            "lambdaT3ExtraQ020101I1": _jsonc(scalar_quartic),
        }
        for channel in ("y2", "y2Conj"):
            point_data.append({"index":index,"channel":channel,"couplings":couplings})
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as folder:
        folder = Path(folder)
        inp = folder / "input.json"
        out = folder / "kernels.json"
        # Matchete hbar is a loop counter: the numerical convention is user-set.
        inp.write_text(json.dumps({"points":point_data,"hbar_value":args.hbar_value}),encoding="utf-8")
        command = [args.wolframscript,"-file",str(args.evaluator_wl),
                   str(args.coefficient),str(inp),str(out)]
        completed = subprocess.run(command, text=True, capture_output=True)
        if completed.returncode != 0:
            raise RuntimeError("Wolfram kernel evaluation failed:\n" +
                               completed.stdout[-3500:] + "\n" + completed.stderr[-3500:])
        results = json.loads(out.read_text(encoding="utf-8"))["results"]
    kernels = {(x["channel"],x["index"]):complex(x["real"],x["imag"])
               for x in results if x["ok"]}
    if len(kernels) != 6:
        raise RuntimeError("Expected exactly six evaluated flavour kernels")
    # Extracted channels correspond to conjugated Yukawas. Symmetrize alpha,beta.
    c5 = np.zeros((3,3),dtype=complex)
    for i in range(3):
        c5 += kernels["y2",i] * np.outer(np.conj(y1[:,i]),np.conj(y2[:,i]))
        c5 += kernels["y2Conj",i] * np.outer(np.conj(y1[:,i]),np.conj(y2c[:,i]))
    raw_asymmetry = float(np.linalg.norm(c5-c5.T)/max(np.linalg.norm(c5),1e-60))
    c5 = 0.5*(c5+c5.T)
    sm = SMWeinbergInitialConditions(
        gY=float(eft["gY"]),g2=float(eft["g2"]),g3=float(eft["g3"]),
        lambdaH=float(eft["lambdaH"]),
        ye=_mat(eft["ye"]),yu=_mat(eft["yu"]),yd=_mat(eft["yd"]),K=c5)
    # A singular full-SM trajectory must not discard the already evaluated C5.
    # Retry with a stiff integrator, then provide explicitly labelled
    # threshold-only observables if the low-energy trajectory is still invalid.
    mu_s = float(data["scales_gev"]["S"])
    evolved_k = None
    solver_status = "Success"
    solver_errors = []
    try:
        result = evolve_sm_weinberg(sm, mu_s, args.mu_low)
        evolved_k = result.K
    except (RuntimeError, ValueError, FloatingPointError) as exc:
        solver_errors.append("DOP853: " + str(exc))
        try:
            state0 = sm.validated()
            solution = solve_ivp(sm_rge._beta,
                (math.log(mu_s), math.log(args.mu_low)),
                sm_rge._pack(state0), method="Radau",
                rtol=1e-7, atol=1e-12)
            if not solution.success or not np.all(np.isfinite(solution.y[:, -1])):
                raise RuntimeError(solution.message)
            evolved_k = sm_rge._unpack(solution.y[:, -1])[-1]
            evolved_k = 0.5*(evolved_k + evolved_k.T)
            solver_status = "StiffRadauSuccess"
        except (RuntimeError, ValueError, FloatingPointError) as stiff_exc:
            solver_errors.append("Radau: " + str(stiff_exc))
            evolved_k = c5.copy()
            solver_status = "LowEnergyRGEFailed_ThresholdOnly"
    mnu=numerical_neutrino_mass_matrix(evolved_k)
    obs=calculate_neutrino_observables(mnu,ordering=args.ordering)
    payload={
        "status":"ProvisionalDirectMatching_NotFullGeneralEFT1Prediction",
        "model":"A_-2", "matching_input":str(args.coefficient),
        "running_input":str(args.running_json),
        "scalar_masses_gev":{"S1":math.sqrt(m1sq),"S2Plus":math.sqrt(mplus2),
                             "S2Minus":math.sqrt(mminus2)},
        "physical_mass_mapping":{"kappa":_jsonc(kappa),"phase":_jsonc(phase),
                                  "bilinear_coefficient_factor":2},
        "explicit_matchete_quartic":_jsonc(scalar_quartic),
        "matchete_lambdaT3":_jsonc(lam),"hbar_numeric_value":args.hbar_value,
        "kernels_gev_inverse":[{"index":i,"y1bar_y2bar":_jsonc(kernels["y2",i]),
                                  "y1bar_y2conjbar":_jsonc(kernels["y2Conj",i])} for i in range(3)],
        "C5_at_scalar_threshold":_dump_matrix(c5),
        "C5_at_low_scale":_dump_matrix(evolved_k) if solver_status != "LowEnergyRGEFailed_ThresholdOnly" else None,
        "C5_at_scalar_threshold_only":_dump_matrix(c5),
        "observables_scale_gev":args.mu_low if solver_status != "LowEnergyRGEFailed_ThresholdOnly" else mu_s,
        "low_energy_rge_status":solver_status,
        "low_energy_rge_errors":solver_errors,
        "sm_boundary_diagnostics":{
            "gY":float(eft["gY"]),"g2":float(eft["g2"]),"g3":float(eft["g3"]),
            "lambdaH":float(eft["lambdaH"]),
            "yu_max_abs":float(np.max(abs(sm.yu))),
            "yd_max_abs":float(np.max(abs(sm.yd))),
            "ye_max_abs":float(np.max(abs(sm.ye)))},
        "neutrino_masses_ev":obs.masses_ev.tolist(),
        "delta_m21_sq_ev2":obs.delta_m21_sq_ev2,
        "delta_m31_sq_ev2":obs.delta_m31_sq_ev2,
        "pmns_abs":obs.pmns_abs.tolist(),
        "takagi_residual":obs.takagi_residual,
        "flavour_unsymmetrized_relative_asymmetry":raw_asymmetry,
        "not_yet_certified":["real SU(2) triplet field and CG normalization",
           "Matchete-RGBeta coupling mapping (quartic explicitly supplied)",
           "bilinear SU(2) normalization", "EFT1 Wilson-operator running",
           "threshold matching at distinct heavy masses",
           "hbar normalization and lepton-flavour symmetrization convention"],
    }
    args.output_json.write_text(json.dumps(payload,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print("A_-2 provisional C5 benchmark:",args.output_json)
    print("low-energy RGE status:", solver_status)
    print("observables scale [GeV]:",payload["observables_scale_gev"])
    print("neutrino masses [eV]:",obs.masses_ev.tolist())
    return 0


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--running-json",type=Path,default=Path("output/general_A_m2_running.json"))
    p.add_argument("--coefficient",type=Path,default=Path("output/general_physical_matching/A_m2/c5_coefficient.txt"))
    p.add_argument("--evaluator-wl",type=Path,default=Path("Lagrangian/EvaluateGeneralC5Numerically.wl"))
    p.add_argument("--output-json",type=Path,default=Path("output/general_A_m2_neutrino_provisional.json"))
    p.add_argument("--wolframscript",default="wolframscript")
    p.add_argument("--kappa-real",type=float,default=2e4)
    p.add_argument("--kappa-imag",type=float,default=1e4)
    p.add_argument("--quartic-real",type=float,default=0.02)
    p.add_argument("--quartic-imag",type=float,default=0.0)
    p.add_argument("--lambda-t3-real",type=float,default=None)
    p.add_argument("--lambda-t3-imag",type=float,default=0.0)
    p.add_argument("--hbar-value",type=float,default=1/(16*math.pi**2))
    p.add_argument("--mu-low",type=float,default=91.1876)
    p.add_argument("--ordering",choices=("NO","IO","AUTO"),default="AUTO")
    return run(p.parse_args())

if __name__=="__main__":
    raise SystemExit(main())
