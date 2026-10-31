"""Workstream B: validate LIVE B/C(-1) one-loop 5-real-scalar-mass closure.

This is a renormalisable mass-sector diagnostic. Dimensionless couplings are
held fixed in the optional short RG step; thresholds and physical C5 are NOT
computed. Source is the opt-in RGBeta Wolfram runner, NOT an inferred beta.

From the project root:
  python -X utf8 -m scripts.VerifyT3BCFullScalarMassClosure \
    --directory output/audits/workstream_b_bc_full_mass \
    --output output/audits/workstream_b_bc_full_mass/closure_v1.json
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Callable, Mapping

import numpy as np

MODELS = ("B_-1", "C_-1")
STAGES = ("UV", "EFT1")
NAMES = ("mHSq", "mS1Sq", "mS2Sq", "muT3GeneralS1S2")
FORMAT = "workstream-b-bc-full-scalar-mass-live-rgbeta-v1"
STATUS = "LiveHiggsAndMixedScalarMassBetasExtractedNotPhysicalCertificate"
CONVENTION = "16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g"
LOOP = 16.0 * math.pi**2
SCALE = 1.0e5  # finite-difference mass-squared step, GeV^2


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError("B/C complete scalar mass closure: " + reason)


def read_payload(path: Path) -> dict[str, Any]:
    require(path.is_file(), f"missing live beta export: {path}")
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    require(isinstance(obj, dict), f"invalid report object: {path}")
    return obj


def check_payload(data: Mapping[str, Any], model: str, stage: str) -> dict[str, str]:
    require(data.get("format") == FORMAT and data.get("status") == STATUS,
            f"{model}/{stage}: unknown live RGBeta export format/status")
    require(data.get("model") == model and data.get("stage") == stage,
            f"{model}/{stage}: model or stage differs")
    require(data.get("registered_mass_symbols") == list(NAMES),
            f"{model}/{stage}: mass registry does not contain exactly four couplings")
    require(data.get("physical_mapping_certified") is False and
            data.get("mass_eigenstate_threshold_matching_performed") is False,
            f"{model}/{stage}: false physical certification")
    meta = data.get("metadata")
    require(isinstance(meta, dict), f"{model}/{stage}: missing metadata")
    for k, v in {
        "SharedScalar": False, "dS1": 2, "dS2": 2,
        "dF": 1 if model == "B_-1" else 3, "alpha": -1,
        "GeneralScalarModel": model, "GeneralScalarStage": stage,
        "HiggsScalarMassSymbol": "mHSq", "RegisteredQuarticCount": 7,
        "MassBetaParameters": list(NAMES), "ReportBetaConvention": CONVENTION,
        "PhysicalMatcheteRGBetaMassBasisCertified": False,
        "DimensionlessCouplingsRunInThisReport": False,
        "FullRenormalizableMassSectorClosed": "PendingIndependentVerification",
    }.items():
        require(meta.get(k) == v, f"{model}/{stage}: invalid native metadata {k}")
    betas = data.get("report_betas")
    require(isinstance(betas, dict) and set(betas) == set(NAMES),
            f"{model}/{stage}: incomplete mass beta inventory")
    for name, expr in betas.items():
        require(isinstance(expr, str) and bool(expr.strip()) and
                not any(x in expr for x in ("BetaTerm[", "$Failed", "Missing[", "Indeterminate")),
                f"{model}/{stage}: unresolved {name} beta")
    require("mHSq" in betas["mHSq"] and
            "mS1Sq" in betas["mHSq"] and "mS2Sq" in betas["mHSq"],
            f"{model}/{stage}: Higgs beta is missing its mass or portal contributions")
    require("mHSq" in betas["mS1Sq"] and "mHSq" in betas["mS2Sq"],
            f"{model}/{stage}: scalar beta missing Higgs portal back-reaction")
    return dict(betas)


def benchmark_environment(model: str, stage: str) -> dict[str, Any]:
    """Deterministic illustrative parameters, NOT a fitted benchmark."""
    h = np.diag([1.1e-5, 0.0065, 0.94]).astype(complex)
    d = np.diag([2.4e-5, 0.00055, 0.024]).astype(complex)
    e = np.diag([2.9e-6, 0.00059, 0.0102]).astype(complex)
    z: dict[str, Any] = {
        "gY": 0.36, "g2": 0.65, "g3": 1.12,
        "yu": h, "yd": d, "ye": e,
        "lambdaH": 0.26, "lambdaS1": 0.22, "lambdaS2": 0.24,
        "lambdaH1": 0.07, "lambdaH2": 0.05, "lambda12": 0.025,
        "lambdaT3": 0.012 + 0.003j,
        "lambdaH1Adj": 0.035, "lambdaH2Adj": -0.018,
        "lambda12Adj": 0.015,
        "lambdaHHdagS2S2": 0.009 - 0.004j,
        "lambdaHHdagS1barS1bar": -0.007 + 0.003j,
        "lambdaS1bar2S2bar2": 0.011 + 0.004j,
        "lambdaS1barS2S2bar2": 0.006 - 0.002j,
        "lambdaS1S1bar2S2bar": -0.005 + 0.003j,
        "lambdaHHdagS1barS2barCross": 0.009 + 0.007j,
        "lambdaT3GeneralQ110101I1": -0.013 + 0.005j,
        "mHSq": -7800., "mS1Sq": 1.44e6, "mS2Sq": 1.96e6,
        "muT3GeneralS1S2": 4200. + 3100.j,
    }
    if stage == "UV":
        z.update({
            "MF": np.diag([750., 980., 1200.]).astype(complex),
            "y1": np.array([[0.08, 0.012, 0.0], [0.004, 0.065, 0.013],
                             [0., 0.008, 0.045]], dtype=complex),
            "y2": np.array([[0.045, -0.003, 0.001], [0.007, 0.075, 0.0],
                             [0.001, 0.003, 0.052]], dtype=complex),
        })
    return z


def real_vector(env: Mapping[str, Any]) -> np.ndarray:
    mu = complex(env[NAMES[-1]])
    return np.array([float(env[x]) for x in NAMES[:3]] + [mu.real, mu.imag], dtype=float)


def set_mass_vector(env: Mapping[str, Any], vector: np.ndarray) -> dict[str, Any]:
    x = np.asarray(vector, dtype=float)
    require(x.shape == (5,) and np.all(np.isfinite(x)), "invalid five-real mass vector")
    result = dict(env)
    for i, name in enumerate(NAMES[:3]):
        result[name] = float(x[i])
    result[NAMES[-1]] = complex(x[3], x[4])
    return result


def evaluate_mass_vector(betas: Mapping[str, str], env: Mapping[str, Any],
                         evaluator: Callable[[str, Mapping[str, Any]], Any] | None = None) -> np.ndarray:
    if evaluator is None:
        from Numerical.core.RGBetaEvaluator import evaluate_inputform_expression
        evaluator = evaluate_inputform_expression
    values = {name: evaluator(betas[name], env) for name in NAMES}
    out: list[float] = []
    for name in NAMES[:3]:
        z = complex(values[name])
        require(math.isfinite(z.real) and math.isfinite(z.imag) and
                abs(z.imag) < 1e-8 * max(1.0, abs(z.real)),
                f"nonreal or nonfinite diagonal mass beta {name}")
        out.append(z.real)
    mu = complex(values[NAMES[-1]])
    require(math.isfinite(mu.real) and math.isfinite(mu.imag), "invalid mixed mass beta")
    out.extend([mu.real, mu.imag])
    return np.asarray(out, dtype=float)


def affine_mass_closure(betas: Mapping[str, str], env: Mapping[str, Any],
                        evaluator: Callable[[str, Mapping[str, Any]], Any] | None = None) -> dict[str, Any]:
    """Check beta(x)=A x+b at fixed dimensionless couplings and fermion masses."""
    def f(x: np.ndarray) -> np.ndarray:
        return evaluate_mass_vector(betas, set_mass_vector(env, x), evaluator)
    zero = np.zeros(5)
    b = f(zero)
    a = np.column_stack([(f(SCALE*np.eye(5)[i])-b)/SCALE for i in range(5)])
    rng = np.random.default_rng(32421)
    residuals = []
    for v in (real_vector(env), rng.normal(size=5)*SCALE,
              rng.normal(size=5)*SCALE, -real_vector(env)):
        actual = f(v)
        predicted = a @ v + b
        res = float(np.max(np.abs(actual-predicted)))
        residuals.append(res)
        require(res < 2e-7*max(1.0,float(np.max(np.abs(actual))),
                                 float(np.max(np.abs(predicted)))),
                f"non-affine mass beta; max residual={res}")
    require(np.all(np.isfinite(a)) and np.all(np.isfinite(b)),
            "nonfinite mass anomalous-dimension matrix")
    # Trace-free SU(2) adjoint portals cannot source the singlet Higgs mass
    # from an identity-proportional heavy-doublet mass.
    for i, value in ((a[0,1],4*env["lambdaH1"]),
                     (a[0,2],4*env["lambdaH2"]),
                     (a[1,0],4*env["lambdaH1"]),
                     (a[2,0],4*env["lambdaH2"])):
        require(abs(i-value)<2e-7*max(1.,abs(value)),
                f"Higgs-doublet portal coefficient mismatch: {i} != {value}")

    # Restrict all extra complex B/C quartics, y1/y2 and mixed mass to zero;
    # independently verify the canonical 1-loop SM+two-doublet Higgs beta.
    sm = dict(env)
    for k in ("lambdaHHdagS2S2","lambdaHHdagS1barS1bar",
              "lambdaS1bar2S2bar2","lambdaS1barS2S2bar2",
              "lambdaS1S1bar2S2bar","lambdaHHdagS1barS2barCross",
              "lambdaT3GeneralQ110101I1","lambdaT3"):
        sm[k] = 0j
    sm[NAMES[-1]] = 0j
    if "y1" in sm:
        sm["y1"] = np.zeros((3,3), dtype=complex)
        sm["y2"] = np.zeros((3,3), dtype=complex)
    ref = (
        (6*sm["lambdaH"]-4.5*sm["g2"]**2-1.5*sm["gY"]**2
         + 6*np.real(np.trace(sm["yu"] @ sm["yu"].conj().T))
         + 6*np.real(np.trace(sm["yd"] @ sm["yd"].conj().T))
         + 2*np.real(np.trace(sm["ye"] @ sm["ye"].conj().T)))*sm["mHSq"]
        + 4*sm["lambdaH1"]*sm["mS1Sq"]
        + 4*sm["lambdaH2"]*sm["mS2Sq"])
    obtained = evaluate_mass_vector(betas, sm, evaluator)[0]
    require(abs(ref-obtained) < 2e-7*max(1.,abs(ref)),
            f"SM+two-doublet Higgs beta fails: live={obtained}, reference={ref}")

    # Optional illustrative mass-only flow, with all dimensionless couplings
    # and MF held constant. This is NOT full one-loop coupled RG running.
    from scipy.linalg import expm
    t = math.log(1000./1500.)
    aug = np.zeros((6,6))
    aug[:5,:5] = a/LOOP
    aug[:5,5] = b/LOOP
    start = real_vector(env)
    final = (expm(t*aug) @ np.r_[start,1.])[:5]
    return {
        "beta_mass_vector_gev2_per_16pi2": f(start).tolist(),
        "affine_beta_matrix_5_by_5": a.tolist(),
        "affine_source_gev2": b.tolist(),
        "affine_residual_max": max(residuals),
        "independent_zero_mixing_higgs_beta": {"live": obtained,"reference": float(ref)},
        "mass_only_frozen_dimensionless_step": {
            "from_scale_gev": 1500., "to_scale_gev": 1000.,
            "initial_mass_vector_gev2": start.tolist(),
            "final_mass_vector_gev2": final.tolist(),
            "dimensionless_running_included": False,
            "fermion_threshold_matching_included": False,
        },
    }


def audit(directory: Path,
          evaluator: Callable[[str, Mapping[str, Any]], Any] | None = None) -> dict[str, Any]:
    results = {}
    for model in MODELS:
        stages = {}
        for stage in STAGES:
            stem = model.replace("_-", "_m").replace("_+", "_p")
            path = directory / f"{stem}_{stage.lower()}_full_mass_v1.json"
            payload = read_payload(path)
            betas = check_payload(payload,model,stage)
            stages[stage] = affine_mass_closure(betas,benchmark_environment(model,stage),evaluator)
        results[model] = stages
    return {
        "status": "LiveFullFiveRealMassSectorAffineClosureCheckedFrozenCouplingsNotPhysicalMatching",
        "models":results,
        "renormalizable_scalar_mass_parameters": list(NAMES),
        "five_real_mass_components": True,
        "physical_matchete_rgbeta_phase_certified": False,
        "dimensionless_couplings_run": False,
        "mass_eigenstate_threshold_matching_performed": False,
        "neutrino_mass_prediction_performed": False,
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--directory",type=Path,default=Path("output/audits/workstream_b_bc_full_mass"))
    p.add_argument("--output",type=Path,default=Path("output/audits/workstream_b_bc_full_mass/closure_v1.json"))
    args = p.parse_args()
    require(not args.output.exists(),"refusing to overwrite an existing closure report")
    result = audit(args.directory)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    for model in MODELS:
        for stage in STAGES:
            r = result["models"][model][stage]
            print(model,stage,"Higgs beta (16pi²) =",r["beta_mass_vector_gev2_per_16pi2"][0],
                  "mass closure residual =",r["affine_residual_max"])
    print(result["status"])
    return 0


if __name__=="__main__":
    raise SystemExit(main())
