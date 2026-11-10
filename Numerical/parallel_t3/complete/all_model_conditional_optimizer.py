"""One independent fit per T3 model, with a strict chi2 goal, using B11+B13+B15+B16
conditional source pipeline for A(-4) and D(-2).

No old analytic-triangle, fixed-background or rough-all16 fallback is allowed.
The other fourteen models have no equivalent certified source implementation:
those model jobs are explicitly BLOCKED, not silently approximated.

For each A/D trial:
    UV y1,y2 -> fresh UV RGBeta -> fresh EFT1 RGBeta,
    F native hard + simultaneous scalar-soft finite source sector (B11),
    source-exact F-tree -> full 45x45x12x12 EFT1 master tensor (B13+B15),
    source/reference quartic coordinate correction (B16),
    seven-term O(hbar) EFT1 Wilson beta at muF, with its log SUBSTITUTED
    for (not added to) the B11 log, -> SM+Weinberg RG -> Takagi / PMNS.

NOT full physical matching: Matchete GammaCC/Weyl normalization, complete
finite MSbar counterterms, complete EFT1 operator ADM, and some RGBeta quartic
basis directions remain uncertified.  The fit is CONDITIONAL.

This module does not modify existing source files or GitHub.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any, Callable, Mapping
import json
import math

import numpy as np
from scipy.optimize import least_squares

FORMAT = "t3-all16-one-fit-target-chi2-conditional-v2"
DEFAULT_TARGET_CHI2 = 0.01
SUPPORTED = ("A_-4", "D_-2")
SEEDS = ("A_-4", "D_-2")
BENCHMARK = Path("configs/thesis_numerics/t3_thesis_A4_Dm2_pilot.json")
B13 = Path("output/parallel_t3/thesis_master_wilson_batch13_cgfix_v3")
B16 = Path("output/parallel_t3/thesis_scalar_coordinate_batch16_v1")
OBS = ("delta_m21_sq_ev2", "delta_m3l_sq_ev2", "sin2_theta12",
       "sin2_theta13", "sin2_theta23")


class SourceBlocked(ValueError):
    """Cannot evaluate requested full-general source rather than approximating."""


def all_models() -> tuple[str, ...]:
    """Same 16-model inventory as the existing project, with no guessed labels."""
    from Numerical.parallel_t3.threshold.models import ALL_MODELS
    models = tuple(ALL_MODELS)
    if len(models) != 16 or len(set(models)) != 16 or not set(SUPPORTED) <= set(models):
        raise SourceBlocked("Unexpected project model inventory")
    return models


def jobs(models: tuple[str, ...]) -> list[dict]:
    """Exactly one independent model optimisation per requested model.

    The two stored Yukawa pairs are only alternative INITIAL CONDITIONS.
    They do not define two separate optimisation jobs.
    """
    full = set(all_models())
    if not models or any(m not in full for m in models) or len(set(models)) != len(models):
        raise SourceBlocked("Unknown, repeated or empty model set")
    return [{"model": m,
             "supported": m in SUPPORTED,
             "backend": "B11+B13+B15+B16_dynamic_Ohbar" if m in SUPPORTED else None,
             "reason": None if m in SUPPORTED else
             "No verified B11+B13+B15+B16 per-trial matching evaluator for this model"}
            for m in models]


def _hash_file(path: Path) -> str:
    if not path.is_file():
        raise SourceBlocked("Missing source file: " + str(path))
    return sha256(path.read_bytes()).hexdigest()


def _json(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    if not isinstance(obj, dict):
        raise SourceBlocked("Expected JSON object: " + str(path))
    return obj


def _encode_complex(matrix: np.ndarray) -> dict:
    a = np.asarray(matrix, complex)
    return {"real": a.real.tolist(), "imag": a.imag.tolist()}


def _complex_scalar(value: Any, label: str) -> complex:
    from Numerical.parallel_t3.complete.two_model_master_wilson_batch13 import decode_complex
    z = np.asarray(decode_complex(value), complex)
    if z.shape != () or not np.isfinite(z).all():
        raise SourceBlocked(label + " is not a finite scalar")
    return complex(z)


def _running_snapshot(model: str, svalues: Mapping[str, Any]):
    """Supply the ACTUAL trial muS masses/mu, not the baseline's masses."""
    from Numerical.joint_scalar_matching.numerical_benchmark_v5 import RunningSnapshot
    m1 = _complex_scalar(svalues["mS1Sq"], "mS1Sq")
    m2 = _complex_scalar(svalues["mS2Sq"], "mS2Sq")
    if (abs(m1.imag) > 1e-12*max(1.,abs(m1.real)) or
        abs(m2.imag) > 1e-12*max(1.,abs(m2.real)) or
        min(m1.real, m2.real) <= 0):
        raise SourceBlocked("Nonreal or nonpositive running scalar masses")
    c = {}
    for name, value in svalues.items():
        if name.startswith("lambda") or name.startswith("muT3"):
            z = _complex_scalar(value, name)
            c[name] = (float(z.real) if z.imag == 0 else
                       {"real": float(z.real), "imag": float(z.imag)})
    return RunningSnapshot(model=model, source=Path("trial-in-memory"),
                           sha256="trial-not-saved", mu_f=1e5, mu_s=1e3,
                           scalar_mass1_sq=float(m1.real),
                           scalar_mass2_sq=float(m2.real), couplings=c,
                           evidence={"source": "fresh UV+EFT1 RG trial",
                                     "native_matchete_scalar_coupling_aliases_certified": False})


@dataclass
class ConditionalPrediction:
    observables: dict
    chi2: float
    c5_source_muS: np.ndarray
    matching: dict
    rg_nfev: dict
    result: dict | None


class LatestTwoModelEvaluator:
    """Per-trial latest B11/B13/B15/B16 *conditional* evaluator.

    Reuses the Batch12 complete renormalizable RG integration, but REPLACES
    its obsolete B11-only C5 callback by recomputing the full tensor beta.
    """

    def __init__(self, root: Path, model: str, *, signed_factor: float = 1.0,
                 max_baseline_relative: float = 2e-6):
        if model not in SUPPORTED:
            raise SourceBlocked(model + ": latest full tensor backend unavailable")
        self.root = Path(root).resolve()
        self.model = model
        self.signed_factor = float(signed_factor)
        if not np.isfinite(self.signed_factor) or self.signed_factor == 0:
            raise SourceBlocked("Explicit nonzero finite signed factor required")
        from Numerical.parallel_t3.complete.two_model_conditional_uvfit_batch12 import (
            SourceSectorEvaluator)
        self._native_engine = SourceSectorEvaluator(self.root, model,
                                                     signed_factor=signed_factor)
        self.cfg = self._native_engine.cfg
        self.target = self._native_engine.target
        self.b13_path = self.root / B13 / model / "result.json"
        self.b16_path = self.root / B16 / model / "result.json"
        self.b13 = _json(self.b13_path)
        self.b16 = _json(self.b16_path)
        if (self.b13.get("status") != "FullTensorSevenChannelStrictOneLoopConditional_NotPhysicalC5"
                or self.b16.get("status") != "ScalarCoordinateSpanReconstructed_ConditionalNotPhysicalC5"
                or self.b13.get("model") != model or self.b16.get("model") != model):
            raise SourceBlocked("Need successful model-specific B13 and B16 source audit")
        self.beta_mode = self.b13["candidate_ADM"]
        self.master_factor = float(self.b13["master_to_native_factor_assumed"])
        if not np.isfinite(self.master_factor) or self.master_factor == 0:
            raise SourceBlocked("Bad formal master-to-native normalization")
        from Numerical.parallel_t3.complete.two_model_joint_scalar_batch9 import load_snapshot
        from Numerical.parallel_t3.complete.two_model_scalar_coordinate_batch16 import build_map
        from Numerical.parallel_t3.complete.native_triplet_basis_batch15 import heavy_source_from_reference
        self.snapshot = load_snapshot(self.root, model)
        self.seed = _json(self.snapshot.quartic_seed)
        self.coordinate_audit, self.coordinate_sectors, _ = build_map(self.seed, model)
        self.native_R, self.triplet_meta = heavy_source_from_reference(model, self.seed)
        self.source_hashes = {
            "B13": _hash_file(self.b13_path), "B16": _hash_file(self.b16_path),
            "scalar_quartic_seed": _hash_file(self.snapshot.quartic_seed),
            "B11": _hash_file(self._native_engine.b11_file),
        }
        if self.b13.get("input_sha256", {}).get("source_scalar_quartic") != self.source_hashes["scalar_quartic_seed"]:
            raise SourceBlocked("B13 and current Matchete scalar seed differ")
        self.count = 0
        self._cache: dict[bytes, ConditionalPrediction] = {}
        self.baseline_relative = self._check_baseline(max_baseline_relative)

    def _run_rg(self, x: np.ndarray):
        """Reuse Batch12's verified fresh UV/EFT1 integrators; no SM running here."""
        from Numerical.parallel_t3.complete.two_model_conditional_uvfit_batch12 import unflatten_uv
        from Numerical.orchestration.GeneralRunningBenchmark import initial_state
        from Numerical.running.UVRunner import run_uv_segment
        from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
        from Numerical.running.IntermediateScalarRunner import run_intermediate_scalar_segment
        from Numerical.core.RGBetaEvaluator import state_environment
        from Numerical.running.IntermediateScalarRGBetaEvaluator import intermediate_scalar_state_environment

        eng = self._native_engine
        uv_payload, eft_payload, baseline = eng._prepare_full_uv()
        y1,y2 = unflatten_uv(np.asarray(x,float))
        inputs = dict(baseline)
        from Numerical.parallel_t3.background.general_batch1_api import encode_native
        inputs["y1"] = encode_native(np.asarray(y1, complex))
        inputs["y2"] = encode_native(np.asarray(y2, complex))
        actual_model, uv0 = initial_state(uv_payload, mu_uv=1e7, mu_f=1e5,
                                          mu_s=1e3, overrides=inputs)
        if actual_model != self.model:
            raise SourceBlocked("Trial UV model differs from requested model")
        uv = run_uv_segment(uv0, uv_payload, 1e5,
                            save_scales_gev=np.geomspace(1e7, 1e5, 3),
                            rtol=eng.rtol, atol=eng.atol)
        if not uv.solver_success:
            raise SourceBlocked("Fresh UV RGBeta integration failed")
        boundary = build_intermediate_scalar_boundary(uv.final_state, matching_scale_gev=1e5)
        eft = run_intermediate_scalar_segment(boundary, eft_payload, 1e3,
                           save_scales_gev=np.geomspace(1e5, 1e3, 3),
                           rtol=eng.rtol, atol=eng.atol)
        if not eft.solver_success:
            raise SourceBlocked("Fresh EFT1 RGBeta integration failed")
        return (state_environment(uv.final_state),
                intermediate_scalar_state_environment(eft.final_state),
                {"UV": int(uv.nfev), "EFT1": int(eft.nfev)})

    def _master_log(self, f: Mapping[str,Any], s: Mapping[str,Any]) -> tuple[np.ndarray,dict]:
        from Numerical.parallel_t3.complete.two_model_master_wilson_batch13 import (
            source_to_master, sm_yukawa_tensor, gauge_sectors,
            scalar_quartic_from_matchete, add_d2_extra_hhss,
            add_d2_extra_scalar_only, projection)
        from Numerical.parallel_t3.complete.two_model_scalar_coordinate_batch16 import (
            rebuilt_quartic_delta, transform_tensor)
        from Numerical.running.GeneralT3DynamicFermionWilson import uv_from_json,nonmixed_boundary
        from Numerical.parallel_t3.background.general_batch1_api import encode_native
        from Numerical.running.Aminus4EFT1Couplings import F_LABELS
        from Numerical.parallel_t3.wilson.adm import beta, SYMMETRIC_TERMS

        fjson = {k:encode_native(v) for k,v in f.items()}
        uv = uv_from_json(self.model, fjson, mu_f=1e5)
        source, tree_info = nonmixed_boundary(self.model, uv, repo_root=self.root, mu_f=1e5)
        if source.raw.shape != (3,3,2,2,8,8):
            raise SourceBlocked("Wrong fresh 3-family F-tree tensor shape")
        # Check both native seeds refer to the same SU2 triplet frame.
        from Numerical.parallel_t3.complete.native_triplet_basis_batch15 import (
            triplet_cg_from_seed, generators_from_cg)
        seed_file = getattr(source, "seed_path", None) or source.path
        tree_seed = _json(Path(seed_file))
        g1=generators_from_cg(triplet_cg_from_seed(tree_seed))
        g2=generators_from_cg(triplet_cg_from_seed(self.seed))
        if max(float(np.linalg.norm(a-b)) for a,b in zip(g1,g2)) > 1e-10:
            raise SourceBlocked("F-tree and scalar quartic SU2 CG frames differ")
        R = np.eye(8)
        if self.model == "D_-2":
            from Numerical.joint_scalar_matching.special_mass_v6 import build_special_basis
            R = np.asarray(build_special_basis(self.model,
                        _running_snapshot(self.model,s),source)["R_original_from_boundary"],float)
        fermions=tuple(F_LABELS)
        c0=source_to_master(source.raw, fermions, R@self.native_R,self.model)
        quartic,_=scalar_quartic_from_matchete(
                    self.seed,f,self.model,R,native_source_from_reference=self.native_R)
        if self.model == "D_-2":
            add_d2_extra_hhss(quartic,self.root,f,R)
            add_d2_extra_scalar_only(quartic,f,R)
        correction,info=rebuilt_quartic_delta(self.coordinate_sectors,f)
        M=np.eye(12);M[4:,4:]=R
        quartic+=transform_tensor(correction,M)
        channels,backend=beta(c0,sm_yukawa_tensor(f,fermions),quartic,
                              gauge_sectors(f,self.model,R),variant=self.beta_mode)
        if set(channels)!=set(SYMMETRIC_TERMS):
            raise SourceBlocked("Lost a full-master one-loop RGE channel")
        dt=math.log(1e3/1e5)/(16*math.pi**2)
        delta=sum((dt*v for v in channels.values()),np.zeros_like(c0))
        log=projection(delta,fermions)
        if log.shape!=(3,3) or not np.isfinite(log).all():
            raise SourceBlocked("Invalid full-tensor EFT1 log")
        return log,{"tree_source":tree_info,
                    "coordinate_sector_modes":{k:v["mode"] for k,v in info.items()},
                    "candidate_ADM_backend":backend,
                    "master_log_norm_GeVinv":float(np.linalg.norm(log))}

    def _evaluate_states(self, f: Mapping[str,Any], s: Mapping[str,Any],
                         rg_nfev: dict, *, detail: bool) -> ConditionalPrediction:
        from Numerical.parallel_t3.complete.two_model_master_wilson_batch13 import replace_eft1_log
        from Numerical.parallel_t3.complete.two_model_conditional_uvfit_batch12 import candidate_c5_as_smeft
        from Numerical.parallel_t3.complete.full_general_completion import finish_smeft
        from types import SimpleNamespace
        native, hard, scalar = self._native_engine.source_matrix(f,s)
        master_log, audit = self._master_log(f,s)
        candidate = replace_eft1_log(native["native_finite_source_sector"],
                                      native["source_derived_eft1_scale_log"],
                                      master_log,self.master_factor)
        rg_state = SimpleNamespace(model=self.model,rg_at_ms=s,ms_gev=1e3,
                                    sources_sha256=self.source_hashes)
        result=finish_smeft(rg_state,
                            candidate_c5_as_smeft(candidate,self.signed_factor),
                            vev=float(self.cfg["vev_GeV"]),target=self.target)
        chi2=result["nufit_gaussian_approximation"]["chi2_gaussian_uncorrelated"]
        return ConditionalPrediction(result["observables"],float(chi2),candidate,
                                     {"native_hard_soft_pole_relative":native["source_pole_closure_relative"],
                                      "native_log_norm_GeVinv":float(np.linalg.norm(
                                          native["source_derived_eft1_scale_log"])),
                                      **audit},rg_nfev,result if detail else None)

    def _check_baseline(self, tol: float) -> float:
        """Fail closed unless this NEW evaluator reproduces saved Batch16 C5."""
        from Numerical.parallel_t3.complete.two_model_master_wilson_batch13 import decode_complex
        traj=self._native_engine.traj
        f=traj.fermion_matching_input()["native_uv_rgbeta_couplings_at_muF"]
        s=traj.scalar_matching_input()["native_eft1_rgbeta_couplings_at_muS"]
        test=self._evaluate_states(f,s,{},detail=False)
        stored=np.asarray(decode_complex(self.b16["conditional_C5_at_muS_GeVinv"]),complex)
        diff=float(np.linalg.norm(test.c5_source_muS-stored)/max(np.linalg.norm(stored),1e-100))
        if not np.isfinite(diff) or diff>tol:
            raise SourceBlocked("Latest trial evaluator does not reproduce B16 saved source: "
                                f"relative difference {diff:.6g}, limit {tol:.6g}")
        return diff

    def verify_fresh_baseline(self, x: np.ndarray, *, tol: float=5e-3) -> float:
        """Ensure a fresh RG integration reproduces the saved same-model B16 point."""
        from Numerical.parallel_t3.complete.two_model_master_wilson_batch13 import decode_complex
        point=self.evaluate(x,full=False)
        old=np.asarray(decode_complex(self.b16["conditional_C5_at_muS_GeVinv"]),complex)
        relative=float(np.linalg.norm(point.c5_source_muS-old)/max(np.linalg.norm(old),1e-100))
        if not np.isfinite(relative) or relative>tol:
            raise SourceBlocked(f"Fresh UV/EFT1 trial differs from saved B16 baseline: {relative:.6g} > {tol:.6g}")
        return relative

    def evaluate(self, x: np.ndarray, *, full: bool=False) -> ConditionalPrediction:
        x=np.asarray(x,float)
        if x.shape!=(18,) or not np.isfinite(x).all():
            raise SourceBlocked("Expected exactly 18 finite real Yukawa entries")
        key=x.tobytes()
        if key in self._cache and (not full or self._cache[key].result is not None):
            return self._cache[key]
        f,s,nfev=self._run_rg(x)
        pred=self._evaluate_states(f,s,nfev,detail=full)
        self.count+=1
        self._cache[key]=pred
        return pred


def vector_from_config(cfg: Mapping[str,Any], seed: str) -> np.ndarray:
    if seed not in SEEDS:
        raise SourceBlocked("Unsupported thesis seed")
    model=cfg["models"][seed]
    a=np.asarray(model["y1_UV_times_1e4"],float)
    b=np.asarray(model["y2_UV_times_1e4"],float)
    if a.shape!=(3,3) or b.shape!=(3,3) or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise SourceBlocked("Malformed original Yukawa benchmark")
    return np.r_[a.ravel(),b.ravel()]*float(cfg["fitting"]["input_Yukawa_scale"])


def choose_start(evaluate: Callable[[np.ndarray, bool], ConditionalPrediction],
                 cfg: Mapping[str, Any], *, policy: str = "best-of-two") -> tuple[str,np.ndarray,dict]:
    """Choose one starting point, never start two independent optimisations.

    For best-of-two, execute each of the two stored benchmark pairs exactly once
    through the requested MODEL evaluator, and start from lower unpenalised chi2.
    This selection is conditional; source physics is unchanged.
    """
    if policy not in ("best-of-two",)+SEEDS:
        raise SourceBlocked("Unknown start policy")
    choices = SEEDS if policy == "best-of-two" else (policy,)
    values = []
    for name in choices:
        x = vector_from_config(cfg, name)
        point = evaluate(x, full=False)
        chi2 = float(point.chi2)
        if not np.isfinite(chi2) or chi2 < 0:
            raise SourceBlocked("Nonfinite or negative benchmark chi2 for " + name)
        values.append((name,x,chi2))
    winner = min(values,key=lambda v:v[2])
    return winner[0],winner[1].copy(),{name:chi for name,_,chi in values}


class _TargetReached(Exception):
    """Early stop: a VALID physics evaluation passed the unpenalised target."""


def optimize_job(evaluate:Callable[[np.ndarray,bool],ConditionalPrediction],
                 x0:np.ndarray, *, max_nfev:int, bound:float,
                 regularization:float=0.,target_chi2:float=DEFAULT_TARGET_CHI2) -> dict:
    """One 18-real UV Yukawa fit, stopping immediately when chi2 < target_chi2.

    chi2 is the five-observable NuFIT uncorrelated Gaussian chi-square, NOT
    sqrt(chi2), nor a least-squares objective that includes regularisation.
    Even an early-stopped point is re-evaluated through the COMPLETE available
    conditional model pipeline at the end before a target PASS is reported.

    In the current project only A_-4/D_-2 have a B11+B13+B15+B16 backend;
    no use of this generic minimiser certifies the physical C5.
    """
    if not 1<=max_nfev<=10000:
        raise SourceBlocked("max_nfev must be 1..10000")
    if not np.isfinite(bound) or not np.isfinite(regularization) or bound<=0 or regularization<0:
        raise SourceBlocked("Invalid bound/regularization")
    if not np.isfinite(target_chi2) or not 0<target_chi2<1e10:
        raise SourceBlocked("target_chi2 must be strictly positive and finite")
    x0=np.asarray(x0,float)
    if x0.shape!=(18,) or not np.isfinite(x0).all() or np.any(abs(x0)>=bound):
        raise SourceBlocked("Seed outside strict 18-real Yukawa bounds")
    from Numerical.parallel_t3.complete.two_model_conditional_uvfit_batch12 import residual5
    parent=getattr(evaluate,"__self__",None)
    target=getattr(parent,"target",None)
    if not isinstance(target,Mapping):
        raise SourceBlocked("Expected bound evaluator method with a NuFIT target")
    memo:dict[bytes,np.ndarray]={}
    ledger:list[dict]=[]
    best=None

    def f(x, *, halt:bool=True):
        nonlocal best
        x=np.asarray(x,float)
        key=x.tobytes()
        if key in memo:
            out=memo[key].copy()
            if halt and best is not None and best["chi2"] < target_chi2:
                raise _TargetReached()
            return out
        try:
            p=evaluate(x,full=False)
            resid=np.asarray(residual5(p.observables,target),float)
            if resid.shape!=(5,) or not np.isfinite(resid).all():
                raise SourceBlocked("Nonfinite five-observable residual")
            chi2=float(resid@resid)
            if not np.isfinite(chi2):
                raise SourceBlocked("Nonfinite chi2")
            if not np.isclose(chi2,p.chi2,rtol=1e-6,atol=1e-9):
                raise SourceBlocked("Observable residual chi2 differs from benchmark chi2")
            status="valid"
            if best is None or chi2 < best["chi2"]:
                best={"chi2":chi2,"x":x.tolist(),"observable_residuals":resid.tolist()}
        except Exception as exc:
            resid=np.full(5,1e4,float)
            chi2=float(resid@resid)
            status="invalid:"+type(exc).__name__+":"+str(exc)[:160]
        reg=(np.sqrt(regularization)*(x-x0)/1e-4 if regularization>0
             else np.array([],float))
        out=np.r_[resid,reg]
        memo[key]=out.copy()
        ledger.append({"number":len(ledger)+1,"status":status,"chi2":chi2})
        # Do not allow scipy's next finite-difference/Jacobian evaluation once
        # the requested model has an acceptable point. No dummy 'convergence'.
        if halt and status=="valid" and best is not None and best["chi2"]<target_chi2:
            raise _TargetReached()
        return out

    initial=f(x0,halt=False)
    if not ledger or ledger[0]["status"]!="valid":
        raise SourceBlocked("Initial full-RG/full-tensor candidate failed: "+
                            (ledger[0]["status"] if ledger else "unknown"))
    start_chi2=float(initial[:5]@initial[:5])
    solver=None
    stopped_by_target=bool(best is not None and best["chi2"]<target_chi2)
    if not stopped_by_target:
        try:
            solver=least_squares(f,x0,method="trf",bounds=(-bound,bound),
                      max_nfev=max_nfev,x_scale=np.full(18,1e-4),
                      ftol=1e-7,xtol=1e-7,gtol=1e-7)
        except _TargetReached:
            stopped_by_target=True
    if best is None:
        raise SourceBlocked("No valid minimizer point")
    xbest=np.asarray(best["x"],float)
    # This is a FRESH full-detail callback for the best UV Yukawa matrices.
    # The production evaluator re-runs the actual UV/EFT1 RG integration.
    pbest=evaluate(xbest,full=True)
    verified_resid=np.asarray(residual5(pbest.observables,target),float)
    verified_chi2=float(verified_resid@verified_resid)
    if (verified_resid.shape!=(5,) or not np.isfinite(verified_resid).all()
        or not np.isclose(verified_chi2,best["chi2"],rtol=1e-5,atol=1e-8)
        or not np.isclose(verified_chi2,pbest.chi2,rtol=1e-6,atol=1e-9)):
        raise SourceBlocked("Best Yukawa candidate did not reproduce on final full benchmark run")
    target_met=bool(verified_chi2 < target_chi2)
    if target_met:
        status="ConditionalTargetChi2Reached_NotPhysicalC5"
        stop_reason="verified_chi2_below_target"
    elif solver is not None and solver.success:
        status="ConditionalSolverConverged_AboveTarget_NotPhysicalC5"
        stop_reason="solver_converged_but_target_not_reached"
    else:
        status="ConditionalBudgetExhausted_AboveTarget_NotPhysicalC5"
        stop_reason="max_nfev_or_solver_stop_without_target"
    return {
        "status":status,"stop_reason":stop_reason,
        "target_chi2_strictly_less_than":float(target_chi2),
        "target_met_after_fresh_benchmark_validation":target_met,
        "early_stop_triggered_inside_solver":stopped_by_target,
        "baseline_chi2":start_chi2,
        "best_chi2":verified_chi2,"verified_benchmark_chi2":verified_chi2,
        "best_18_real_UV_yukawas":xbest.tolist(),
        "best_y1_UV":xbest[:9].reshape(3,3).tolist(),
        "best_y2_UV":xbest[9:].reshape(3,3).tolist(),
        "best_five_observables":pbest.observables,
        "best_C5_muS_source_GeVinv":_encode_complex(pbest.c5_source_muS),
        "matching_audit":pbest.matching,"best_rg_nfev":pbest.rg_nfev,
        "optimizer_success":bool(solver.success) if solver is not None else False,
        "optimizer_message":(str(solver.message) if solver is not None else
                             "Stopped immediately at requested chi2 target"),
        "optimizer_nfev":int(solver.nfev) if solver is not None else None,
        "distinct_objective_evaluations":len(ledger),
        "valid_trials":sum(v["status"]=="valid" for v in ledger),
        "invalid_trials":sum(v["status"]!="valid" for v in ledger),
        "trials":ledger,"regularization_weight":float(regularization),
        "fit_parameter_count":18,"observable_count":5,
        "full_uv_and_eft1_rgbeta_reintegrated_each_trial":True,
        "latest_candidate_full_tensor_eft1_log_recomputed_each_trial":True,
        "complete_physical_MSbar_matching":False,
        "full_EFT1_ADM_certified":False,
        "GammaCC_to_Weyl_certified":False}
