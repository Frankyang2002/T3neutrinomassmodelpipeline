"""A(-4), D(-2): full RG + *conditional source-sector* neutrino objective.

UV(1e7) -> all three F(1e5) -> EFT1(S1,S2) -> both scalars(1e3)
-> SM+Weinberg(100 GeV).  The only C5 fed to SMEFT is the native LLHH
pole/log-closed *source sector* from Batch11. A signed overall multiplier
must be supplied explicitly; unity is NOT an independently verified value.

The UV Yukawa fit integrates the existing complete renormalizable RGBeta UV
and EFT1 ODEs on EVERY trial, with exactly 18 real UV y1/y2 entries. It is
NOT a physical fit: finite operator matching, full dim-5 anomalous dimension,
Matchete/Weyl sign and D-extra quartic normalization are incomplete.

Do not add Batch9's scalar log/bubble or Batch10 reference HHSS on top of
Batch11: those structures overlap with the native source-level calculation.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Mapping
import json
import math

import numpy as np
import sympy as sp
from scipy.optimize import least_squares

from Numerical.parallel_t3.complete.two_model_native_laurent_batch11 import (
    B10, FORMAT as B11_FORMAT, MODELS, NATIVE_STATUS, SOFT_PATHS,
    ONELOOP, ALL_SYMBOLS, EPS, MUBAR2, parse_laurent,
    combine_native, evaluate_soft, _candidate_soft_values,
)
from Numerical.parallel_t3.complete.two_model_hard_general_batch10 import (
    matrix, matrix_json, symmetric,
)
from Numerical.parallel_t3.background.general_batch1_api import (
    decode_native, encode_native, load_batch1, read_json,
)

FORMAT = 't3-thesis-a4-d2-conditional-uv-rg-neutrino-fit-batch12-v1'
B11 = Path('output/parallel_t3/thesis_native_laurent_batch11_v1')
B1 = Path('output/parallel_t3/full_general_batch1/thesis_full_general_batch1_v1')
B2 = Path('output/parallel_t3/full_general_batch2/thesis_full_general_batch2_raw_v1')
LAURENT = Path('output/parallel_t3/full_general_batch2/thesis_full_general_batch2_laurent_v2')
BENCHMARK = Path('configs/thesis_numerics/t3_thesis_A4_Dm2_pilot.json')
PARAMS = ('W12','W12bar','WExtra','WExtrabar','lam','lambar',
          'extra','extrabar','m1','m2','mminus','mplus','hbar','MUBAR2')
SOURCE_STATUS = 'NativeHardSoftLaurentAndScaleLogClosed_ConditionalOperatorConvention'


def require(ok: bool, reason: str) -> None:
    if not ok:
        raise ValueError('Batch12: '+reason)


def digest(path: Path) -> str:
    require(path.is_file(), 'Missing '+str(path))
    return sha256(path.read_bytes()).hexdigest()


def flatten_uv(y1: np.ndarray, y2: np.ndarray) -> np.ndarray:
    """Row-major 9+9 REAL entries, input defined at mu_UV=1e7 GeV."""
    mats = [np.asarray(z) for z in (y1,y2)]
    for x in mats:
        require(x.shape==(3,3) and np.isrealobj(x) and np.isfinite(x).all(),
                'UV fit expects two finite REAL 3x3 matrices')
    return np.concatenate([x.reshape(-1) for x in mats]).astype(float)


def unflatten_uv(x: np.ndarray) -> tuple[np.ndarray,np.ndarray]:
    v=np.asarray(x)
    require(v.shape==(18,) and np.isrealobj(v) and np.isfinite(v).all(),
            'UV fit requires exactly 18 finite REAL numbers')
    return np.asarray(v[:9],float).reshape(3,3), np.asarray(v[9:],float).reshape(3,3)


def residual5(obs: Mapping[str,Any],target: Mapping[str,Any]) -> np.ndarray:
    """NuFIT 6.1's existing five-quantity UNCORRELATED Gaussian approximation."""
    from Numerical.parallel_t3.complete.final_one_batch_fit import residual5 as source
    return source(dict(obs),dict(target))


def candidate_c5_as_smeft(c5: np.ndarray, signed_factor: float) -> np.ndarray:
    require(not isinstance(signed_factor,bool) and
            math.isfinite(float(signed_factor)) and float(signed_factor)!=0,
            'Signed GammaCC->Weyl convention factor must be explicit, finite, nonzero')
    c=np.asarray(c5,complex)
    require(c.shape==(3,3) and np.isfinite(c).all() and
            np.allclose(c,c.T,rtol=1e-10,atol=1e-24),
            'Candidate 3x3 source sector must be complex symmetric')
    return float(signed_factor)*c


class NativeSoftFast:
    """Vectorized exact SymPy Laurent coefficient, compiled once per model.

    Unlike repeating nine SymPy subs per optimizer trial, evaluate the original
    ordered native expression over a 3x3 flavour array using NumPy broadcasts.
    Sympy's parser still validates each source token before compiling.
    """
    def __init__(self, text:str, model:str):
        require(model in MODELS, 'Unsupported model')
        self.model=model
        reg,pole=parse_laurent(text,model)
        args=[ALL_SYMBOLS[p] for p in PARAMS]
        for label,expr in (('regular',reg),('pole',pole)):
            require(expr.free_symbols.issubset(set(args)),
                    'Unknown '+label+' native symbol')
        self._reg=sp.lambdify(args,reg,modules='numpy',cse=True)
        self._pole=sp.lambdify(args,pole,modules='numpy',cse=True)
        self.exact_pole_expression=str(pole)
        self.symbol_count=len(set(reg.free_symbols)|set(pole.free_symbols))

    def evaluate(self, *, w12:np.ndarray,wextra:np.ndarray,lam:complex,
                 extra:complex,m1:float,m2:float,mminus:float,mplus:float,
                 mu_bar:float) -> tuple[np.ndarray,np.ndarray]:
        W=np.asarray(w12,complex); X=np.asarray(wextra,complex)
        require(W.shape==X.shape==(3,3) and np.isfinite(W).all() and
                np.isfinite(X).all(), 'Bad three-flavour Wilson matrices')
        require(all(math.isfinite(float(s)) and float(s)>0 for s in
                    (m1,m2,mminus,mplus,mu_bar)),
                'Nonpositive scalar mass or matching scale')
        a=[W,np.conj(W),X,np.conj(X),complex(lam),complex(lam).conjugate(),
           complex(extra),complex(extra).conjugate(),
           m1,m2,mminus,mplus,ONELOOP,mu_bar**2]
        outs=[]
        for fn in (self._reg,self._pole):
            val=np.asarray(fn(*a),dtype=complex)
            if val.shape==(): val=np.full((3,3),val.item(),complex)
            require(val.shape==(3,3) and np.isfinite(val).all(),
                    'Vectorized source coefficient is not finite 3x3')
            outs.append(val)
        return outs[0],outs[1]


def _soft_params(fvalues:Mapping[str,Any],svalues:Mapping[str,Any],
                 model:str, *,mu_f:float) -> dict:
    """Same strict-one-loop frozen F couplings and S mass convention as B11."""
    m=np.asarray(decode_native(fvalues['MF']),complex)
    require(m.shape==(3,3),'MF must be 3x3')
    diag=np.diag(m)
    rel=float(np.linalg.norm(m-np.diag(diag))/max(np.linalg.norm(m),1e-300))
    require(rel<1e-8 and np.min(diag.real)>0 and
            np.max(np.abs(diag.imag))<1e-8*np.max(diag.real),
            'Heavy MF is not in the source-supported nearly-diagonal positive basis'
            ' (relative offdiag %.3g); do not infer an SVD flavour convention' % rel)
    y1=np.asarray(decode_native(fvalues['y1']),complex)
    y2=np.asarray(decode_native(fvalues['y2']),complex)
    require(y1.shape==y2.shape==(3,3) and np.isfinite(y1).all() and
            np.isfinite(y2).all(),'Need actual 3x3 threshold Yukawas')
    masses=diag.real
    w12=(y1/masses[None,:])@y2.T
    wextra=np.zeros_like(w12)
    extra=0j
    if model=='D_-2':
        y2conj=np.asarray(decode_native(fvalues['y2Conj']),complex)
        require(y2conj.shape==(3,3) and np.count_nonzero(y2conj)==0,
                'D2 has a nonzero y2Conj: extra F matching not implemented')
        wextra=(y1/masses[None,:])@y2conj.T
        extra=complex(np.asarray(decode_native(fvalues['lambdaT3GeneralQ020101I1'])))
    a=float(np.real(np.asarray(decode_native(svalues['mS1Sq']))))
    b=float(np.real(np.asarray(decode_native(svalues['mS2Sq']))))
    require(a>0 and b>0,'Scalar mass² is nonpositive')
    m1,m2=math.sqrt(a),math.sqrt(b)
    minus=plus=m2
    if model=='D_-2':
        delta=complex(np.asarray(decode_native(svalues['muT3GeneralS2Sq'])))
        require(abs(delta.imag)<1e-10*max(1.,abs(delta.real)),
                'D2 native split-real threshold does not cover a complex bilinear')
        require(b-2*delta.real>0 and b+2*delta.real>0,
                'D2 has a nonpositive split-real scalar eigenmass²')
        minus,plus=math.sqrt(b-2*delta.real),math.sqrt(b+2*delta.real)
    return {'w12':w12,'wextra':wextra,
            'lam':complex(np.asarray(decode_native(fvalues['lambdaT3']))),
            'extra':extra,'m1':m1,'m2':m2,
            'mminus':minus,'mplus':plus}


@dataclass
class Evaluation:
    c5_raw:np.ndarray
    pole_residual:float
    source_log:np.ndarray
    hard:np.ndarray
    scalar:np.ndarray
    observables:dict
    rg_nfev:dict
    full_smeft:dict|None


class SourceSectorEvaluator:
    def __init__(self, root:Path, model:str, *, signed_factor:float,
                 batch11_dir:Path=B11, rtol:float=1e-8,atol:float=1e-11):
        require(model in MODELS,'Unknown two-model label')
        self.root=Path(root).resolve()
        self.model=model
        self.signed_factor=float(signed_factor)
        self.rtol=rtol
        self.atol=atol
        self.b11_file=self.root/batch11_dir/model/'result.json'
        self.b11=read_json(self.b11_file)
        require(self.b11.get('status')==SOURCE_STATUS and
                self.b11.get('model')==model and
                self.b11.get('physical_MSbar_C5_claimed') is False,
                'Must use successful, explicitly conditional Batch11 result')
        soft_path=self.root/SOFT_PATHS[model]
        require(self.b11.get('native_source_hash_sha256')==digest(soft_path),
                'Native Matchete scalar source changed since Batch11')
        self.soft=NativeSoftFast(soft_path.read_text(encoding='utf-8-sig'),model)
        self.b1_folder=self.root/B1/model
        self.traj=load_batch1(self.b1_folder)
        require(self.traj.model==model and
                self.traj.manifest['scales_gev']=={'UV':1e7,'F':1e5,'S':1e3,'low':100},
                'Unexpected hierarchy in saved full general Batch1')
        self.cfg=read_json(self.root/BENCHMARK)
        self.target=read_json(self.root/self.cfg['fitting']['target'])
        self.scale=self.traj.manifest['scales_gev']
        self.ledger_path=(self.root/B2/model/'one_loop_operator_ledger.json')
        self.ledger=read_json(self.ledger_path)
        from Numerical.parallel_t3.complete.hard_llhh_source import extract_inventory
        self.native_hard_inventory=extract_inventory(model,self.ledger)
        self.evidence=None
        self._rg_payload=None

    def source_matrix(self, fvalues:Mapping[str,Any],svalues:Mapping[str,Any]):
        """Recompute the actual Matchete native source sector at trial thresholds."""
        from Numerical.parallel_t3.complete.hard_llhh_source import numerical_unsubtracted_source
        encoded={name:encode_native(value) for name,value in fvalues.items()}
        formal=numerical_unsubtracted_source(self.model,self.ledger,encoded,
                                            mu_f=float(self.scale['F']))
        hreg=symmetric(matrix(formal['baseline_ordered_holomorphic_regular_GeVinv'],
                              symmetric=False))
        hpol=symmetric(matrix(formal['baseline_ordered_holomorphic_pole_GeVinv'],
                              symmetric=False))
        params=_soft_params(encoded,{k:encode_native(v) for k,v in svalues.items()},
                            self.model,mu_f=float(self.scale['F']))
        sreg0,spol0=self.soft.evaluate(**params,mu_bar=float(self.scale['S']))
        sreg=symmetric(sreg0)
        spol=symmetric(spol0)
        joined=combine_native(hreg,hpol,sreg,spol,float(self.scale['F']),
                              float(self.scale['S']))
        require(joined['source_pole_closure_relative'] < 2e-7,
                'Native hard/soft poles no longer close for Yukawa trial')
        return joined,hreg,sreg

    def reference(self)->Evaluation:
        """Use precisely the saved, source-matched Batch1 point (NO rerun)."""
        f=self.traj.fermion_matching_input()['native_uv_rgbeta_couplings_at_muF']
        s=self.traj.scalar_matching_input()['native_eft1_rgbeta_couplings_at_muS']
        calc,h,soft=self.source_matrix(f,s)
        c5=calc['native_finite_source_sector']
        old=matrix(self.b11['source_finite_sector_candidate_GeVinv'])
        difference=float(np.linalg.norm(c5-old)/max(np.linalg.norm(old),1e-300))
        require(difference<2e-7,
                'Compiled source kernel does not reproduce completed Batch11 coefficient'
                ' (relative %.4g)'%difference)
        # The existing evidence interface contains strong native file hashes.
        from Numerical.parallel_t3.complete.full_general_completion import (
            existing_evidence,
        )
        evidence=existing_evidence(self.model,self.root/B1/self.model,
                                   self.root/B2/self.model,self.root/LAURENT/self.model)
        return self._observe(calc,h,soft,evidence,rg_nfev={},full=True)

    def _observe(self,calc:dict,h:np.ndarray,soft:np.ndarray,evidence:Any,
                 *,rg_nfev:dict,full:bool) ->Evaluation:
        from Numerical.parallel_t3.complete.full_general_completion import finish_smeft
        c5=candidate_c5_as_smeft(calc['native_finite_source_sector'],self.signed_factor)
        # Reuse actual SM + Weinberg RGE / charged-lepton diagonalisation.
        out=finish_smeft(evidence,c5,vev=float(self.cfg['vev_GeV']),target=self.target)
        return Evaluation(calc['native_finite_source_sector'],
                calc['source_pole_closure_relative'],
                calc['source_derived_eft1_scale_log'],h,soft,
                out['observables'],rg_nfev,out if full else None)

    def _prepare_full_uv(self):
        """Cache actual symbolic beta payloads; never regenerate them in a fit."""
        if self._rg_payload is not None:return self._rg_payload
        from Numerical.orchestration.GeneralRunningBenchmark import load_canonical_general_exports
        from Numerical.parallel_t3.background.runner import Options
        from scripts.RunT3ThesisGeneralBatch1 import thesis_overrides,read_benchmark
        symid=str(self.traj.manifest['symbolic_run_id'])
        uv,eft,provenance=load_canonical_general_exports(self.model,symid,generate=False)
        opts=Options(root=self.root,mu_uv=1e7,mu_f=1e5,mu_s=1e3)
        cfg=read_benchmark(self.root/BENCHMARK)
        names=list(self.traj.manifest['coupling_order']['uv'])
        overrides=thesis_overrides(cfg,self.model,names,opts)
        require(self.traj.manifest.get('eft1_boundary_deltas',{})=={},
                'UV trial needs independently supplied one-loop renormalizable F shifts')
        self._rg_payload=(uv,eft,overrides)
        return self._rg_payload

    def trial(self,x:np.ndarray,*,full:bool=False)->Evaluation:
        """Reintegrate complete UV and EFT1 renormalizable ODEs at every UV trial."""
        uv_payload,eft_payload,base=self._prepare_full_uv()
        from Numerical.orchestration.GeneralRunningBenchmark import initial_state
        from Numerical.running.UVRunner import run_uv_segment
        from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
        from Numerical.running.IntermediateScalarRunner import run_intermediate_scalar_segment
        from Numerical.core.RGBetaEvaluator import state_environment
        from Numerical.running.IntermediateScalarRGBetaEvaluator import intermediate_scalar_state_environment
        y1,y2=unflatten_uv(x)
        overrides=dict(base)
        overrides['y1']=encode_native(np.asarray(y1,complex))
        overrides['y2']=encode_native(np.asarray(y2,complex))
        source_model,uvstate=initial_state(uv_payload,mu_uv=1e7,mu_f=1e5,
                                         mu_s=1e3,overrides=overrides)
        require(source_model==self.model,'UV state has changed model')
        uv=run_uv_segment(uvstate,uv_payload,1e5,
                          save_scales_gev=np.geomspace(1e7,1e5,3),
                          rtol=self.rtol,atol=self.atol)
        require(uv.solver_success,'UV RGBeta integration failed')
        eftbound=build_intermediate_scalar_boundary(uv.final_state,
                       matching_scale_gev=1e5)
        eft=run_intermediate_scalar_segment(eftbound,eft_payload,1e3,
                          save_scales_gev=np.geomspace(1e5,1e3,3),
                          rtol=self.rtol,atol=self.atol)
        require(eft.solver_success,'EFT1 RGBeta integration failed')
        f=state_environment(uv.final_state)
        s=intermediate_scalar_state_environment(eft.final_state)
        joined,h,scalar=self.source_matrix(f,s)
        ev=SimpleNamespace(model=self.model,rg_at_ms=s,
                           ms_gev=1e3, sources_sha256={
                              'batch11':digest(self.b11_file),
                              'beta_source':self.traj.manifest.get('source_sha256',{}),
                           })
        return self._observe(joined,h,scalar,ev,
                             rg_nfev={'uv':int(uv.nfev),'eft1':int(eft.nfev)},
                             full=full)


def _report(e:Evaluation, *,signed_factor:float,detail:bool)->dict:
    ob=e.observables
    out={'status':'ConditionalNumericalSMEFTAndNeutrinosComputed_NotPhysical',
         'source_sector_C5_at_muS_GeVinv':matrix_json(e.c5_raw),
         'assumed_signed_GammaCC_to_Weyl_factor':signed_factor,
         'conditional_C5_at_muS_GeVinv':matrix_json(signed_factor*e.c5_raw),
         'source_pole_closure_relative':e.pole_residual,
         'observables':ob,'nufit_chi2_gaussian':None,
         'UV_EFT1_RG_nfev':e.rg_nfev}
    if e.full_smeft is not None:
        out['smeft_and_neutrino']=e.full_smeft
        out['nufit_chi2_gaussian']=e.full_smeft['nufit_gaussian_approximation']['chi2_gaussian_uncorrelated']
    if detail:
        out['source_pole_partition_log_GeVinv']=matrix_json(e.source_log)
        out['native_hard_regular_GeVinv']=matrix_json(e.hard)
        out['native_scalar_regular_GeVinv']=matrix_json(e.scalar)
    return out


def _fit(engine:SourceSectorEvaluator, *,max_nfev:int,
         yukawa_bound:float,penalty:float=1e4)->dict:
    from scripts.RunT3ThesisGeneralBatch1 import real_yukawas
    y1,y2=real_yukawas(engine.cfg,engine.model)
    x0=flatten_uv(y1,y2)
    require(max_nfev>=1,'max_nfev must be positive')
    require(math.isfinite(yukawa_bound) and
            1e-4<yukawa_bound<0.2 and np.max(abs(x0))<yukawa_bound,
            'Yukawa bound must exceed initial UV Yukawas and remain perturbative')
    # Check a newly integrated initial UV point before optimizer activation;
    # prevents fitting a different numerical model than saved Batch1.
    startup=engine.trial(x0,full=True)
    reference=engine.reference()
    compare=float(np.linalg.norm(startup.c5_raw-reference.c5_raw)/max(
                  np.linalg.norm(reference.c5_raw),1e-300))
    require(compare<5e-3,
            'Rerun UV/EFT1 baseline disagrees with saved Batch1 source-sector'
            ' (relative %.4g)'%compare)
    memo={};trials=[];best=None
    def residual(x):
        nonlocal best
        key=np.asarray(x,float).tobytes()
        if key in memo:return memo[key].copy()
        try:
            point=engine.trial(x,full=False)
            out=residual5(point.observables,engine.target)
            require(out.shape==(5,) and np.isfinite(out).all(),
                    'Nonfinite NuFIT 5-vector')
            stat='ConditionalTrialEvaluated'
        except Exception as exc:
            out=np.full(5,penalty,float)
            stat='InvalidTrial:'+type(exc).__name__+':'+str(exc)[:170]
        chi2=float(out@out)
        if best is None or (stat=='ConditionalTrialEvaluated' and chi2<best['chi2']):
            if stat=='ConditionalTrialEvaluated':
                best={'chi2':chi2,'UV_y1_y2_18_real':np.asarray(x,float).tolist(),
                      'residuals':out.tolist(),'trial_number':len(trials)+1}
        trials.append({'number':len(trials)+1,'status':stat,'chi2':chi2})
        memo[key]=out.copy()
        return out
    sol=least_squares(residual,x0,bounds=(-yukawa_bound,yukawa_bound),
                      max_nfev=max_nfev,x_scale=np.full(18,1e-4),
                      method='trf',ftol=1e-7,xtol=1e-7,gtol=1e-7)
    require(best is not None,'No valid conditional trial evaluated')
    best_eval=engine.trial(np.asarray(best['UV_y1_y2_18_real']),full=True)
    best_mats=unflatten_uv(np.asarray(best['UV_y1_y2_18_real']))
    return {'status':'Conditional18RealUVYukawaFit_NotPhysical',
            'best':{**best,'UV_y1':best_mats[0].tolist(),
                    'UV_y2':best_mats[1].tolist(),
                    'prediction':_report(best_eval,signed_factor=engine.signed_factor,
                                         detail=True)},
            'baseline_rerun_relative_C5':compare,
            'optimizer_success':bool(sol.success),
            'optimizer_message':str(sol.message),
            'iterations':int(sol.nfev),'unique_trial_evaluations':len(trials),
            'max_nfev_requested':max_nfev,
            'trial_status_count':{
                'valid':sum(x['status']=='ConditionalTrialEvaluated' for x in trials),
                'invalid':sum(x['status']!='ConditionalTrialEvaluated' for x in trials)},
            'trials':trials,
            'parameterization':'18 real entries of y1/y2 at mu_UV=1e7 GeV, row-major',
            'RG_reintegrated_at_every_unique_trial':True,
            'one_loop_matching_complete':False,
            'physical_prediction_claimed':False}


def run_one(root:Path,model:str,*,signed_factor:float,
            fit:bool=False,max_nfev:int=25,yukawa_bound:float=0.005,
            batch11_dir:Path=B11)->dict:
    engine=SourceSectorEvaluator(root,model,signed_factor=signed_factor,
                                 batch11_dir=batch11_dir)
    baseline=engine.reference()
    results={'format':FORMAT,'model':model,
             'status':'Conditional100GeVNeutrinoBenchmarkComputed_NotPhysicalC5',
             'hierarchy':[['UV_1e7','F_all_1e5'],['EFT1_S1_S2','S_joint_1e3'],
                          ['SMEFT','100GeV']],
             'source_sector_type':'native holomorphic LLHH F-hard + native joint-scalar-soft + one source-derived log; no additive B9/B10 proxies',
             'signed_normalization_assumption':signed_factor,
             'baseline':_report(baseline,signed_factor=signed_factor,detail=True),
             'input_source_sha256':{
                 'batch11':digest(engine.b11_file),
                 'soft':digest(root/SOFT_PATHS[model]),
                 'hard_ledger':digest(engine.ledger_path),
                 'batch1_manifest':digest(engine.b1_folder/'manifest.json')},
             'nufit_note':'Five uncorrelated Gaussian residuals; not official full NuFIT likelihood',
             'full_finite_matching_certified':False,
             'GammaCC_to_Weyl_normalization_certified':False,
             'D2_extra_matchete_rgbeta_map_certified':False,
             'full_EFT1_dim5_ADM_included':False,
             'UV_and_EFT1_reference_RG_from_existing_saved_batch1':True,
             'all_F_then_all_S_sequence_preserved':True,
             'physical_neutrino_prediction_claimed':False}
    if fit:
        results['conditional_18_real_UV_Yukawa_fit']=_fit(
             engine,max_nfev=max_nfev,yukawa_bound=yukawa_bound)
        results['status']='ConditionalUVRGRefitAnd100GeVNeutrinosComputed_NotPhysicalC5'
    return results
