"""Batch 21: B(+1) conditional general-scalar neutrino fits for three studies.

Read-only historical scenario seeds. Reuses Batch20's independently exported
18+1 native scalar EFT1 equations. Frozen-UV uses no UV integration; full-UV
uses Stage31's coupled UV equations for EVERY trial. C5 comes from B(+1)'s
existing hierarchical Stage36 matching artifact, not the pole-bearing raw
one-scale Matchete expression. The extra quartic's *direct* one-loop C5
response was checked in Batch19, but its RGBeta/Matchete map and the omitted
Wilson-EFT evolution remain conditional.

  python -m scripts.FitT3BPlus1GeneralThreeStudies --all --run-id bp1_three_v1
  python -m scripts.FitT3BPlus1GeneralThreeStudies --case small_uv_frozen_5e-9 --run-id bp1_small_v1
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from types import SimpleNamespace
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
from scipy.stats import qmc

from scripts import RunT3BPlus1GeneralEFT1 as b20

ROOT = Path(__file__).resolve().parents[1]
MODEL = 'B_+1'
EXTRA = b20.NEW
STEM = 'T3_dS1_2_dS2_2_dF_1_alpha_p1.json'
CASES = {
    'small_uv_frozen_5e-9': ('output/studies/fermion_start_lambda_comparison/5e-9', False, 5e-9),
    'large_uv_frozen_0p5': ('output/studies/fermion_start_lambda_comparison/0p5', False, 0.5),
    'small_uv_running_5e-9': ('output/studies/lambda_5e-9_refit', True, 5e-9),
}
OUT = ROOT / 'output/studies/general_scalar_extended'
C5_REL = 'output/full/hypercharge/T3_B_alpha_p1/data/final_weinberg_coefficient.json'
BETA_RUN = 'output/studies/general_scalar_extended/b_p1_eft1_running/b_p1_eft1_general_v1'
UV_BETAS = 'output/general_t3_stage22/stage31_B_p1_full_uv_betas.json'
TARGET_PASS = 5.86
SCALAR_MIN_TEV = .5
SCALAR_MAX_TEV = 5.
RHO_MAX = .5
LOOP = 16. * math.pi**2
NAMES = tuple(f'{y}_{i}{j}' for y in ('y1', 'y2') for i in range(3) for j in range(3))
FORMAT = 'batch21-b-p1-three-scenario-conditional-general-scalar-fit-v1'


def require(condition: bool, reason: str) -> None:
    if not condition:
        raise ValueError('Batch21 B(+1) three scenarios: '+reason)


def read(path: Path) -> dict:
    require(path.is_file(), 'missing '+str(path))
    obj=json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(obj,dict), 'expected object in '+str(path))
    return obj


def sha(path: Path) -> str:
    require(path.is_file(),'missing provenance file '+str(path))
    return sha256(path.read_bytes()).hexdigest()


def scalar_complex(raw: Any) -> complex:
    require(isinstance(raw,dict) and set(raw)=={'real','imag'}, 'complex scalar requires real/imag')
    z=complex(float(raw['real']),float(raw['imag']))
    require(math.isfinite(z.real) and math.isfinite(z.imag),'nonfinite complex scalar')
    return z


def checked_run_id(value: str) -> str:
    require(isinstance(value,str) and bool(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',value)),
            'invalid run-id (1–64 alphanumeric/_/-)')
    return value


def input_files(case: str, root: Path) -> tuple[Path,Path,bool,float]:
    require(case in CASES, 'unknown study '+str(case))
    sub,active,lam=CASES[case]
    base=root/sub
    return base/'fitted'/STEM,base/'fits'/STEM,active,lam


def load_historical(case: str, root: Path) -> dict:
    cfgpath, fitpath, uv_on, expected_lambda=input_files(case,root)
    cfg,fit=read(cfgpath),read(fitpath)
    rep={'d_s1':2,'d_s2':2,'d_f':1,'alpha':1,'shared_scalar':False}
    require(cfg.get('representation')==rep,'historical fitted representation is not B(+1)')
    require(cfg.get('base_state',{}).get('mu_gev')==1e7,'historical input UV scale changed')
    require(cfg.get('scales')=={'mu_fermion_threshold_gev':1e5,
       'mu_scalar_threshold_gev':1e3,'mu_low_gev':1e2},'historical threshold plan changed')
    ordinary=cfg['base_state']['ordinary']
    lam=scalar_complex(ordinary['lambdaT3'])
    require(lam.imag==0 and math.isclose(lam.real,expected_lambda,rel_tol=1e-13),
            'historical topology boundary differs from original study')
    require(fit.get('model')==STEM[:-5] and
            math.isclose(float(fit['lambdaT3']),expected_lambda,rel_tol=1e-13),
            'historical fit record does not match original lambda')
    if uv_on:
        require('fermion_start_study' not in cfg and
                isinstance(cfg.get('lambda_5e9_study'),dict),
                'full UV study provenance changed')
    else:
        meta=cfg.get('fermion_start_study',{})
        require(meta.get('case')==('0p5' if expected_lambda==.5 else '5e-9')
                and 'UV RGEs disabled' in meta.get('boundary_definition',''),
                'frozen UV study provenance changed')
    from Numerical.orchestration import RunSmallLambdaRefit as base
    y=[base.matr(ordinary[name]) for name in ('y1','y2')]
    require(all(m.shape==(3,3) and np.max(abs(m.imag))<1e-12 for m in y),
            'fitter supports 18 real historical Yukawa entries')
    seed=np.concatenate([m.real.reshape(-1) for m in y])
    require(np.all(np.isfinite(seed)), 'nonfinite historical Yukawa input')
    return {'cfg':cfg,'fit':fit,'cfg_path':cfgpath,'fit_path':fitpath,
            'uv_on':uv_on,'lambda_fixed':expected_lambda,'seed':seed}


def options(*,sobol_points:int,local_nfev:int,seed_count:int,local_fraction:float,
            floor_fraction:float,abs_yukawa_max:float,extra_fraction:float|None,
            extra_real:float,
            rho_max:float,scalar_min_tev:float,scalar_max_tev:float) -> None:
    require(type(sobol_points) is int and 2<=sobol_points<=4096
            and sobol_points & (sobol_points-1)==0,
            '--sobol-points must be a power of two >=2')
    require(type(local_nfev) is int and 0<=local_nfev<=1000,
            '--max-nfev must be 0–1000 (0=pilot only)')
    require(type(seed_count) is int and 1<=seed_count<=8,
            '--seed-count must be 1–8')
    require(0<local_fraction<=3 and 0<floor_fraction<=1,
            'invalid local search fractions')
    require(0<abs_yukawa_max<=math.sqrt(4*math.pi),
            'Yukawa search absolute cap out of range')
    if extra_fraction is not None:
        require(0<extra_fraction<=.5,'extra fraction must be in (0,0.5]')
    require(math.isfinite(extra_real) and abs(extra_real)<=.5,
            'extra real boundary must be finite and |value| <= 0.5')
    require(0<rho_max<1 and 0<scalar_min_tev<scalar_max_tev,
            'invalid scalar/fermion screens')


def search_box(seed: np.ndarray, *, local_fraction:float, floor_fraction:float,
               abs_yukawa_max:float) -> tuple[np.ndarray,np.ndarray]:
    x=np.asarray(seed,dtype=float)
    require(x.shape==(18,) and np.isfinite(x).all(), 'expected 18 real seed Yukawas')
    require(np.max(np.abs(x))<abs_yukawa_max, 'historical Yukawa seed exceeds search bound')
    half=np.maximum(local_fraction*np.abs(x),floor_fraction*max(1e-6,float(np.max(abs(x)))))
    lo=np.maximum(-abs_yukawa_max,x-half)
    hi=np.minimum(abs_yukawa_max,x+half)
    require(np.all(lo<hi),'collapsed Yukawa search bounds; adjust --yukawa-abs-max')
    return lo,hi


def candidate_config(historical: dict, point: np.ndarray) -> dict:
    from Numerical.orchestration import RunSmallLambdaRefit as refit
    p=np.asarray(point,dtype=float)
    require(p.shape==(18,) and np.isfinite(p).all(),'expected 18 finite candidate parameters')
    cfg=deepcopy(historical['cfg'])
    ordinary=cfg['base_state']['ordinary']
    for i,key in enumerate(('y1','y2')):
        refit.setmat(ordinary[key],p[9*i:9*i+9].reshape(3,3).astype(complex))
    require(scalar_complex(ordinary['lambdaT3'])==complex(historical['lambda_fixed']),
            'trial altered original fixed topology coupling')
    cfg['benchmark_search']={'enabled':False}
    cfg['sensitivity']={'enabled':False}
    return cfg


def verify_beta_source(root: Path, run_id: str, state: Any):
    base=root/'output/studies/general_scalar_extended/b_p1_eft1_running'/run_id
    rec=read(base/'report.json')
    require(rec.get('format')==b20.RESULT_FORMAT and
            rec.get('status')=='SuccessCoupledBPlus1GeneralEFT1Running' and
            rec.get('model')==MODEL and rec.get('new_stage29_quartic_evolved_in_eft1') is True and
            rec.get('zero_control_vs_original_stage34',{}).get('status')==
            'ZeroExtraAgreesWithSavedStage34',
            'Batch20 live B(+1) report not validated or zero-control failed')
    path=base/'native_eft1_betas.json'
    hashes=rec.get('input_sha256',{})
    keys=[k for k in hashes if isinstance(k,str) and k.replace('\\','/').endswith('/native_eft1_betas.json')]
    require(len(keys)==1 and hashes[keys[0]]==sha(path),
            'Batch20 native beta source digest does not match report')
    payload=read(path)
    betas,layout=b20.validate_export(payload,state)
    return payload,layout,path,base/'report.json'


def prepare_uv(root: Path, config: dict, *, uv_on: bool):
    from scripts import Stage31CoupledUVRunning as s31
    from Numerical.core.StateVector import unpack_uv_state
    registry=read(root/'RGE/general/generated/general_rgbeta_quartic_registry.json')
    entry,new,reused=s31.inventories(registry,MODEL)
    require(new==(EXTRA,) and reused==(), 'B(+1) general quartic registry changed')
    vec,layout=s31.prepare_state(config,entry,(),{})
    if not uv_on:
        return unpack_uv_state(vec,layout,mu_gev=b20.MU_F),None
    beta_path=root/UV_BETAS
    beta=s31.validate_export(read(beta_path),entry,new,reused,tuple(layout.names))
    return (vec,layout,beta),beta_path


def run_uv(boundary: Any, q: complex):
    from scripts import Stage30CoupledA0UVRunning as s30
    from Numerical.core.StateVector import unpack_uv_state
    vec,layout,betas=boundary
    y0=np.concatenate((vec,[q.real,q.imag]))
    rhs=lambda t,y:s30.combined_rhs(t,y,layout=layout,names=(EXTRA,),betas=betas)
    sol=solve_ivp(rhs,(math.log(1e7),math.log(b20.MU_F)),y0,
                  method='DOP853',rtol=1e-8,atol=1e-11)
    require(sol.success and np.isfinite(sol.y).all() and sol.y.shape[0]==y0.size,
            'live coupled UV integration failed: '+sol.message)
    end=sol.y[:,-1]
    uv=unpack_uv_state(end[:layout.size],layout,mu_gev=b20.MU_F)
    return uv,complex(end[-2],end[-1]),int(sol.nfev)


def run_eft1(uv_state: Any, extra: complex, payload: dict, *, points:int=3):
    from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
    from Numerical.running.IntermediateScalarStateVector import pack_intermediate_scalar_state
    from Numerical.running.ScalarThresholdBoundary import scalar_threshold_masses
    intermediate=build_intermediate_scalar_boundary(uv_state,matching_scale_gev=b20.MU_F)
    betas,layout=b20.validate_export(payload,intermediate)
    initial=np.r_[pack_intermediate_scalar_state(intermediate)[0],extra.real,extra.imag]
    sol,_=b20.integrate(initial,layout,betas,points=points)
    _,end=b20.state_rows(sol,layout)
    masses=np.asarray(scalar_threshold_masses(end).masses_gev,dtype=float)
    require(masses.shape==(2,) and np.all(np.isfinite(masses)), 'invalid EFT1 scalar masses')
    return end,masses,complex(sol.y[-2,-1],sol.y[-1,-1]),int(sol.nfev),sol


class MatchedKernel:
    """Fast evaluation of the exact Stage36 B(+1) hard + direct kernel.

    Build from the verified model-specific symbolic matching source. Compare
    numerically against Stage36 at the historical seed (once per study).
    """
    def __init__(self,path:Path):
        import sympy as sp
        from RGE.matching.FinalWeinbergAdapter import load_hierarchical_majorana_c5
        from scripts.Stage36FinalWeinbergRunning import verify_hierarchical_artifact
        verify_hierarchical_artifact(path)
        symbolic=load_hierarchical_majorana_c5(path,n_lepton=3,n_heavy=3)
        expr=sp.simplify(symbolic['hard_kernel']+symbolic['running_kernel'])
        self.names=('MF','MS1','MS2','MS','lambdaT3','hbar')
        symbols=[sp.Symbol(x) for x in self.names]
        require(expr.free_symbols<=set(symbols),
                'C5 matching kernel contains an unexpected symbolic coupling')
        self.fun=sp.lambdify(symbols,expr,modules='numpy',cse=True)
        self.path=path

    def evaluate(self,uv:Any,eft:Any) -> np.ndarray:
        from Numerical.running.FinalC5TrajectoryAdapter import extract_final_c5_inputs
        view=SimpleNamespace(uv_threshold_state=uv,eft1_threshold_state=eft,
                             mu_scalar_threshold_gev=b20.MU_S)
        inp=extract_final_c5_inputs(view)
        args=inp.as_config()['t3']
        mf=np.asarray(args['MF'],dtype=float)
        weights=np.asarray([complex(self.fun(float(m),args['MS1'],args['MS2'],
                    args['MS'],args['lambdaT3'],args['hbar'])) for m in mf],dtype=complex)
        require(weights.shape==(3,) and np.isfinite(weights.real).all() and
                np.isfinite(weights.imag).all(), 'nonfinite matched hierarchical C5 kernel')
        y1=np.asarray(inp.y1,dtype=complex).conj()
        y2=np.asarray(inp.y2,dtype=complex).conj()
        c5=(y1*weights)@y2.T+(y2*weights)@y1.T
        require(c5.shape==(3,3) and np.isfinite(c5.real).all() and
                np.allclose(c5,c5.T,rtol=1e-10,atol=1e-17),
                'matched C5 nonfinite or not symmetric')
        return c5

    def parity_check(self,uv:Any,eft:Any) -> float:
        from Numerical.running.FinalC5TrajectoryAdapter import extract_final_c5_inputs
        from Numerical.running.SMWeinbergStage import evaluate_final_weinberg_json
        view=SimpleNamespace(uv_threshold_state=uv,eft1_threshold_state=eft,
                             mu_scalar_threshold_gev=b20.MU_S)
        official=evaluate_final_weinberg_json(self.path,
                         extract_final_c5_inputs(view).as_config())
        fast=self.evaluate(uv,eft)
        error=float(np.max(abs(fast-official))/max(1e-30,float(np.max(abs(official)))))
        require(error<3e-8,'fast model-specific kernel differs from Stage36 evaluator')
        return error


def diagnostics(uv:Any,eft:Any,masses:np.ndarray,*,rho_max:float,
                scalar_min_tev:float,scalar_max_tev:float,extra_at_muF:complex,
                extra_at_muS:complex) -> dict:
    from Numerical.running.IntermediateScalarRGBetaEvaluator import intermediate_scalar_state_environment
    from Numerical.running.FermionThresholdBoundary import fermion_threshold_masses
    mass_f=np.asarray(fermion_threshold_masses(uv).masses_gev,dtype=float)
    scalar_f=np.sqrt(np.array([uv.mS1Sq,uv.mS2Sq],dtype=float))
    require(np.all(mass_f>0) and np.all(scalar_f>0) and np.all(masses>0),
            'nonpositive fermion/scalar running mass')
    rho=float(np.max(scalar_f)/np.min(mass_f))
    ym=max(float(np.linalg.svd(np.asarray(y,dtype=complex),compute_uv=False)[0])
           for y in (uv.y1,uv.y2))
    env=intermediate_scalar_state_environment(eft)
    quartic=max([abs(complex(env[k])) for k in env if k.startswith('lambda')] +
                [abs(extra_at_muF),abs(extra_at_muS)])
    return {'rhoF':rho,'scalar_masses_muS_tev':np.sort(masses/1e3).tolist(),
            'fermion_masses_muF_tev':np.sort(mass_f/1e3).tolist(),
            'max_yukawa_singular_muF':ym,'max_scalar_quartic_abs':quartic,
            'mass_screen_pass':bool(rho<=rho_max and np.min(masses)>=1e3*scalar_min_tev
                                 and np.max(masses)<=1e3*scalar_max_tev),
            'perturbativity_screen_pass':bool(ym<=math.sqrt(4*math.pi) and quartic<=4*math.pi)}


def residual_from_fit(fit:Any) -> np.ndarray:
    v=np.asarray(fit.whitened_residual,dtype=float)
    require(v.shape==(5,) and np.all(np.isfinite(v)),
            'NuFIT target residual must have five finite entries')
    require(math.isclose(float(fit.chi2),float(v@v),rel_tol=1e-8,abs_tol=1e-10),
            'reported chi2 differs from target residual')
    return v


class Evaluator:
    def __init__(self,historical:dict,betas:dict,uv_betas:Any,kernel:MatchedKernel,
                 target:Any,*,root:Path,extra_fraction:float|None,extra_real:float,
                 rho_max:float,
                 scalar_min_tev:float,scalar_max_tev:float):
        self.historical=historical;self.betas=betas;self.uv_betas=uv_betas
        self.kernel=kernel;self.target=target;self.root=Path(root);self.extra_fraction=extra_fraction;self.extra_real=extra_real
        self.rho_max=rho_max;self.scalar_min_tev=scalar_min_tev
        self.scalar_max_tev=scalar_max_tev;self.cache={};self.best=None
        self.calls=0;self.failed=0;self.mass_safe=0;self.total_eft_nfev=0
        self.frozen=None;self.parity=None
        self.extra_boundary=complex(extra_real if extra_fraction is None else
                                   extra_fraction*historical['lambda_fixed'])
        if not historical['uv_on']:
            config=candidate_config(historical,historical['seed'])
            uv0,_=prepare_uv(self.root,config,uv_on=False)
            eft,masses,qS,nfev,_=run_eft1(uv0,self.extra_boundary,self.betas)
            self.frozen=(eft,masses,qS,nfev)
            self.total_eft_nfev+=nfev

    def evaluate(self,x:np.ndarray,*,label:str)->dict:
        p=np.asarray(x,dtype=float)
        key=p.tobytes()
        if key in self.cache:return self.cache[key]
        self.calls+=1
        row={'label':label,'chi2':None,'screen_pass':False,
             'parameters':{k:float(v) for k,v in zip(NAMES,p)}}
        try:
            cfg=candidate_config(self.historical,p)
            uv_in,_=prepare_uv(self.root,cfg,uv_on=self.historical['uv_on'])
            if self.historical['uv_on']:
                uv,qF,uv_nfev=run_uv(uv_in,self.extra_boundary)
                eft,masses,qS,eft_nfev,_=run_eft1(uv,qF,self.betas)
                self.total_eft_nfev+=eft_nfev
            else:
                uv=uv_in
                eft,masses,qS,eft_nfev=self.frozen
                qF=self.extra_boundary;uv_nfev=0
            checks=diagnostics(uv,eft,masses,rho_max=self.rho_max,
                      scalar_min_tev=self.scalar_min_tev,scalar_max_tev=self.scalar_max_tev,
                      extra_at_muF=qF,extra_at_muS=qS)
            row.update(checks)
            row.update({'lambdaT3_fixed':self.historical['lambda_fixed'],
                        'extra_quartic_at_start':{'real':self.extra_boundary.real,
                                                  'imag':self.extra_boundary.imag},
                        'extra_quartic_at_muS':{'real':qS.real,'imag':qS.imag},
                        'uv_solver_nfev':uv_nfev,'eft1_solver_nfev':eft_nfev})
            if checks['mass_screen_pass'] and checks['perturbativity_screen_pass']:
                if self.parity is None:
                    self.parity=self.kernel.parity_check(uv,eft)
                c5=self.kernel.evaluate(uv,eft)
                from scripts.RunT3GeneralScalarBCExploratoryFit import low_energy_fit
                fit,low=low_energy_fit(c5,eft,self.target,
                    vev=float(cfg.get('numerical',{}).get('vev_gev',246.22)))
                residual=residual_from_fit(fit)
                row.update({'status':'MassSafeMatchedNeutrinoTrial','screen_pass':True,
                    'chi2':float(fit.chi2),
                    'whitened_residual':residual.tolist(),
                    'prediction':dict(zip(fit.observable_names,map(float,fit.prediction))),
                    'c5_at_muS':{'real':c5.real.tolist(),'imag':c5.imag.tolist()},
                    'c5_at_mu100':{'real':low.K.real.tolist(),'imag':low.K.imag.tolist()}})
                self.mass_safe+=1
                if self.best is None or row['chi2']<self.best['chi2']:
                    self.best=deepcopy(row)
            else:
                row['status']='MassOrCouplingScreenRejected'
                row['whitened_residual']=[1e5]*5
        except (ValueError,TypeError,RuntimeError,OverflowError,ArithmeticError,
                FloatingPointError, np.linalg.LinAlgError) as exc:
            self.failed+=1
            row.update({'status':'NumericalTrialFailed','reason':str(exc)[:550],
                        'whitened_residual':[1e6]*5})
        self.cache[key]=row
        return row


def fit_case(case:str,run_id:str,*,root:Path=ROOT,beta_run_id:str='b_p1_eft1_general_v1',
             sobol_points:int=8,local_nfev:int=25,seed_count:int=1,
             random_seed:int=20261028,local_fraction:float=.8,
             floor_fraction:float=.03,abs_yukawa_max:float=2.,
             extra_fraction:float|None=None,extra_real:float=-.01,rho_max:float=RHO_MAX,
             scalar_min_tev:float=SCALAR_MIN_TEV,
             scalar_max_tev:float=SCALAR_MAX_TEV) -> dict:
    from scripts import Stage31CoupledUVRunning as s31
    from Numerical.orchestration import RunSmallLambdaRefit as base
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    options(sobol_points=sobol_points,local_nfev=local_nfev,
            seed_count=seed_count,local_fraction=local_fraction,
            floor_fraction=floor_fraction,abs_yukawa_max=abs_yukawa_max,
            extra_fraction=extra_fraction,extra_real=extra_real,rho_max=rho_max,
            scalar_min_tev=scalar_min_tev,scalar_max_tev=scalar_max_tev)
    run_id=checked_run_id(run_id)
    dest=(root/'output/studies/general_scalar_extended'/case/run_id).resolve()
    require(dest.is_relative_to((root/'output/studies/general_scalar_extended').resolve()) and
            not dest.exists(),'target already exists or escapes output root: '+str(dest))
    historical=load_historical(case,root)
    registry=read(root/'RGE/general/generated/general_rgbeta_quartic_registry.json')
    b20.validate_registry(registry)
    config=candidate_config(historical,historical['seed'])
    uv_init,_=prepare_uv(root,config,uv_on=False)
    from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
    initial=build_intermediate_scalar_boundary(uv_init,matching_scale_gev=b20.MU_F)
    payload,layout,beta_path,beta_report=verify_beta_source(root,beta_run_id,initial)
    uv_source=root/UV_BETAS
    if historical['uv_on']:
        s31.validate_export(read(uv_source),s31.inventories(registry,MODEL)[0],
                          (EXTRA,),(),tuple(prepare_uv(root,config,uv_on=True)[0][1].names))
    kernel_path=root/C5_REL
    kernel=MatchedKernel(kernel_path)
    target=OscillationFitTarget.from_json(base.TARGET)
    require(target.ordering=='NO' and len(target.observable_names)==5,
            'expected the same five-observable NO NuFIT target')
    evaluator=Evaluator(historical,payload,uv_source if historical['uv_on'] else None,
                       kernel,target,root=root,extra_fraction=extra_fraction,extra_real=extra_real,rho_max=rho_max,
                       scalar_min_tev=scalar_min_tev,scalar_max_tev=scalar_max_tev)
    lo,hi=search_box(historical['seed'],local_fraction=local_fraction,
                     floor_fraction=floor_fraction,abs_yukawa_max=abs_yukawa_max)
    trials=[(historical['seed'],'historical_seed')]
    generator=qmc.Sobol(d=18,scramble=True,seed=random_seed)
    for i,z in enumerate(generator.random_base2(int(math.log2(sobol_points)))):
        trials.append((lo+(hi-lo)*z,f'sobol_{i+1}'))
    valid=[];problems=[]
    for j,(p,label) in enumerate(trials):
        row=evaluator.evaluate(p,label=label)
        if row['chi2'] is not None:valid.append(row)
        elif len(problems)<15:problems.append({'trial':label,'status':row['status'],
                                            'reason':row.get('reason','screen failed')})
        print(f'B(+1) {case} {j+1}/{len(trials)}: {row["status"]} chi2={row["chi2"]}',flush=True)
    valid.sort(key=lambda r:r['chi2'])
    local=[]
    if valid and valid[0]['chi2']>TARGET_PASS and local_nfev>0:
        for start in valid[:seed_count]:
            x0=np.asarray([start['parameters'][k] for k in NAMES])
            def fun(z):
                x=lo+(hi-lo)*np.asarray(z)
                return np.asarray(evaluator.evaluate(x,label='least_squares')['whitened_residual'])
            z0=np.clip((x0-lo)/(hi-lo),1e-8,1-1e-8)
            sol=least_squares(fun,z0,bounds=(np.zeros(18),np.ones(18)),
               max_nfev=local_nfev,xtol=1e-6,ftol=1e-6,gtol=1e-6)
            local.append({'nfev':int(sol.nfev),'message':str(sol.message),
                          'success':bool(sol.success),'best_chi2':
                          evaluator.best['chi2'] if evaluator.best else None})
            print(f'B(+1) {case} local {len(local)}: best chi2='
                  f'{evaluator.best["chi2"] if evaluator.best else None}',flush=True)
            if evaluator.best and evaluator.best['chi2']<=TARGET_PASS:break
    best=evaluator.best
    status='CandidateConstrainedAccepted' if best and best['chi2']<=TARGET_PASS else (
           'CandidateConstrainedPoorFit' if best else 'NoMassSafeNeutrinoFit')
    files=[historical['cfg_path'],historical['fit_path'],beta_path,beta_report,
           kernel_path,base.TARGET,root/'RGE/general/generated/general_rgbeta_quartic_registry.json',
           root/'scripts/RunT3BPlus1GeneralEFT1.py',
           root/'scripts/FitT3BPlus1GeneralThreeStudies.py']
    if historical['uv_on']:files.append(uv_source)
    report={'format':FORMAT,'status':status,'case':case,'model':MODEL,
        'run_id':run_id,'uv_running_per_trial':historical['uv_on'],
        'historical_lambdaT3_fixed':historical['lambda_fixed'],
        'extra_quartic_boundary_rule':('fixed real RGBeta boundary' if extra_fraction is None
                                       else 'extra_fraction * abs(lambdaT3)'),
        'extra_quartic_boundary_real':float(evaluator.extra_boundary.real),
        'extra_fraction':extra_fraction,'conditional_one_scale_matching':True,
        'hierarchical_existing_c5_artifact_reused':True,
        'matched_wilson_eft1_evolved':False,'physical_quartic_normalization_certified':False,
        'screens':{'rhoF_max':rho_max,'scalar_min_tev':scalar_min_tev,
                   'scalar_max_tev':scalar_max_tev,'chi2_max':TARGET_PASS},
        'search':{'sobol_points':sobol_points,'max_nfev':local_nfev,
                  'seed_count':seed_count,'seed':random_seed,
                  'local_fraction':local_fraction,'floor_fraction':floor_fraction,
                  'yukawa_abs_max':abs_yukawa_max},
        'counts':{'candidate_evaluations':evaluator.calls,'mass_safe_c5_evaluations':evaluator.mass_safe,
                  'numerical_failures':evaluator.failed,
                  'total_eft1_rhs_evaluations':evaluator.total_eft_nfev},
        'matching_kernel_stage36_parity_relative_error':evaluator.parity,
        'best_mass_safe':best,'local_optimisation':local,'trial_problems':problems,
        'source_sha256':{str(p.relative_to(root)) if p.is_relative_to(root) else str(p):sha(p)
                         for p in files},
        'historical_results_overwritten':False,
        'limitations':['One-scale grouped C5 evaluation from existing model-specific hierarchical artifact',
                       'Native extra quartic included in UV/EFT1 renormalisable running only',
                       'No new Wilson coefficient beta functions or sequential general-scalar matching',
                       'Matchete-RGBeta physical coupling map not certified',
                       '18 real Yukawa fit variables to five observed values: no predictive claim']}
    dest.mkdir(parents=True,exist_ok=False)
    (dest/'fit.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    if best:
        fitted=candidate_config(historical,np.array([best['parameters'][k] for k in NAMES]))
        fitted['general_scalar_candidate_fit']={'case':case,'model':MODEL,
          'status':status,'chi2':best['chi2'],'run_id':run_id,
          'physical_matchete_rgbeta_map_certified':False,
          'sequential_matching_performed':False}
        (dest/'fitted_config.json').write_text(json.dumps(fitted,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    (dest/'summary.md').write_text(summary(report),encoding='utf-8')
    print('Saved:',dest/'summary.md',flush=True)
    return report


def summary(report:dict)->str:
    best=report['best_mass_safe']
    s=['# Batch 21 — B(+1) general-scalar conditional neutrino fit','',
       f"Study: `{report['case']}`; status: `{report['status']}`.",
       f"UV prescription: {'full coupled UV per trial' if report['uv_running_per_trial'] else 'frozen UV at muF'}.",
       f"Original fixed lambdaT3: `{report['historical_lambdaT3_fixed']:.8g}`; "
       f"new scalar quartic boundary: `{report['extra_quartic_boundary_real']:.8g}` "
       f"({report['extra_quartic_boundary_rule']}).",
       f"Candidate evaluations: `{report['counts']['candidate_evaluations']}`; "
       f"mass-safe neutrino evaluations: `{report['counts']['mass_safe_c5_evaluations']}`; "
       f"failed numerical evaluations: `{report['counts']['numerical_failures']}`.",
       f"Hierarchical Stage36 C5 fast-kernel parity relative error: "
       f"`{report['matching_kernel_stage36_parity_relative_error']}`.",'']
    if best:
        s.extend([f"Best screened chi2: **{best['chi2']:.9g}**",
                  f"rhoF: `{best['rhoF']:.8g}`",
                  'Scalar running masses at muS (TeV): '+', '.join(f'{v:.9g}' for v in best['scalar_masses_muS_tev']),
                  f"Max Yukawa singular muF: `{best['max_yukawa_singular_muF']:.8g}`",
                  f"Max scalar quartic coordinate at muS or extra-boundary endpoints: `{best['max_scalar_quartic_abs']:.8g}`",
                  '', '| Observable | Prediction |','|---|---:|',
                  *[f'| `{k}` | {v:.11g} |' for k,v in best['prediction'].items()],''])
    else:s.extend(['No mass-safe fit reached the neutrino evaluator.',''])
    s+=['## Interpretation','',
        'The new complex B(+1) quartic is integrated with the full native RGEs in UV (when enabled) and EFT1.',
        'The existing B(+1) hard-plus-direct C5 matching source is reused; no new sequential general-scalar matching is claimed.',
        'The Matchete-RGBeta extra-quartic normalization and physical vacuum stability are not certified.',
        'The five-variable NuFIT chi2 is a fit diagnostic, not a predictivity statement.',
        'Earlier studies and Batch20 artifacts are read-only.','']
    return '\n'.join(s)+'\n'


def main()->None:
    ap=argparse.ArgumentParser(description=__doc__)
    group=ap.add_mutually_exclusive_group(required=True)
    group.add_argument('--case',choices=tuple(CASES))
    group.add_argument('--all',action='store_true')
    ap.add_argument('--run-id',required=True)
    ap.add_argument('--beta-run-id',default='b_p1_eft1_general_v1')
    ap.add_argument('--sobol-points',type=int,default=8)
    ap.add_argument('--max-nfev',type=int,default=25)
    ap.add_argument('--seed-count',type=int,default=1)
    ap.add_argument('--seed',type=int,default=20261028)
    ap.add_argument('--local-fraction',type=float,default=.8)
    ap.add_argument('--floor-fraction',type=float,default=.03)
    ap.add_argument('--yukawa-abs-max',type=float,default=2.)
    ap.add_argument('--extra-real',type=float,default=-.01,
        help='Illustrative fixed B(+1) RGBeta extra-quartic boundary; default -0.01')
    ap.add_argument('--extra-fraction',type=float,
        help='Optional alternative: extra quartic = this factor * abs(lambdaT3)')
    ap.add_argument('--rho-max',type=float,default=.5)
    ap.add_argument('--scalar-min-tev',type=float,default=.5)
    ap.add_argument('--scalar-max-tev',type=float,default=5.)
    a=ap.parse_args()
    cases=tuple(CASES) if a.all else (a.case,)
    try:
        for case in cases:
            fit_case(case,a.run_id,beta_run_id=a.beta_run_id,
              sobol_points=a.sobol_points,local_nfev=a.max_nfev,
              seed_count=a.seed_count,random_seed=a.seed,
              local_fraction=a.local_fraction,floor_fraction=a.floor_fraction,
              abs_yukawa_max=a.yukawa_abs_max,extra_fraction=a.extra_fraction,
              extra_real=a.extra_real,
              rho_max=a.rho_max,scalar_min_tev=a.scalar_min_tev,
              scalar_max_tev=a.scalar_max_tev)
    except (ValueError,TypeError,RuntimeError,OSError,KeyError,OverflowError,
            json.JSONDecodeError) as exc:
        ap.exit(2,str(exc)+'\n')

if __name__=='__main__': main()
