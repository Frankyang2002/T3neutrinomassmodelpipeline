"""v28: B/C(-1) *joint* two-scalar matching-scale placement from candidate EFT1 RGEs.

UV -> integrate out ALL heavy fermions -> keep BOTH scalars in EFT1 ->
    integrate out ALL heavy scalars at ONE scale -> SMEFT.

Only renormalisable candidate RGE transport is performed, never scalar or Wilson
matching. The optional 5-real-mass backend uses the already-source-verified
RGBeta betas in RunT3BCFullScalarMassEFT1. The 4-mass-coordinate quadratic
backend is an explicitly downgraded fallback when full Higgs beta files are
unavailable. Both use the SAME live saved source, not invented beta functions.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

MODELS = ("B_-1", "C_-1")
FORMAT = 't3-parallel-a-bc-joint-scalar-eft1-threshold-state-v28'
STATUS = 'CandidateEFT1BothScalarsRunningAtOneCommonScale_NoPhysicalScalarMatching'
ORDER = 'UV -> integrate out ALL fermions -> integrate out ALL scalars -> SMEFT'

class JointScaleBlocked(ValueError):
    pass


def need(cond: Any, message: str) -> None:
    if not bool(cond):
        raise JointScaleBlocked('Workstream A v28: ' + message)


def sha(path: Path) -> str:
    need(path.is_file(), 'missing source file: '+str(path))
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path: Path) -> dict:
    need(path.is_file(), 'missing input file: '+str(path))
    result = json.loads(path.read_text(encoding='utf-8-sig'), parse_constant=lambda s: (_ for _ in ()).throw(JointScaleBlocked('nonfinite JSON '+s)))
    need(isinstance(result, dict), 'JSON object required')
    return result


def checked_complex(x: Any, name: str) -> complex:
    need(isinstance(x, dict) and set(x)=={'real','imag'}, 'bad complex pair '+name)
    z=complex(float(x['real']),float(x['imag']))
    need(math.isfinite(z.real) and math.isfinite(z.imag),'nonfinite '+name)
    return z


def joint_scale_from_running(t_low: float, t_high: float,
                             eigenvalues: Callable[[float], np.ndarray],
                             *, scheme: str = 'geometric', samples: int = 81) -> dict:
    """Locate mu_S = sqrt(mX(mu_S)*mY(mu_S)) (or mY, when requested).

    Require one isolated, bracketed sign crossing. The "upper" prescription
    matches both scalars simultaneously at the Y mass; no intermediate EFT2.
    """
    need(scheme in ('geometric','upper'),'unknown joint scalar scale prescription')
    need(math.isfinite(t_low) and math.isfinite(t_high) and t_high > t_low,
         'invalid EFT1 logarithmic scale interval')
    need(11<=samples<=1001,'invalid root bracket sampling')
    evaluation_count=0
    def evaluate(t:float)->tuple[float,np.ndarray]:
        nonlocal evaluation_count
        evaluation_count+=1
        eig=np.asarray(eigenvalues(t),dtype=float)
        need(eig.shape==(2,) and np.isfinite(eig).all() and np.all(eig>0)
             and eig[0]<=eig[1],'tachyonic, nonfinite or unordered mixed-scalar eigenvalues')
        log_masses=.5*np.log(eig)
        target=.5*(log_masses[0]+log_masses[1]) if scheme=='geometric' else log_masses[1]
        return t-target,eig
    grid=np.linspace(t_low,t_high,samples)
    values=np.asarray([evaluate(float(t))[0] for t in grid])
    crosses=np.flatnonzero(values[:-1]*values[1:]<0)
    zeros=np.flatnonzero(np.abs(values)<1e-12)
    need(len(crosses)+len(zeros)==1, 'matching-scale root missing or non-unique in EFT1 interval')
    if len(zeros)==1:
        tstar=float(grid[zeros[0]])
    else:
        i=int(crosses[0]);tstar=brentq(lambda t:evaluate(float(t))[0],float(grid[i]),float(grid[i+1]),xtol=2e-13,rtol=1e-13)
    residual,eig=evaluate(tstar)
    mu=math.exp(tstar)
    mx,my=np.sqrt(eig)
    need(abs(residual)<1e-10 and t_low<=tstar<=t_high,'joint root residual/interval invalid')
    return {'matching_scale_gev':mu,'scheme':scheme,
            'self_consistent_root_log_residual':abs(residual),
            'mX_gev':float(mx),'mY_gev':float(my),
            'scalar_mass_squared_eigenvalues_gev2':eig.tolist(),
            'log_mu_over_mX':math.log(mu/mx),
            'log_mu_over_mY':math.log(mu/my),
            'log_mY_over_mX':math.log(my/mx),
            'max_absolute_scalar_matching_log':max(abs(math.log(mu/mx)),abs(math.log(mu/my))),
            'sample_count':samples,'function_evaluations':evaluation_count,
            'two_scalar_eigenstates_both_retained_in_EFT1_until_matching':True,
            'intermediate_single_scalar_EFT2_used':False}


def compare_endpoint(actual: np.ndarray, reference: np.ndarray, *, tol:float=2e-6)->float:
    a=np.asarray(actual,dtype=float);b=np.asarray(reference,dtype=float)
    need(a.shape==b.shape and a.ndim==1 and np.isfinite(a).all() and np.isfinite(b).all(),
         'saved EFT1 endpoint array invalid')
    err=float(np.linalg.norm(a-b)/max(1.,np.linalg.norm(b)))
    need(err<tol, 'replayed EFT1 1 TeV endpoint does not match saved candidate result')
    return err


def _imports():
    from scripts import RunGeneralT3QuadraticEFT1 as qeft
    from scripts import RunT3BCFullScalarMassEFT1 as full
    from scripts import IntermediateEFTRunning as s34
    from scripts import FermionThresholdHandoff as s33
    from scripts import VerifyT3BCFullScalarMassClosure as closure
    from Numerical.experimental.GeneralT3ScalarRGState import (
        GeneralT3ScalarRGState,load_registry,pack_general_scalar_state,
        registry_model,unpack_general_scalar_state)
    from Numerical.experimental.GeneralT3ScalarQuadraticMasses import load_manifest,quadratic_kind,quadratic_mass_hessian
    from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis,verify_su2_covariance
    return (qeft,full,s34,s33,closure,GeneralT3ScalarRGState,load_registry,
            pack_general_scalar_state,registry_model,unpack_general_scalar_state,
            load_manifest,quadratic_kind,quadratic_mass_hessian,MixedDoubletMassBasis,verify_su2_covariance)


def verify_saved_beta_sources(saved: dict, files: dict[str,Path]) -> None:
    """Recheck live beta/registry/manifest hashes recorded by older Windows/POSIX qeft.

    Compare full normalized paths; unlike a basename lookup, collisions and
    unrelated source with the same filename cannot pass this check.
    """
    recorded=saved.get('input_sha256')
    need(isinstance(recorded,dict),'missing saved quadratic EFT1 source hashes')
    normalized={}
    for key,value in recorded.items():
        need(isinstance(key,str) and isinstance(value,str) and len(value)==64 and
             all(c in '0123456789abcdef' for c in value),
             'invalid saved EFT1 source digest')
        canonical=key.replace('\\','/').casefold()
        need(canonical not in normalized,'duplicate normalized source hash paths')
        normalized[canonical]=value
    for name,path in files.items():
        key=str(path).replace('\\','/').casefold()
        need(normalized.get(key)==sha(Path(path)),
             'saved quadratic EFT1 source changed: '+name)


def _stage33_or_saved_qeft(model: str, backend: str, root: Path, *,
                           stage33_dir: Path, qeft_dir: Path,
                           mass_dir: Path, higgs_mass_squared:float,
                           mu2_fraction:float, mu2_phase:float):
    (qeft,full,s34,s33,closure,GeneralState,load_registry,pack,reg_model,unpack,
     load_manifest,quadratic_kind,hessian,Basis,gauge_check)=_imports()
    registry=load_registry(qeft.REGISTRY)
    manifest=load_manifest(qeft.MANIFEST)
    sources={'registry':Path(qeft.REGISTRY),'manifest':Path(qeft.MANIFEST)}
    tag=s34.safe_model(model)
    beta_q=Path(qeft_dir)/f'{tag}_eft1_quadratic_live_rgbeta.json'
    sources['native_full_quartic_quadratic_beta']=beta_q
    if backend=='five_mass':
        _,baseline,sidecar,reused,hashes=s34.verify_boundary(model,stage33_dir)
        info=reg_model(registry,baseline.representation)
        need(quadratic_kind(manifest,model)=='mixed', 'five_mass B/C quadratic kind changed')
        need(set(sidecar)=={x['name'] for x in info['new_rgbeta_symbols']}
             and set(reused)=={x['existing_name'] for x in info['preexisting_rgbeta_symbols']},
             'incomplete Stage33 quartic couplings')
        state=GeneralState(baseline,{k:qeft.decode_complex(z,k) for k,z in sidecar.items()}).validated(registry)
        ybase,layout=pack(state,registry)
        need(min(baseline.mS1Sq,baseline.mS2Sq)>0,'nonpositive Stage33 scalar masses')
        mu2=mu2_fraction*min(baseline.mS1Sq,baseline.mS2Sq)*complex(math.cos(mu2_phase),math.sin(mu2_phase))
        y0=np.r_[ybase,mu2.real,mu2.imag,higgs_mass_squared]
        report=qeft.checked_payload(model,baseline.representation,'mixed',layout,registry,beta_q)
        mass_file=Path(mass_dir)/f'{tag}_eft1_full_mass_v1.json'
        sources['native_five_scalar_mass_betas']=mass_file
        mass_payload=closure.read_payload(mass_file)
        betas=full.merge_beta_reports(report,closure.check_payload(mass_payload,model,'EFT1'),layout.names)
        rhs=lambda t,y:full.beta_derivative(np.asarray(y),mu_gev=math.exp(t),layout=layout,registry=registry,betas=betas)
        saved=None
        mu_start=float(baseline.mu_gev)
    elif backend=='quadratic_only':
        savedfile=Path(qeft_dir)/f'{tag}_eft1_quadratic_running.json'
        sources['saved_candidate_quadratic_run']=savedfile
        saved=read(savedfile)
        need(saved.get('format')==qeft.FORMAT and
             saved.get('status')=='CandidateQuadraticAndFullQuarticEFT1Integrated' and
             saved.get('model')==model and saved.get('quadratic_kind')=='mixed' and
             saved.get('physical_matchete_rgbeta_mapping_certified') is False and
             saved.get('scalar_mass_threshold_matching_performed') is False and
             saved.get('wilson_coefficient_matching_performed') is False,
             'saved quadratic EFT1 evidence incompatible or overclaimed')
        baseline=s34.decode_intermediate_state(saved['initial_baseline'])
        state=GeneralState(baseline,{k:qeft.decode_complex(z,k) for k,z in saved['initial_new_quartics'].items()}).validated(registry)
        ybase,layout=pack(state,registry)
        mu2=checked_complex(saved['mu2_initial_gev2'],'mu2_initial_gev2')
        y0=np.r_[ybase,mu2.real,mu2.imag]
        report=qeft.checked_payload(model,baseline.representation,'mixed',layout,registry,beta_q)
        need(set(saved.get('beta_keys',[]))==set(report),'saved beta inventory differs from current native RGBeta')
        verify_saved_beta_sources(saved,{'native_beta':beta_q, 'registry':Path(qeft.REGISTRY),
                                         'manifest':Path(qeft.MANIFEST)})
        rhs=lambda t,y:qeft.derivative(np.asarray(y),scale_gev=math.exp(t),layout=layout,kind='mixed',registry=registry,report=report)
        mu_start=float(baseline.mu_gev)
    else:
        raise JointScaleBlocked('unknown backend')
    need(abs(mu_start-1e5)/1e5<1e-12, 'fermion threshold must be 1e5 GeV for this Stage33 candidate')
    return locals()


def run_model(root:Path, model:str, *,backend:str='auto',scheme:str='geometric',
              qeft_dir:Path|None=None,mass_dir:Path|None=None,stage33_dir:Path|None=None,
              higgs_mass_squared:float=-7800.,mu2_fraction:float=.01,mu2_phase:float=.6,
              rtol:float=1e-8,atol:float=1e-11)->dict:
    root=Path(root).resolve()
    need(model in MODELS,'only B_-1 and C_-1 supported')
    need(backend in ('auto','five_mass','quadratic_only'),'unknown RGE backend')
    need(math.isfinite(higgs_mass_squared) and math.isfinite(mu2_fraction)
         and 0<=mu2_fraction<=.45 and math.isfinite(mu2_phase),'invalid mass benchmark')
    need(0<rtol<=1e-6 and 0<atol<=1e-7,'unsupported integration tolerances')
    qeft_dir=Path(qeft_dir or root/'output/general_t3_quadratic_eft1')
    mass_dir=Path(mass_dir or root/'output/audits/workstream_b_bc_full_mass')
    # The Stage33 DATA directory is defined by the installed repository,
    # rather than guessed from old and potentially stale report filenames.
    if stage33_dir is None:
        from scripts import IntermediateEFTRunning as s34
        stage33_dir=Path(s34.DATA)
    else:stage33_dir=Path(stage33_dir)
    tag=model.replace('-','m').replace('+','p')
    mass_path=mass_dir/f'{tag}_eft1_full_mass_v1.json'
    qrun_path=qeft_dir/f'{tag}_eft1_quadratic_running.json'
    if backend=='auto':
        backend='five_mass' if mass_path.is_file() else 'quadratic_only'
    need(backend!='five_mass' or mass_path.is_file(),
         f'missing native full-mass beta {mass_path}; run existing Workstream B export or use --backend quadratic_only')
    need(backend!='quadratic_only' or qrun_path.is_file(),
         f'missing saved candidate quadratic EFT1 {qrun_path}; generate via scripts.RunGeneralT3QuadraticEFT1')
    ctx=_stage33_or_saved_qeft(model,backend,root,stage33_dir=stage33_dir,
                               qeft_dir=qeft_dir,mass_dir=mass_dir,
                               higgs_mass_squared=higgs_mass_squared,
                               mu2_fraction=mu2_fraction,mu2_phase=mu2_phase)
    rhs=ctx['rhs'];y0=ctx['y0'];layout=ctx['layout'];registry=ctx['registry']
    unpack=ctx['unpack'];Basis=ctx['Basis'];gauge_check=ctx['gauge_check'];hessian=ctx['hessian']
    s34=ctx['s34'];manifest=ctx['manifest'];sources=ctx['sources'];qeft=ctx['qeft']
    stop=float(s34.STOP)
    mu0=ctx['mu_start']
    need(0<stop<mu0,'invalid candidate EFT1 integration scale interval')
    sol=solve_ivp(rhs,(math.log(mu0),math.log(stop)),y0,method='DOP853',
                  rtol=rtol,atol=atol,dense_output=True)
    need(sol.success and sol.sol is not None and np.isfinite(sol.y).all(),
         'candidate EFT1 RG integration failed: '+str(sol.message))
    def state_at(t:float):
        z=np.asarray(sol.sol(t),dtype=float)
        need(z.shape==y0.shape and np.isfinite(z).all(),'invalid dense EFT1 state')
        st=unpack(z[:layout.size],layout,registry,uv=False,mu_gev=math.exp(t))
        mix=complex(z[layout.size],z[layout.size+1])
        basis=Basis.build(model,float(st.baseline.mS1Sq),float(st.baseline.mS2Sq),mix)
        need(not basis.tachyonic,'tachyon at sampled EFT1 scale')
        return z,st,basis
    # Existing candidate endpoint is verified only in the backend that shares
    # its four-mass-coordinate beta system. The 5-mass benchmark is independent.
    qsaved=ctx['saved']
    endpoint_residual=None
    if qsaved is not None:
        z,st,be=state_at(math.log(stop))
        endpoint_residual=compare_endpoint(np.array([st.baseline.mS1Sq,st.baseline.mS2Sq,
                              be.mu_sq.real,be.mu_sq.imag]),
                       np.array([qsaved['final_baseline']['mS1Sq'],
                                 qsaved['final_baseline']['mS2Sq'],
                                 checked_complex(qsaved['mu2_final_gev2'],'final_mu').real,
                                 checked_complex(qsaved['mu2_final_gev2'],'final_mu').imag]))
        for name,zval in st.extra.items():
            need(name in qsaved['final_new_quartics'],'missing saved new quartic: '+name)
            ref=checked_complex(qsaved['final_new_quartics'][name],name)
            need(abs(zval-ref)<2e-6*max(1.,abs(ref)),
                 'saved new-quartic endpoint mismatch: '+name)
    joint=joint_scale_from_running(math.log(stop),math.log(mu0),
           lambda t:state_at(t)[2].eigenvalues_sq,scheme=scheme)
    other_scheme='upper' if scheme=='geometric' else 'geometric'
    other_joint=joint_scale_from_running(math.log(stop),math.log(mu0),
           lambda t:state_at(t)[2].eigenvalues_sq,scheme=other_scheme)
    need(other_joint['matching_scale_gev'] >= joint['matching_scale_gev'] if scheme=='geometric'
         else other_joint['matching_scale_gev'] <= joint['matching_scale_gev'],
         'common-scale scheme ordering inconsistent')
    mu=joint['matching_scale_gev']
    _,st,basis=state_at(math.log(mu))
    need(mu<mu0 and mu>stop, 'common scale outside usable interval')
    gauge_error=float(gauge_check(basis,trials=20,seed=2812))
    independent=hessian(manifest,model,basis.m_s1_sq,basis.m_s2_sq,basis.mu_sq)
    reference=np.repeat(basis.eigenvalues_sq,4)
    eig=np.asarray(independent['mass_squared_eigenvalues_gev2'],dtype=float)
    need(eig.shape==(8,) and np.allclose(eig,reference,rtol=3e-11,atol=4e-6),
         'independent Cartesian-real scalar Hessian disagrees with complex mass basis')
    U=np.asarray(basis.U,dtype=complex)
    projectors=[np.outer(U[:,i],U[:,i].conj()) for i in range(2)]
    need(np.allclose(sum(projectors),np.eye(2),atol=2e-12) and
         np.allclose(projectors[0]@projectors[1],np.zeros((2,2)),atol=2e-12),
         'joint scalar eigenprojectors incomplete')
    from scripts import FermionThresholdHandoff as s33
    quartics={k:s33.json_state(z) for k,z in st.extra.items()}
    result={'format':FORMAT,'status':STATUS,'model':model,'backend':backend,
       'backend_full_higgs_mass_beta':backend=='five_mass',
       'source_mass_parameters_at_start':{'mu2_fraction':mu2_fraction,'mu2_phase':mu2_phase,
                                          'mHSq_gev2':higgs_mass_squared} if backend=='five_mass'
                      else {'taken_from_saved_candidate_quadratic_eft1':True},
       'different_benchmark_from_v21_v27_possible':True,
       'EFT_order':ORDER,'mu_fermion_first_gev':mu0,
       'mu_saved_candidate_endpoint_gev':stop,
       'common_scalar_matching_scale':joint,
       'other_joint_scale_scheme':other_joint,
       'state_at_common_scale':{'mu_gev':mu,'base':s33.json_state(st.baseline),
                    'extra_quartics':quartics,
                    'muT3GeneralS1S2_gev2':s33.json_state(basis.mu_sq),
                    'mHSq_gev2':float(sol.sol(math.log(mu))[layout.size+2]) if backend=='five_mass' else None,
                    'mass_basis':basis.info()},
       'native_scalar_component_tensor_transport_repeated_here':False,
       'independent_real_hessian_spectrum_verified':True,
       'epsilon_SU2_covariance_max_absolute_error':gauge_error,
       'scalar_eigenprojectors_complete':True,
       'saved_quadratic_endpoint_relative_residual':endpoint_residual,
       'one_loop_beta_convention':'16*pi^2*dX/dln(mu); gauge BetaTerm divided by 2*g',
       'rge_solver':'scipy solve_ivp DOP853 dense_output',
       'rge_nfev':int(sol.nfev),'rge_rtol':rtol,'rge_atol':atol,
       'source_sha256':{k:sha(v) for k,v in sources.items()},
       'source_files':{k:Path(v).resolve().relative_to(root).as_posix() for k,v in sources.items()},
       'certification':{'all_two_scalar_eigenstates_retained_to_common_matching_scale':True,
         'physical_Matchete_RGBeta_coupling_dictionary':False,
         'fermion_threshold_Wilson_boundary_computed':False,
         'signed_GammaCC_to_Weyl':False,'full_one_loop_dim5_ADM':False,
         'joint_two_scalar_threshold_Wilson_matching_performed':False,
         'physical_SMEFT_C5_generated':False,'physical_Wilson_evolution_authorized':False}}
    return result
