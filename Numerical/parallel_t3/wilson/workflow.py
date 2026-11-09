"""Full 45-Weyl, (4+2dS1+2dS2)-real-component coupled Wilson evolution.

No separate C5 ODE: C5 is an injection/projection of the SAME C tensor.
All fermions F were integrated out before the supplied EFT1 background.
"""
from __future__ import annotations
from pathlib import Path
from dataclasses import dataclass
import hashlib,json,math
import numpy as np
from scipy.integrate import solve_ivp
from .basis import ModelBasis,build_basis
from .coordinates import embed_c5,project_c5,symmetrize,symmetry_residual
from .adm import beta,total_beta,gauge_ward_residual,VARIANTS
from .interchange import Background,BlockedInput,sha256,load_background,load_boundary,FORMAT,BASIS_ID

LOOP=1/(16*math.pi**2)


def _save(path:Path,**payload)->None:
    if path.exists():raise BlockedInput('Blocked_OutputExists: '+str(path))
    path.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path,**payload)


def _summary(x:np.ndarray)->float:return float(np.linalg.norm(x))


def run_wilson(basis:ModelBasis,initial:np.ndarray,background:Background,
               mf:float,ms:float,*,adm_variant:str='fonseca_eq22_six_term',
               solver:str='RK2',steps:int=4,rtol:float=1e-6,atol:float=1e-21):
    """Return endpoints and diagnostics, numerically operational but not physical.

    'RK2': midpoint with steps in ln(mu), cheap qualitative survey.
    'DOP853': adaptive SciPy solve_ivp, recommended once B data are available.
    Both run ALL independent full-tensor coordinates at once, including
    generated LLHH, reverse mixing, self-running and non-LL components.
    """
    shape=(basis.nf,basis.nf,basis.ns,basis.ns)
    c=np.asarray(initial,np.complex128).copy()
    if c.shape!=shape or not np.isfinite(c).all():raise BlockedInput('Blocked_InitialTensor')
    if not 0<ms<mf:raise BlockedInput('Blocked_ThresholdOrderingMFgtMS')
    if adm_variant not in VARIANTS:raise BlockedInput('Blocked_ADMVariantNotSupported')
    if background.basis!=basis:raise BlockedInput('Blocked_BackgroundBasisMismatch')
    if not (steps>0 and 0<rtol<1 and 0<atol<1):raise ValueError('Invalid ODE controls')
    if abs(background.t[0]-math.log(ms))>1e-8 or abs(background.t[-1]-math.log(mf))>1e-8:
        raise BlockedInput('Blocked_InterpolationEndpoints')
    if symmetry_residual(c)>max(1e-24,1e-10*np.max(np.abs(c))):
        raise BlockedInput('Blocked_InitialPairSymmetry')
    t0,t1=math.log(mf),math.log(ms)
    counts={'calls':0,'backends':set()}
    def deriv(t,x):
        y,lam,gauge=background.at(t)
        channels,tag=beta(x,y,lam,gauge,variant=adm_variant)
        counts['calls']+=1;counts['backends'].add(tag)
        out=LOOP*total_beta(channels)
        if not np.isfinite(out).all():raise FloatingPointError('Nonfinite Wilson derivative')
        return out
    times=np.linspace(t0,t1,steps+1)
    projection=[]
    def record(t,x):
        c5,res=project_c5(x,basis)
        projection.append({'mu_gev':float(math.exp(t)),
                           'C5_frobenius_gev_inverse':_summary(c5),
                           'C5_projection_absolute_residual':res,
                           'C_full_frobenius_gev_inverse':_summary(x)})
    record(t0,c)
    if solver=='RK2':
        for lo,hi in zip(times[:-1],times[1:]):
            h=hi-lo;k1=deriv(lo,c);k2=deriv((lo+hi)/2,c+h*k1/2)
            c=symmetrize(c+h*k2)
            record(hi,c)
        nfev=counts['calls']
    elif solver=='DOP853':
        size=c.size
        def func(t,v):return deriv(t,v.reshape(shape)).reshape(-1)
        sol=solve_ivp(func,(t0,t1),c.reshape(-1),method='DOP853',t_eval=times,
                      rtol=rtol,atol=atol, max_step=max(abs(t1-t0)/steps,1e-6))
        if not sol.success:raise RuntimeError('Wilson integration failed: '+sol.message)
        for ix in range(1,len(sol.t)):
            record(float(sol.t[ix]),sol.y[:,ix].reshape(shape))
        c=symmetrize(sol.y[:,-1].reshape(shape));nfev=int(sol.nfev)
    else:raise ValueError('Unknown ODE solver')
    c5_initial,res_initial=project_c5(initial,basis)
    c5_final,res_final=project_c5(c,basis)
    b0,_=beta(initial,*background.at(t0),variant=adm_variant)
    b_end,_=beta(c,*background.at(t1),variant=adm_variant)
    beta_start=total_beta(b0)
    est=c5_initial+LOOP*(t1-t0)*project_c5(beta_start,basis)[0]
    return {
        'C_initial':np.asarray(initial,np.complex128), 'C_final':c,
        'C5_initial':c5_initial,'C5_final':c5_final,
        'C_LLSS_at_MS':c[np.ix_([f for p in basis.leptons for f in p],
                                  [f for p in basis.leptons for f in p],
                                  range(4,basis.ns),range(4,basis.ns))],
        'beta_at_MF':beta_start,'beta_at_MS':total_beta(b_end),
        'beta_channels_at_MF':b0,'beta_channels_at_MS':b_end,
        'leading_log_C5_linearized':est,
        'metadata':{
           'status':'Candidate_Qualitative', 'physical_certified':False,
           'adm_variant':adm_variant,
           'source_v1_unstarred_known_covariance_failure':adm_variant=='source_v1_unstarred',
           'loop_order':'conditional_one_loop_RGE_improved_not_physical_matching',
           'MF_gev':float(mf),'MS_gev':float(ms),
           'source_ODE':'16pi2*dC/dln(mu)=beta(C,C_conjugate)',
           'solver':solver,'steps':steps,'rhs_evaluations':nfev,
           'execution_backends':sorted(counts['backends']),
           'C5_at_MF_norm':_summary(c5_initial),'C5_at_MS_norm':_summary(c5_final),
           'C5_boundary_projection_residual':res_initial,
           'C5_final_projection_residual':res_final,
           'fermion_pair_symmetry_residual':symmetry_residual(c),
           'linearized_leading_log_C5_relative_difference':float(
              np.linalg.norm(est-c5_final)/max(np.linalg.norm(c5_final),1e-300)),
           'C5_mutual_mixing_present':True,
           'no_posthoc_fixed_order_leading_log':True,
           'full_psi2phi2_operators_active':True,
           'all_F_integrated_out_at_MF':True,
           'both_S1_S2_remain_until_MS':True,
           'other_dim5_classes':'symmetry_forbidden_in_T3_Z2_EFT1_under_3family_SM_Weyl_assumption',
           'other_dim5_exclusion_basis':{
               'phi5':'Z2 even forces even heavy fields; odd H doublet count, non-singlet SU2 center',
               'phi_F2':'no Z2-even gauge-singlet scalar (H is a doublet)',
               'psi2_F':'no two-SM-left-Weyl gauge-singlet bilinear without H'
           },
           'other_operator_closure_automated_native_certified':False,
           'warnings': ['six-vs-seven formula and gauge completion not certified',
                        'GammaCC/Weyl signed map not certified',
                        'complete Matchete/RGBeta scalar tensor map not certified',
                        'MS simultaneous scalar threshold belongs to C',
                        'one-loop hard fermion matching belongs to B'],
           'projections':projection}}


def write_results(result:dict,basis:ModelBasis,output:Path,provenance:dict):
    if output.exists() and any(output.iterdir()):
        raise BlockedInput('Blocked_OutputDirectoryNonempty: '+str(output))
    output.mkdir(parents=True,exist_ok=True)
    ini=output/'canonical_wilson_MF.npz'; end=output/'wilson_at_MS.npz'
    _save(ini,C_full=result['C_initial'],C5_EFT1_MF=result['C5_initial'],
          beta_C_at_MF_16pi2=result['beta_at_MF'])
    _save(end,C_full=result['C_final'],C5_EFT1_MS=result['C5_final'],
          C_LLSS_at_MS=result['C_LLSS_at_MS'],
          C5_EFT1_MF=result['C5_initial'],
          C5_linearized_leading_log_diagnostic=result['leading_log_C5_linearized'],
          beta_C_at_MS_16pi2=result['beta_at_MS'])
    meta={'format':'t3-parallel-a-wilson-output-v1',**result['metadata'],
          'basis':basis.to_json(), 'provenance':provenance,
          'outputs':{'canonical_wilson_MF.npz':sha256(ini),
                     'wilson_at_MS.npz':sha256(end)},
          'C5_projection_definition':'C5[p,q]=4*C[L[p,1],L[q,1],H[2,R],H[2,R]]',
          'C5_not_smeft_after_scalar_matching':True,
          'blocked_physics':['Blocked_ADM_PhysicalNormalization',
                             'Blocked_GammaCCWeylSignedNormalization',
                             'Blocked_MatcheteRGBetaCompleteQuarticMap',
                             'Blocked_NativeOneLoopFHardMatching',
                             'Blocked_FullMassAndScalarEigenbasisForSpecialModels']}
    (output/'manifest.json').write_text(json.dumps(meta,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    return meta


def run_from_files(root:Path,model:str,background_npz:Path,background_json:Path,
                   boundary_npz:Path,boundary_json:Path,output:Path,*,
                   adm_variant='fonseca_eq22_six_term',solver='RK2',steps=4,
                   assume_formal_identity=False):
    root=Path(root);basis=build_basis(model)
    back=load_background(root,basis,Path(background_npz),Path(background_json))
    ini,src=load_boundary(root,basis,Path(boundary_npz),Path(boundary_json),
                          assume_formal_identity=assume_formal_identity)
    bm=json.loads(Path(background_json).read_text(encoding='utf-8-sig'))
    sm=json.loads(Path(boundary_json).read_text(encoding='utf-8-sig'))
    if abs(float(bm['MF_gev'])-float(sm['MF_gev']))>1e-8:
        raise BlockedInput('Blocked_BoundaryMFBackgroundMFDisagree')
    if sm.get('model')!=bm.get('model'):raise BlockedInput('Blocked_MultiModelInputs')
    out=run_wilson(basis,ini,back,float(bm['MF_gev']),float(bm['MS_gev']),
                   adm_variant=adm_variant,solver=solver,steps=steps)
    return write_results(out,basis,Path(output),{'input_boundary':src, 'input_background':back.provenance,
        'assume_formal_identity_used':assume_formal_identity,
        'source_repository':'Frankyang2002/T3neutrinomassmodelpipeline',
        'verified_public_main':'a944e0030b8ef24b23342d99e0d8f33cef3cffd5'})


def synthetic_case(model='A_-4',mf=1e5,ms=1e3,*, portal=0.02):
    """Self-contained FIxTURE, not a physical UV/RGBeta or Matchete boundary.

    All 16 carry an exact SU2xU1 gauge-invariant synthetic LLSS, a
    conjugate portal and an independently nonzero hard C5. These are
    algebraic fixtures, NOT source-normalized matched couplings.
    """
    basis=build_basis(model)
    c5=np.array([[1.5e-16+0.4e-16j,0.3e-16-0.2e-16j,0.1e-16],
                 [0.3e-16-0.2e-16j,1.7e-16+0.1e-16j,0.2e-16+0.1e-16j],
                 [0.1e-16,0.2e-16+0.1e-16j,1.3e-16-0.2e-16j]],complex)
    c=embed_c5(c5,basis)
    lam=np.zeros((basis.ns,)*4,float)
    from .synthetic_invariants import synthetic_llss,synthetic_portal
    x=np.array([[2e-10+1e-10j,1e-10,0],[1e-10,3e-10,0],[0,0,2e-10-1e-10j]],complex)
    c+=synthetic_llss(basis,x)
    lam=synthetic_portal(basis,portal)
    # In the synthetic case gauge strengths are zero by choice, to isolate the
    # explicit quartic mixing; it is NOT an RGBeta running background.
    t=np.log([ms,mf]);y=np.zeros((2,basis.nf,basis.nf,basis.ns),complex)
    L=np.stack([lam,lam]);g=np.zeros((2,3))
    back=Background(basis,t,y,L,g,basis.generators(),
                    {'status':'SYNTHETIC_NoUVBInput','synthetic':True})
    return basis,c,back,{'MF_gev':mf,'MS_gev':ms,'C5_hard_nonzero':True,
                          'S1_S2_both_active':True,'not_physical':True}
