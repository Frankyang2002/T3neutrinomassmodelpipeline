"""Experimental read-only B(-1)/C(-1) EFT1(H,X,Y) -> EFT2(H,X) tree projection.

Reintegrates the SAME continuous_v1 scalar background as the existing coupled
Wilson calculation, retains the complex mixed quadratic mu2(t), finds the
running Y crossing on THAT trajectory, and integrates the original 12-real-
scalar coupled Wilson system only down to muY. No EFT2 RGE, no finite scalar
threshold correction, no physical neutrino spectrum. Never changes source data.

Run from repository root, after the prerequisite direct/selected Wilson reports:
 python -u -X utf8 -m scripts.RunT3BCYThresholdProjection --all --run-id continuous_v1
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from hashlib import sha256

import numpy as np
from scipy.optimize import brentq

from scripts import RunT3BCCoupledWilsonSequential as coupled
from Numerical.experimental import T3BCCoupledWilsonTensor as tensor
from Numerical.experimental import T3BCScalarYThresholdProjection as projection
from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis

continuous=coupled.continuous
neutrino=coupled.neutrino
direct=coupled.direct
selected=coupled.selected
ROOT=coupled.ROOT
MODELS=('B_-1','C_-1')
FORMAT='t3-bc-y-threshold-projected-wilson-conditional-v1'
STATUS='EFT1YThresholdLocatedEFT2TreeProjectedFiniteMatchingMissing'
DEFAULT_OUTPUT_ID='bc_y_threshold_projection_v1'


def require(ok:bool,message:str)->None:
    if not ok:raise ValueError('B/C Y-threshold projection: '+message)


def scalar_background_with_mu2(model:str, record:dict, sidecar:dict, *,samples:int):
    """Reuses source continuous coupled RGBeta solver; retain its dropped mu2(t)."""
    from Numerical.experimental.GeneralT3ScalarRGState import (
        GeneralT3ScalarRGState,pack_general_scalar_state,unpack_general_scalar_state)
    from scripts import IntermediateEFTRunning as s34
    require(model in MODELS and 17<=samples<=129 and samples%2==1,
            'samples must be odd 17..129')
    require(sidecar.get('status')=='CandidateQuadraticAndFullQuarticEFT1Integrated'
            and sidecar.get('format')==continuous.qeft.FORMAT
            and sidecar.get('model')==model
            and record.get('quadratic_kind')=='mixed',
            'saved full-scalar EFT1 source incompatible')
    registry=continuous.load_registry(continuous.REGISTRY)
    base=s34.decode_intermediate_state(sidecar['initial_baseline'])
    extra={name:neutrino.complex_pair(raw,name)
           for name,raw in sidecar['initial_new_quartics'].items()}
    initial=GeneralT3ScalarRGState(base,extra).validated(registry)
    vector,layout=pack_general_scalar_state(initial,registry)
    init_mu2=neutrino.complex_pair(sidecar['mu2_initial_gev2'],'initial mu2')
    vector=np.concatenate((vector,[init_mu2.real,init_mu2.imag]))
    beta_path=coupled.sources(model)['eft1_beta']
    betas=continuous.beta_inventory(continuous.read(beta_path),model=model,
        stage='EFT1',representation=base.representation,layout=layout,
        registry=registry,kind='mixed')
    result=continuous.integrate_segment(vector,mu_start=s34.START,mu_end=s34.STOP,
        layout=layout,registry=registry,stage='EFT1',kind='mixed',
        betas=betas,points=samples)
    mu=np.exp(np.asarray(result.t,dtype=float))
    states=tuple(unpack_general_scalar_state(result.y[:layout.size,i],layout,
        registry,uv=False,mu_gev=float(mu[i])) for i in range(samples))
    mu2=np.asarray(result.y[layout.size]+1j*result.y[layout.size+1],dtype=complex)
    # Reproduce the existing coupled-runner provenance/endpoint validation.
    previous=neutrino.decode_eft(record['eft1_scalar_endpoint_state'])
    errors={}
    for name in ('lambdaT3','mS1Sq','mS2Sq',*selected.ACTIVE_QUARTICS):
        now=complex(getattr(states[-1].baseline,name))
        reference=complex(getattr(previous,name))
        err=abs(now-reference)/max(abs(reference),1e-12)
        require(err<2e-6,'continuous endpoint mismatch: '+name)
        errors[name]=float(err)
    require(set(states[-1].extra)==set(sidecar['final_new_quartics']),
            'new quartic inventory mismatch')
    for name,value in states[-1].extra.items():
        reference=neutrino.complex_pair(sidecar['final_new_quartics'][name],name)
        require(abs(value-reference)<2e-6*max(abs(reference),1e-10),
                'new quartic endpoint mismatch: '+name)
    reference_mu2=neutrino.complex_pair(sidecar['mu2_final_gev2'],'final mu2')
    err=abs(mu2[-1]-reference_mu2)/max(abs(reference_mu2),1e-12)
    require(err<2e-6,'mu2 endpoint mismatch')
    errors['mixed_mu2']=float(err)
    require(np.all(np.isfinite(mu)) and np.all(np.diff(mu)<0)
            and np.isfinite(mu2).all(), 'invalid quadratic trajectory')
    return states,mu,mu2,{'endpoint_parity':errors,'eft1_nfev':int(result.nfev)}


def locate_y_crossing(model:str,states:tuple,mu:np.ndarray,mu2:np.ndarray):
    """Find upper eigenvalue of A(t) crossing exp(2t) on ONE EFT1 trajectory.

    Interpolate quadratic entries linearly in ln(mu) *between* RGBeta output
    nodes; do not splice in standalone quadratic-study threshold numbers.
    """
    log_mu=np.log(mu)
    m1=np.array([float(np.real(st.baseline.mS1Sq)) for st in states])
    m2=np.array([float(np.real(st.baseline.mS2Sq)) for st in states])
    require(np.max(np.abs([complex(st.baseline.mS1Sq).imag for st in states]))<1e-6
            and np.max(np.abs([complex(st.baseline.mS2Sq).imag for st in states]))<1e-6,
            'non-Hermitian diagonal quadratic masses')
    def values(t:float,i:int)->tuple[float,float,complex]:
        w=(t-log_mu[i])/(log_mu[i+1]-log_mu[i])
        return (float((1-w)*m1[i]+w*m1[i+1]),
                float((1-w)*m2[i]+w*m2[i+1]),
                complex((1-w)*mu2[i]+w*mu2[i+1]))
    def crossing(t:float,i:int)->float:
        a,b,z=values(t,i)
        upper=float(np.linalg.eigvalsh(np.array([[a,z],[z.conjugate(),b]]) )[-1])
        return upper - math.exp(2*t)
    roots=[]
    for i in range(len(log_mu)-1):
        f0=crossing(float(log_mu[i]),i)
        f1=crossing(float(log_mu[i+1]),i)
        if f0<0 and f1>0:
            root=brentq(lambda t:crossing(t,i),float(log_mu[i+1]),
                        float(log_mu[i]),xtol=2e-13)
            roots.append((i,root))
    require(len(roots)==1,'expected exactly one descending EFT1 Y mass crossing, got '+str(len(roots)))
    index,t_y=roots[0]
    a,b,z=values(t_y,index)
    basis=MixedDoubletMassBasis.build(model,a,b,z)
    require(not basis.tachyonic and not basis.degenerate,
            'tachyonic/degenerate Y-threshold mass basis')
    mu_y=math.exp(t_y)
    require(abs(basis.eigenvalues_sq[1]-mu_y**2) < 1e-8*mu_y**2,
            'running Y crossing not self-consistent')
    require(mu[index]>mu_y>mu[index+1],'Y scale outside bracket')
    return basis,mu_y,index,{'interpolation':'linear quadratic entries in ln(mu)',
        'bracket_gev':[float(mu[index]),float(mu[index+1])],
        'high_eigenmass_gev':float(np.sqrt(basis.eigenvalues_sq[1])),
        'low_eigenmass_at_muY_gev':float(np.sqrt(basis.eigenvalues_sq[0])),
        'crossing_relative_squared_residual':float(abs(basis.eigenvalues_sq[1]/mu_y**2-1))}


def source_paths(model:str,folder:Path,record:dict,source_id:str,extra_id:str)->dict[str,Path]:
    tag=continuous.model_tag(model)
    d=coupled.sources(model)
    d.update({'continuous':folder/tag/'continuous_trajectory.json',
        'sidecar':folder/tag/record['saved_candidate_eft1_endpoint'],
        'prior_direct':folder/source_id/f'{tag}_direct_wilson.json',
        'prior_extra':folder/extra_id/f'{tag}_extra_wilson.json',
        'runner':Path(__file__).resolve(),
        'projection_module':Path(projection.__file__).resolve(),
        'mass_basis_module':ROOT/'Numerical/experimental/GeneralT3MixedScalarMassBasis.py',
        'coupled_runner':Path(coupled.__file__).resolve(),
        'coupled_tensor_module':Path(tensor.__file__).resolve(),
        'direct_runner':Path(direct.__file__).resolve(),
        'selected_runner':Path(selected.__file__).resolve(),
        'neutrino_runner':Path(neutrino.__file__).resolve(),
        'continuous_runner':Path(continuous.__file__).resolve(),
        'registry':continuous.REGISTRY})
    return d


def calculate(model:str,folder:Path,out:Path,*,samples:int=65,
              rtol:float=5e-7,direct_output_id:str=direct.DEFAULT_OUTPUT_ID,
              selected_output_id:str=selected.DEFAULT_OUTPUT_ID)->dict:
    require(model in MODELS,'only B_-1/C_-1')
    tag=continuous.model_tag(model)
    rec_path=folder/tag/'continuous_trajectory.json'
    record=neutrino.verify_continuous(model,rec_path)
    sidecar=neutrino.load(folder/tag/record['saved_candidate_eft1_endpoint'])
    paths=source_paths(model,folder,record,direct_output_id,selected_output_id)
    prior=neutrino.load(paths['prior_direct'])
    prior_extra=neutrino.load(paths['prior_extra'])
    require(prior.get('format')==direct.FORMAT and prior.get('status')==direct.STATUS
        and prior.get('model')==model and prior.get('saved_points')==samples,
        'missing/incompatible direct-source parity report')
    require(prior_extra.get('format')==selected.FORMAT and prior_extra.get('status')==selected.STATUS
        and prior_extra.get('model')==model and prior_extra.get('eft1_source_samples')==samples,
        'missing/incompatible selected-source parity report')
    hashes=coupled.fingerprints(paths)
    json_path=out/f'{tag}_y_threshold_wilson.json'
    npz_path=out/f'{tag}_y_threshold_wilson.npz'
    require(not json_path.exists() and not npz_path.exists(),
            'output exists; choose another --output-id (source outputs are immutable)')

    stage17=neutrino.load(paths['stage17'])
    stage19=neutrino.load(paths['stage19'])
    seed=neutrino.load(paths['seed'])
    wilson=neutrino.load(paths['wilson_rge'])
    matched=neutrino.load(paths['matching'])
    selected.certify_saved_selected_source(model,stage17,stage19,seed,wilson,matched)
    uv=neutrino.decode_uv(record['uv_fermion_threshold_state'])
    mf,y1,y2,_heavyres=direct.heavy_inputs(uv)
    source=direct.compile_direct_source(paths['matching'],paths['wilson_rge'])

    states,mu,mu2,background=scalar_background_with_mu2(model,record,sidecar,samples=samples)
    basis,mu_y,idx,crossing=locate_y_crossing(model,states,mu,mu2)
    from RGE.running.intermediate.ScalarOnlyWilsonTensorRGE import build_scalar_only_rge_context
    context=build_scalar_only_rge_context(paths['seed'],paths['quartic_seed'],paths['old_eft1_rgbeta'])
    static=tensor.build_static_kinetics(context)
    require(static['ns']==12,'unexpected real scalar component count')
    mat22=neutrino.load(paths['stage22'])
    map29=neutrino.load(paths['stage29'])
    quartet=tensor.QuarticEvaluator.build(context.quartic,mat22,map29,model)
    template=coupled.template_at_context(context)
    ids=static['indices']
    flavors=tensor.flavor_pairs(mf,y1,y2)
    tree=tensor.independent_tensor(context.coefficient,flavors)
    initial=tree.copy()
    hard,relocation,grouped=direct.grouped_hard_boundary(
        source,mf,y1,y2,complex(states[0].baseline.lambdaT3),mu_f=float(mu[0]))
    coupled.embed_weinberg(initial,template,ids,grouped)
    projected_boundary,boundary_res=tensor.weinberg_projection(initial,template,ids)
    require(boundary_res<1e-9 and np.allclose(
        tensor.flavor_to_matrix(projected_boundary),grouped,rtol=1e-10,atol=1e-20),
        'grouped hard Weinberg boundary embedding failed')
    direct_s=neutrino.complex_matrix(prior['c5_at_scalar_threshold_gev_inverse'],'prior direct C5')
    extra_s=neutrino.complex_matrix(prior_extra['selected_extra_delta_C5_at_scalar_threshold_gev_inverse'],
                                    'prior selected C5')
    parity=coupled.check_source_parity(model,tree,states,static,quartet,
        template,ids,source,mf,y1,y2,mu=mu,hard=grouped,
        original=direct_s,selected_original=extra_s)

    kernels=[tensor.kernels(st,static,quartet.evaluate(st)) for st in states]
    log_y=float(math.log(mu_y))
    y_kernel=tensor.interpolate_kernels(np.log(mu),kernels,log_y)
    mu_segment=np.r_[mu[:idx+1],mu_y]
    k_segment=kernels[:idx+1]+[y_kernel]
    require(len(mu_segment)>=9,'too few nodes above Y crossing for Wilson integration')
    trajectory,solution=tensor.integrate_coupled(initial,mu_segment,k_segment,rtol=rtol)
    c_original=np.asarray(trajectory[...,-1],dtype=complex)

    T=projection.mixing_real_matrix(basis.U)
    checks={'field_identity':projection.check_field_identity(basis.U,T),
            'gauge_covariance':projection.check_gauge_covariance(T)}
    c_mass=projection.rotate_wilson(c_original,T)
    c_eft2=projection.light_projection(c_mass)
    checks['full_polynomial_identity']=projection.check_tensor_identity(c_original,c_mass,T)
    checks['light_restriction_identity']=projection.check_light_identity(c_mass,c_eft2)
    checks['old_symmetry']=projection.symmetry_residual(c_original)
    checks['mass_symmetry']=projection.symmetry_residual(c_mass)
    checks['light_symmetry']=projection.symmetry_residual(c_eft2)
    require(max(checks['light_symmetry'].values())<1e-9,'EFT2 Wilson symmetries broken')
    # H is untouched by scalar diagonalization. Its existing Weinberg projection
    # is the diagnostic HARD+RUNNING C5 contribution, not final EFT2 C5.
    c5_old,old_res=tensor.weinberg_projection(c_original,template,ids)
    c5_mass,mass_res=tensor.weinberg_projection(c_mass,template,ids)
    require(np.allclose(c5_old,c5_mass,rtol=3e-10,atol=1e-24),
            'pure-HH Weinberg invariant changed under S1/S2 rotation')
    checks['weinberg_residual_before']=float(old_res)
    checks['weinberg_residual_after']=float(mass_res)
    checks['higgs_higgs_preserved_max_abs']=float(np.max(abs(c_original[:,:,:4,:4,:]
                                                          -c_eft2[:,:,:4,:4,:])))
    require(checks['higgs_higgs_preserved_max_abs']<1e-19,
            'Higgs-only block altered by scalar field rotation')
    survivors=projection.classify_survivors(c_eft2)
    require(any(v['max_abs']>0 for v in survivors.values()),'empty projected Wilson output')

    output={
        'format':FORMAT,'status':STATUS,'model':model,
        'input_sha256':hashes,'eft1_samples':samples,'wilson_rtol':float(rtol),
        'source_trajectory_id':folder.name,
        'source_scalar_trajectory':'continuous full-scalar EFT1; NOT standalone quadratic threshold study',
        'Y_threshold_gev':float(mu_y),'Y_mass_basis':basis.info(),
        'Y_crossing':crossing,'source_scalar_endpoint_checks':background,
        'source_normalization_checks':parity,'tensor_checks':checks,
        'eft1_wilson_initial_shape':list(initial.shape),
        'eft1_wilson_Y_shape':list(c_original.shape),
        'eft2_tree_wilson_shape':list(c_eft2.shape),
        'physical_lepton_pairs':[list(p) for p in tensor.PAIRS],
        'active_eft2_scalars':['H','X'],
        'real_component_order_eft1':['H_R1','H_I1','H_R2','H_I2',
            'S1_R1','S1_I1','S1_R2','S1_I2','S2_R1','S2_I1','S2_R2','S2_I2'],
        'real_component_order_mass':['H_R1','H_I1','H_R2','H_I2',
            'X_R1','X_I1','X_R2','X_I2','Y_R1','Y_I1','Y_R2','Y_I2'],
        'surviving_blocks':survivors,
        'grouped_fermion_hard_C5':neutrino.matrix_payload(grouped),
        'diagnostic_HH_C5_at_Y':neutrino.matrix_payload(
            tensor.flavor_to_matrix(c5_mass)),
        'wilson_ode_nfev':int(solution.nfev),
        'wilson_ode_status':solution.message,
        'tensor_archive_file':npz_path.name,
        'tensor_archive_convention':'real_old=T_real@real_new; C_mass=T_real^T C_original T_real on scalar legs',
        'scalar_threshold_tree_projection_computed':True,
        'finite_Y_threshold_delta_C_computed':False,
        'finite_Y_threshold_delta_C_status':'MISSING_NOT_SET_TO_PHYSICAL_ZERO',
        'EFT2_wilson_running_performed':False,
        'X_threshold_located_in_EFT2':False,
        'finite_X_threshold_matching_performed':False,
        'complete_neutrino_mass_spectrum_predicted':False,
        'full_three_generation_SM_Yukawa_wilson_anomalous_dimension_included':False,
        'additional_quartic_mapping_physically_certified':False,
        'full_general_scalar_matching_certified':False,
        'full_sequential_matching_performed':False,
        'notes':('This output is a computed EFT1 Wilson evolution and an exact field-basis '
            'tree projection at the Y crossing on the SAME continuous_v1 candidate '
            'trajectory. The finite scalar threshold correction is unknown, not zero. '
            'The EFT2 beta functions/operator closure and physical C5 matching remain undone.')}
    out.mkdir(parents=True,exist_ok=True)
    # Atomic archive write before publishing metadata; refuse prior output above.
    tmp=npz_path.with_suffix('.npz.tmp')
    with tmp.open('wb') as handle:
        np.savez_compressed(handle,
            C_EFT1_at_Y_original_real=c_original,
            C_EFT1_at_Y_XY_real=c_mass,
            C_EFT2_tree_at_Y_HX_real=c_eft2,
            real_old_from_new=T,U_complex=basis.U,
            Y_threshold_gev=np.array([mu_y]),
            Y_eigenmasses_sq_gev2=basis.eigenvalues_sq,
            eft1_mu_above_threshold=mu_segment)
    tmp.replace(npz_path)
    continuous.write_json_atomic(json_path,output)
    return output


def main()->int:
    ap=argparse.ArgumentParser(description=__doc__)
    choose=ap.add_mutually_exclusive_group(required=True)
    choose.add_argument('--all',action='store_true')
    choose.add_argument('--model',choices=MODELS)
    ap.add_argument('--run-id',required=True,type=continuous.run_id)
    ap.add_argument('--source-root',type=Path,default=continuous.OUTPUT)
    ap.add_argument('--output-id',default=DEFAULT_OUTPUT_ID)
    ap.add_argument('--direct-output-id',default=direct.DEFAULT_OUTPUT_ID)
    ap.add_argument('--selected-output-id',default=selected.DEFAULT_OUTPUT_ID)
    ap.add_argument('--samples',type=int,default=65)
    ap.add_argument('--rtol',type=float,default=5e-7)
    args=ap.parse_args()
    for name in (args.output_id,args.direct_output_id,args.selected_output_id):
        continuous.run_id(name)
    require(len({args.output_id,args.direct_output_id,args.selected_output_id})==3,
            'output IDs must differ')
    require(17<=args.samples<=129 and args.samples%2==1,'samples must be odd 17..129')
    require(1e-10<args.rtol<1e-3 and math.isfinite(args.rtol),'invalid Wilson rtol')
    folder=args.source_root/args.run_id
    require(folder.is_dir(),'missing continuous scalar study: '+str(folder))
    output=folder/args.output_id
    errors=[]
    for model in (MODELS if args.all else (args.model,)):
        try:
            result=calculate(model,folder,output,samples=args.samples,rtol=args.rtol,
               direct_output_id=args.direct_output_id,
               selected_output_id=args.selected_output_id)
            print(f'{model}: Y crossing={result["Y_threshold_gev"]:.8f} GeV; '
                  f'EFT2 tree tensor={result["eft2_tree_wilson_shape"]}; '
                  'finite matching=MISSING',flush=True)
        except (OSError,ValueError,RuntimeError,KeyError,TypeError,
                ArithmeticError,np.linalg.LinAlgError) as exc:
            errors.append(model)
            print(f'{model}: BLOCKED {type(exc).__name__}: {exc}',flush=True)
    return 2 if errors else 0


if __name__=='__main__':raise SystemExit(main())
