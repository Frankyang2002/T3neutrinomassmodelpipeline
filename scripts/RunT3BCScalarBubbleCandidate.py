"""Experimental B/C scalar-bubble finite Y matching sensitivity on completed EFT2.

Reads immutable Y-threshold and conditional EFT2 outputs, reconstructs their
verified full EFT1 quartic tensor, computes the renormalised scalar-pair
(LL phi phi)x(phi^4) one-loop bubble contribution at Y in two regimes, and
propagates each *candidate* Wilson boundary through the EXISTING EFT2 and
SM+Weinberg solvers. This is NOT the full finite, physical matching.

Run from project root:
  python -u -X utf8 -m scripts.RunT3BCScalarBubbleCandidate --all \
    --run-id continuous_v1 --output-id bc_scalar_bubble_candidate_v1

This stage does NOT modify existing outputs or any production source.
"""
from __future__ import annotations
import argparse
import csv
import math
from hashlib import sha256
from pathlib import Path
import numpy as np

from Numerical.experimental import T3BCScalarBubbleThreshold as bubble
from Numerical.experimental import T3BCEFT2HX as hx
from Numerical.experimental import T3BCCoupledWilsonTensor as tensor
from Numerical.experimental import T3BCScalarYThresholdProjection as projection
from scripts import RunT3BCEFT2Sequential as sequential
from scripts import RunT3BCYThresholdProjection as ystage

continuous=ystage.continuous
neutrino=ystage.neutrino
coupled=ystage.coupled
MODELS=ystage.MODELS
FORMAT='t3-bc-scalar-pair-finite-bubble-diagnostic-v1'
STATUS='ConditionalScalarBubbleFiniteThresholdSensitivityNotPhysicalMatching'
DEFAULT_OUTPUT_ID='bc_scalar_bubble_candidate_v1'


def require(condition,message):
    if not condition:raise ValueError('B/C bubble stage: '+message)


def digest(path):
    path=Path(path)
    require(path.is_file(),'missing file '+str(path))
    return sha256(path.read_bytes()).hexdigest()


def mass_response(C0,mu_y,mu_x,scalar_solution,static,template,ids,*,samples,rtol,vev=246.22,low_scale=100.):
    """Propagate corrected EFT2 boundary using EXISTING project routines."""
    from Numerical.running.SMWeinbergEvolution import SMWeinbergInitialConditions,evolve_sm_weinberg
    from Numerical.fitting.ParameterScan import (
        neutrino_mass_from_c5,charged_lepton_mass_basis_matrix)
    from Numerical.fitting.OscillationFit import prediction_from_mass_matrix
    cx,ode=hx.evolve_wilson(C0,mu_y,mu_x,scalar_solution,static,
                              samples=samples,rtol=rtol)
    c5,weinberg_res=tensor.weinberg_projection(cx,template,ids,threshold=6e-5)
    mat=neutrino.validated_c5(tensor.flavor_to_matrix(c5))
    scalar_x=hx.unpack_running(scalar_solution.sol(float(np.log(mu_x))))
    sm=SMWeinbergInitialConditions(
        gY=scalar_x['gY'],g2=scalar_x['g2'],g3=scalar_x['g3'],
        lambdaH=scalar_x['quartics'][0],
        ye=scalar_x['ye'],yu=scalar_x['yu'],yd=scalar_x['yd'],K=mat).validated()
    low=evolve_sm_weinberg(sm,mu_x,low_scale,rtol=1e-8,atol=1e-11)
    mass=neutrino_mass_from_c5(low.K,vev_gev=vev)
    flavour=charged_lepton_mass_basis_matrix(mass,low.ye)
    masses=np.linalg.svd(flavour,compute_uv=False)*1e9
    prediction=prediction_from_mass_matrix(flavour,ordering='NO')
    require(np.isfinite(masses).all(), 'nonfinite mass response')
    return dict(C5_at_X=mat,C5_at_100gev=low.K,
      masses_ev_descending=masses.tolist(),
      oscillation_normal_ordering=prediction.as_dict(),
      weinberg_projection_residual=float(weinberg_res),
      eft2_wilson_nfev=int(ode.nfev),sm_weinberg_nfev=int(low.nfev))


def calculate(model,folder,out,*,y_id,base_id,samples=65,
              eft2_samples=33,rtol=5e-7):
    tag=continuous.model_tag(model)
    yj=folder/y_id/(tag+'_y_threshold_wilson.json')
    yz=folder/y_id/(tag+'_y_threshold_wilson.npz')
    bj=folder/base_id/(tag+'_eft2_sequential.json')
    bz=folder/base_id/(tag+'_eft2_sequential.npz')
    output_json=out/(tag+'_scalar_bubble_candidate.json')
    output_npz=out/(tag+'_scalar_bubble_candidate.npz')
    require(not output_json.exists() and not output_npz.exists(),
            'outputs already exist, choose new --output-id')
    info,Cy,T,U,my,m2,record,source_hashes=sequential.read_input(
         model,folder,y_id,samples,rtol)
    original=neutrino.load(bj)
    require(original.get('format')==sequential.FORMAT and
            original.get('status')==sequential.STATUS and
            original.get('model')==model and
            original.get('finite_Y_threshold_delta_C_computed') is False and
            original.get('finite_X_threshold_delta_C5_computed') is False and
            original.get('physical_neutrino_mass_prediction_certified') is False,
            'unexpected baseline EFT2 provenance/certification')
    require(abs(original['mu_Y_gev']/my-1)<1e-10 and
            original['eft2_wilson_samples']==eft2_samples and
            original['wilson_rtol']==rtol,
            'baseline EFT2 settings inconsistent with saved Y threshold')
    with np.load(yz,allow_pickle=False) as z:
        Cfull=np.asarray(z['C_EFT1_at_Y_XY_real'],complex)
        Ctree=np.asarray(z['C_EFT2_tree_at_Y_HX_real'],complex)
    with np.load(bz,allow_pickle=False) as z:
        initial=np.asarray(z['scalar_Y'],float)
        c_prior=np.asarray(z['C_EFT2_at_Y'],complex)
    require(np.allclose(Ctree,Cy,atol=1e-25,rtol=1e-12) and
            np.allclose(Ctree,c_prior,atol=1e-25,rtol=1e-12),
            'Y/EFT2 Wilson boundary mismatch')
    require(np.array_equal(Cfull[:,:,:8,:8,:],Cy),
            'full Y tensor and projected tensor differ at tree boundary')
    sidecar=neutrino.load(folder/tag/record['saved_candidate_eft1_endpoint'])
    states,mu,mu2,_background=ystage.scalar_background_with_mu2(
        model,record,sidecar,samples=samples)
    basis,mu_check,idx,_cross=ystage.locate_y_crossing(model,states,mu,mu2)
    require(abs(mu_check/my-1)<2e-11 and np.linalg.norm(basis.U-U)<2e-10,
            'reconstructed Y crossing/rotation inconsistent with source')
    paths=coupled.sources(model)
    from RGE.running.intermediate.ScalarOnlyWilsonTensorRGE import build_scalar_only_rge_context
    context=build_scalar_only_rge_context(paths['seed'],paths['quartic_seed'],
                                          paths['old_eft1_rgbeta'])
    static_original=tensor.build_static_kinetics(context)
    static_hx=hx.rotate_static(static_original,T)
    quartic_evaluator=tensor.QuarticEvaluator.build(
        context.quartic,neutrino.load(paths['stage22']),
        neutrino.load(paths['stage29']),model)
    t=np.log(my)
    w=float((t-np.log(mu[idx]))/(np.log(mu[idx+1])-np.log(mu[idx])))
    require(0<w<1,'Y crossing outside interpolation nodes')
    qsource=(1-w)*quartic_evaluator.evaluate(states[idx]) +\
            w*quartic_evaluator.evaluate(states[idx+1])
    qmass=np.einsum('abcd,ap,bq,cr,ds->pqrs',qsource,T,T,T,T,optimize=True)
    fit,fit_res=hx.fit_quartics(qmass[:8,:8,:8,:8])
    require(np.allclose(fit,original['eft2_lambda_Y'],rtol=3e-9,atol=1e-10),
            'renormalisable EFT2 quartic boundary changed')
    require(np.allclose(initial[57:63],fit,rtol=3e-9,atol=1e-10),
            'saved EFT2 packed scalar state inconsistent')
    mu_x,scalar_x,sol=hx.solve_threshold(my,initial)
    require(abs(mu_x/original['mu_X_corrected_gev']-1)<2e-9,
            'EFT2 crossing has shifted from saved baseline')
    template=coupled.template_at_context(context)
    ids=static_hx['indices']
    deltas={};responses={};diagnostics={}
    for mode in ('full-mass','hard-massless-X'):
        delta,diag=bubble.bubble_boundary(Cfull,qmass,my,float(m2[0]),mode=mode)
        shifted=Cy+delta
        require(max(projection.symmetry_residual(shifted).values())<2e-6,
                'candidate boundary violates tensor symmetry')
        response=mass_response(shifted,my,mu_x,sol,static_hx,template,ids,
                               samples=eft2_samples,rtol=rtol)
        deltas[mode]=delta
        diagnostics[mode]=diag
        responses[mode]=response
    # Re-evolve unshifted boundary as a strict no-correction parity control.
    control=mass_response(Cy,my,mu_x,sol,static_hx,template,ids,
                          samples=eft2_samples,rtol=rtol)
    reference=np.asarray(original['takagi_masses_ev_descending'],float)
    control_error=float(np.linalg.norm(np.asarray(control['masses_ev_descending'])-reference)/
                        max(np.linalg.norm(reference),1e-30))
    require(control_error<3e-8,'unchanged one-loop evolution does not reproduce baseline')
    report={
      'format':FORMAT,'status':STATUS,'model':model,
      'source_trajectory_id':folder.name,
      'source_input_SHA256':source_hashes,
      'source_output_SHA256':{str(x.relative_to(folder)):digest(x) for x in (yj,yz,bj,bz)},
      'runner_sha256':digest(Path(__file__).resolve()),
      'bubble_module_sha256':digest(Path(bubble.__file__).resolve()),
      'source_EFT2_runner_sha256':digest(Path(sequential.__file__).resolve()),
      'mu_Y_gev':my,'mu_X_gev':mu_x,
      'mass_ratio_X_over_Y':float(math.sqrt(m2[0])/my),
      'quartic_basis_fit_residual':fit_res,
      'unshifted_replay_relative_mass_error':control_error,
      'scalar_pair_YY_bubble_zero_at_muY_MSbar':True,
      'scalar_pair_XX_bubble_zero_at_muX_MSbar':
          abs(bubble.x_equal_mass_bubble_at_threshold(mu_x))<1e-14,
      'Y_scalar_bubble_LLphi_phi_matching_candidate_computed':True,
      'full_finite_Y_matching_computed':False,
      'full_finite_X_matching_computed':False,
      'all_finite_scalar_matching_certified':False,
      'physical_neutrino_mass_prediction_certified':False,
      'full_three_generation_SM_Yukawa_Wilson_anomalous_dimension':False,
      'all_quartic_coupling_maps_physically_certified':False,
      'finite_scalar_renormalisable_thresholds_computed':False,
      'near_degenerate_matching_expansion_validated':False,
      'sequential_higher_dim_operator_tower_included':False,
      'candidate_bubble_modes':diagnostics,
      'baseline_masses_ev_descending':reference.tolist(),
      'candidate_masses_ev_descending':
          {m:responses[m]['masses_ev_descending'] for m in responses},
      'candidate_weinberg_residual_at_X':
          {m:responses[m]['weinberg_projection_residual'] for m in responses},
      'candidate_relative_C5_shift_at_X':{
          m:float(np.linalg.norm(responses[m]['C5_at_X']-control['C5_at_X'])/
               max(np.linalg.norm(control['C5_at_X']),1e-100)) for m in responses},
      'bubble_schemes_NOT_equivalent':'full-mass retains power corrections; hard-massless-X is leading dimension-five hard-region approximation',
      'status_notes':('This report is a controlled single-(Wilson x quartic) scalar-bubble '
          'finite sensitivity study, NOT complete scalar-threshold matching. The full-mass '
          'variant includes unresummed higher-dimensional mass-ratio effects; the massless-X '
          'variant has no parametrically small expansion parameter for these benchmark masses. '
          'Independent Matchete amplitudes and complete EFT operator/renormalisation matching '
          'remain mandatory before any physical certification.'),
      'saved_tensor_archive':output_npz.name,
    }
    out.mkdir(parents=True,exist_ok=True)
    tmp=output_npz.with_suffix('.npz.tmp')
    with tmp.open('wb') as handle:
        np.savez_compressed(handle,
          C_at_Y_tree=Cy,
          delta_Y_scalar_bubble_full_mass=deltas['full-mass'],
          delta_Y_scalar_bubble_hard_massless_X=deltas['hard-massless-X'],
          C5_at_X_baseline=control['C5_at_X'],
          C5_at_X_full_mass=responses['full-mass']['C5_at_X'],
          C5_at_X_hard_massless_X=responses['hard-massless-X']['C5_at_X'],
          C5_at_100gev_baseline=control['C5_at_100gev'],
          C5_at_100gev_full_mass=responses['full-mass']['C5_at_100gev'],
          C5_at_100gev_hard_massless_X=responses['hard-massless-X']['C5_at_100gev'],
          q_mass_basis_at_Y=qmass)
    tmp.replace(output_npz)
    continuous.write_json_atomic(output_json,report)
    return report


def summary(out,models):
    rows=[]
    for model in models:
        tag=continuous.model_tag(model)
        r=neutrino.load(out/(tag+'_scalar_bubble_candidate.json'))
        require(r['status']==STATUS and not r['physical_neutrino_mass_prediction_certified'],
                'invalid candidate record')
        rows.append(dict(model=model,
          mu_Y_GeV=r['mu_Y_gev'],mu_X_GeV=r['mu_X_gev'],
          mass_ratio_XY=r['mass_ratio_X_over_Y'],
          B0_full_mass=r['candidate_bubble_modes']['full-mass']['B0_YX'],
          B0_hard_massless_X=r['candidate_bubble_modes']['hard-massless-X']['B0_YX'],
          baseline_largest_eV=r['baseline_masses_ev_descending'][0],
          full_mass_candidate_largest_eV=r['candidate_masses_ev_descending']['full-mass'][0],
          hard_massless_X_candidate_largest_eV=r['candidate_masses_ev_descending']['hard-massless-X'][0],
          full_physical_matching='NOT_CERTIFIED'))
    continuous.write_json_atomic(out/'summary.json',{
        'format':FORMAT,'status':STATUS,'physical_certification':False,'rows':rows})
    with (out/'summary.csv').open('w',encoding='utf8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    return rows


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    choice=ap.add_mutually_exclusive_group(required=True)
    choice.add_argument('--all',action='store_true')
    choice.add_argument('--model',choices=MODELS)
    ap.add_argument('--run-id',default='continuous_v1',type=continuous.run_id)
    ap.add_argument('--source-root',type=Path,default=continuous.OUTPUT)
    ap.add_argument('--y-output-id',default=sequential.DEFAULT_Y_ID)
    ap.add_argument('--base-output-id',default=sequential.DEFAULT_OUTPUT_ID)
    ap.add_argument('--output-id',default=DEFAULT_OUTPUT_ID)
    ap.add_argument('--samples',type=int,default=65)
    ap.add_argument('--eft2-samples',type=int,default=33)
    ap.add_argument('--rtol',type=float,default=5e-7)
    args=ap.parse_args()
    for name in (args.run_id,args.y_output_id,args.base_output_id,args.output_id):
        continuous.run_id(name)
    require(len(set((args.y_output_id,args.base_output_id,args.output_id)))==3,
            'source/destination study IDs must differ')
    require(args.samples>=17 and args.samples<=129 and args.samples%2==1,
            'samples must be odd 17..129')
    require(args.eft2_samples>=9 and args.eft2_samples%2==1,
            'eft2 samples must be odd >=9')
    require(math.isfinite(args.rtol) and 1e-10<args.rtol<1e-3,'invalid Wilson tolerance')
    folder=args.source_root/args.run_id
    require(folder.is_dir(),'missing continuous run '+str(folder))
    out=folder/args.output_id
    errors=[]
    for model in (MODELS if args.all else (args.model,)):
        try:
            r=calculate(model,folder,out,y_id=args.y_output_id,base_id=args.base_output_id,
                samples=args.samples,eft2_samples=args.eft2_samples,rtol=args.rtol)
            print(f'{model}: Y/X={r["mass_ratio_X_over_Y"]:.7f} '
                  f'B0full={r["candidate_bubble_modes"]["full-mass"]["B0_YX"]:.7g} '
                  f'masses(full)={r["candidate_masses_ev_descending"]["full-mass"]} '
                  'physical matching=NOT CERTIFIED',flush=True)
        except (ValueError,OSError,KeyError,RuntimeError,TypeError,
                ArithmeticError,np.linalg.LinAlgError) as exc:
            errors.append(model)
            print(f'{model}: BLOCKED {type(exc).__name__}: {exc}',flush=True)
    if not errors and args.all:
        rows=summary(out,MODELS)
        print('Scalar-bubble candidate summary:',out/'summary.csv','models=',len(rows))
    return 2 if errors else 0

if __name__=='__main__':raise SystemExit(main())
