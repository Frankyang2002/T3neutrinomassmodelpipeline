"""Large-batch B/C(-1) sequential *conditional* numerical EFT2(H,X) continuation.

Reads the previously saved Y-threshold projection, verifies that its complete
source inventory is unchanged, replays the SAME EFT1 scalar background to Y,
extracts the closed inert-doublet quartics, integrates the genuine X-only EFT2
renormalisable RGEs and the complete 8-real-scalar project seven-term Wilson
system, solves mu_X=m_X(mu_X), tree-projects to SM C5, evolves SM+C5 to 100 GeV,
and reports provisional neutrino mass singular values.

This does NOT calculate the finite Y/X one-loop matching; physical completion is
NOT certified. No remote/local production source files are changed.

  python -u -X utf8 -m scripts.RunT3BCEFT2Sequential --all --run-id continuous_v1
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from hashlib import sha256
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from Numerical.experimental import T3BCEFT2HX as eft2
from Numerical.experimental import T3BCCoupledWilsonTensor as wilson
from Numerical.experimental import T3BCScalarYThresholdProjection as projection
from Numerical.experimental.GeneralT3MixedScalarMassBasis import MixedDoubletMassBasis
from scripts import RunT3BCYThresholdProjection as ystage

coupled=ystage.coupled
neutrino=ystage.neutrino
continuous=ystage.continuous
MODELS=ystage.MODELS
FORMAT='t3-bc-conditional-sequential-hx-to-sm-v1'
STATUS='EFT2HXRunningXCrossingSMWeinbergConditionalFiniteScalarMatchingMissing'
DEFAULT_OUTPUT_ID='bc_eft2_sequential_v1'
DEFAULT_Y_ID=ystage.DEFAULT_OUTPUT_ID


def require(test, message):
    if not test:raise ValueError('B/C EFT2 sequential: '+message)


def digest(path):
    require(path.is_file(),'missing input '+str(path))
    return sha256(path.read_bytes()).hexdigest()


def interpolate_sm(states,index,w):
    a=states[index].baseline.sm;b=states[index+1].baseline.sm
    names=('gY','g2','g3','lambdaH')
    fields={n:float((1-w)*getattr(a,n)+w*getattr(b,n)) for n in names}
    for n in ('ye','yu','yd'):
        fields[n]=(1-w)*np.asarray(getattr(a,n),complex)+w*np.asarray(getattr(b,n),complex)
    from Numerical.core.State import SMNumericalState
    return SMNumericalState(**fields).validated()


def read_input(model,folder,y_id,samples,rtol):
    tag=continuous.model_tag(model)
    yp=folder/y_id/f'{tag}_y_threshold_wilson.json'
    zp=folder/y_id/f'{tag}_y_threshold_wilson.npz'
    data=neutrino.load(yp)
    require(data.get('format')==ystage.FORMAT and data.get('status')==ystage.STATUS,
            'incompatible Y-threshold report')
    require(data.get('model')==model and data.get('source_trajectory_id')==folder.name,
            'incorrect Y-threshold model/trajectory')
    require(data.get('finite_Y_threshold_delta_C_computed') is False and
            data.get('full_general_scalar_matching_certified') is False,
            'upstream physical-certification semantics changed')
    require(data.get('eft1_samples')==samples and data.get('wilson_rtol')==rtol,
            'upstream Y run used different sample count / Wilson tolerance')
    require(zp.is_file(),'missing Y tensor archive')
    record=neutrino.verify_continuous(model,folder/tag/'continuous_trajectory.json')
    paths=ystage.source_paths(model,folder,record,continuous.run_id('direct_wilson_sequential_v1'),
                              continuous.run_id('bc_extra_wilson_sequential_v1'))
    now=coupled.fingerprints(paths)
    require(now==data['input_sha256'],
            'upstream EFT1 code/source SHA256 changed; rerun Y stage with fresh output-id')
    with np.load(zp,allow_pickle=False) as z:
        C=np.array(z['C_EFT2_tree_at_Y_HX_real'],dtype=complex)
        T=np.array(z['real_old_from_new'],dtype=float)
        U=np.array(z['U_complex'],dtype=complex)
        mu_y=float(np.asarray(z['Y_threshold_gev']).ravel()[0])
        m2=np.array(z['Y_eigenmasses_sq_gev2'],dtype=float)
    require(C.shape==(3,3,8,8,6) and T.shape==(12,12) and U.shape==(2,2),
            'unexpected source tensor shape')
    require(max(projection.symmetry_residual(C).values())<1e-8,
            'Y tensor violates exchange symmetries')
    require(np.isclose(mu_y,float(data['Y_threshold_gev']),rtol=1e-11) and
            np.linalg.norm(T-projection.mixing_real_matrix(U))<1e-10,
            'Y tensor and report disagree')
    return data,C,T,U,mu_y,m2,record,{'Y_report':digest(yp),'Y_archive':digest(zp),
                                     **now}


def calculate(model,folder,out,*,y_id=DEFAULT_Y_ID,samples=65,
              rtol=5e-7,eft2_samples=33,low_scale=100.,vev=246.22,resume=False):
    require(model in MODELS,'unknown B/C model')
    tag=continuous.model_tag(model)
    jp=out/f'{tag}_eft2_sequential.json';zp=out/f'{tag}_eft2_sequential.npz'
    if jp.exists() or zp.exists():
        require(resume and jp.is_file() and zp.is_file(),
                'output exists: choose a NEW --output-id, or use --resume')
        earlier=neutrino.load(jp)
        _source,_cy,_T,_U,_my,_sq,_rec,hashes=read_input(
            model,folder,y_id,samples,rtol)
        require(earlier.get('format')==FORMAT and earlier.get('status')==STATUS
                and earlier.get('source_sha256')==hashes
                and earlier.get('runner_sha256')==digest(Path(__file__).resolve())
                and earlier.get('eft2_module_sha256')==digest(Path(eft2.__file__).resolve())
                and earlier.get('eft2_wilson_samples')==eft2_samples,
                'saved output provenance or settings changed; choose NEW --output-id')
        return earlier
    data,Cy,T,U,muy,masses_sq,record,source_hashes=read_input(
        model,folder,y_id,samples,rtol)
    sidecar=neutrino.load(folder/tag/record['saved_candidate_eft1_endpoint'])
    states,mu,mu2,background=ystage.scalar_background_with_mu2(
        model,record,sidecar,samples=samples)
    basis,again,idx,crossing=ystage.locate_y_crossing(model,states,mu,mu2)
    require(abs(again/muy-1)<2e-11 and
            np.linalg.norm(basis.U-U)<2e-10 and
            np.allclose(basis.eigenvalues_sq,masses_sq,rtol=2e-11),
            'Y crossing/basis cannot be reproduced from the unchanged EFT1 source')

    from RGE.running.intermediate.ScalarOnlyWilsonTensorRGE import build_scalar_only_rge_context
    paths=coupled.sources(model)
    context=build_scalar_only_rge_context(paths['seed'],paths['quartic_seed'],
                                           paths['old_eft1_rgbeta'])
    static_original=wilson.build_static_kinetics(context)
    static_hx=eft2.rotate_static(static_original,T)
    quartic_evaluator=wilson.QuarticEvaluator.build(
        context.quartic,neutrino.load(paths['stage22']),
        neutrino.load(paths['stage29']),model)
    t=np.log(muy);tl=np.log(mu[idx]);tr=np.log(mu[idx+1])
    w=float((t-tl)/(tr-tl))
    require(0<w<1,'Y crossing interpolation weight outside adjacent RGBeta nodes')
    # Interpolate the *potential tensor*, not a set of guessed quartic names.
    qa=quartic_evaluator.evaluate(states[idx]);qb=quartic_evaluator.evaluate(states[idx+1])
    qy=(1-w)*qa+w*qb
    quartic_hx=eft2.rotate_potential(qy,T)
    fitted,fit_res=eft2.fit_quartics(quartic_hx)
    sm=interpolate_sm(states,idx,w)
    lambda_h_error=float(abs(sm.lambdaH-fitted[0])/max(1e-9,abs(sm.lambdaH)))
    require(lambda_h_error<1e-6,
            'Higgs quartic tensor disagrees with source SM coupling: '+str(lambda_h_error))
    # No independent EFT1 Higgs quadratic state is stored in the source.
    # Tree-level proxy is therefore explicitly not a physical threshold match.
    mh2=-0.5*fitted[0]*vev**2
    mx2=float(masses_sq[0])
    initial=eft2.pack_running(sm,fitted,mx2,mh2)
    mux,y_x,scalar_solution=eft2.solve_threshold(muy,initial)
    d_x=eft2.unpack_running(y_x)
    require(muy>mux>low_scale and d_x['mx2']>0,'unphysical ordered scalar thresholds')
    Cx,wsol=eft2.evolve_wilson(Cy,muy,mux,scalar_solution,static_hx,
                               samples=eft2_samples,rtol=rtol)
    require(max(projection.symmetry_residual(Cx).values())<2e-6,
            'EFT2 Wilson violated exchange symmetry')
    template=coupled.template_at_context(context)
    ids=static_hx['indices']
    c5_y,ry=wilson.weinberg_projection(Cy,template,ids,threshold=5e-5)
    c5_x,rx=wilson.weinberg_projection(Cx,template,ids,threshold=5e-5)
    c5matrix=neutrino.validated_c5(wilson.flavor_to_matrix(c5_x))
    # Full EFT2->SM tree restriction: the H-only Weinberg coefficient.
    # All finite X-threshold corrections remain unknown.
    from Numerical.running.SMWeinbergEvolution import (
        SMWeinbergInitialConditions,evolve_sm_weinberg)
    final_state=SMWeinbergInitialConditions(
        gY=d_x['gY'],g2=d_x['g2'],g3=d_x['g3'],lambdaH=d_x['quartics'][0],
        ye=d_x['ye'],yu=d_x['yu'],yd=d_x['yd'],K=c5matrix).validated()
    low=evolve_sm_weinberg(final_state,mux,low_scale,rtol=1e-8,atol=1e-11)
    from Numerical.fitting.ParameterScan import (
        neutrino_mass_from_c5,charged_lepton_mass_basis_matrix)
    from Numerical.fitting.OscillationFit import prediction_from_mass_matrix
    m=neutrino_mass_from_c5(low.K,vev_gev=vev)
    flavour=charged_lepton_mass_basis_matrix(m,low.ye)
    eigen=np.linalg.svd(flavour,compute_uv=False)*1e9
    prediction=prediction_from_mass_matrix(flavour,ordering='NO')
    require(np.isfinite(eigen).all() and (eigen>=0).all(),
            'nonfinite conditional neutrino singular values')
    report={
        'format':FORMAT,'status':STATUS,'model':model,
        'source_trajectory_id':folder.name,'source_Y_output_id':y_id,
        'source_sha256':source_hashes,'runner_sha256':digest(Path(__file__).resolve()),
        'eft2_module_sha256':digest(Path(eft2.__file__).resolve()),
        'mu_Y_gev':muy,'mu_X_corrected_gev':mux,
        'mu_X_relative_crossing_residual':abs(np.sqrt(d_x['mx2'])/mux-1),
        'mu_X_provisional_full_EFT1_not_used':True,
        'eft2_integrated_with_Y_removed':True,
        'eft2_quartic_source':'rotated and fully fitted original EFT1 scalar quartic tensor',
        'eft2_lambda_basis':list(eft2.NAMES),
        'eft2_lambda_Y':fitted.tolist(),
        'eft2_lambda_X':d_x['quartics'].tolist(),
        'eft2_quartic_fit_relative_residual':fit_res,
        'eft2_higgs_quartic_source_parity_error':lambda_h_error,
        'eft2_scalar_beta':'Z2 inert-doublet one-loop 16pi^2 d/dln(mu)',
        'eft2_scalar_nfev':int(scalar_solution.nfev),
        'eft2_wilson_master':'project seven-term one-generation L/eC gauge-index closure',
        'eft2_wilson_nfev':int(wsol.nfev),
        'eft2_wilson_samples':eft2_samples,'wilson_rtol':rtol,
        'C5_weinberg_projection_residual_Y':ry,
        'C5_weinberg_projection_residual_X':rx,
        'C5_at_X_tree_gev_inverse':neutrino.matrix_payload(c5matrix),
        'C5_at_100gev_gev_inverse':neutrino.matrix_payload(low.K),
        'mnu_at_100gev_ev':neutrino.matrix_payload(flavour*1e9),
        'takagi_masses_ev_descending':eigen.tolist(),
        'oscillation_observables_normal_ordering':prediction.as_dict(),
        'takagi_residual':float(prediction.takagi_residual),
        'sm_weinberg_rge_nfev':int(low.nfev),
        'H_quadratic_muY_source':'TREE_LEVEL_PROXY_FROM_LAMBDA_H_AND_VEV_NOT_EFT1_RUNNING',
        'H_quadratic_muY_gev2':float(mh2),
        'finite_Y_threshold_delta_C_computed':False,
        'finite_X_threshold_delta_C5_computed':False,
        'finite_threshold_corrections_status':'BOTH_MISSING_NOT_PHYSICAL_ZERO',
        'physical_neutrino_mass_prediction_certified':False,
        'full_three_generation_SM_Yukawa_Wilson_anomalous_dimension':False,
        'additional_quartic_mapping_physically_certified':False,
        'full_sequential_matching_performed':False,
        'leading_log_sequential_candidate_completed':True,
        'tensor_archive_file':zp.name,
        'source_endpoint_checks':background,
        'notes':('Numerical outputs are conditional leading-log/tree scalar-threshold benchmarks. '
            'The Y and X finite one-loop thresholds were not derived and must NOT be '
            'interpreted as zero. One-generation Wilson gauge-index master and some '
            'remaining additional quartic maps are not fully certified. The Higgs '
            'quadratic input is a tree-level proxy because the source omits its running.')
    }
    out.mkdir(parents=True,exist_ok=True)
    tmp=zp.with_suffix('.npz.tmp')
    with tmp.open('wb') as f:
        np.savez_compressed(f,C_EFT2_at_Y=Cy,C_EFT2_at_X=Cx,
            C_EFT3_HH_at_X=Cx[:,:,:4,:4,:],
            C5_at_X=c5matrix,C5_at_100gev=low.K,
            lambda_inert_Y=fitted,lambda_inert_X=d_x['quartics'],
            Y_threshold_gev=np.array([muy]),X_threshold_gev=np.array([mux]),
            scalar_Y=initial,scalar_X=y_x)
    tmp.replace(zp)
    continuous.write_json_atomic(jp,report)
    return report


def summary(out,models):
    rows=[]
    for model in models:
        rec=neutrino.load(out/f'{continuous.model_tag(model)}_eft2_sequential.json')
        require(rec.get('status')==STATUS and rec.get('physical_neutrino_mass_prediction_certified') is False,
                'bad candidate output')
        rows.append({'model':model,'mu_Y_GeV':rec['mu_Y_gev'],
          'mu_X_GeV':rec['mu_X_corrected_gev'],
          'mnu_largest_eV':rec['takagi_masses_ev_descending'][0],
          'mnu_middle_eV':rec['takagi_masses_ev_descending'][1],
          'mnu_smallest_eV':rec['takagi_masses_ev_descending'][2],
          'quartic_fit_residual':rec['eft2_quartic_fit_relative_residual'],
          'finite_scalar_matching':'MISSING'})
    continuous.write_json_atomic(out/'summary.json',dict(format=FORMAT,
                physically_certified=False,rows=rows))
    with (out/'summary.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    return rows


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    which=ap.add_mutually_exclusive_group(required=True)
    which.add_argument('--all',action='store_true')
    which.add_argument('--model',choices=MODELS)
    ap.add_argument('--run-id',default='continuous_v1',type=continuous.run_id)
    ap.add_argument('--source-root',type=Path,default=continuous.OUTPUT)
    ap.add_argument('--y-output-id',default=DEFAULT_Y_ID)
    ap.add_argument('--output-id',default=DEFAULT_OUTPUT_ID)
    ap.add_argument('--samples',type=int,default=65)
    ap.add_argument('--rtol',type=float,default=5e-7)
    ap.add_argument('--eft2-samples',type=int,default=33)
    ap.add_argument('--resume',action='store_true')
    args=ap.parse_args()
    for item in (args.run_id,args.y_output_id,args.output_id):continuous.run_id(item)
    require(args.output_id!=args.y_output_id and
            args.output_id not in ('neutrino_v1','bc_coupled_wilson_v1'),
            'output-id would overwrite existing study')
    require(17<=args.samples<=129 and args.samples%2==1,'samples 17..129 odd')
    require(args.eft2_samples>=9 and args.eft2_samples%2==1,
            'eft2-samples >=9 odd')
    require(math.isfinite(args.rtol) and 1e-10<args.rtol<1e-3,'invalid rtol')
    folder=args.source_root/args.run_id
    require(folder.is_dir(),'continuous run folder missing: '+str(folder))
    out=folder/args.output_id;models=MODELS if args.all else (args.model,)
    errors=[]
    for m in models:
        try:
            r=calculate(m,folder,out,y_id=args.y_output_id,samples=args.samples,
                        rtol=args.rtol,eft2_samples=args.eft2_samples,resume=args.resume)
            print(f'{m}: muY={r["mu_Y_gev"]:.8f} muX={r["mu_X_corrected_gev"]:.8f} GeV; '
                  f'conditional masses(eV)={r["takagi_masses_ev_descending"]}; '
                  f'finite Y/X=MISSING',flush=True)
        except (ValueError,OSError,KeyError,RuntimeError,TypeError,
                ArithmeticError,np.linalg.LinAlgError) as exc:
            errors.append(m)
            print(f'{m}: BLOCKED {type(exc).__name__}: {exc}',flush=True)
    if not errors and args.all:
        rows=summary(out,models)
        print('Conditional EFT2+SM summary:',out/'summary.csv',
              'models=',len(rows),flush=True)
    return 2 if errors else 0


if __name__=='__main__':raise SystemExit(main())
