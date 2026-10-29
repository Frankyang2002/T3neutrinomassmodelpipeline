"""B(-1)/C(-1) coupled L/eC EFT1 Wilson tensor candidate, on saved scalar RG flow.

Evolves the complete seven-term one-generation L/eC psi^2 phi^2 closure,
retaining six physical lepton flavor pairs. Includes the published Matchete
LLSS seed, grouped fermion hard boundary, baseline scalar quartics and Stage22
additional quartics (latter with provisional Stage29 map). It does not compute
finite scalar matching nor the full three-family SM Yukawa Wilson anomalous
matrix. Fail-closed direct/Stage17 source parity checks are mandatory.
"""
from __future__ import annotations
import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping
import math

import numpy as np
from scipy.integrate import simpson

from scripts import RunT3BCExtraWilsonSequential as selected
from scripts import RunT3DirectWilsonSequential as direct
from scripts import RunT3FullScalarContinuous as continuous
from scripts import RunT3FullScalarNeutrino as neutrino
from Numerical.experimental import T3BCCoupledWilsonTensor as tensor

ROOT=continuous.ROOT
MODELS=selected.MODELS
FORMAT='t3-bc-coupled-lec-wilson-conditional-v1'
STATUS='CoupledLEcWilsonCandidateOnFullScalarEFT1'
DEFAULT_OUTPUT_ID='bc_coupled_wilson_v1'


def require(ok:bool, message:str)->None:
    if not ok:raise ValueError('B/C coupled sequential: '+message)


def sources(model:str)->dict[str,Path]:
    base=selected.sources(model)
    assert model in MODELS
    stage='T3_'+model[0]+'_alpha_m1'
    base.update({
       'quartic_seed':ROOT/f'output/full/hypercharge/{stage}/data/eft1_after_F_scalar_quartic_seed.json',
       'old_eft1_rgbeta':ROOT/f'output/full/hypercharge/{stage}/data/eft1_rgbeta_rge.json',
       'stage22':ROOT/'output/general_t3_stage22/matchete_tensors.json',
       'stage29':ROOT/'output/general_t3_stage22/stage29_convention_uv_handoff.json',
    })
    return base


def fingerprints(paths:Mapping[str,Path])->dict[str,str]:
    return {name:selected.digest(path) for name,path in paths.items()}


def full_scalar_background(model:str,record:Mapping[str,Any],sidecar:Mapping[str,Any],*,samples:int):
    """Re-run EFT1 only, returning the entire general-scalar state at each node."""
    from Numerical.experimental.GeneralT3ScalarRGState import (
      GeneralT3ScalarRGState,pack_general_scalar_state,unpack_general_scalar_state,
    )
    from scripts import IntermediateEFTRunning as s34
    require(model in MODELS and 17<=samples<=129 and samples%2==1,'--samples odd 17..129')
    require(sidecar.get('status')=='CandidateQuadraticAndFullQuarticEFT1Integrated' and
            sidecar.get('format')==continuous.qeft.FORMAT and sidecar.get('model')==model and
            record.get('quadratic_kind')=='mixed', model+': saved full-scalar EFT1 sidecar mismatch')
    registry=continuous.load_registry(continuous.REGISTRY)
    base=s34.decode_intermediate_state(sidecar['initial_baseline'])
    extra={name:neutrino.complex_pair(raw,name) for name,raw in sidecar['initial_new_quartics'].items()}
    initial=GeneralT3ScalarRGState(base,extra).validated(registry)
    vector,layout=pack_general_scalar_state(initial,registry)
    mu2=neutrino.complex_pair(sidecar['mu2_initial_gev2'],'initial mu2')
    vector=np.concatenate((vector,[mu2.real,mu2.imag]))
    path=sources(model)['eft1_beta']
    betas=continuous.beta_inventory(continuous.read(path), model=model,
             stage='EFT1',representation=base.representation,layout=layout,
             registry=registry,kind='mixed')
    result=continuous.integrate_segment(vector,mu_start=s34.START,
              mu_end=s34.STOP,layout=layout,registry=registry,stage='EFT1',
              kind='mixed',betas=betas,points=samples)
    states=tuple(unpack_general_scalar_state(result.y[:layout.size,k],layout,
                  registry,uv=False,mu_gev=float(math.exp(result.t[k])))
                  for k in range(samples))
    prior=neutrino.decode_eft(record['eft1_scalar_endpoint_state'])
    errors={}
    for name in ('lambdaT3','mS1Sq','mS2Sq',*selected.ACTIVE_QUARTICS):
        got=complex(getattr(states[-1].baseline,name));expected=complex(getattr(prior,name))
        err=abs(got-expected)/max(abs(expected),1e-12)
        require(err<2e-6,'saved endpoint parity failed for '+name)
        errors[name]=float(err)
    require(set(sidecar['final_new_quartics'])==set(states[-1].extra),
            'new quartic endpoint inventory mismatch')
    for name,value in states[-1].extra.items():
        target=neutrino.complex_pair(sidecar['final_new_quartics'][name],name)
        require(abs(value-target)<2e-6*max(abs(target),1e-10),'new quartic endpoint mismatch '+name)
    final_mu2=complex(result.y[layout.size,-1],result.y[layout.size+1,-1])
    target_mu2=neutrino.complex_pair(sidecar['mu2_final_gev2'],'mixed mu2 endpoint')
    require(abs(final_mu2-target_mu2)<2e-6*max(abs(target_mu2),1e-10),'mixed quadratic endpoint mismatch')
    return states,np.exp(np.asarray(result.t,dtype=float)),{'endpoint_parity':errors,'eft1_nfev':int(result.nfev)}


def template_at_context(context:Any)->dict:
    from RGE.matching.WeinbergTensorAdapter import build_weinberg_wilson_tensor
    return build_weinberg_wilson_tensor(context.model,context.fermion_basis,1)


def embed_weinberg(C:np.ndarray,template:Mapping, ids:np.ndarray, value:np.ndarray)->None:
    """Add complete hard C5 at the grouped F threshold in real-field basis."""
    for (i,j,a,b),w in template.items():
        i=i-1;j=j-1;a=a-1;b=b-1
        if i not in ids or j not in ids:continue
        ii=int(np.flatnonzero(ids==i)[0]);jj=int(np.flatnonzero(ids==j)[0])
        for f,(p,q) in enumerate(tensor.PAIRS):
            C[ii,jj,a,b,f]+=complex(w)*value[p,q]


def check_source_parity(model:str, C_tree:np.ndarray, backgrounds:list,
                        static:dict, quartet:tensor.QuarticEvaluator,
                        template:Mapping,ids:np.ndarray,
                        direct_source:Any,masses:np.ndarray,y1:np.ndarray,y2:np.ndarray,
                        *, mu:np.ndarray,hard:np.ndarray,original:np.ndarray,
                        selected_original:np.ndarray):
    """Check direct and Stage17 source normalization *before* solving Wilson ODE."""
    direct_betas=[]; extra_betas=[]
    first_active=None
    for j,st in enumerate(backgrounds):
        q_base=quartet.evaluate(st,include_extra=False)
        q_selected=quartet.evaluate(st,include_extra=True,omit_active=False)
        q_other=quartet.evaluate(st,include_extra=True,omit_active=True)
        kb=tensor.kernels(st,static,q_base)
        ks=tensor.kernels(st,static,q_selected-q_other+q_base)
        be=tensor.rhs_master(C_tree,ks)-tensor.rhs_master(C_tree,kb)
        bd=tensor.rhs_master(C_tree,kb)
        fdir,res1=tensor.weinberg_projection(bd,template,ids)
        fextra,res2=tensor.weinberg_projection(be,template,ids)
        require(max(res1,res2)<2e-5,'Weinberg tensor source projection failed')
        direct_betas.append(tensor.flavor_to_matrix(fdir))
        extra_betas.append(tensor.flavor_to_matrix(fextra))
        if j==0:first_active=(fextra,q_selected-q_other)
    # The previous direct runner integrated lambdaT3(t) and a full-flavor
    # coefficient, and the earlier selected source used precisely Stage17.
    dflow=simpson(np.array(direct_betas),x=np.log(mu),axis=0)/tensor.LOOP
    aflow=simpson(np.array(extra_betas),x=np.log(mu),axis=0)/tensor.LOOP
    expected_direct=neutrino.validated_c5(original-hard)
    expected_extra=neutrino.validated_c5(selected_original)
    er_direct=np.linalg.norm(dflow-expected_direct)/max(np.linalg.norm(expected_direct),1e-300)
    er_extra=np.linalg.norm(aflow-expected_extra)/max(np.linalg.norm(expected_extra),1e-300)
    require(er_direct<3e-4,model+': master direct source disagrees with saved hierarchical beta: '+str(er_direct))
    require(er_extra<3e-4,model+': master extra source disagrees with certified Stage17: '+str(er_extra))
    return {'direct_source_parity':float(er_direct),'stage17_extra_source_parity':float(er_extra),
            'first_extra_beta_norm':float(np.linalg.norm(first_active[0]))}


def calculate(model:str,*,folder:Path,out:Path,samples:int=65,resume:bool=False,
              direct_output_id:str=direct.DEFAULT_OUTPUT_ID,
              selected_output_id:str=selected.DEFAULT_OUTPUT_ID,
              rtol:float=5e-7)->dict:
    require(model in MODELS,'only B_-1 and C_-1 supported')
    tag=continuous.model_tag(model)
    cpath=folder/tag/'continuous_trajectory.json'
    record=neutrino.verify_continuous(model,cpath)
    sidepath=folder/tag/record['saved_candidate_eft1_endpoint']
    sidecar=neutrino.load(sidepath)
    prev_direct_path=folder/direct_output_id/f'{tag}_direct_wilson.json'
    prev_extra_path=folder/selected_output_id/f'{tag}_extra_wilson.json'
    prior=neutrino.load(prev_direct_path)
    prior_extra=neutrino.load(prev_extra_path)
    require(prior.get('format')==direct.FORMAT and prior.get('status')==direct.STATUS and
        prior.get('model')==model and prior.get('saved_points')==samples,
        model+': incompatible saved direct sequential baseline')
    require(prior_extra.get('format')==selected.FORMAT and prior_extra.get('status')==selected.STATUS and
        prior_extra.get('model')==model and prior_extra.get('eft1_source_samples')==samples,
        model+': incompatible saved B/C selected extra source')
    paths=sources(model)
    paths.update({'continuous':cpath,'sidecar':sidepath,'prior_direct':prev_direct_path,
                  'prior_extra':prev_extra_path,'runner':Path(__file__).resolve(),
                  'engine':Path(tensor.__file__).resolve(),
                  'direct_runner':ROOT/'scripts/RunT3DirectWilsonSequential.py',
                  'extra_runner':ROOT/'scripts/RunT3BCExtraWilsonSequential.py'})
    hashes=fingerprints(paths)
    output_path=out/f'{tag}_coupled_wilson.json'
    if output_path.is_file():
        old=neutrino.load(output_path)
        require(resume and old.get('format')==FORMAT and old.get('status')==STATUS and
                old.get('input_sha256')==hashes and old.get('eft1_samples')==samples and
                old.get('wilson_rtol')==rtol,
                model+': output already exists or source hash changed: choose new --output-id')
        return old
    stage17=neutrino.load(paths['stage17']);stage19=neutrino.load(paths['stage19'])
    seed=neutrino.load(paths['seed']);wilson=neutrino.load(paths['wilson_rge'])
    matched=neutrino.load(paths['matching'])
    cert,norm=selected.certify_saved_selected_source(model,stage17,stage19,seed,wilson,matched)
    uv=neutrino.decode_uv(record['uv_fermion_threshold_state'])
    eft=neutrino.decode_eft(record['eft1_scalar_endpoint_state'])
    mf,y1,y2,heavyres=direct.heavy_inputs(uv)
    source=direct.compile_direct_source(paths['matching'],paths['wilson_rge'])
    states,mu,background=full_scalar_background(model,record,sidecar,samples=samples)
    from RGE.running.intermediate.ScalarOnlyWilsonTensorRGE import build_scalar_only_rge_context
    context=build_scalar_only_rge_context(paths['seed'],paths['quartic_seed'],paths['old_eft1_rgbeta'])
    static=tensor.build_static_kinetics(context)
    require(static['nf']==context.inputs.fermion_dimension and static['ns']==12,
            'unexpected master gauge/Yukawa input dimensions')
    mat22=neutrino.load(paths['stage22']);map29=neutrino.load(paths['stage29'])
    quartet=tensor.QuarticEvaluator.build(context.quartic,mat22,map29,model)
    template=template_at_context(context)
    ids=static['indices']
    flavors=tensor.flavor_pairs(mf,y1,y2)
    tree=tensor.independent_tensor(context.coefficient,flavors)
    initial=tree.copy()
    hard,reloc,grouped=direct.grouped_hard_boundary(source,mf,y1,y2,
                           complex(states[0].baseline.lambdaT3),mu_f=float(mu[0]))
    embed_weinberg(initial,template,ids,grouped)
    projector_initial,res=tensor.weinberg_projection(initial,template,ids)
    require(res<1e-9 and np.allclose(tensor.flavor_to_matrix(projector_initial),grouped,rtol=1e-10,atol=1e-20),
            model+': grouped fermion boundary embedding parity failed')
    direct_s=neutrino.complex_matrix(prior['c5_at_scalar_threshold_gev_inverse'],'direct C5')
    extra_s=neutrino.complex_matrix(prior_extra['selected_extra_delta_C5_at_scalar_threshold_gev_inverse'],'extra C5')
    parity=check_source_parity(model,tree,states,static,quartet,template,ids,
                    source,mf,y1,y2,mu=mu,hard=grouped,original=direct_s,
                    selected_original=extra_s)
    # Generate the complete one-generation L/eC closure: 6 fermion pairs x
    # 78 scalar pairs x 6 physical lepton flavor entries.
    computed=[tensor.kernels(st,static,quartet.evaluate(st)) for st in states]
    y,solution=tensor.integrate_coupled(initial,mu,computed,rtol=rtol)
    projected,residual=tensor.weinberg_projection(y[..., -1],template,ids)
    c5=neutrino.validated_c5(tensor.flavor_to_matrix(projected))
    old_c5=neutrino.validated_c5(direct_s+extra_s)
    low=neutrino.low_energy(model,c5,eft)
    ratio=float(np.linalg.norm(c5-old_c5)/max(np.linalg.norm(old_c5),1e-300))
    back={'format':FORMAT,'status':STATUS,'model':model,'eft1_samples':samples,
          'input_sha256':hashes,'wilson_rtol':float(rtol),'coupled_master_rge':'RGE.general.WilsonTensorRGE seven terms',
          'fermion_index_closure':'L1,L2,eC one-generation, 468 independent complex gauge-index components per flavor pair',
          'lepton_flavor_pairs':6,'sm_yukawa_projection':'one-generation SM master with sqrt(Tr[Y Ydag]) effective Yukawas',
          'source_normalization_checks':parity,
          'quartic_potential_source':'existing Matchete baseline seed + Stage17 certified active two + conditional Stage22/29 remaining five',
          'additional_quartic_mapping_physically_certified':False,
          'finite_scalar_threshold_matching_included':False,
          'full_three_generation_wilson_anomalous_dimension_included':False,
          'full_general_scalar_matching_certified':False,
          'one_generation_gauge_index_wilson_evolution_completed':True,
          'hard_threshold_grouped_scale_gev':float(mu[0]),
          'scalar_threshold_gev':float(mu[-1]),
          'scalar_ode_nfev':background['eft1_nfev'],
          'wilson_ode_nfev':int(solution.nfev),
          'wilson_ode_status':solution.message,
          'wilson_boundary_weinberg_residual':res,
          'scalar_threshold_weinberg_residual':residual,
          'saved_scalar_endpoint_parity':background['endpoint_parity'],
          'prior_direct_plus_selected_C5':neutrino.matrix_payload(old_c5),
          'coupled_C5_at_scalar_threshold_gev_inverse':neutrino.matrix_payload(c5),
          'relative_coupled_C5_vs_prior_selected':ratio,
          'prior_selected_masses_ev_descending':prior_extra['updated_masses_ev_descending'],
          'hard_grouped_C5_at_fermion_scale_gev_inverse':neutrino.matrix_payload(grouped),
          'uv_running_repeated':False,'eft1_scalar_running_repeated':True,
          'full_sequential_matching_performed':False,
          'notes':('Conditional coupled L/eC Wilson one-generation gauge-index closure. '
             'General scalar additional quartic Stage29 maps are not physical coupling '
             'certificates; three-family SM Yukawa anomalous mixing and finite one-loop '
             'scalar threshold corrections remain absent. Not a complete C5 prediction.'),
          **low}
    neutrino.atomic_json(output_path,back)
    return back


def summary(out:Path, chosen:tuple[str,...])->dict:
    rows=[]
    for model in chosen:
        path=out/f'{continuous.model_tag(model)}_coupled_wilson.json'
        x=neutrino.load(path)
        require(x.get('status')==STATUS and x.get('format')==FORMAT and
                x.get('model')==model and x.get('full_general_scalar_matching_certified') is False,
                model+': incompatible summary source')
        rows.append({'model':model,'status':x['status'],
              'relative_C5_change_vs_previous_selected':x['relative_coupled_C5_vs_prior_selected'],
              'direct_source_parity':x['source_normalization_checks']['direct_source_parity'],
              'extra_source_parity':x['source_normalization_checks']['stage17_extra_source_parity'],
              'weinberg_projection_residual':x['scalar_threshold_weinberg_residual'],
              'old_largest_mnu_ev':x['prior_selected_masses_ev_descending'][0],
              'new_largest_mnu_ev':x['takagi_masses_ev_descending'][0],
              'new_middle_mnu_ev':x['takagi_masses_ev_descending'][1],
              'new_smallest_mnu_ev':x['takagi_masses_ev_descending'][2]})
    payload={'format':FORMAT,'status':'CoupledLEcWilsonSummary',
             'model_count':len(rows),'full_general_scalar_matching_certified':False,
             'rows':rows}
    neutrino.atomic_json(out/'summary.json',payload)
    path=out/'summary.csv';path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.tmp')
    with temp.open('w',newline='',encoding='utf-8') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)
    temp.replace(path)
    return payload


def main()->int:
    parser=argparse.ArgumentParser(description=__doc__)
    s=parser.add_mutually_exclusive_group(required=True)
    s.add_argument('--all',action='store_true')
    s.add_argument('--model',choices=MODELS)
    parser.add_argument('--run-id',required=True,type=continuous.run_id)
    parser.add_argument('--source-root',type=Path,default=continuous.OUTPUT)
    parser.add_argument('--output-id',default=DEFAULT_OUTPUT_ID)
    parser.add_argument('--direct-output-id',default=direct.DEFAULT_OUTPUT_ID)
    parser.add_argument('--selected-output-id',default=selected.DEFAULT_OUTPUT_ID)
    parser.add_argument('--samples',type=int,default=65)
    parser.add_argument('--rtol',type=float,default=5e-7)
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--combine-only',action='store_true')
    args=parser.parse_args()
    for val in (args.output_id,args.direct_output_id,args.selected_output_id):continuous.run_id(val)
    require(len({args.output_id,args.direct_output_id,args.selected_output_id})==3,
            'output IDs must be different')
    require(args.samples>=17 and args.samples<=129 and args.samples%2==1,'--samples odd 17..129')
    require(math.isfinite(args.rtol) and 1e-10<args.rtol<1e-3,'--rtol outside allowed range')
    require(not args.combine_only or args.all,'--combine-only requires --all')
    folder=args.source_root/args.run_id
    require(folder.is_dir(),'missing continuous run folder '+str(folder))
    out=folder/args.output_id
    models=MODELS if args.all else (args.model,)
    errors=[]
    if not args.combine_only:
        for model in models:
            try:
                x=calculate(model,folder=folder,out=out,samples=args.samples,
                    rtol=args.rtol,resume=args.resume,
                    direct_output_id=args.direct_output_id,
                    selected_output_id=args.selected_output_id)
                print(f'{model}: coupled L/eC Wilson candidate; '
                      f'change vs selected={x["relative_coupled_C5_vs_prior_selected"]:.7g}; '
                      f'direct parity={x["source_normalization_checks"]["direct_source_parity"]:.3g}; '
                      f'extra parity={x["source_normalization_checks"]["stage17_extra_source_parity"]:.3g}; '
                      f'masses={x["takagi_masses_ev_descending"]}',flush=True)
            except (OSError,ValueError,RuntimeError,KeyError,TypeError,
                    ArithmeticError,np.linalg.LinAlgError) as err:
                errors.append(model)
                print(f'{model}: BLOCKED {type(err).__name__}: {err}',flush=True)
    if errors:
        print('Blocked:',', '.join(errors))
        return 2
    if args.all:
        payload=summary(out,models)
        print('Coupled candidate table:',out/'summary.csv','models=',payload['model_count'])
    return 0


if __name__=='__main__':raise SystemExit(main())
