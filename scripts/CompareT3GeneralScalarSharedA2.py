"""Use saved 16x3 general-scalar Sobol fits; compare same-case A(-2) Yukawas.

This is a conditional candidate study reusing the exact existing exploratory
fitting/evaluation implementations. It does NOT replace scalar Lagrangians,
matching sources, or existing fitted results. Shared comparisons evaluate
only UV Yukawa substitutions; other model-specific initial data are unchanged.

From the repository root:
 python -m scripts.CompareT3GeneralScalarSharedA2 --campaign bulk48_v1 --check
 python -m scripts.CompareT3GeneralScalarSharedA2 --campaign bulk48_v1 --run-id shared_a2_v1
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'output/studies/general_scalar_extended/bulk_48'
OUT=ROOT/'output/studies/general_scalar_extended/six_benchmarks'
CASES=('small_uv_frozen_5e-9','large_uv_frozen_0p5','small_uv_running_5e-9')
MODELS=('A_+0','A_-2','A_-4','B_+1','B_-1','B_-3',
        'C_+1','C_-1','C_-3','D_+0','D_+2','D_-2',
        'E_+0','E_+2','E_-2','E_-4')
OBS=('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23')
ACCEPTED={'CandidateConstrainedAccepted','ExistingConditionalAccepted'}
FORMAT='t3-six-conditional-general-scalar-benchmarks-v1'


def require(ok:bool,why:str)->None:
    if not ok:raise ValueError('General-scalar six benchmarks: '+why)

def read(path:Path)->dict:
    require(path.is_file(),'missing input '+str(path))
    d=json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(d,dict),'bad JSON object '+str(path))
    return d

def sha(path:Path)->str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def safe(model:str)->str:
    require(model in MODELS,'unrecognized model '+str(model))
    return model.replace('+','p').replace('-','m')

def atomic(path:Path,value:Any):
    path.parent.mkdir(parents=True,exist_ok=True)
    p=path.with_name(path.name+'.tmp')
    p.write_text(json.dumps(value,indent=2,allow_nan=False,ensure_ascii=False)+'\n',encoding='utf-8')
    p.replace(path)


def campaign_rows(source:Path)->dict[tuple[str,str],dict]:
    raw=read(source/'campaign.json')
    require(raw.get('format')=='t3-48-model-general-quartic-conditional-campaign-v1',
            'unsupported bulk 48 source')
    require(raw.get('conditional_only') is True and
            raw.get('no_sequential_wilson_matching') is True,
            'bulk campaign physics flags changed')
    rows=raw.get('rows')
    require(isinstance(rows,list) and len(rows)==48,'expected exactly 48 saved study cells')
    mapping={}
    for row in rows:
        key=(row.get('study'),row.get('model'))
        require(key[0] in CASES and key[1] in MODELS and key not in mapping,
                'duplicate/unrecognized case: '+str(key))
        mapping[key]=row
    require(set(mapping)=={(c,m) for c in CASES for m in MODELS},
            'incomplete six-benchmark source')
    return mapping


def donor(source:Path,case:str,*,allow_poor:bool)->dict:
    require(case in CASES,'unknown case')
    fitpath=source/case/'A_m2'/'fit.json'
    if not fitpath.is_file():
        return {'status':'MissingDonorFit','reason':str(fitpath)}
    fit=read(fitpath)
    best=fit.get('best_mass_safe')
    if not isinstance(best,dict):
        return {'status':'NoMassSafeA2Donor','reason':fit.get('status')}
    status=fit.get('status')
    if status!='CandidateConstrainedAccepted' and not allow_poor:
        return {'status':'A2DonorNotAccepted','reason':status,'chi2':best.get('chi2')}
    vals=best.get('parameters')
    require(isinstance(vals,dict),'A(-2) donor best Yukawa vector absent')
    # Existing RunT3GeneralScalar48Campaign.NAMES is exactly 18 canonical coords.
    from scripts import RunT3GeneralScalar48Campaign as bulk
    from numpy import asarray, isfinite
    v=asarray([float(vals[n]) for n in bulk.NAMES],dtype=float)
    require(v.shape==(18,) and isfinite(v).all(),'nonfinite A(-2) Yukawas')
    cfgpath=source/case/'A_m2'/'fitted_config.json'
    require(cfgpath.is_file(),'A(-2) donor fitted_config.json not saved')
    cfg=read(cfgpath)
    from Numerical.orchestration import RunSmallLambdaRefit as reference
    ordinary=cfg['base_state']['ordinary']
    stored=asarray([reference.matr(ordinary[k]).real.ravel() for k in ('y1','y2')]).ravel()
    require(stored.shape==(18,) and max(abs(stored-v))<1e-9,
            'donor best values differ from saved fitted configuration')
    return {'status':'AcceptedDonor' if status=='CandidateConstrainedAccepted' else 'PoorFitDonor',
            'chi2':float(best['chi2']),'parameters':{k:float(x) for k,x in zip(bulk.NAMES,v)},
            'source_fit':fitpath.relative_to(ROOT).as_posix(),'source_fit_sha256':sha(fitpath),
            'source_fitted_config':cfgpath.relative_to(ROOT).as_posix(),
            'source_fitted_config_sha256':sha(cfgpath)}


def evaluate_bulk(case:str,model:str,point:dict,source:Path,*,extra_real:float)->dict:
    """Run the ORIGINAL 48-campaign evaluator, not the reduced Sobol evaluator."""
    import numpy as np
    from scripts import RunT3GeneralScalar48Campaign as bulk
    from Numerical.orchestration import RunSmallLambdaRefit as reference
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    registry=bulk.read(bulk.REGISTRY)
    entry,new,reused=bulk.inventory(registry,model)
    require(not reused,'bulk 48 shared comparator does not support reused special quartics')
    history=bulk.historical(case,model,entry)
    beta_path=source/'betas'/(bulk.safe(model)+'.json')
    payload=bulk.read(beta_path)
    kernel,matching_path,kind=bulk.matching(model)
    target=OscillationFitTarget.from_json(reference.TARGET)
    evaluator=bulk.Evaluator(case,model,entry,new,history,payload,kernel,target,extra_real)
    x=np.array([point[k] for k in bulk.NAMES])
    result=evaluator.evaluate(x,'shared_a_m2')
    return {'status':'CandidateSharedNeutrinoEvaluated' if result.get('chi2') is not None
            else result.get('status','SharedEvaluationRejected'),
            'chi2':result.get('chi2'),'prediction':result.get('prediction'),
            'physics_source':'bulk_48_native_extra_quartic_candidate',
            'source_matching':matching_path.relative_to(ROOT).as_posix(),
            'source_matching_sha256':sha(matching_path),
            'source_beta_sha256':sha(beta_path),
            'source_c5_kind':kind,
            'screen':result.get('screen'),
            'reason':result.get('reason'),
            'rhoF':result.get('rhoF'),
            'scalar_masses_muS_tev':result.get('scalar_masses_muS_tev'),
            'max_yukawa_singular_muF':result.get('max_yukawa_singular_muF'),
            'max_scalar_quartic_abs':result.get('max_scalar_quartic_abs')}


def evaluate_bp1(case:str,point:dict,*,extra_real:float,beta_id:str)->dict:
    """Original B(+1) Batch21 evaluator, preserving its mass/perturbativity screens."""
    import numpy as np
    from scripts import FitT3BPlus1GeneralThreeStudies as b21
    from scripts import RunT3BPlus1GeneralEFT1 as b20
    from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from Numerical.orchestration import RunSmallLambdaRefit as reference
    historical=b21.load_historical(case,ROOT)
    registry=b21.read(ROOT/'RGE/general/generated/general_rgbeta_quartic_registry.json')
    b20.validate_registry(registry)
    config=b21.candidate_config(historical,historical['seed'])
    uv,_=b21.prepare_uv(ROOT,config,uv_on=False)
    initial=build_intermediate_scalar_boundary(uv,matching_scale_gev=b20.MU_F)
    payload,layout,beta_path,beta_report=b21.verify_beta_source(ROOT,beta_id,initial)
    kernel=b21.MatchedKernel(ROOT/b21.C5_REL)
    target=OscillationFitTarget.from_json(reference.TARGET)
    uv_path=ROOT/b21.UV_BETAS if historical['uv_on'] else None
    evaluator=b21.Evaluator(historical,payload,uv_path,kernel,target,root=ROOT,
          extra_fraction=None,extra_real=extra_real,rho_max=b21.RHO_MAX,
          scalar_min_tev=b21.SCALAR_MIN_TEV,scalar_max_tev=b21.SCALAR_MAX_TEV)
    row=evaluator.evaluate(np.asarray([point[k] for k in b21.NAMES]),label='shared_a_m2')
    return {'status':'CandidateSharedNeutrinoEvaluated' if row.get('chi2') is not None
            else row.get('status','SharedEvaluationRejected'),
            'chi2':row.get('chi2'),'prediction':row.get('prediction'),
            'physics_source':'BPlus1_Batch21_candidate',
            'source_beta_sha256':sha(beta_path),'source_matching_sha256':sha(ROOT/b21.C5_REL),
            'screen':row.get('screen_pass'),'reason':row.get('reason'),
            'rhoF':row.get('rhoF')}


def evaluate_bc(case:str,model:str,point:dict,*,quartic_fraction:float,mu2_fraction:float)->dict:
    """Existing B/C mixed-scalar candidate with its own conditional matching."""
    from scripts import RunT3GeneralScalarBCExploratoryFit as bc
    from scripts import RunT3GeneralScalarBCPilot as pilot
    from scripts import RunGeneralT3NumericalC5 as c5src
    from Numerical.experimental.GeneralT3ScalarRGState import load_registry
    from Numerical.orchestration import RunSmallLambdaRefit as reference
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    # This case's historical source is model-dependent and read-only.
    source=pilot.study_source(case,model,load_registry(bc.b7.REGISTRY),root=ROOT)
    registry=load_registry(bc.b7.REGISTRY)
    bundle=bc.beta_bundle(model,source,registry)
    name={'B_-1':'B_m1_final','C_-1':'C_m1_final'}[model]
    coefficient=c5src.DEFAULT_SOURCE/(name+'_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')
    compiled=bc.compile_source(coefficient)
    target=OscillationFitTarget.from_json(reference.TARGET)
    ev=bc.CandidateEvaluator(case,model,source,registry,bundle,compiled,target,
          quartic_fraction=quartic_fraction,mu2_fraction=mu2_fraction)
    fit,details=ev.evaluate(point,detail=True)
    prediction={k:float(v) for k,v in zip(fit.observable_names,fit.prediction)}
    return {'status':'CandidateSharedNeutrinoEvaluated',
            'chi2':float(fit.chi2),'prediction':prediction,
            'physics_source':'BC_mixed_scalar_one_scale_candidate',
            'source_matching_sha256':sha(coefficient),
            'rhoF':None,'screen':None,
            'reason':'B/C special one-scale evaluation; coordinate screens require independent validation'}


def run_shared(args, source:Path, campaign:dict, out:Path)->list[dict]:
    donors={c:donor(source,c,allow_poor=args.allow_poor_donor) for c in CASES}
    old={}
    results=[]
    saved_file=out/'shared_results.json'
    if args.resume and saved_file.is_file():
        previous=read(saved_file)
        require(previous.get('source_sha256')==sha(source/'campaign.json') and
                previous.get('donors')==donors and
                previous.get('settings')==settings(args),
                'existing checkpoint differs from current settings/inputs')
        old={(r['case'],r['model']):r for r in previous['rows']}
    atomic(out/'a_m2_donors.json',donors)
    for case in CASES:
        for model in MODELS:
            ident=(case,model)
            if ident in old:
                results.append(old[ident]);continue
            d=donors[case]
            result={'case':case,'model':model,'yukawa_mode':'shared_a_m2',
                    'donor_status':d['status'],'donor_chi2':d.get('chi2'),
                    'physics_sequential_matching_certified':False,
                    'scalar_quadratic_terms_fully_included':False}
            if d['status'] not in ('AcceptedDonor','PoorFitDonor'):
                result.update(status='UnavailableDonor',reason=d.get('reason'))
            elif args.dry_run:
                result.update(status='ReadyNotEvaluated')
            else:
                try:
                    # The original physics runner switches by *actual model-specific* scheme.
                    if model in ('B_-1','C_-1'):
                        x=evaluate_bc(case,model,d['parameters'],
                          quartic_fraction=args.bc_quartic_fraction,
                          mu2_fraction=args.bc_mu2_fraction)
                    elif model=='B_+1':
                        x=evaluate_bp1(case,d['parameters'],extra_real=args.extra_real,
                          beta_id=args.bp1_beta_id)
                    else:
                        x=evaluate_bulk(case,model,d['parameters'],source,
                          extra_real=args.extra_real)
                    result.update(x)
                except Exception as error:
                    result.update(status='SharedEvaluationFailed',reason=str(error)[:600])
            results.append(result)
            atomic(saved_file,{'format':FORMAT,'source_sha256':sha(source/'campaign.json'),
                      'donors':donors,'settings':settings(args),'rows':results})
            print(f'{case} {model}: {result["status"]} {result.get("chi2")}',flush=True)
    return results


def settings(args)->dict:
    return {'allow_poor_donor':args.allow_poor_donor,'extra_real':args.extra_real,
            'bc_quartic_fraction':args.bc_quartic_fraction,
            'bc_mu2_fraction':args.bc_mu2_fraction,'bp1_beta_id':args.bp1_beta_id,
            'dry_run':args.dry_run}


def summarize(campaign:dict,shared:list[dict])->list[dict]:
    rows=[]
    for case in CASES:
        for model in MODELS:
            a=campaign[(case,model)]
            b=next(z for z in shared if z['case']==case and z['model']==model)
            for mode,item in [('optimised',a),('shared_a_m2',b)]:
                pred=item.get('prediction') or (item if mode=='optimised' else {})
                rows.append({'case':case,'model':model,'yukawa_mode':mode,
                   'status':item.get('status'),'chi2':item.get('chi2'),
                   'donor_status':b.get('donor_status') if mode=='shared_a_m2' else None,
                   'donor_chi2':b.get('donor_chi2') if mode=='shared_a_m2' else None,
                   **{k:pred.get(k) for k in OBS},
                   'conditional_matching_only':True,
                   'general_scalar_quadratic_extension_incomplete':True})
    require(len(rows)==96,'expected exactly 96 model-scenario-mode rows')
    return rows


def main(argv=None)->int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--campaign',default='bulk48_v1')
    p.add_argument('--run-id',default='shared_a2_v1')
    p.add_argument('--check',action='store_true')
    p.add_argument('--dry-run',action='store_true')
    p.add_argument('--resume',action='store_true')
    p.add_argument('--allow-poor-donor',action='store_true')
    p.add_argument('--extra-real',type=float,default=-.01)
    p.add_argument('--bc-quartic-fraction',type=float,default=.1)
    p.add_argument('--bc-mu2-fraction',type=float,default=0.)
    p.add_argument('--bp1-beta-id',default='b_p1_eft1_general_v1')
    args=p.parse_args(argv)
    require(re.fullmatch('[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.run_id) is not None,'bad run id')
    require(0<abs(args.extra_real)<=.5 and 0<=args.bc_quartic_fraction<=1 and
            0<=args.bc_mu2_fraction<=.2,'bad numerical settings')
    source=BASE/args.campaign
    rows=campaign_rows(source)
    if args.check:
        print('VALIDATED: 48 existing candidate Sobol scenario cells.')
        print('Campaign status counts:',{s:sum(r.get('status')==s for r in rows.values()) for s in sorted(set(r.get('status') for r in rows.values()))})
        for case in CASES:
            d=donor(source,case,allow_poor=args.allow_poor_donor)
            print(f'A(-2) donor {case}: {d["status"]}, chi2={d.get("chi2")}')
        return 0
    dest=OUT/args.run_id
    require(not dest.exists() or args.resume,
            'Output exists: choose a new --run-id or use --resume')
    results=run_shared(args,source,rows,dest)
    table=summarize(rows,results)
    atomic(dest/'comparison.json',{'format':FORMAT,'source_campaign_sha256':sha(source/'campaign.json'),
        'study_count':3,'model_count':16,'modes':2,'rows':table,
        'limitations':['All 96 rows are conditional exploratory results or documented failures; not verified general-scalar matching.',
        'Reduced/no-mixing C5 for six self-bilinear models; B/C use conditional one-scale matching.',
        'Missing/poor A(-2) donors are not replaced by historical reduced-potential Yukawas.',
        'Full general-scalar EFT Wilson operator matching and running not certified.']})
    csvpath=dest/'comparison.csv'
    with csvpath.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(table[0]))
        writer.writeheader();writer.writerows(table)
    n=sum(x['status']=='CandidateSharedNeutrinoEvaluated' for x in results)
    print('Wrote '+str(csvpath)+f'; {n}/48 shared cases evaluated; 48 saved original Sobol cases included')
    return 0

if __name__=='__main__':raise SystemExit(main())
