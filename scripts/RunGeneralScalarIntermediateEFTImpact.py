"""General-scalar, conditional benchmark: intermediate EFT versus SMEFT K running.

Uses the exact bulk48 model-specific UV/EFT1 beta sources and C5 kernels.
Special B(-1), C(-1), B(+1) evaluators are deliberately not substituted.
No new Yukawa fitting; matching approximations are those of the source campaign.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from scripts import RunT3GeneralScalar48Campaign as bulk
from scripts import CompareT3GeneralScalarSharedA2 as shared
from scripts.RunT3GeneralScalarBCExploratoryFit import low_energy_fit
from Numerical.fitting.OscillationFit import OscillationFitTarget, evaluate_mass_matrix_fit
from Numerical.fitting.ParameterScan import neutrino_mass_from_c5, charged_lepton_mass_basis_matrix
from Numerical.orchestration import RunSmallLambdaRefit as refit
from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
from Numerical.running.IntermediateScalarStateVector import pack_intermediate_scalar_state, unpack_intermediate_scalar_state
from Numerical.running.ScalarThresholdBoundary import scalar_threshold_masses

OBS = ('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23')
CASE_NAMES = ('00','10','01','11')
SPECIAL = {'B_-1','C_-1','B_+1'}


def finite_fit(fit):
    return {'chi2': float(fit.chi2), 'prediction': {k: float(v) for k,v in zip(fit.observable_names,fit.prediction)}}


def frozen_eft(uv, qf):
    init = build_intermediate_scalar_boundary(uv, matching_scale_gev=1e5)
    vector, layout = pack_intermediate_scalar_state(init)
    end = unpack_intermediate_scalar_state(vector, layout, mu_gev=1e3)
    masses = np.asarray(scalar_threshold_masses(end).masses_gev, dtype=float)
    return end, masses, qf, 0


def predictions(c5, eft, target, vev):
    # Real running of SM couplings / Ye for a common low-energy charged-lepton basis.
    fit, low = low_energy_fit(c5, eft, target, vev=vev)
    run = finite_fit(fit)
    m = neutrino_mass_from_c5(np.asarray(c5, dtype=complex), vev_gev=vev)
    frozen = evaluate_mass_matrix_fit(charged_lepton_mass_basis_matrix(m, low.ye), target)
    return finite_fit(frozen), run


def angle(s):
    if not -1e-9 <= s <= 1+1e-9: raise ValueError('Unphysical sin² angle')
    return math.degrees(math.asin(math.sqrt(max(0.,min(1.,s)))))


def compute(model, source, case, target):
    registry=bulk.read(bulk.REGISTRY)
    entry,new,reused=bulk.inventory(registry,model)
    if reused: raise ValueError('Special/reused quartics require model-specific evaluator')
    fitted=source/case/bulk.safe(model)/'fitted_config.json'
    fitpath=source/case/bulk.safe(model)/'fit.json'
    if not fitted.is_file(): raise FileNotFoundError('Missing individually fitted general-scalar configuration: '+str(fitted))
    saved=bulk.read(fitpath)
    if not isinstance(saved.get('best_mass_safe'),dict): raise ValueError('No saved mass-safe fit')
    cfg=bulk.read(fitted)
    hist=bulk.historical(case,model,entry)
    beta_path=source/'betas'/(bulk.safe(model)+'.json')
    payload=bulk.read(beta_path)
    kernel,match_path,kind=bulk.matching(model)
    extra0=complex(float(saved['extra_quartics_rgbeta_boundary_real']))
    q0=tuple(extra0 for _ in new)
    vec,layout,unpack=bulk.uv_base(cfg,entry)
    if hist['uv_on']:
        uvbeta,_=bulk.uv_betas(model,entry,new,layout)
        uv,qf,_=bulk.run_uv(vec,layout,uvbeta,new,q0)
    else:
        uv=unpack(vec,layout,mu_gev=1e5)
        qf=q0
    running_eft=bulk.run_eft(uv,qf,payload,model,new,())
    frozen=frozen_eft(uv,qf)
    results={}
    screen={}
    for enabled,data in ((False,frozen),(True,running_eft)):
        eft,masses,qs,_=data
        # Record screening but do not eliminate the comparison when a frozen branch
        # violates a fit-era screen; this is an ablation diagnostic.
        screen['running' if enabled else 'frozen']=bulk.testscreen(uv,eft,masses,qf,qs)
        c5=kernel.evaluate(uv,eft)
        vev=float(cfg.get('numerical',{}).get('vev_gev',246.22))
        no_sm,yes_sm=predictions(c5,eft,target,vev)
        results['10' if enabled else '00']=no_sm
        results['11' if enabled else '01']=yes_sm
    effects={}
    for k in ('10','01','11'):
        effects[k]={}
        for obs in OBS:
            base=results['00']['prediction'][obs]
            value=results[k]['prediction'][obs]
            effects[k][obs+'_percent']=100*(value/base-1) if base else None
            if obs.startswith('sin2_theta'):
                effects[k]['theta_'+obs[-2:]+'_shift_deg']=angle(value)-angle(base)
    return {'cases':results,'relative_to_00':effects,'screens':screen,
            'source_fit':str(fitpath.relative_to(bulk.ROOT)),
            'source_config':str(fitted.relative_to(bulk.ROOT)),
            'source_beta':str(beta_path.relative_to(bulk.ROOT)),
            'source_matching':str(match_path.relative_to(bulk.ROOT)),
            'matching_kind':kind,'uv_running':bool(hist['uv_on']),
            'notes':'Conditional general-scalar matching; not certified full sequential Wilson matching'}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--campaign',default='bulk48_v1')
    ap.add_argument('--case',choices=shared.CASES,default='large_uv_frozen_0p5')
    ap.add_argument('--model',choices=shared.MODELS,help='One representation; omitted means all')
    ap.add_argument('--run',action='store_true')
    ap.add_argument('--output-id',default='eft_impact_v1')
    args=ap.parse_args()
    source=bulk.OUT/args.campaign
    if not (source/'campaign.json').is_file():
        ap.error('General-scalar campaign not found: '+str(source))
    selection=[args.model] if args.model else shared.MODELS
    output=bulk.ROOT/'output/studies/general_scalar_extended/intermediate_eft_impact'/args.output_id
    if not args.run:
        print('General-scalar source:',source)
        for name in selection:
            fitfile=source/args.case/bulk.safe(name)/'fitted_config.json'
            print(name, 'SpecialModelUnsupported' if name in SPECIAL else ('Ready' if fitfile.is_file() else 'MissingFit'))
        print('Dry run. Add --run to evaluate.')
        return 0
    target=OscillationFitTarget.from_json(refit.TARGET)
    rows=[]
    for name in selection:
        r={'model':name,'case':args.case}
        try:
            if name in SPECIAL: raise NotImplementedError('Special B(-1), C(-1), B(+1) require separate matching/evaluator adaptation')
            r.update(status='Success',**compute(name,source,args.case,target))
        except Exception as exc:
            r.update(status='Unsupported' if isinstance(exc,NotImplementedError) else 'Failed',error=f'{type(exc).__name__}: {exc}')
        rows.append(r)
        print(name,r['status'],r.get('error',''),flush=True)
    output.mkdir(parents=True,exist_ok=True)
    (output/'results.json').write_text(json.dumps({'case':args.case,'rows':rows},indent=2,allow_nan=False)+'\n',encoding='utf-8')
    with (output/'results.csv').open('w',newline='',encoding='utf-8') as f:
        fields=['model','case','status','error']+[f'{k}_{o}' for k in CASE_NAMES for o in OBS]+[f'{k}_{o}_pct' for k in ('10','01','11') for o in OBS]+[f'{k}_theta_{v}_shift_deg' for k in ('10','01','11') for v in ('12','13','23')]
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for r in rows:
            flat={k:r.get(k,'') for k in ('model','case','status','error')}
            if r['status']=='Success':
                for c in CASE_NAMES:
                    for obs in OBS: flat[f'{c}_{obs}']=r['cases'][c]['prediction'][obs]
                for c in ('10','01','11'):
                    for obs in OBS: flat[f'{c}_{obs}_pct']=r['relative_to_00'][c][obs+'_percent']
                    for v in ('12','13','23'): flat[f'{c}_theta_{v}_shift_deg']=r['relative_to_00'][c]['theta_'+v+'_shift_deg']
            w.writerow(flat)
    print('Saved:',output)
    return 0 if any(r['status']=='Success' for r in rows) else 1

if __name__=='__main__':raise SystemExit(main())
