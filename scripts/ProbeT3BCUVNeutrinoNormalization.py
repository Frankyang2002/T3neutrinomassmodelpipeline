"""Batch 16: B/C UV-running mass-safe candidate C5/chi2 normalization probe.

Numerical PROBES only, not oscillation refits. Rescales saved fitted or historical
Yukawa textures by small factors, and lambdaT3 = lambda_historical/f**2*k.
Re-runs candidate native RGBeta UV + EFT1 with the revised inputs, checks
fermion-first running-mass conditions, then evaluates the existing *conditional*
one-scale mixed-mass Matchete C5 and 100 GeV SM+Weinberg observables.

No physical sequential C5 or Matchete/RGBeta map is certified. Independent
immutable output directory. Original historical and general-scalar results are
never written.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
import csv
import json
import math
from pathlib import Path

import numpy as np

from scripts import RunT3GeneralScalarBCExploratoryFit as b11
from scripts import RunT3GeneralScalarBCPilot as b9
from scripts import RunT3ZeroExtraBoundaryBenchmark as b7
from scripts import DiagnoseT3BCUVScalarHierarchy as b14
from scripts import StudyT3BCUVScalarMassFeasibility as b15

ROOT = b11.ROOT
OUT = ROOT / 'output/studies/general_scalar_extended/diagnostics'
FORMAT = 't3-bc-uv-neutrino-normalization-probe-v1'
MODELS = ('B_-1', 'C_-1')
FACTORS = {'B_-1': (0.005, 0.01, 0.015, 0.02),
           'C_-1': (0.005, 0.01, 0.015)}
SOURCES = ('saved_poor_fit', 'historical_fit')
CSV_COLUMNS = ('model','yukawa_source','yukawa_factor','lambda_multiplier',
    'lambdaT3_uv','rhoF','scalar_light_muS_tev','scalar_heavy_muS_tev',
    'status','chi2','delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12',
    'sin2_theta13','sin2_theta23','suggested_lambda_multiplier_linear_C5_only','reason')


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError('T3 BC normalization probe: '+message)


def parse_positive_csv(raw: str, label: str, *, maximum: float = 100.) -> tuple[float,...]:
    try:
        result=tuple(float(v.strip()) for v in raw.split(','))
    except (AttributeError,ValueError) as exc:
        raise ValueError(label+': expected comma-separated positive numbers') from exc
    require(1<=len(result)<=16 and len(set(result))==len(result)
            and all(math.isfinite(v) and 0<v<=maximum for v in result),
            label+': expected 1..16 unique finite positive factors')
    return result


def scaled_lambda_source(source: dict, *, yukawa_factor: float,
                         lambda_multiplier: float) -> tuple[dict,float]:
    """Update BOTH duplicated topology-coupling representations in source."""
    require(math.isfinite(yukawa_factor) and 0<yukawa_factor<=1.,
            'invalid Yukawa scaling factor')
    require(math.isfinite(lambda_multiplier) and 0<lambda_multiplier<=100.,
            'invalid lambda multiplier')
    original=complex(source['lambdaT3'])
    require(original.imag==0 and original.real>0,'expected real positive topology coupling')
    value=float(original.real/yukawa_factor**2*lambda_multiplier)
    require(math.isfinite(value) and 0<value<=1.,
            'new lambdaT3 not finite or exceeds exploratory bound 1')
    out=deepcopy(source)
    out['lambdaT3']=complex(value)
    out['cfg']['base_state']['ordinary']['lambdaT3']=b7.complex_json(complex(value))
    return out,value


def seeded_yukawa(saved: dict, historical: dict, *, source_name: str,
                  factor: float) -> np.ndarray:
    require(source_name in SOURCES, 'invalid Yukawa seed source')
    require(math.isfinite(factor) and 0<factor<=1.,'invalid seed Yukawa factor')
    if source_name=='saved_poor_fit':
        params=saved['best_parameters']
        require(set(params)==set(b11.NAMES),'saved fit Yukawa inventory mismatch')
        seed=np.asarray([params[k] for k in b11.NAMES],dtype=float)
    else:
        seed=b11.seed_vector(historical)
    require(seed.shape==(18,) and np.isfinite(seed).all(), 'invalid seed vector')
    return seed*factor


def check_mass_screen(uv_spectrum:dict, eft_spectrum:dict|None,
                      *,ratio_max:float,scalar_min_tev:float,
                      scalar_max_tev:float)->dict:
    return b15.classify({'final':uv_spectrum},
             None if eft_spectrum is None else {'final':eft_spectrum},
             ratio_max=ratio_max,scalar_min_tev=scalar_min_tev,
             scalar_max_tev=scalar_max_tev)


def run_one(model: str,record:dict,registry:dict,bundle:dict,compiled,target,
            *,quartic_fraction:float,mu2_fraction:float,
            ratio_max:float,scalar_min_tev:float,scalar_max_tev:float)->dict:
    """Full UV/EFT1 candidate replay; NO C5 for an invalid mass screen."""
    from Numerical.experimental.GeneralT3ScalarRGState import (
        unpack_general_scalar_state,project_general_after_fermion,pack_general_scalar_state)
    from scripts import AnalyzeT3BCMatchingConventions as b10
    y0,layout,_,_=b9.prepare_initial(record,registry,
        'general_quartics_and_mixing' if mu2_fraction>0 else 'general_quartics',
        quartic_fraction,mu2_fraction)
    uv=b7.integrate_segment(y0,
        lambda t,y:b7.rhs_at(t,y,layout=layout,registry=registry,
                              betas=bundle['UV'],stage='UV'),
        1e7,1e5,n_points=3)
    yf=np.asarray(uv.y[:,-1],dtype=float)
    uv_state=unpack_general_scalar_state(yf[:layout.size],layout,registry,
                                         uv=True,mu_gev=1e5)
    us=b15.scalar_spectrum(yf,layout,registry,scale=1e5,stage='UV')
    check=check_mass_screen(us,None,ratio_max=ratio_max,
        scalar_min_tev=scalar_min_tev,scalar_max_tev=scalar_max_tev)
    row={'rhoF':check['ratio_muF'],'uv_nfev':int(uv.nfev),
         'scalar_light_muF_tev':(us['scalar_masses_gev'][0]/1000
                   if us['scalar_masses_gev'][0] is not None else None),
         'scalar_heavy_muF_tev':(us['scalar_masses_gev'][1]/1000
                   if us['scalar_masses_gev'][1] is not None else None),
         'min_fermion_muF_tev':us['min_fermion_mass_gev']/1000,
         'quadratic_mu2_at_muF_gev2':us['mu2_abs_gev2'],
         'screen_pass':False,'chi2':None,'prediction':None,
         'match_performed':False,'candidate_one_scale_only':True,
         'physical_map_certified':False}
    if not check['strict_hierarchy_at_muF']:
        row.update(status='MassScreenFailAtMuF',reason=check['reason'],
                   **{'screen':check})
        return row
    after=project_general_after_fermion(uv_state,registry)
    ye,el=pack_general_scalar_state(after,registry)
    ye=np.concatenate((ye,yf[layout.size:]))
    new,reused=record['new'],record['reused']
    betas=b7.check_payload(bundle['payloads']['EFT1'],model=model,stage='EFT1',
        rep=el.representation,base_names=tuple(el.names[:-len(new)]),
        new_names=new,reused_names=reused)
    eft=b7.integrate_segment(ye,
        lambda t,y:b7.rhs_at(t,y,layout=el,registry=registry,
                              betas=betas,stage='EFT1'),
        1e5,1e3,n_points=3)
    ys=np.asarray(eft.y[:,-1],dtype=float)
    es=b15.scalar_spectrum(ys,el,registry,scale=1e3,stage='EFT1',
                           heavy_at_F=tuple(np.linalg.svd(
                               uv_state.baseline.MF,compute_uv=False)))
    check=check_mass_screen(us,es,ratio_max=ratio_max,
        scalar_min_tev=scalar_min_tev,scalar_max_tev=scalar_max_tev)
    row.update({'eft1_nfev':int(eft.nfev), 'screen':check,
         'screen_pass':check['screen_pass'],
         'scalar_light_muS_tev':(es['scalar_masses_gev'][0]/1000
                   if es['scalar_masses_gev'][0] is not None else None),
         'scalar_heavy_muS_tev':(es['scalar_masses_gev'][1]/1000
                   if es['scalar_masses_gev'][1] is not None else None)})
    if not check['screen_pass']:
        row.update(status='MassScreenFailAtMuS',reason=check['reason'])
        return row
    endpoint=b7.snapshot(ys,el,registry,stage='EFT1',mu_gev=1e3,
                         new=new,reused=reused)
    eft_state=unpack_general_scalar_state(ys[:el.size],el,registry,
                                          uv=False,mu_gev=1e3).baseline
    heavy=b11.heavy_flavour_at_muF(uv_state.baseline)
    c5=b11.c5_at_endpoint(model,endpoint,heavy,compiled)
    fit,low=b11.low_energy_fit(c5,eft_state,target,
        vev=float(record['cfg'].get('numerical',{}).get('vev_gev',246.22)))
    prediction=dict(zip(fit.observable_names,map(float,fit.prediction)))
    require('delta_m3l_sq_ev2' in prediction,'missing atmospheric splitting')
    row.update(status='CandidateOscillationEvaluated',reason='ScreenPassOnly',
        chi2=float(fit.chi2),prediction=prediction,match_performed=True,
        candidate_C5_muS_gev_inv=b10.matrix_json(c5),
        candidate_C5_muLow_gev_inv=b10.matrix_json(low.K),
        heavy_mass_basis_residual=float(heavy[3]),
        absolute_neutrino_fit_pass=bool(float(fit.chi2)<=b9.FIT_CUTOFF))
    return row


def estimate_lambda_multiplier(row:dict,*, target_dm3l:float=0.002511)->float|None:
    """Linear-C5 scaling estimate, NOT actual tuned solution or RGE prediction."""
    pred=row.get('prediction')
    if not isinstance(pred,dict):return None
    dm=float(pred.get('delta_m3l_sq_ev2',float('nan')))
    if not math.isfinite(dm) or dm<=0:return None
    return math.sqrt(target_dm3l/dm)


def destination(root:Path,run_id:str)->Path:
    b9.parse_run_id(run_id)
    dest=root/'output/studies/general_scalar_extended/diagnostics'/run_id
    require(not dest.exists(),'run ID already exists; choose a new run ID')
    return dest


def summary_markdown(report:dict)->str:
    lines=['# Batch 16 — hierarchy-safe neutrino normalization probes','',
       '**Not a fitted or certified full-general-scalar neutrino calculation.**',
       'Candidate native RGBeta UV/EFT1, existing conditional one-scale Matchete C5,',
       'and SM+Weinberg running to 100 GeV. No UV threshold matching is changed.','',
       f"Mass screen: rhoF <= {report['screen']['rho_max']:.3g}; "
       f"scalar masses at 1 TeV in [{report['screen']['scalar_min_tev']:.3g}, "
       f"{report['screen']['scalar_max_tev']:.3g}] TeV.",'',
       '| Model | Yukawa seed | Y factor | Lambda multiplier | Lambda UV | rhoF | Scalars at 1 TeV (TeV) | chi2 | Status |',
       '|---|---|---:|---:|---:|---:|---|---:|---|']
    for model in MODELS:
        if model not in report['models']:continue
        for r in report['models'][model]['trials']:
            mass=(f"{r['scalar_light_muS_tev']:.3g}–{r['scalar_heavy_muS_tev']:.3g}"
                  if r.get('scalar_light_muS_tev') is not None else '—')
            chi=(f"{r['chi2']:.7g}" if r['chi2'] is not None else '—')
            rho=(f"{r['rhoF']:.5g}" if r.get('rhoF') is not None else '—')
            lines.append(f"| {model} | {r['yukawa_source']} | {r['yukawa_factor']:.4g} | "
                         f"{r['lambda_multiplier']:.4g} | {r['lambdaT3_uv']:.4g} | "
                         f"{rho} | {mass} | {chi} | {r['status']} |")
    lines+=['','## Interpretation','',
       '- The 1/f² lambda scaling is a starting hypothesis, NOT an RG or C5 identity.',
       '- Any chi2 <=5.86 is a seed-level candidate fit conditional on the existing one-scale C5.',
       '- No Yukawa or lambda optimisation is performed; numeric trials are distinct from fits.',
       '- The beta functions and physical threshold ordering have not been newly certified.',
       '- Inputs are source-hash-checked and historical outputs remain untouched.','']
    for model in MODELS:
        if model not in report['models']:continue
        best=report['models'][model]['best']
        if best:
            lines+=[f"### {model}: lowest computed chi2",'',
              f"Yukawa source `{best['yukawa_source']}`, factor {best['yukawa_factor']:g}, "
              f"lambda multiplier {best['lambda_multiplier']:g}, chi2 = {best['chi2']:.8g}.",
              f"rhoF = {best['rhoF']:.5g}, scalar masses = "
              f"{best['scalar_light_muS_tev']:.4g}–{best['scalar_heavy_muS_tev']:.4g} TeV.",'']
    return '\n'.join(lines)+'\n'


def execute(*,run_id:str,models:tuple[str,...]=MODELS,
            yukawa_factors:dict[str,tuple[float,...]]|None=None,
            lambda_multipliers:tuple[float,...]=(0.1,1.,10.),
            seed_sources:tuple[str,...]=SOURCES,
            ratio_max:float=0.5,scalar_min_tev:float=0.5,
            scalar_max_tev:float=5.,root:Path=ROOT,
            run_ids:dict[str,str]|None=None)->dict:
    require(models and len(set(models))==len(models) and set(models)<=set(MODELS),
            'models must be unique B/C')
    require(seed_sources and len(set(seed_sources))==len(seed_sources)
            and set(seed_sources)<=set(SOURCES), 'invalid Yukawa seed sources')
    require(0<ratio_max<1 and 0<scalar_min_tev<scalar_max_tev,
            'invalid hierarchy/mass screen')
    require(1<=len(lambda_multipliers)<=10 and len(set(lambda_multipliers))==len(lambda_multipliers)
            and all(math.isfinite(v) and 0<v<=100 for v in lambda_multipliers),
            'invalid lambda multipliers')
    yf=FACTORS if yukawa_factors is None else yukawa_factors
    dest=destination(root,run_id)
    from Numerical.experimental.GeneralT3ScalarRGState import load_registry
    from Numerical.orchestration import RunSmallLambdaRefit as base
    from Numerical.fitting.OscillationFit import OscillationFitTarget
    from scripts.RunGeneralT3NumericalC5 import DEFAULT_SOURCE
    from scripts import AnalyzeT3BCMatchingConventions as b10

    registry=load_registry(root/'RGE/general/generated/general_rgbeta_quartic_registry.json')
    target=OscillationFitTarget.from_json(base.TARGET)
    result={'format':FORMAT,'run_id':run_id,'status':'CandidateOneScaleNormalizationProbes',
       'input_case':b14.CASE,'models':{},'screen':{'rho_max':ratio_max,
                'scalar_min_tev':scalar_min_tev,'scalar_max_tev':scalar_max_tev},
       'lambda_scaling_formula':'lambdaT3 = 5e-9 / (Y_factor)^2 * lambda_multiplier',
       'physical_matching_certified':False,'sequential_C5_matched':False,
       'optimisation_performed':False,'historical_outputs_overwritten':False}
    csvrows=[]
    for model in models:
        factors=tuple(yf[model]); require(factors and len(set(factors))==len(factors)
                       and all(0<x<=0.2 and math.isfinite(x) for x in factors),
                       'invalid model Yukawa factors')
        saved,historical,settings=b14.fit_source(root,model,
             (run_ids or b14.DEFAULT_RUNS)[model],registry)
        require(settings['quartic_fraction']==0.1 and settings['mu2_fraction']==0,
                'Batch16 expects Batch15 general-quartic=0.1, mu2=0 baseline')
        # Source and beta expressions are checked by the established Batch11 loader.
        bundle=b11.beta_bundle(model,historical,registry)
        coeff=Path(DEFAULT_SOURCE)/({'B_-1':'B_m1_final',
                  'C_-1':'C_m1_final'}[model]+ '_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')
        compiled=b11.compile_source(coeff)
        b10.check_convention_evidence(b7.read_json(b10.EVIDENCE23),
            b7.read_json(b10.EVIDENCE26),models=(model,))
        provenance={'saved_fit':str(settings['saved_fit_path'].relative_to(root)),
                    'saved_fit_sha256':sha256(settings['saved_fit_path'].read_bytes()).hexdigest(),
                    'historical_config_sha256':sha256(historical['cfg_path'].read_bytes()).hexdigest(),
                    'ordered_C5_sha256':sha256(coeff.read_bytes()).hexdigest(),
                    'target_sha256':sha256(Path(base.TARGET).read_bytes()).hexdigest(),
                    'registry_sha256':sha256((root/'RGE/general/generated/general_rgbeta_quartic_registry.json').read_bytes()).hexdigest(),
                    'convention_stage23_sha256':sha256(Path(b10.EVIDENCE23).read_bytes()).hexdigest(),
                    'convention_stage26_sha256':sha256(Path(b10.EVIDENCE26).read_bytes()).hexdigest(),
                    'uv_beta_sha256':sha256(bundle['files']['UV'].read_bytes()).hexdigest(),
                    'eft1_beta_sha256':sha256(bundle['files']['EFT1'].read_bytes()).hexdigest()}
        trials=[]
        for source_name in seed_sources:
            for f in factors:
                for k in lambda_multipliers:
                    modified,lam=scaled_lambda_source(historical,
                         yukawa_factor=f,lambda_multiplier=k)
                    seed=seeded_yukawa(saved,historical,
                          source_name=source_name,factor=f)
                    rec=b11.trial_record(modified,seed)
                    row={'model':model,'yukawa_source':source_name,
                         'yukawa_factor':float(f),'lambda_multiplier':float(k),
                         'lambdaT3_uv':float(lam),'run_id':run_id}
                    try:
                        got=run_one(model,rec,registry,bundle,compiled,target,
                             quartic_fraction=settings['quartic_fraction'],
                             mu2_fraction=settings['mu2_fraction'],ratio_max=ratio_max,
                             scalar_min_tev=scalar_min_tev,scalar_max_tev=scalar_max_tev)
                        row.update(got)
                        row['suggested_lambda_multiplier_linear_C5_only']=estimate_lambda_multiplier(row)
                    except Exception as exc:
                        row.update(status='NumericalProbeFailed',reason=str(exc)[:400],
                                   screen_pass=False,chi2=None,prediction=None,
                                   suggested_lambda_multiplier_linear_C5_only=None)
                    trials.append(row)
                    csvrow={c:row.get(c) for c in CSV_COLUMNS}
                    if row.get('prediction'):
                        csvrow.update({key:val for key,val in row['prediction'].items()
                                       if key in CSV_COLUMNS})
                    csvrows.append(csvrow)
                    print(f"{model}: {source_name} Y×{f:g} lambda×{k:g}: "
                          f"{row['status']}, rhoF={row.get('rhoF')}, chi2={row.get('chi2')}",flush=True)
        valid=[r for r in trials if r['status']=='CandidateOscillationEvaluated'
               and r['screen_pass'] and r['chi2'] is not None]
        best=min(valid,key=lambda r:r['chi2']) if valid else None
        require(not all(t['status']=='NumericalProbeFailed' for t in trials),
                model+': every numerical probe failed; no results written')
        result['models'][model]={'model':model,'n_trials':len(trials),
                'n_screen_pass_with_chi2':len(valid),
                'best':best,'provenance':provenance,'trials':trials}
    dest.mkdir(parents=True,exist_ok=False)
    b7.save_json(dest/'probe.json',result)
    with (dest/'probe.csv').open('w',newline='',encoding='utf-8') as fp:
        w=csv.DictWriter(fp,fieldnames=CSV_COLUMNS)
        w.writeheader();w.writerows(csvrows)
    (dest/'summary.md').write_text(summary_markdown(result),encoding='utf-8')
    print('Saved:',dest/'summary.md',flush=True)
    return result


def main()->None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id',required=True)
    p.add_argument('--model',choices=MODELS)
    p.add_argument('--yukawa-factors',help='Optional single-model override: CSV positive factors <=0.2')
    p.add_argument('--lambda-multipliers',default='0.1,1,10')
    p.add_argument('--seed-source',choices=('both',)+SOURCES,default='both')
    p.add_argument('--rho-max',type=float,default=0.5)
    p.add_argument('--scalar-min-tev',type=float,default=0.5)
    p.add_argument('--scalar-max-tev',type=float,default=5.)
    p.add_argument('--b-run-id',default=b14.DEFAULT_RUNS['B_-1'])
    p.add_argument('--c-run-id',default=b14.DEFAULT_RUNS['C_-1'])
    a=p.parse_args()
    require(not(a.yukawa_factors and not a.model),
            '--yukawa-factors requires --model')
    factors=(dict(FACTORS,**{a.model:parse_positive_csv(a.yukawa_factors,
             'Y factors',maximum=0.2)}) if a.yukawa_factors else None)
    sources=SOURCES if a.seed_source=='both' else (a.seed_source,)
    execute(run_id=a.run_id,models=(a.model,) if a.model else MODELS,
            yukawa_factors=factors,
            lambda_multipliers=parse_positive_csv(a.lambda_multipliers,
                       'lambda multipliers'),seed_sources=sources,
            ratio_max=a.rho_max,scalar_min_tev=a.scalar_min_tev,
            scalar_max_tev=a.scalar_max_tev,
            run_ids={'B_-1':a.b_run_id,'C_-1':a.c_run_id})

if __name__=='__main__':main()
