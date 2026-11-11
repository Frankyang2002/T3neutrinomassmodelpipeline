"""Consolidate existing coupled-IVP T3 calculations into a 16-model rough-results campaign.

Source selection is explicit: latest experimental cached non-HHSS context for A_-4
and D_-2; historical source-preserving coupled IVPs for the other 14 models.
The selected data are *conditional* numerical approximations and must not be
labelled fully renormalised or physically normalised EFT matching.

Writes exclusively into a new output directory. Standard library only.
Project path: scripts/ConsolidateT3Rough16Results.py
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
from datetime import datetime, timezone
from pathlib import Path

SCHEMA = 't3-conditional-coupled-rough16-consolidation-v1'
REGISTRY = 'RGE/general/generated/general_rgbeta_quartic_registry.json'
BASE = Path('output/studies/general_scalar_extended')
HISTORICAL = {
    'A_-4': BASE / 'experimental_simultaneous_eft1_seven_term/simultaneous_A_m4_v2',
    'mixed': BASE / 'experimental_simultaneous_eft1_mixed_BC/simultaneous_mixed_BC_v1',
    'nonmixed': BASE / 'experimental_simultaneous_eft1_general_nonmixed/simultaneous_nonmixed13_v1',
}
EXPERIMENT_DIR = BASE / 'experimental_coupled_cached_scalar_a4d2'
RESULT_DIR = BASE / 'rough16_consolidated'
CASE_CODES = tuple(f'{i:04b}' for i in range(16))
STAGE = ('UV','EFT1','Wilson','SMEFT')
EXPERIMENTAL = {'A_-4', 'D_-2'}
MIXED = {'B_-1','C_-1'}
MODEL_FIELDS = ('model','scalar_context','chi2_1111','m1_eV','m2_eV','m3_eV',
                'dm21_sq_eV2','dm3l_sq_eV2','sin2_theta12','sin2_theta13','sin2_theta23',
                'chi2_0111','chi2_1011','chi2_1101','chi2_1110',
                'c5_norm_at_scalar_threshold_GeVinv','source_cases_sha256')
CASE_FIELDS = ('model','scalar_context','case','UV','EFT1','Wilson','SMEFT',
               'chi2','m1_eV','m2_eV','m3_eV',
               'dm21_sq_eV2','dm3l_sq_eV2','sin2_theta12','sin2_theta13','sin2_theta23',
               'c5_norm_at_scalar_threshold_GeVinv','c5_norm_at_low_scale_GeVinv',
               'coupled_ode_nfev')
DELTA_FIELDS = ('model','scalar_context','stage_frozen','case', 'chi2_full','chi2_frozen',
                'delta_chi2_frozen_minus_full','m3_full_eV','m3_frozen_eV',
                'delta_m3_eV','c5_scalar_threshold_relative_difference')
COMPARE_FIELDS = ('model','case','chi2_prior','chi2_selected','delta_chi2',
                  'm3_prior_eV','m3_selected_eV','delta_m3_eV',
                  'c5_scalar_threshold_relative_difference')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def load(path):
    with path.open(encoding='utf-8-sig') as file:
        return json.load(file, parse_constant=lambda x: (_ for _ in ()).throw(
            ValueError('Nonfinite JSON constant in '+str(path)+': '+x)))


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for part in iter(lambda: f.read(1 << 20), b''):
            h.update(part)
    return h.hexdigest()


def finite_walk(value, location):
    if isinstance(value, float):
        require(math.isfinite(value), 'nonfinite numeric '+location)
    elif isinstance(value, dict):
        for k,v in value.items():
            finite_walk(v, location+'.'+str(k))
    elif isinstance(value, list):
        for i,v in enumerate(value):
            finite_walk(v, location+'['+str(i)+']')


def float_at(obj, key, location, *, nullable=False):
    value = obj.get(key)
    if value is None and nullable:
        return None
    require(isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(float(value)),
            f'{location}.{key} must be finite numeric')
    return float(value)


def decode_c5(obj, location):
    require(isinstance(obj,dict),'invalid C5: '+location)
    a,b=obj.get('real'),obj.get('imag')
    require(isinstance(a,list) and isinstance(b,list) and len(a)==len(b)==3,
            'C5 must have 3x3 real/imag: '+location)
    result=[]
    for i in range(3):
        require(isinstance(a[i],list) and isinstance(b[i],list) and len(a[i])==len(b[i])==3,
                'C5 row dimensions invalid: '+location)
        row=[]
        for j in range(3):
            x,y=a[i][j],b[i][j]
            require(isinstance(x,(int,float)) and isinstance(y,(int,float)) and
                    math.isfinite(float(x)) and math.isfinite(float(y)),
                    'nonfinite C5 entry: '+location)
            row.append(complex(x,y))
        result.append(row)
    norm=math.sqrt(sum(abs(z)**2 for row in result for z in row))
    asym=math.sqrt(sum(abs(result[i][j]-result[j][i])**2
                       for i in range(3) for j in range(3)))
    require(asym <= 1e-8*max(norm,1e-35),
            f'{location}: C5 transpose symmetry failed ({asym/max(norm,1e-35):.3g})')
    return result,norm


def relative_c5(left,right):
    a,na=decode_c5(left,'selected C5')
    b,nb=decode_c5(right,'reference C5')
    err=math.sqrt(sum(abs(a[i][j]-b[i][j])**2
                      for i in range(3) for j in range(3)))
    return err/max(nb,1e-100)


def tag(model):
    return model.replace('-', 'm').replace('+','p')


def historical_path(project,model):
    which='A_-4' if model=='A_-4' else 'mixed' if model in MIXED else 'nonmixed'
    return project/HISTORICAL[which]/tag(model)/'cases.json'


def candidate_path(project,run_id,model):
    return project/EXPERIMENT_DIR/run_id/tag(model)/'cases.json'


def verify_model(path,model):
    require(path.is_file(),f'{model}: cases file does not exist: {path}')
    raw=load(path)
    require(isinstance(raw,dict) and raw.get('model')==model and raw.get('status')=='Computed',
            f'{model}: wrong model or incomplete result at {path}')
    cases=raw.get('cases')
    require(isinstance(cases,dict) and set(cases)==set(CASE_CODES),
            f'{model}: expected exactly 16 U/E/W/S cases; got {sorted(cases) if isinstance(cases,dict) else type(cases)}')
    require(raw.get('full_sequential_general_C5_certified') is not True,
            f'{model}: unexpected claim of fully certified result; review provenance')
    finite_walk(raw,model)
    rows={}
    for code in CASE_CODES:
        x=cases[code]
        require(isinstance(x,dict),f'{model}/{code}: malformed case')
        require(x.get('case',code)==code,f'{model}/{code}: case label mismatch')
        for key,index in (('uv_on',0),('eft1_on',1),('wilson_log_on',2),('smeft_on',3)):
            if key in x:
                require(x[key] == (code[index]=='1'),f'{model}/{code}: bit flag {key} mismatch')
        chi=float_at(x,'chi2',f'{model}/{code}')
        masses=x.get('masses_ev_ascending')
        require(isinstance(masses,list) and len(masses)==3 and
                all(isinstance(m,(int,float)) and not isinstance(m,bool) and math.isfinite(float(m)) and m>=0 for m in masses),
                f'{model}/{code}: invalid neutrino masses')
        require(masses==sorted(masses),f'{model}/{code}: masses not ascending')
        pred=x.get('prediction')
        require(isinstance(pred,dict),f'{model}/{code}: missing oscillation observables')
        threshold,norm=decode_c5(x.get('c5_at_scalar_threshold'),f'{model}/{code} scalar threshold')
        _,low_norm=decode_c5(x.get('c5_at_low_scale'),f'{model}/{code} low scale')
        rows[code]={
            'case':code,'chi2':chi,'m1_eV':float(masses[0]),'m2_eV':float(masses[1]),'m3_eV':float(masses[2]),
            'dm21_sq_eV2':float_at(pred,'delta_m21_sq_ev2',f'{model}/{code}/prediction'),
            'dm3l_sq_eV2':float_at(pred,'delta_m3l_sq_ev2',f'{model}/{code}/prediction'),
            'sin2_theta12':float_at(pred,'sin2_theta12',f'{model}/{code}/prediction'),
            'sin2_theta13':float_at(pred,'sin2_theta13',f'{model}/{code}/prediction'),
            'sin2_theta23':float_at(pred,'sin2_theta23',f'{model}/{code}/prediction'),
            'c5_norm_at_scalar_threshold_GeVinv':norm,
            'c5_norm_at_low_scale_GeVinv':low_norm,
            'coupled_ode_nfev':x.get('combined_ode_nfev'),
        }
    return raw,rows


def write_csv(path,fieldnames,records):
    with path.open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fieldnames,extrasaction='ignore')
        writer.writeheader()
        writer.writerows(records)


def svg_heatmap(path,models,rows):
    """Self-contained chart: delta chi2 = chi2(case) - chi2(1111)."""
    cellw,cellh=47,25
    left,top=95,105
    w=left+16*cellw+40
    h=top+16*cellh+110
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
           '<rect width="100%" height="100%" fill="#fff"/>',
           '<text x="15" y="25" font-family="Arial" font-size="17" font-weight="bold">T3 rough 16-model stage-freeze sensitivity</text>',
           '<text x="15" y="45" font-family="Arial" font-size="12">Cells show chi2(case) - chi2(1111); saturated at +/- 1000 for colour only.</text>',
           '<text x="15" y="63" font-family="Arial" font-size="11">Selected A(-4), D(-2) cached non-HHSS; other 14 historical candidate contexts.</text>']
    for j,code in enumerate(CASE_CODES):
        x=left+j*cellw+cellw/2
        parts.append(f'<text x="{x}" y="{top-12}" text-anchor="middle" font-family="Arial" font-size="10">{code}</text>')
    for i,model in enumerate(models):
        y=top+i*cellh
        parts.append(f'<text x="{left-9}" y="{y+17}" text-anchor="end" font-family="Arial" font-size="12">{html.escape(model)}</text>')
        baseline=rows[model]['1111']['chi2']
        for j,code in enumerate(CASE_CODES):
            delta=rows[model][code]['chi2']-baseline
            mag=min(abs(delta),1000.)/1000.
            if delta>0:
                r,g,b=round(245-50*mag),round(245-180*mag),round(245-180*mag)
            else:
                r,g,b=round(245-180*mag),round(245-95*mag),round(245-40*mag)
            colour=f'#{r:02x}{g:02x}{b:02x}'
            x=left+j*cellw
            tooltip=f'{model}: {code}; delta chi2 = {delta:+.6g}'
            parts.append(f'<rect x="{x}" y="{y}" width="{cellw-2}" height="{cellh-2}" fill="{colour}" stroke="#dddddd" stroke-width=".5"><title>{html.escape(tooltip)}</title></rect>')
    parts.append(f'<text x="15" y="{h-38}" font-family="Arial" font-size="12">White = no effect; red = worsened fit; blue-green = improved fit. Comparison is conditional, not physically certified.</text>')
    parts.append(f'<text x="15" y="{h-18}" font-family="Arial" font-size="11">Columns are UV / EFT1 / Wilson / SMEFT on-off bits; 1111 is the all-on baseline.</text>')
    parts.append('</svg>')
    path.write_text('\n'.join(parts),encoding='utf-8')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--experimental-run',required=True,help='Successful Batch 3C run ID, not a path')
    parser.add_argument('--output',type=Path,required=True,help='Fresh output directory; must not exist')
    args=parser.parse_args(argv)
    project=args.project.resolve()
    output=args.output.resolve()
    require(project.is_dir(),f'missing project {project}')
    require(not output.exists(),f'refusing to overwrite existing output: {output}')
    require(args.experimental_run and args.experimental_run.replace('_','').replace('-','').isalnum(),
            'invalid experimental run ID')
    reg=load(project/REGISTRY)
    models=sorted(reg['models'])
    require(len(models)==16 and EXPERIMENTAL <= set(models) and MIXED <= set(models),
            'expected 16 distinct models in general quartic registry')
    diag_path=project/EXPERIMENT_DIR/args.experimental_run/'quartic_context_diagnostics.json'
    summary_path=project/EXPERIMENT_DIR/args.experimental_run/'summary.json'
    require(diag_path.is_file() and summary_path.is_file(),'missing successful Batch 3C diagnostics or summary')
    diag=load(diag_path)
    experiment_summary=load(summary_path)
    require(diag.get('status')=='ConditionalExperimentCompleted' and
            diag.get('source_HHSS_preserved') is True and
            diag.get('cached_nonHHSS_used') is True,
            'Batch 3C run not complete or not preserving source HHSS')
    require(experiment_summary.get('computed')==2 and
            set(experiment_summary.get('requested_models',[]))==EXPERIMENTAL,
            'Batch 3C expected exactly A_-4 and D_-2 completed')
    for model in EXPERIMENTAL:
        item=diag.get('scalar_response',{}).get(model,{})
        require(item.get('number_of_context_calls',0)>0,
                f'{model}: Batch 3C cached tensor was never evaluated')
    selected={}; records={}; selected_rows={}; old_rows={}; all_cases=[]; baseline=[]; freezes=[]; changes=[]
    for model in models:
        chosen=candidate_path(project,args.experimental_run,model) if model in EXPERIMENTAL else historical_path(project,model)
        raw,rows=verify_model(chosen,model)
        kind=('cached_candidate_nonHHSS_source_HHSS' if model in EXPERIMENTAL else
              'historical_mixed_cartesian' if model in MIXED else
              'historical_provisional_nonHHSS')
        selected[model]={'source':str(chosen.relative_to(project)), 'sha256':sha(chosen),
                         'scalar_context':kind,'cases':16,
                         'physical_matching_certified':False}
        records[model]=raw['cases'];selected_rows[model]=rows
        main=rows['1111']
        block={'model':model,'scalar_context':kind,'chi2_1111':main['chi2'],
               'm1_eV':main['m1_eV'],'m2_eV':main['m2_eV'],'m3_eV':main['m3_eV'],
               **{k:main[k] for k in ('dm21_sq_eV2','dm3l_sq_eV2','sin2_theta12','sin2_theta13','sin2_theta23')},
               **{'chi2_'+c:rows[c]['chi2'] for c in ('0111','1011','1101','1110')},
               'c5_norm_at_scalar_threshold_GeVinv':main['c5_norm_at_scalar_threshold_GeVinv'],
               'source_cases_sha256':sha(chosen)}
        baseline.append(block)
        for code in CASE_CODES:
            rr={'model':model,'scalar_context':kind,
                **{STAGE[i]:int(code[i]) for i in range(4)},**rows[code]}
            all_cases.append(rr)
        for i,stage in enumerate(STAGE):
            code=''.join('0' if j==i else '1' for j in range(4))
            freezes.append({'model':model,'scalar_context':kind,'stage_frozen':stage,'case':code,
                            'chi2_full':main['chi2'],'chi2_frozen':rows[code]['chi2'],
                            'delta_chi2_frozen_minus_full':rows[code]['chi2']-main['chi2'],
                            'm3_full_eV':main['m3_eV'],'m3_frozen_eV':rows[code]['m3_eV'],
                            'delta_m3_eV':rows[code]['m3_eV']-main['m3_eV'],
                            'c5_scalar_threshold_relative_difference':relative_c5(
                                records[model][code]['c5_at_scalar_threshold'],
                                records[model]['1111']['c5_at_scalar_threshold'])})
        if model in EXPERIMENTAL:
            hist=historical_path(project,model)
            prior,prev_rows=verify_model(hist,model)
            old_rows[model]=prev_rows
            selected[model]['previous_source']=str(hist.relative_to(project))
            selected[model]['previous_sha256']=sha(hist)
            for code in CASE_CODES:
                a,b=prev_rows[code],rows[code]
                changes.append({'model':model,'case':code,'chi2_prior':a['chi2'],
                                'chi2_selected':b['chi2'],'delta_chi2':b['chi2']-a['chi2'],
                                'm3_prior_eV':a['m3_eV'],'m3_selected_eV':b['m3_eV'],
                                'delta_m3_eV':b['m3_eV']-a['m3_eV'],
                                'c5_scalar_threshold_relative_difference':relative_c5(
                                    records[model][code]['c5_at_scalar_threshold'],
                                    prior['cases'][code]['c5_at_scalar_threshold'])})
    require(len(baseline)==16 and len(all_cases)==256 and len(freezes)==64 and len(changes)==32,
            'internal consolidated row-count error')
    output.mkdir(parents=True)
    write_csv(output/'models_1111.csv',MODEL_FIELDS,baseline)
    write_csv(output/'all_256_stage_freeze_cases.csv',CASE_FIELDS,all_cases)
    write_csv(output/'single_stage_freeze_effects.csv',DELTA_FIELDS,freezes)
    write_csv(output/'A4_D2_cached_vs_previous.csv',COMPARE_FIELDS,changes)
    svg_heatmap(output/'stage_freeze_delta_chi2.svg',models,selected_rows)
    report={'format':SCHEMA,'status':'CompleteConditionalRoughResults',
            'generated_utc':datetime.now(timezone.utc).isoformat(),
            'models_computed':len(models),'selected_cases':len(all_cases),
            'single_stage_freeze_comparisons':len(freezes),
            'A4_D2_before_after_comparisons':len(changes),
            'scales_gev':{'UV':1e7,'fermions_out':1e5,'both_scalars_out':1e3,'low':100.},
            'EFT_order':'UV -> integrate all heavy fermions -> EFT1 with both heavy scalars -> jointly integrate both scalars -> SMEFT',
            'selection':selected,
            'model_1111':[{'model':r['model'],'scalar_context':r['scalar_context'],
                            'chi2':r['chi2_1111'],'masses_eV':[r['m1_eV'],r['m2_eV'],r['m3_eV']]}
                           for r in baseline],
            'experimental_run':args.experimental_run,
            'experimental_summary_sha256':sha(summary_path),
            'experimental_diagnostics_sha256':sha(diag_path),
            'limitations':[
                'Absolute Matchete-to-Weyl normalization and signs are not physically certified.',
                'Fermion hard finite one-loop UV-minus-EFT1 matching has not been completed.',
                'Full native Matchete/RGBeta scalar quartic coordinate dictionary is not certified.',
                'D_-2 extra quadratic mass direction and full Hessian remain uncertified.',
                'Joint two-scalar threshold uses saved conditional finite-bubble prescription, not complete one-loop matching.',
                'Only A_-4 and D_-2 have Batch 3C cached non-HHSS; other models retain historical provisional quartic contexts.',
                'Neutrino mass and chi2 numbers are benchmark outputs, not global fits with verified physical normalization.'
            ],
            'files':['models_1111.csv','all_256_stage_freeze_cases.csv',
                     'single_stage_freeze_effects.csv','A4_D2_cached_vs_previous.csv',
                     'stage_freeze_delta_chi2.svg','summary.json'],
            'historical_results_modified':False}
    (output/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('PASS: 16 models, 256 stage-freeze cases, 64 single-stage effects, 32 before/after comparisons')
    print('Conditional results:',output)
    for row in baseline:
        print(f"{row['model']:5s}  chi2={row['chi2_1111']:.6g}  m3={row['m3_eV']:.7g} eV  {row['scalar_context']}")
    return 0


if __name__=='__main__':
    raise SystemExit(main())
