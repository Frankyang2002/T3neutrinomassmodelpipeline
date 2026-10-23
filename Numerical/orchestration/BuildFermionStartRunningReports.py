"""Generate true fermion-start (011) numerical trajectories and offline reports.

Reads the completed comparison configs, reruns ONLY EFT1 and Weinberg RGEs.
No UV integration, fitting, symbolic generation or modifications to saved fits.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from html import escape

import numpy as np

from Numerical.orchestration import RunSmallLambdaRefit as base
from Numerical.orchestration.RunFermionStartLambdaStudy import (
    ROOT, FIELDS, CASES, FermionStartEvaluator, artifacts_for, vector_values,
)
from Numerical.orchestration.ModelNames import shorten
from Numerical.running.SMWeinbergEvolution import evolve_sm_weinberg
from Numerical.running.ScalarThresholdBoundary import build_sm_weinberg_initial_conditions
from Numerical.fitting.ParameterScan import neutrino_mass_from_c5, charged_lepton_mass_basis_matrix
from Numerical.fitting.OscillationFit import evaluate_mass_matrix_fit

OUT = base.ROOT / 'Reports' / 'output' / 'display_fermion_start'
MODES = ('optimised', 'shared_a_m2')


def atomic_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix('.json.tmp')
    temp.write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    temp.replace(path)


def prediction(K, ye, evaluator):
    mnu = neutrino_mass_from_c5(np.asarray(K, dtype=complex), vev_gev=evaluator.vev)
    mnu = charged_lepton_mass_basis_matrix(mnu, ye)
    fit = evaluate_mass_matrix_fit(mnu, evaluator.target)
    return {'chi2': float(fit.chi2), **{
        str(k): float(v) for k, v in zip(fit.observable_names, fit.prediction)
    }}


def generate_one(case, mode, model, artifact_map, n_samples):
    source = ROOT / case / 'comparison' / mode / 'models' / model
    cfg_path = source / 'fermion_start_config.json'
    endpoint_path = source / 'evaluation.json'
    cfg = base.load(cfg_path)
    saved = base.load(endpoint_path)
    if saved.get('status') != 'Success':
        raise ValueError('Endpoint comparison did not succeed')
    ev = FermionStartEvaluator(cfg, artifacts_for(model, artifact_map))
    result, trajectory, low = ev.evaluate(vector_values(cfg), return_trajectory=True)
    if not math.isclose(float(result.chi2), float(saved['chi2']), rel_tol=1e-5, abs_tol=1e-5):
        raise ValueError(f'Endpoint chi2 mismatch: {result.chi2} vs {saved["chi2"]}')
    for key in FIELDS:
        i = list(result.observable_names).index(key)
        if not math.isclose(float(result.prediction[i]), float(saved['prediction'][key]), rel_tol=2e-5, abs_tol=2e-10):
            raise ValueError(f'Endpoint observable mismatch: {key}')

    # EFT1 is the actual scalar-only ODE solution, not interpolation or a UV reconstruction.
    eft = trajectory.eft1
    eft_keys = ('gY', 'g2', 'g3', 'lambdaH', 'lambdaT3',
                'lambdaS1', 'lambdaS2', 'lambdaH1', 'lambdaH2', 'lambda12',
                'lambda3', 'lambda4', 'lambda5', 'lambdaS')
    eft_series = {}
    for index in range(eft.n_points):
        state = eft.state_at_index(index)
        scalars = {'gY': float(state.sm.gY), 'g2': float(state.sm.g2),
                   'g3': float(state.sm.g3), 'lambdaH': float(state.sm.lambdaH)}
        for key in eft_keys:
            if hasattr(state, key):
                z = getattr(state, key)
                scalars[key] = float(abs(z)) if isinstance(z, complex) else float(z)
        for key, value in scalars.items():
            eft_series.setdefault(key, []).append(value)
    eft_mu = [float(mu) for mu in eft.mu_gev]
    if not (math.isclose(eft_mu[0], ev.mu_f, rel_tol=1e-10) and
            math.isclose(eft_mu[-1], ev.mu_s, rel_tol=1e-10)):
        raise ValueError('EFT1 scale range mismatch')

    # Evolve Weinberg EFT from the SAME matched C5 initial state to each sample scale.
    # Each sample is a true RGE integration from mu_S; no matrix interpolation.
    initial = build_sm_weinberg_initial_conditions(
        trajectory.final_sm_boundary, low.c5_at_scalar_threshold)
    mu_grid = np.geomspace(ev.mu_s, ev.mu_low, n_samples)
    mu_grid[0], mu_grid[-1] = ev.mu_s, ev.mu_low
    final_series = {'chi2': [], 'c5_frobenius': [], 'g2': [], 'gY': []}
    for key in FIELDS: final_series[key] = []
    for idx, mu in enumerate(mu_grid):
        if idx == 0:
            K, ye = initial.K, initial.ye
            gY, g2 = initial.gY, initial.g2
        elif idx == len(mu_grid) - 1:
            K, ye = low.weinberg.K, low.weinberg.ye
            gY, g2 = low.weinberg.gY, low.weinberg.g2
        else:
            item = evolve_sm_weinberg(initial, ev.mu_s, float(mu),
                                      rtol=ev.rtol, atol=ev.atol)
            K, ye = item.K, item.ye
            gY, g2 = item.gY, item.g2
        obs = prediction(K, ye, ev)
        final_series['chi2'].append(obs['chi2'])
        final_series['c5_frobenius'].append(float(np.linalg.norm(K)))
        final_series['gY'].append(float(gY))
        final_series['g2'].append(float(g2))
        for key in FIELDS: final_series[key].append(obs[key])
    if not math.isclose(final_series['chi2'][-1], float(result.chi2),
                        abs_tol=2e-5, rel_tol=1e-5):
        raise ValueError('Weinberg endpoint cross-check failed')

    data = {'study': 'Fermion Start — UV Off', 'case': case, 'mode': mode,
            'model': model, 'label': shorten(model), 'lambdaT3':
            float(cfg['base_state']['ordinary']['lambdaT3']['real']),
            'UV_running': False, 'fermion_matching': True,
            'EFT1_running': True, 'scalar_matching': True, 'Weinberg_running': True,
            'scales_gev': {'mu_F': ev.mu_f, 'mu_S': ev.mu_s, 'mu_low': ev.mu_low},
            'eft1': {'mu_gev': eft_mu, 'series': eft_series},
            'weinberg': {'mu_gev': [float(x) for x in mu_grid], 'series': final_series},
            'endpoint': saved}
    dest = OUT / case / 'running' / mode / 'models' / model / 'trajectory.json'
    atomic_json(dest, data)
    return data


def write_csv(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.writer(fh)
        writer.writerow(['model', 'mu_gev', *payload['series'].keys()])
        for j, mu in enumerate(payload['mu_gev']):
            writer.writerow([payload.get('label', ''), mu] +
                            [vals[j] for vals in payload['series'].values()])


HTML = r'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title><style>
body{font:14px system-ui,Arial,sans-serif;max-width:1200px;margin:20px auto;padding:0 16px;color:#17243b}
header{padding-bottom:12px;border-bottom:1px solid #ccd5df} select{padding:7px;margin:4px 16px 8px 0}
.chart{width:100%;height:430px;border:1px solid #ccd5df;border-radius:8px;margin-top:12px}
.legend{display:flex;gap:15px;flex-wrap:wrap;margin:12px 0}.legend label{cursor:pointer}
small{color:#555}.error{color:#a12} a{color:#245ca8}
</style></head><body><header><h2>__TITLE__</h2>
<p>Fermion Start — UV Off. Physical interval: μF → μS → μlow. Fermion matching, EFT1 running, scalar matching and Weinberg running enabled. No UV trajectory or artificial continuation.</p>
<p><a href="../../index.html">Case overview</a> · <a href="../../optimised.html">Optimised endpoints</a> · <a href="../../shared_a_m2.html">Shared A(-2) endpoints</a></p></header>
<label>Region <select id="region"><option value="weinberg">SM + Weinberg (μS → μlow)</option><option value="eft1">Scalar-only EFT1 (μF → μS)</option></select></label>
<label>Quantity <select id="quantity"></select></label><label><input id="relative" type="checkbox"> Normalise to first scale</label>
<div id="legend" class="legend"></div><svg id="chart" class="chart" viewBox="0 0 1100 430" role="img" aria-label="Running trajectories"></svg>
<p><small id="note"></small></p><script>const INPUT=__DATA__;
const region=document.getElementById('region'), quantity=document.getElementById('quantity'), norm=document.getElementById('relative');
const NS='http://www.w3.org/2000/svg';let active={};
function tag(name,attrs,text){let e=document.createElementNS(NS,name);Object.entries(attrs||{}).forEach(([k,v])=>e.setAttribute(k,v));if(text!==undefined)e.textContent=text;return e}
const colors=['#2062ab','#cf652c','#32824c','#843c99','#a36b21','#c13c67','#248889','#6f7c2a','#4a596b','#a64743','#7c47bd','#4a829f','#ac922c','#644c36','#d05f98','#49937c'];
function changeRegion(){let q=region.value,keys=[...new Set(INPUT.filter(a=>a[q]).flatMap(a=>Object.keys(a[q].series)))];
let preference=q==='weinberg'?'delta_m3l_sq_ev2':'lambdaT3'; quantity.innerHTML='';keys.sort().forEach(k=>{let o=document.createElement('option');o.value=k;o.textContent=k;quantity.appendChild(o)});quantity.value=keys.includes(preference)?preference:keys[0];render()}
function render(){let q=region.value,k=quantity.value,svg=document.getElementById('chart');svg.replaceChildren();let models=INPUT.filter(a=>a[q]&&a[q].series[k]);
let legend=document.getElementById('legend');legend.replaceChildren();models.forEach((a,i)=>{let key=a.label;if(!(key in active))active[key]=true;let l=document.createElement('label'),cb=document.createElement('input');cb.type='checkbox';cb.checked=active[key];cb.onchange=()=>{active[key]=cb.checked;render()};l.append(cb,document.createTextNode(' '+key));l.style.color=colors[i%colors.length];legend.appendChild(l)});
let curves=models.filter(a=>active[a.label]).map(a=>{let x=a[q].mu_gev.map(Math.log10),v=a[q].series[k];let div=norm.checked?v[0]:1;let y=v.map(z=>div!==0?z/div:z);return {a,x,y}}); if(!curves.length)return;
let xx=curves.flatMap(v=>v.x),yy=curves.flatMap(v=>v.y);let xmin=Math.min(...xx),xmax=Math.max(...xx),ymin=Math.min(...yy),ymax=Math.max(...yy); if(ymin===ymax){ymin-=1;ymax+=1}let pad=.05*(ymax-ymin);ymin-=pad;ymax+=pad;
let X=x=>80+950*(x-xmin)/(xmax-xmin||1),Y=y=>360-305*(y-ymin)/(ymax-ymin);for(let t=0;t<=5;t++){let y=55+305*t/5;svg.appendChild(tag('line',{x1:80,x2:1030,y1:y,y2:y,stroke:'#d9e0e7'}));svg.appendChild(tag('text',{x:73,y:y+4,'font-size':12,'text-anchor':'end'},(ymax-(ymax-ymin)*t/5).toPrecision(4)))}
for(let t=0;t<=5;t++){let x=80+950*t/5;svg.appendChild(tag('line',{x1:x,x2:x,y1:55,y2:360,stroke:'#edf0f5'}));svg.appendChild(tag('text',{x:x,y:381,'font-size':12,'text-anchor':'middle'},(xmin+(xmax-xmin)*t/5).toFixed(2)))}
svg.appendChild(tag('text',{x:550,y:412,'text-anchor':'middle','font-size':13},'log₁₀(μ / GeV)'));
svg.appendChild(tag('text',{x:80,y:24,'font-size':15},k+(norm.checked?' / initial value':'')));
curves.forEach(c=>{let i=models.findIndex(a=>a.label===c.a.label);let d=c.x.map((x,j)=>(j?'L':'M')+X(x).toFixed(2)+','+Y(c.y[j]).toFixed(2)).join(' ');svg.appendChild(tag('path',{d,stroke:colors[i%colors.length],'stroke-width':2.2,fill:'none'}))});
document.getElementById('note').textContent='These are separately solved physical EFT trajectories. Mixing observables and C5 shown only BELOW μS; no invented C5 curve above μS. χ² is computed at each final-EFT scale for illustration, not a fitted observable at that scale.'}
region.onchange=changeRegion;quantity.onchange=render;norm.onchange=render;changeRegion();</script></body></html>'''


def build_report(case, mode, records, errors):
    root = OUT / case / 'running' / mode
    root.mkdir(parents=True, exist_ok=True)
    title = f'Fermion Start — UV Off | {case} | {mode}'
    content = HTML.replace('__TITLE__', escape(title)).replace('__DATA__',json.dumps(records,separators=(',', ':')))
    (root / 'index.html').write_text(content, encoding='utf-8')
    atomic_json(root / 'summary.json', {'case':case,'mode':mode,'models':len(records),'errors':errors,
                                      'boundary':'fermion threshold','uv_running':False})
    print('Report:',root/'index.html',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case', choices=(*CASES,'all'),default='all')
    p.add_argument('--mode',choices=(*MODES,'all'),default='all')
    p.add_argument('--model',help='Optional exact or short model name')
    p.add_argument('--samples',type=int,default=15,help='Final-EFT sample count, including endpoints')
    p.add_argument('--artifact-map',type=Path)
    a=p.parse_args()
    if a.samples<2: p.error('--samples must be >= 2')
    mapping=base.load(a.artifact_map) if a.artifact_map else None
    cases=tuple(CASES) if a.case=='all' else (a.case,)
    modes=MODES if a.mode=='all' else (a.mode,)
    failed=0
    for case in cases:
        for mode in modes:
            source=ROOT/case/'comparison'/mode/'summary.json'
            summary=base.load(source)
            rows=[];errors=[]
            for record in summary['models']:
                model=record['model']
                if a.model and a.model not in (model,shorten(model)):continue
                if record['status']!='Success':
                    errors.append({'model':model,'reason':'endpoint status '+record['status']});continue
                try:
                    target=OUT/case/'running'/mode/'models'/model
                    data=generate_one(case,mode,model,mapping,a.samples)
                    for key in ('eft1','weinberg'):
                        write_csv(target/(key+'.csv'),dict(data[key],label=data['label']))
                    rows.append(data)
                    print(case,mode,shorten(model),'OK',flush=True)
                except Exception as exc:
                    failed+=1
                    errors.append({'model':model,'reason':f'{type(exc).__name__}: {exc}'})
                    print(case,mode,shorten(model),'FAILED',exc,flush=True)
            build_report(case,mode,rows,errors)
    if failed: raise SystemExit(f'{failed} trajectories failed; inspect running/*/summary.json')

if __name__=='__main__':main()
