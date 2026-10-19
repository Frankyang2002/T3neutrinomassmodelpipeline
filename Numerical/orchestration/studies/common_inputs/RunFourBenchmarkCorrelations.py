"""Evaluate four fixed Yukawa/scalar benchmarks for 16 T3 models and build an interactive dashboard.

No optimisation is performed. Existing fitted references, common-input files and
previous comparison outputs are read-only. The new output folder is separate.

Usage:
 python -m Numerical.orchestration.studies.common_inputs.RunFourBenchmarkCorrelations
 python -m Numerical.orchestration.studies.common_inputs.RunFourBenchmarkCorrelations --run
 python -m Numerical.orchestration.studies.common_inputs.RunFourBenchmarkCorrelations --display-only
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from copy import deepcopy
from pathlib import Path

from Numerical.orchestration.studies.common_inputs.RunCommonParameterComparisonV2 import (
    PROJECT_ROOT, SCENARIO, candidate_for_target, evaluate, load, models_and_cases,
)
from Numerical.fitting.OscillationFit import OscillationFitTarget

OUT = PROJECT_ROOT / 'Reports/output/full/four_benchmark_correlations' / SCENARIO
COMMON = PROJECT_ROOT / 'Reports/output/full/common_yukawa_least_squares' / SCENARIO / 'optimized_common_input.json'
CSV_PATH = OUT / 'four_benchmark_predictions.csv'
HTML_PATH = OUT / 'interactive_comparison' / 'index.html'
OBS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12', 'sin2_theta13', 'sin2_theta23')
BENCHMARKS = {'LL': (1., 1.), 'HL': (100., 1.), 'LH': (1., 100.), 'HH': (100., 100.)}
FIELDS = ('model', 'benchmark', 'yukawa_factor', 'scalar_factor', 'status', 'error', 'chi2', *OBS)


def _multiply_complex_matrix(matrix, factor):
    for part in ('real', 'imag'):
        matrix[part] = [[float(v) * factor for v in row] for row in matrix[part]]


def benchmark_payload(source, target, benchmark):
    result = candidate_for_target(source, target)
    ordinary = result['base_state']['ordinary']
    yf, sf = BENCHMARKS[benchmark]
    for key in ('y1', 'y2'):
        _multiply_complex_matrix(ordinary[key], yf)
    lam = ordinary['lambdaT3']
    if isinstance(lam, dict):
        # Some versions encode scalar couplings as {real, imag}; preserve shape.
        for part in ('real', 'imag'):
            if part in lam:
                val = lam[part]
                if isinstance(val, list):
                    raise ValueError('lambdaT3 is matrix-valued; scalar benchmark definition needs confirmation')
                lam[part] = float(val) * sf
    elif isinstance(lam, (int, float)):
        ordinary['lambdaT3'] = float(lam) * sf
    else:
        raise TypeError(f'Unexpected lambdaT3 type: {type(lam).__name__}')
    return result


def read_rows():
    if not CSV_PATH.exists():
        return []
    with CSV_PATH.open(newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        if not set(FIELDS).issubset(reader.fieldnames or []):
            raise ValueError(f'Unexpected CSV header in {CSV_PATH}')
        return list(reader)


def save_rows(rows):
    OUT.mkdir(parents=True, exist_ok=True)
    tmp = CSV_PATH.with_suffix('.csv.tmp')
    with tmp.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(CSV_PATH)


def run(limit_models=0, retry_failed=False):
    configs, cases = models_and_cases()
    names = sorted(configs)
    if limit_models:
        names = names[:limit_models]
    if not COMMON.is_file():
        raise FileNotFoundError(f'Optimised common input not found: {COMMON}')
    shared = load(COMMON)['shared_ordinary']
    source_name = load(COMMON)['source_model']
    source = load(configs[source_name])
    for key in ('y1', 'y2', 'lambdaT3', 'MF', 'mS1Sq', 'mS2Sq'):
        source['base_state']['ordinary'][key] = deepcopy(shared[key])
    target = OscillationFitTarget.from_json(PROJECT_ROOT / 'configs/nufit_6_1_ic24_no.json')
    rows = read_rows()
    seen = {(r['model'], r['benchmark']) for r in rows if not (retry_failed and r['status'] != 'Success')}
    if retry_failed:
        rows = [r for r in rows if r['status'] == 'Success' or r['model'] not in names]
    for name in names:
        original = load(configs[name])
        for bench in BENCHMARKS:
            if (name, bench) in seen:
                continue
            yf, sf = BENCHMARKS[bench]
            row = dict(model=name, benchmark=bench, yukawa_factor=yf, scalar_factor=sf,
                       status='Failed', error='', chi2='', **{k: '' for k in OBS})
            try:
                payload = benchmark_payload(source, original, bench)
                chi2, predictions = evaluate(payload, cases[name], target)
                if not math.isfinite(chi2) or any(not math.isfinite(predictions[k]) for k in OBS):
                    raise ValueError('Nonfinite prediction')
                row.update(status='Success', chi2=chi2, **predictions)
            except Exception as exc:
                row['error'] = f'{type(exc).__name__}: {exc}'
            rows.append(row)
            save_rows(rows)  # checkpoint after every evaluation, including failures
            print(f'{name} {bench}: {row["status"]} chi2={row["chi2"]} {row["error"]}', flush=True)
    return rows


def build_display(rows):
    if not rows:
        raise ValueError('No benchmark predictions available; run with --run first')
    try:
        from plotly.offline.offline import get_plotlyjs
    except ImportError as exc:
        raise RuntimeError('Install plotly: python -m pip install plotly') from exc
    data = []
    for row in rows:
        item = {k: row.get(k, '') for k in ('model', 'benchmark', 'status', 'error')}
        for k in ('chi2', *OBS):
            v = row.get(k, '')
            item[k] = float(v) if v != '' and math.isfinite(float(v)) else None
        data.append(item)
    payload = json.dumps(data, allow_nan=False, separators=(',', ':')).replace('</', '<\\/')
    page = r'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>T3 four-benchmark correlations</title>
<style>
:root{color-scheme:light;--ink:#17243b;--muted:#556375;--line:#d7dee8;--paper:#fff;--surface:#f7f9fc;--accent:#225ba8}
*{box-sizing:border-box}body{font-family:system-ui,-apple-system,Segoe UI,Arial,sans-serif;margin:0;background:var(--surface);color:var(--ink)}
header{padding:20px 28px;background:#fff;border-bottom:1px solid var(--line)}header h1{font-size:20px;margin:0 0 5px;font-weight:650}header p{font-size:12px;color:var(--muted);margin:0;line-height:1.5}
.layout{display:grid;grid-template-columns:275px minmax(0,1fr);min-height:calc(100vh - 100px)}aside{padding:20px 17px;background:#fff;border-right:1px solid var(--line);overflow:auto}main{padding:20px 22px;min-width:0}
label.caption{display:block;font-size:12px;font-weight:650;margin:0 0 8px}select{width:100%;padding:10px 8px;border:1px solid #bac7d8;border-radius:7px;background:#fff;color:var(--ink);font-size:13px;margin:0 0 16px}
button{padding:8px 12px;border:1px solid #ccd5e2;border-radius:7px;background:white;color:var(--ink);font-size:12px;cursor:pointer;white-space:nowrap}button:hover{background:#eef4ff;border-color:#acc6e9}button:focus-visible,select:focus-visible,input:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
.panel{border:1px solid var(--line);border-radius:10px;background:var(--paper);padding:16px;box-shadow:0 1px 2px #182b480b}.zoom{display:flex;gap:7px;margin-bottom:12px;flex-wrap:wrap}.graph{width:100%;height:min(75vh,740px);min-height:430px}.small{font-size:12px;color:var(--muted);line-height:1.5}.note{margin-top:12px;padding:11px 13px;background:#f5f8fc;border-left:3px solid #b6cce6;border-radius:4px;font-size:12px;line-height:1.5}#status{font-size:12px;color:var(--muted);margin-top:9px}.toggle{font-size:12px;display:block;margin:10px 0;cursor:pointer}.toggle input{accent-color:var(--accent)}
@media(max-width:830px){.layout{grid-template-columns:1fr}aside{border-right:0;border-bottom:1px solid var(--line)}main{padding:12px}.panel{padding:8px}header{padding:15px}.graph{min-height:360px;height:65vh}}
</style></head><body>
<header><h1>T3 neutrino correlations · four benchmark comparison</h1><p>Previously optimised LL input · Yukawa and λT3 high values = 100 × LL · fixed mass parameters</p></header>
<div class="layout"><aside>
<label class="caption" for="mode">Comparison</label><select id="mode"><option value="models">Between models · fixed benchmark</option><option value="benchmarks">Within model · four benchmarks</option></select>
<div id="benchwrap"><label class="caption" for="bench">Benchmark</label><select id="bench"><option>LL</option><option>HL</option><option>LH</option><option>HH</option></select></div>
<div id="modelwrap" style="display:none"><label class="caption" for="model">Model</label><select id="model"></select></div>
<label class="caption" for="graph">Correlation</label><select id="graph"><option value="13">Δm²₃₁ vs sin²θ₁₃</option><option value="23">Δm²₃₁ vs sin²θ₂₃</option><option value="cp">δCP vs sin²θ₁₃</option></select>
<label class="caption">Experimental regions</label><label class="toggle"><input type="checkbox" id="one" checked> Illustrative 1σ ellipse</label><label class="toggle"><input type="checkbox" id="three" checked> Illustrative 3σ ellipse</label>
<div id="modelchecks"><label class="caption">Visible models</label><div style="display:flex;gap:6px;margin-bottom:9px"><button type="button" id="allmodels">Select all</button><button type="button" id="nomodels">Deselect all</button></div><div id="checklist" style="max-height:34vh;overflow:auto;border:1px solid var(--line);padding:8px;border-radius:6px"></div></div><p class="small">Select a benchmark to compare representations, or select a model to compare LL, HL, LH and HH.</p>
</aside><main><div class="panel"><div class="zoom"><button id="reset">Reset zoom</button></div><div id="plot" class="graph"></div><div id="status"></div></div>
<div class="note">NuFIT ellipses use independent marginal 1σ uncertainties: (x−x₀)²/σₓ² + (y−y₀)²/σᵧ² = n². These are <b>illustrative diagonal-Gaussian ellipses, not official NuFIT joint confidence regions</b>. No δCP values are present in the existing five-observable evaluator; that panel reports unavailable data rather than inventing predictions. Lines connect discrete benchmark points in LL → HL → LH → HH order; they are visual guides, not physical trajectories. Different models are not connected.</div>
</main></div>
<script>/*PLOTLY*/</script><script>
const DATA=/*DATA*/;const BENCH=['LL','HL','LH','HH'];const COLORS={LL:'#2273b9',HL:'#e67e22',LH:'#258a62',HH:'#9b59b6'};
const C={sin2_theta13:0.02248,sin2_theta23:0.47,delta_m3l_sq_ev2:0.002511};const S={sin2_theta13:0.00057,sin2_theta23:0.0155,delta_m3l_sq_ev2:0.0000205};
const $=id=>document.getElementById(id);const models=[...new Set(DATA.map(x=>x.model))].sort();models.forEach(m=>{let o=document.createElement('option');o.textContent=m;o.value=m;$('model').appendChild(o);let lab=document.createElement('label');lab.className='toggle';lab.style.overflowWrap='anywhere';let cb=document.createElement('input');cb.type='checkbox';cb.checked=true;cb.value=m;cb.className='model-toggle';cb.addEventListener('change',()=>draw(true));lab.append(cb,document.createTextNode(' '+m));$('checklist').appendChild(lab)});$('allmodels').onclick=()=>{document.querySelectorAll('.model-toggle').forEach(c=>c.checked=true);draw(true)};$('nomodels').onclick=()=>{document.querySelectorAll('.model-toggle').forEach(c=>c.checked=false);draw(true)};
function ellipse(x,y,n){let xs=[],ys=[];for(let i=0;i<=180;i++){let a=2*Math.PI*i/180;xs.push(C[x]+n*S[x]*Math.cos(a));ys.push(C[y]+n*S[y]*Math.sin(a))}return {type:'scatter',mode:'lines',x:xs,y:ys,fill:'toself',line:{color:n===1?'#2588c7':'#86bce0',width:1},fillcolor:n===1?'rgba(30,132,209,.22)':'rgba(110,187,231,.12)',name:n+'σ (illustrative)',hoverinfo:'skip'}}
function draw(reset=false){const mode=$('mode').value, bench=$('bench').value, model=$('model').value, graph=$('graph').value; $('benchwrap').style.display=mode==='models'?'':'none';$('modelwrap').style.display=mode==='benchmarks'?'':'none';$('modelchecks').style.display=mode==='models'?'':'none';const x=graph==='23'?'sin2_theta23':'sin2_theta13',y=graph==='cp'?'delta_cp_deg':'delta_m3l_sq_ev2';let traces=[];if(graph!=='cp'){if($('three').checked)traces.push(ellipse(x,y,3));if($('one').checked)traces.push(ellipse(x,y,1));traces.push({type:'scatter',mode:'markers',x:[C[x]],y:[C[y]],name:'NuFIT centre',marker:{color:'#25364a',symbol:'cross',size:12}})}
const visible=new Set([...document.querySelectorAll('.model-toggle:checked')].map(c=>c.value));let selected=DATA.filter(r=>mode==='models'?(r.benchmark===bench&&visible.has(r.model)):r.model===model);let valid=selected.filter(r=>r.status==='Success'&&r[x]!=null&&r[y]!=null);if(mode==='models'){valid.sort((a,b)=>a.model.localeCompare(b.model));valid.forEach((r,i)=>{traces.push({type:'scatter',mode:'lines+markers',x:[r[x]],y:[r[y]],name:r.model,marker:{size:12,color:`hsl(${Math.round(i*360/Math.max(models.length,1))},65%,44%)`,symbol:i%2?'diamond':'circle',line:{width:1,color:'#fff'}},customdata:[[r.benchmark,r.chi2]],hovertemplate:r.model+'<br>Benchmark: %{customdata[0]}<br>x=%{x:.8g}<br>y=%{y:.8g}<br>χ²=%{customdata[1]:.3g}<extra></extra>'})})}else{valid.sort((a,b)=>BENCH.indexOf(a.benchmark)-BENCH.indexOf(b.benchmark));if(valid.length){traces.push({type:'scatter',mode:'lines+markers',x:valid.map(r=>r[x]),y:valid.map(r=>r[y]),name:model,connectgaps:false,line:{color:'#225ba8',width:2},marker:{size:11,color:valid.map(r=>COLORS[r.benchmark]),line:{width:1,color:'#fff'}},text:valid.map(r=>r.benchmark),customdata:valid.map(r=>[r.chi2]),hovertemplate:'%{text}<br>x=%{x:.8g}<br>y=%{y:.8g}<br>χ²=%{customdata[0]:.3g}<extra></extra>'})}}

const xlabel=x==='sin2_theta13'?'sin²θ₁₃':'sin²θ₂₃',ylabel=graph==='cp'?'δCP [degrees]':'Δm²₃₁ [eV²]';let layout={title:{text:ylabel+' versus '+xlabel,font:{size:20}},xaxis:{title:xlabel,gridcolor:'#e9edf2',zeroline:false},yaxis:{title:ylabel,gridcolor:'#e9edf2',zeroline:false},legend:{orientation:'v',x:1.02,y:1,font:{size:11}},margin:{l:95,r:mode==='models'?245:145,t:75,b:85},paper_bgcolor:'white',plot_bgcolor:'white',hovermode:'closest',uirevision:reset?String(Date.now()):'stable'};Plotly.react('plot',traces,layout,{responsive:true,displaylogo:false,toImageButtonOptions:{format:'png',scale:2}});let failed=selected.filter(r=>r.status!=='Success');$('status').textContent=`${mode==='models'?bench:model}: ${valid.length} plotted, ${failed.length} failed, ${selected.length} rows available.`+(graph==='cp'?' δCP unavailable: no CP observable was saved.':'')+(failed.length?' Failed cases: '+failed.map(r=>r.model+' '+r.benchmark+': '+(r.error||'unknown error')).join(' | '):'');}
for(const id of ['mode','bench','model','graph','one','three'])$(id).addEventListener('change',()=>draw(true));$('reset').addEventListener('click',()=>draw(true));draw(true);
</script></body></html>'''
    page = page.replace('/*PLOTLY*/', get_plotlyjs()).replace('/*DATA*/', payload)
    HTML_PATH.parent.mkdir(parents=True, exist_ok=True)
    HTML_PATH.write_text(page, encoding='utf-8')
    print('Interactive display:', HTML_PATH)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', action='store_true', help='Evaluate missing model/benchmark pairs and build display')
    p.add_argument('--display-only', action='store_true', help='Build HTML from existing CSV')
    p.add_argument('--retry-failed', action='store_true', help='Re-evaluate failed CSV rows (use with --run)')
    p.add_argument('--limit-models', type=int, default=0, help='Pilot first N models')
    args = p.parse_args()
    if args.limit_models < 0:
        p.error('--limit-models must be nonnegative')
    print('Optimised input:', COMMON)
    print('Output CSV:', CSV_PATH)
    print('Interactive HTML:', HTML_PATH)
    print('Four benchmarks:', BENCHMARKS)
    if args.run:
        rows = run(args.limit_models, args.retry_failed)
        build_display(rows)
    elif args.display_only:
        build_display(read_rows())
    else:
        print('Dry run only. Use --run to evaluate; --display-only for an existing CSV.')


if __name__ == '__main__':
    main()
