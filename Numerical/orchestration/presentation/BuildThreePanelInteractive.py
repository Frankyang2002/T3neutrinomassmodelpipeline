"""Build a self-contained interactive three-panel display from completed scan_points.csv.

Usage:
 python -m Numerical.orchestration.presentation.BuildThreePanelInteractive
 python -m Numerical.orchestration.presentation.BuildThreePanelInteractive --csv PATH --output PATH
No numerical re-evaluation; input CSV is read-only. Uses bundled Plotly JS.
"""
from __future__ import annotations
import argparse
import csv
import html
import json
from pathlib import Path
import numpy as np
from plotly.offline.offline import get_plotlyjs

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_CSV = ROOT / 'Reports/output/full/common_yukawa_three_panel_scans/smallY_smallL/scan_points.csv'
DEFAULT_OUTPUT = ROOT / 'display/T3_common_yukawa_three_panel_interactive.html'
CENTRE = {'sin2_theta13': .02248, 'sin2_theta23': .47, 'delta_m3l_sq_ev2': .002511}
SIGMA = {'sin2_theta13': .00057, 'sin2_theta23': .0155, 'delta_m3l_sq_ev2': .0000205}


def build(source: Path, destination: Path) -> None:
    with source.open(newline='', encoding='utf-8-sig') as stream:
        rows = list(csv.DictReader(stream))
    if not rows or not {'model','direction','t','status','chi2','sin2_theta13','sin2_theta23','delta_m3l_sq_ev2','delta_cp_deg'} <= set(rows[0]):
        raise ValueError('CSV missing required three-panel scan columns')
    data = {}
    failures = []
    for row in rows:
        model = row['model']
        data.setdefault(model, {'S1': [], 'S2': [], 'S3': []})
        if row['status'] != 'Success':
            failures.append({'model': model, 'direction': row['direction'], 't': row['t'], 'error': row.get('error','')})
            continue
        record = {k: float(row[k]) for k in ('t','chi2','sin2_theta13','sin2_theta23','delta_m3l_sq_ev2','delta_cp_deg')}
        if not all(np.isfinite(v) for v in record.values()):
            raise ValueError(f'Nonfinite successful row: {model}')
        data[model][row['direction']].append(record)
    for model in data:
        for direction in data[model]:
            data[model][direction].sort(key=lambda r:r['t'])
    payload = json.dumps(data, allow_nan=False, separators=(',',':')).replace('</','<\\/')
    failed = json.dumps(failures, allow_nan=False).replace('</','<\\/')
    js = get_plotlyjs()
    page = r'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>T3 | Neutrino correlation explorer</title>
<style>body{font:14px system-ui,-apple-system,Segoe UI,sans-serif;margin:0;background:#f7f9fc;color:#1b2939}header{padding:22px 28px;background:#13283d;color:white}h1{font-size:21px;margin:0 0 6px}header p{margin:0;color:#cbd9e7}.controls{display:flex;gap:16px;align-items:center;flex-wrap:wrap;padding:17px 26px;background:white;border-bottom:1px solid #dbe2eb}select{max-width:100%;padding:8px;border:1px solid #abbaca;border-radius:5px}label{white-space:nowrap}.plots{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:9px;padding:12px} .plot{height:440px;background:white;border:1px solid #dce3eb;border-radius:8px;min-width:0}.foot{padding:12px 26px 25px;color:#475569;line-height:1.55}.note{background:#fff6e6;padding:10px;border-radius:6px;margin-top:8px}@media(max-width:1000px){.plots{grid-template-columns:1fr}.plot{height:410px}}</style></head><body>
<header><h1>T3 neutrino-oscillation correlation explorer</h1><p>Optimised common Yukawas · 16 representations · S1/S2/S3 flavour scans</p></header>
<section class="controls"><label>Model <select id="model"></select></label><label><input id="s1" type="checkbox" checked>S1</label><label><input id="s2" type="checkbox" checked>S2</label><label><input id="s3" type="checkbox" checked>S3</label><label><input id="one" type="checkbox" checked>Gaussian 1σ ellipse</label><label><input id="three" type="checkbox" checked>Gaussian 3σ ellipse</label><label><input id="fit" type="checkbox" checked>Fit centre</label><button id="reset">Reset zoom</button></section>
<div class="plots"><div id="p1" class="plot"></div><div id="p2" class="plot"></div><div id="p3" class="plot"></div></div>
<div class="foot"><div id="status"></div><div class="note">Experimental regions are <b>diagonal-Gaussian approximations</b> from marginal NuFIT 6.1 (NO) uncertainties; they are <b>not official joint Δχ² contours</b>. Ellipses use (x/σx)²+(y/σy)²=1 or 9. No experimental δCP contour is assumed. Dashed segments connect available scan points; missing evaluations are not interpolated. The CP angle is periodic (0° ≡ 360°).</div><div id="fail"></div></div>
<script>/* PLOTLY */</script><script>
const DATA=/* DATA */;const FAIL=/* FAIL */;
const C={sin2_theta13:0.02248,sin2_theta23:0.47,delta_m3l_sq_ev2:0.002511};const S={sin2_theta13:0.00057,sin2_theta23:0.0155,delta_m3l_sq_ev2:0.0000205};
const palette={S1:'#087fc2',S2:'#ed6c30',S3:'#209452'};
const configs=[{id:'p1',x:'sin2_theta13',y:'delta_m3l_sq_ev2',xlabel:'sin² θ₁₃',ylabel:'Δm²₃₁ [eV²]'},{id:'p2',x:'sin2_theta23',y:'delta_m3l_sq_ev2',xlabel:'sin² θ₂₃',ylabel:'Δm²₃₁ [eV²]'},{id:'p3',x:'sin2_theta13',y:'delta_cp_deg',xlabel:'sin² θ₁₃',ylabel:'δCP [degrees]'}];
function ellipse(xkey,ykey,n){let xs=[],ys=[];for(let k=0;k<=120;k++){let a=2*Math.PI*k/120;xs.push(C[xkey]+n*S[xkey]*Math.cos(a));ys.push(C[ykey]+n*S[ykey]*Math.sin(a));}return {x:xs,y:ys,type:'scatter',mode:'lines',fill:'toself',line:{color:n===1?'#1479c3':'#82c7ec',width:1},fillcolor:n===1?'rgba(30,147,225,.25)':'rgba(95,188,240,.12)',hoverinfo:'skip',name:`Gaussian ${n}σ`,showlegend:false};}
function plot(reset=false){let model=document.getElementById('model').value;let all=DATA[model];configs.forEach(cfg=>{let traces=[];if(cfg.id!=='p3'){if(document.getElementById('three').checked)traces.push(ellipse(cfg.x,cfg.y,3));if(document.getElementById('one').checked)traces.push(ellipse(cfg.x,cfg.y,1));if(document.getElementById('fit').checked)traces.push({x:[C[cfg.x]],y:[C[cfg.y]],mode:'markers',type:'scatter',marker:{symbol:'cross',color:'#101827',size:10},name:'NuFIT centre',hovertemplate:'NuFIT centre<extra></extra>'});}
for(let d of ['S1','S2','S3']){if(!document.getElementById(d.toLowerCase()).checked)continue;let a=all[d];if(!a.length)continue;traces.push({x:a.map(r=>r[cfg.x]),y:a.map(r=>r[cfg.y]),customdata:a.map(r=>[r.t,r.chi2]),type:'scatter',mode:'lines+markers',name:d,line:{color:palette[d],width:2},marker:{color:palette[d],size:8},hovertemplate:d+'<br>t=%{customdata[0]:+.3f}<br>x=%{x:.6g}<br>y=%{y:.6g}<br>χ²=%{customdata[1]:.2f}<extra></extra>'});}
let layout={margin:{l:68,r:15,t:38,b:68},title:{text:cfg.ylabel+' vs '+cfg.xlabel,font:{size:15}},xaxis:{title:cfg.xlabel,zeroline:false,showgrid:true,gridcolor:'#e9edf2'},yaxis:{title:cfg.ylabel,zeroline:false,showgrid:true,gridcolor:'#e9edf2',...(cfg.id==='p3'?{range:[0,360]}:{})},legend:{orientation:'h',y:1.15},paper_bgcolor:'white',plot_bgcolor:'white',hovermode:'closest',uirevision:reset?String(Date.now()):'keep'};Plotly.react(cfg.id,traces,layout,{responsive:true,displaylogo:false});});
let available=Object.values(all).reduce((n,a)=>n+a.length,0);document.getElementById('status').textContent=`${model} — ${available}/9 successful scan points. Baseline at t=0 is shared by S1, S2 and S3.`;let f=FAIL.filter(r=>r.model===model);document.getElementById('fail').textContent=f.length?`Failed points: ${f.map(r=>r.direction+' t='+r.t+' ('+r.error+')').join('; ')}`:'';}
const select=document.getElementById('model');Object.keys(DATA).sort().forEach(n=>{let option=document.createElement('option');option.value=n;option.textContent=n;select.appendChild(option);});document.querySelectorAll('select,input').forEach(el=>el.addEventListener('change',()=>plot(true)));document.getElementById('reset').addEventListener('click',()=>plot(true));plot(true);
</script></body></html>'''
    page = page.replace('/* PLOTLY */', js).replace('/* DATA */', payload).replace('/* FAIL */', failed)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        raise FileExistsError(f'Refusing to overwrite existing display: {destination}')
    destination.write_text(page, encoding='utf-8')
    print(f'Created {destination} ({len(rows)} rows; {len(failures)} failed)')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--csv',type=Path,default=DEFAULT_CSV)
    p.add_argument('--output',type=Path,default=DEFAULT_OUTPUT)
    args=p.parse_args()
    build(args.csv,args.output)
if __name__=='__main__':main()
