"""Self-contained interactive four-scenario curves from saved diagnostics only.

The intermediate hard+direct line is a checked display reconstruction, NOT an
intermediate-EFT Wilson coefficient. No matching or running is performed.
"""
from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Mapping
import numpy as np

from Numerical.plotting.RunningTrajectoryData import SCENARIOS, QUANTITIES, SavedRun, extract_trajectory
from Numerical.plotting.WeinbergMatchedContinuation import reconstruct_matched_continuation

LABELS = {
    'smallY_smallL': 'Y=0.005, λref=0.1',
    'smallY_largeL': 'Y=0.005, λref=1.0',
    'largeY_smallL': 'Y=0.5, λref=0.1',
    'largeY_largeL': 'Y=0.5, λref=1.0',
}


def build_interactive_data(runs: Mapping[str, SavedRun]) -> dict:
    if not runs:
        raise ValueError('At least one saved run is required')
    result = {'model': next(iter(runs.values())).model_key, 'scenarios': {},
              'quantities': {name: {'label': spec[3], 'traces': {}} for name, spec in QUANTITIES.items()},
              'warnings': []}
    result['quantities']['weinberg_eft'] = {'label': '||C₅||F [GeV⁻¹] — final and checked display continuation', 'traces': {}}
    for scenario in SCENARIOS:
        if scenario not in runs:
            continue
        payload = runs[scenario].payload
        result['scenarios'][scenario] = LABELS[scenario]
        for quantity in QUANTITIES:
            mu, values = extract_trajectory(payload, quantity)
            result['quantities'][quantity]['traces'][scenario] = [
                {'x': mu.tolist(), 'y': values.tolist(), 'part': 'stored'}]
        mu, values = extract_trajectory(payload, 'c5')
        curves = [{'x': mu.tolist(), 'y': values.tolist(), 'part': 'final SM+C5'}]
        try:
            matched = reconstruct_matched_continuation(payload)
        except (KeyError, ValueError, TypeError, IndexError) as exc:
            result['warnings'].append(f'{scenario}: matched continuation unavailable ({exc})')
        else:
            norm = np.linalg.norm(matched.combined.reshape(len(matched.mu_gev), 9), axis=1)
            curves.append({'x': matched.mu_gev.tolist(), 'y': norm.tolist(),
                           'part': 'display reconstruction (not EFT Wilson coefficient)'})
        result['quantities']['weinberg_eft']['traces'][scenario] = curves
    return result


def write_interactive_comparison(runs: Mapping[str, SavedRun], output_path: Path) -> tuple[Path, list[str]]:
    data = build_interactive_data(runs)
    # Escape characters with special meaning in HTML script elements.
    payload = json.dumps(data, ensure_ascii=True, allow_nan=False).replace('<', '\\u003c')
    choices = '\n'.join(f'<option value="{html.escape(k)}">{html.escape(k)} — {html.escape(v["label"])}</option>'
                        for k, v in data['quantities'].items())
    page = r'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Interactive T3 within-model comparison</title>
<style>
body{font-family:system-ui,Arial,sans-serif;color:#20252d;max-width:1200px;margin:2rem auto;padding:0 1rem}
.controls{display:flex;flex-wrap:wrap;align-items:center;gap:1rem;margin:1rem 0}
select,button{padding:.45rem;font:inherit}fieldset{border:1px solid #bbb;border-radius:6px}
#plot{width:100%;height:560px;border:1px solid #d4d6da;touch-action:none;cursor:crosshair}
.legend{display:flex;flex-wrap:wrap;gap:1rem}label{white-space:nowrap}
.note{font-size:.9rem;color:#495466}#tooltip{font-size:.85rem;min-height:1.5rem} .warning{color:#8b4000}
</style></head><body>
<h1 id="heading"></h1><p>Four benchmark scenarios for the same T3 representation. All curves are read from saved numerical diagnostics; no physics is recomputed.</p>
<div class="controls"><label for="quantity">Quantity</label><select id="quantity">__CHOICES__</select>
<label for="yaxis">Vertical scale</label><select id="yaxis"><option value="auto">Auto</option><option value="log">Logarithmic</option><option value="linear">Linear</option></select>
<button id="reset">Reset zoom</button></div><fieldset><legend>Scenarios</legend><div id="legend" class="legend"></div></fieldset>
<canvas id="plot" aria-label="Interactive running trajectories"></canvas><div id="tooltip"></div>
<p class="note">Drag horizontally to zoom the renormalisation-scale range; double-click or use Reset zoom to restore it. Hover for numerical values. Colours distinguish Yukawa magnitudes; solid/dashed lines distinguish reference scalar couplings. Dotted lines in the combined Weinberg view are matched <em>display reconstructions</em>, not intermediate-EFT Wilson coefficients. Intermediate direct LLSS → Weinberg generation is a separate selectable quantity.</p>
<p class="warning" id="warnings"></p>
<script>const DATA=__DATA__;
const colors={smallY_smallL:'#1766bc',smallY_largeL:'#1766bc',largeY_smallL:'#d36d16',largeY_largeL:'#d36d16'};
const canvas=document.getElementById('plot'),ctx=canvas.getContext('2d');
const q=document.getElementById('quantity'),ys=document.getElementById('yaxis'),tip=document.getElementById('tooltip');
const enabled={};Object.keys(DATA.scenarios).forEach(s=>enabled[s]=true);
document.getElementById('heading').textContent=DATA.model+' — interactive comparison';
document.getElementById('warnings').textContent=DATA.warnings.join(' | ');
const legend=document.getElementById('legend');
Object.entries(DATA.scenarios).forEach(([s,label])=>{const el=document.createElement('label');const input=document.createElement('input');input.type='checkbox';input.checked=true;input.addEventListener('change',()=>{enabled[s]=input.checked;draw()});el.append(input,document.createTextNode(' '+label));el.style.color=colors[s];legend.appendChild(el)});
let zoom=null,drag=null,plotRect=null,lastBounds=null;
function traces(){const spec=DATA.quantities[q.value];return Object.entries(spec.traces).filter(([s])=>enabled[s]).flatMap(([s,lines])=>lines.map(line=>({...line,scenario:s})));}
function dimensions(){const box=canvas.getBoundingClientRect();const pixelRatio=window.devicePixelRatio||1;const w=Math.max(300,box.width),h=Math.max(250,box.height);canvas.width=Math.round(w*pixelRatio);canvas.height=Math.round(h*pixelRatio);ctx.setTransform(pixelRatio,0,0,pixelRatio,0,0);return {w,h};}
function boundsOf(series){let X=[],Y=[];for(const c of series){for(let i=0;i<c.x.length;i++){let x=c.x[i],y=c.y[i];if(Number.isFinite(x)&&x>0&&Number.isFinite(y)){X.push(Math.log10(x));Y.push(y)}}}
if(!X.length)return null;const xmin=zoom?zoom[0]:Math.min(...X),xmax=zoom?zoom[1]:Math.max(...X);const visible=[];
for(const c of series)for(let i=0;i<c.x.length;i++){const x=Math.log10(c.x[i]),y=c.y[i];if(x>=xmin&&x<=xmax&&Number.isFinite(y))visible.push(y)}
let requested=ys.value,logY=requested==='log'||(requested==='auto'&&visible.length>0&&visible.every(v=>v>0)&&Math.max(...visible)/Math.max(Math.min(...visible),1e-100)>20);
let val=visible.filter(v=>!logY||v>0).map(v=>logY?Math.log10(v):v);if(!val.length)return null;let lo=Math.min(...val),hi=Math.max(...val);if(lo===hi){lo-=.5;hi+=.5}else{const pad=(hi-lo)*.09;lo-=pad;hi+=pad}return {xmin,xmax,ymin:lo,ymax:hi,logY};}
function draw(){const {w,h}=dimensions();ctx.clearRect(0,0,w,h);ctx.fillStyle='#fff';ctx.fillRect(0,0,w,h);const margin={left:82,right:22,top:24,bottom:64};const rect={x:margin.left,y:margin.top,w:w-margin.left-margin.right,h:h-margin.top-margin.bottom};plotRect=rect;
const series=traces(),b=boundsOf(series);lastBounds=b;
if(!b){ctx.fillStyle='#444';ctx.fillText('No valid values in selected range. Try Reset zoom or Auto scale.',50,100);return}
const px=x=>rect.x+(Math.log10(x)-b.xmin)/(b.xmax-b.xmin||1)*rect.w,py=y=>rect.y+rect.h-((b.logY?Math.log10(y):y)-b.ymin)/(b.ymax-b.ymin)*rect.h;
ctx.strokeStyle='#d9dde3';ctx.lineWidth=1;ctx.fillStyle='#47505d';ctx.font='12px system-ui';
for(let i=0;i<=5;i++){let u=i/5,x=rect.x+rect.w*u,y=rect.y+rect.h*(1-u);ctx.beginPath();ctx.moveTo(x,rect.y);ctx.lineTo(x,rect.y+rect.h);ctx.moveTo(rect.x,y);ctx.lineTo(rect.x+rect.w,y);ctx.stroke();let xt=10**(b.xmin+(b.xmax-b.xmin)*u),yt=b.ymin+(b.ymax-b.ymin)*u;if(b.logY)yt=10**yt;ctx.fillText(xt.toExponential(1),x-20,rect.y+rect.h+19);ctx.fillText(yt.toExponential(1),5,y+4)}
ctx.save();ctx.beginPath();ctx.rect(rect.x,rect.y,rect.w,rect.h);ctx.clip();
series.forEach(c=>{ctx.beginPath();ctx.strokeStyle=colors[c.scenario];ctx.lineWidth=2;ctx.setLineDash(c.part.startsWith('display')?[2,5]:c.scenario.endsWith('largeL')?[8,5]:[]);let started=false;for(let i=0;i<c.x.length;i++){let x=c.x[i],y=c.y[i];if(!Number.isFinite(x)||x<=0||!Number.isFinite(y)||(b.logY&&y<=0)){started=false;continue}const xx=px(x),yy=py(y);if(!started){ctx.moveTo(xx,yy);started=true}else ctx.lineTo(xx,yy)}ctx.stroke()});ctx.setLineDash([]);ctx.restore();
ctx.fillStyle='#222';ctx.textAlign='center';ctx.fillText('Renormalisation scale μ [GeV] (logarithmic)',rect.x+rect.w/2,h-12);ctx.save();ctx.translate(16,rect.y+rect.h/2);ctx.rotate(-Math.PI/2);ctx.fillText(DATA.quantities[q.value].label,0,0);ctx.restore();ctx.textAlign='left';
if(drag){ctx.fillStyle='rgba(50,85,175,.12)';ctx.fillRect(Math.min(drag.start,drag.end),rect.y,Math.abs(drag.end-drag.start),rect.h)}}
function position(e){const r=canvas.getBoundingClientRect();return e.clientX-r.left}
canvas.addEventListener('pointerdown',e=>{if(!plotRect||!lastBounds)return;let x=position(e);if(x<plotRect.x||x>plotRect.x+plotRect.w)return;drag={start:x,end:x};canvas.setPointerCapture(e.pointerId)});
canvas.addEventListener('pointermove',e=>{const x=position(e);if(drag){drag.end=Math.max(plotRect.x,Math.min(plotRect.x+plotRect.w,x));draw();return}if(!lastBounds||!plotRect||x<plotRect.x||x>plotRect.x+plotRect.w){tip.textContent='';return}const mu=10**(lastBounds.xmin+(x-plotRect.x)/plotRect.w*(lastBounds.xmax-lastBounds.xmin));let closest=[];for(const c of traces()){let idx=0,best=Infinity;for(let i=0;i<c.x.length;i++){const d=Math.abs(Math.log10(c.x[i]/mu));if(d<best){best=d;idx=i}}if(best<.09)closest.push(`${DATA.scenarios[c.scenario]} (${c.part}): μ=${c.x[idx].toExponential(3)}, value=${c.y[idx].toExponential(5)}`)}tip.textContent=closest.join(' | ')});
canvas.addEventListener('pointerup',e=>{if(!drag||!lastBounds)return;const d=drag;drag=null;if(Math.abs(d.end-d.start)>8){const x1=lastBounds.xmin+(Math.min(d.start,d.end)-plotRect.x)/plotRect.w*(lastBounds.xmax-lastBounds.xmin);const x2=lastBounds.xmin+(Math.max(d.start,d.end)-plotRect.x)/plotRect.w*(lastBounds.xmax-lastBounds.xmin);zoom=[x1,x2]}draw()});
canvas.addEventListener('dblclick',()=>{zoom=null;draw()});document.getElementById('reset').addEventListener('click',()=>{zoom=null;draw()});q.addEventListener('change',()=>{zoom=null;draw()});ys.addEventListener('change',draw);window.addEventListener('resize',draw);draw();
</script></body></html>'''
    page = page.replace('__CHOICES__', choices).replace('__DATA__', payload)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(page, encoding='utf-8')
    return output_path, data['warnings']
