"""Interactive correlations from saved T3 scan points (presentation only).

Rebuild plots and filter their display without rerunning matching, RGEs, or fits.
Gamma is a sensitivity factor, NOT a Yukawa coupling. A separate Yukawa-range
filter is available only when input CSV contains an explicit Yukawa quantity.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OBS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12', 'sin2_theta13', 'sin2_theta23')
LABELS = {'delta_m21_sq_ev2': 'Δm²₂₁ [eV²]', 'delta_m3l_sq_ev2': 'Δm²₃ℓ [eV²]',
          'sin2_theta12': 'sin²θ₁₂', 'sin2_theta13': 'sin²θ₁₃', 'sin2_theta23': 'sin²θ₂₃'}
DIRS = ('T1', 'T2', 'T3')
ELLIPSE_DELTAS = {'1sigma': 2.30, '3sigma': 11.83}
# Recognise explicit, named saved reference inputs only. Never infer Y from gamma.
YUKAWA_COLUMNS = ('reference_yukawa_target', 'reference_yukawa', 'yukawa_value')


def finite(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError, OverflowError):
        return None


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def gather(source: Path):
    """Preserve original point selection and deduplication, with optional Y metadata."""
    source = Path(source)
    points, target = [], None
    for path in sorted(source.rglob('sensitivity_diagnostics.json')):
        try:
            doc = read_json(path)
            if not isinstance(doc.get('directions'), dict):
                continue
            if target is None and isinstance(doc.get('target'), dict):
                target = doc['target']
            for direction in DIRS:
                for entry in doc['directions'].get(direction, {}).get('points', []):
                    pred = entry.get('prediction', {})
                    gamma = finite(entry.get('gamma'))
                    vals = {key: finite(pred.get(key)) for key in OBS}
                    if gamma is not None and gamma > 0 and all(v is not None for v in vals.values()):
                        points.append(dict(model=path.parent.name, direction=direction,
                                           gamma=gamma, **vals))
        except (OSError, ValueError, TypeError, KeyError):
            continue
    for path in sorted(source.rglob('scan_points.csv')):
        with path.open(newline='', encoding='utf-8-sig') as stream:
            for row in csv.DictReader(stream):
                if row.get('status') != 'Success' or row.get('direction') not in DIRS:
                    continue
                gamma = finite(row.get('gamma'))
                vals = {key: finite(row.get(key)) for key in OBS}
                if gamma is None or gamma <= 0 or any(v is None for v in vals.values()):
                    continue
                p = dict(model=row.get('model') or path.parent.name,
                         direction=row['direction'], gamma=gamma, **vals)
                for field in YUKAWA_COLUMNS:
                    yref = finite(row.get(field))
                    if yref is not None:
                        p['yukawa'] = yref
                        break
                points.append(p)
    unique = {}
    for p in points:
        key = (p['model'], p['direction'], round(p['gamma'], 12))
        if key not in unique or ('yukawa' in p and 'yukawa' not in unique[key]):
            unique[key] = p
    return sorted(unique.values(), key=lambda p: (p['model'], p['direction'], p['gamma'])), target


def target_payload(target):
    if not isinstance(target, dict):
        return {}
    names = list(target.get('observable_names') or [])
    means = target.get('central_values') or []
    covariance = target.get('covariance') or []
    inds = [i for i, k in enumerate(names) if k in OBS and i < len(means)
            and finite(means[i]) is not None]
    if not inds:
        return {}
    mat = []
    for i in inds:
        row = []
        for j in inds:
            try:
                row.append(finite(covariance[i][j]))
            except (IndexError, TypeError):
                row.append(None)
        mat.append(row)
    return {'observable_names': [names[i] for i in inds],
            'central_values': [float(means[i]) for i in inds],
            'covariance': mat, 'source': target.get('source', ''),
            'note': target.get('note', ''), 'ellipse_delta_chi2': ELLIPSE_DELTAS}


HTML = r'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>T3 interactive correlations</title><style>
:root{color-scheme:light;--ink:#17243b;--muted:#556375;--line:#d7dee8;--bg:#f7f9fc}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:13px system-ui,Arial}
header{padding:15px 22px;background:#fff;border-bottom:1px solid var(--line)}h2{margin:0 0 5px}
.small{color:var(--muted);font-size:12px;line-height:1.5}.layout{display:grid;grid-template-columns:270px minmax(0,1fr)}
aside{background:white;border-right:1px solid var(--line);padding:16px;max-height:calc(100vh - 70px);overflow:auto}
main{min-width:0;padding:16px}.panel{background:#fff;border:1px solid var(--line);border-radius:8px;padding:14px}
select,input[type=number]{width:100%;min-width:0;margin:4px 0 9px;padding:6px;border:1px solid #b8c5d5;border-radius:5px}
button{padding:7px 10px;margin:3px;border:1px solid #c4cddb;border-radius:5px;background:white;cursor:pointer}
button:hover{background:#eef4ff}label{font-size:12px}.choice{display:block;margin:5px 0}
.grid{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:8px}.pair{display:grid;grid-template-columns:1fr 1fr;gap:8px}
.toolbar{display:flex;flex-wrap:wrap;align-items:center;gap:3px;margin-bottom:9px}
canvas{width:100%;height:min(67vh,720px);min-height:360px;display:block}
#caption{padding:10px;background:#f4f7fb;font-size:12px;line-height:1.5}
@media(max-width:850px){.layout{grid-template-columns:1fr}aside{max-height:none;border-bottom:1px solid var(--line)}main{padding:8px}.grid{grid-template-columns:repeat(2,minmax(0,1fr))}}
</style></head><body><header><h2>T3 observable correlations</h2>
<div class="small">Saved sensitivity trajectories only. Filters change presentation, not numerical data.</div></header>
<div class="layout"><aside><label>X observable<select id="x"></select></label><label>Y observable<select id="y"></select></label>
<button id="all">Show all</button><button id="none">Hide all</button>
<h4>Gamma (sensitivity factor)</h4><div class="pair"><label>Minimum γ<input type="number" step="any" id="gammaMin"></label>
<label>Maximum γ<input type="number" step="any" id="gammaMax"></label></div>
<button id="gammaReset">Full gamma range</button><p class="small" id="gammaStatus"></p>
<div id="yukawaSection" hidden><h4>Saved reference Yukawa</h4>
<div class="pair"><label>Minimum Y<input type="number" step="any" id="yukawaMin"></label>
<label>Maximum Y<input type="number" step="any" id="yukawaMax"></label></div>
<button id="yukawaReset">Full Yukawa range</button><p class="small">Uses an explicit saved reference Yukawa value, not γ. Points without this metadata are excluded while a custom Y filter is active.</p></div>
<h4>Directions</h4><div id="dirs"></div><h4>Models</h4><div id="models"></div></aside>
<main><div class="panel"><div class="toolbar"><button id="reset">Fit visible points</button>
<button id="xIn">X +</button><button id="xOut">X −</button><button id="yIn">Y +</button><button id="yOut">Y −</button>
<label class="choice"><input type="checkbox" id="regions" checked> Gaussian 1σ / 3σ</label></div>
<div class="grid"><label>X minimum<input type="number" step="any" id="xmin"></label><label>X maximum<input type="number" step="any" id="xmax"></label>
<label>Y minimum<input type="number" step="any" id="ymin"></label><label>Y maximum<input type="number" step="any" id="ymax"></label></div>
<div class="toolbar"><button id="applyAxes">Apply axis limits</button><button id="clearAxes">Reset axis limits</button>
<span class="small" id="axisError" aria-live="polite"></span></div>
<canvas id="plot"></canvas><div id="caption"></div><p class="small" id="status"></p></div></main></div>
<script>'use strict';
const DATA=__DATA__,OBS=__OBS__,LABELS=__LABELS__;
const $=id=>document.getElementById(id),canvas=$('plot'),ctx=canvas.getContext('2d');
const COLORS={T1:'#2376c9',T2:'#ed713d',T3:'#269a57'};
const models=[...new Set(DATA.points.map(p=>p.model))].sort();
let enabled=Object.fromEntries(models.map(m=>[m,true]));let dirs={T1:true,T2:true,T3:true};
let x=OBS[3],y=OBS[1],manual=null,axesError='';
const gammas=DATA.points.map(p=>p.gamma).filter(v=>Number.isFinite(v)&&v>0);
const gFull=[Math.min(...gammas),Math.max(...gammas)];let gMin=gFull[0],gMax=gFull[1];
const yvals=DATA.points.map(p=>p.yukawa).filter(Number.isFinite);
const yFull=yvals.length?[Math.min(...yvals),Math.max(...yvals)]:null;
let yMin=yFull?.[0],yMax=yFull?.[1],yCustom=false;
function fillOptions(){for(let id of ['x','y']){for(let k of OBS){let opt=document.createElement('option');opt.value=k;opt.textContent=LABELS[k];$(id).append(opt)}$(id).value=id==='x'?x:y;$(id).onchange=()=>{x=$('x').value;y=$('y').value;manual=null;draw()}}}
fillOptions();
for(let d of ['T1','T2','T3']){let label=document.createElement('label');label.className='choice';let cb=document.createElement('input');cb.type='checkbox';cb.checked=true;cb.onchange=()=>{dirs[d]=cb.checked;manual=null;draw()};label.append(cb,document.createTextNode(' '+d));$('dirs').append(label)}
for(let m of models){let label=document.createElement('label');label.className='choice';let cb=document.createElement('input');cb.type='checkbox';cb.checked=true;cb.onchange=()=>{enabled[m]=cb.checked;manual=null;draw()};label.append(cb,document.createTextNode(' '+m));$('models').append(label)}
$('all').onclick=()=>{models.forEach(m=>enabled[m]=true);$('models').querySelectorAll('input').forEach(c=>c.checked=true);manual=null;draw()};
$('none').onclick=()=>{models.forEach(m=>enabled[m]=false);$('models').querySelectorAll('input').forEach(c=>c.checked=false);manual=null;draw()};
$('gammaMin').value=gMin;$('gammaMax').value=gMax;
function setGamma(){let a=Number($('gammaMin').value),b=Number($('gammaMax').value);if(!(a>0&&b>0&&a<=b)){ $('gammaStatus').textContent='Require 0 < minimum ≤ maximum.';return}gMin=a;gMax=b;manual=null;draw()}
$('gammaMin').onchange=setGamma;$('gammaMax').onchange=setGamma;
$('gammaReset').onclick=()=>{$('gammaMin').value=gFull[0];$('gammaMax').value=gFull[1];setGamma()};
if(yFull){$('yukawaSection').hidden=false;$('yukawaMin').value=yMin;$('yukawaMax').value=yMax;
  function setY(){const a=Number($('yukawaMin').value),b=Number($('yukawaMax').value);if(!Number.isFinite(a)||!Number.isFinite(b)||a>b)return;yMin=a;yMax=b;yCustom=true;manual=null;draw()}
  $('yukawaMin').onchange=setY;$('yukawaMax').onchange=setY;
  $('yukawaReset').onclick=()=>{$('yukawaMin').value=yFull[0];$('yukawaMax').value=yFull[1];yCustom=false;manual=null;draw()};}
function selected(){return DATA.points.filter(p=>enabled[p.model]&&dirs[p.direction]&&p.gamma>=gMin-1e-12*Math.abs(gMin)&&p.gamma<=gMax+1e-12*Math.abs(gMax)&&(!yCustom||(Number.isFinite(p.yukawa)&&p.yukawa>=yMin&&p.yukawa<=yMax)))}
function pair(){const t=DATA.target||{},names=t.observable_names||[],v=t.central_values||[],c=t.covariance||[];let i=names.indexOf(x),j=names.indexOf(y);if(i<0||j<0)return null;
 const mx=Number(v[i]),my=Number(v[j]),a=Number(c[i]?.[i]),b=Number(c[i]?.[j]),d=Number(c[j]?.[j]);
 if(![mx,my,a,b,d].every(Number.isFinite)||a<=0||d<=0||a*d-b*b<=0)return null;return {mx,my,a,b,d}}
function ellipse(g,delta){const diff=g.a-g.d,disc=Math.sqrt(Math.max(0,diff*diff+4*g.b*g.b)),l1=(g.a+g.d+disc)/2,l2=(g.a+g.d-disc)/2;if(l2<=0)return null;return {cx:g.mx,cy:g.my,rx:Math.sqrt(delta*l1),ry:Math.sqrt(delta*l2),angle:.5*Math.atan2(2*g.b,diff)}}
function ellipseBox(e){const c=Math.cos(e.angle),s=Math.sin(e.angle),dx=Math.hypot(e.rx*c,e.ry*s),dy=Math.hypot(e.rx*s,e.ry*c);return [e.cx-dx,e.cx+dx,e.cy-dy,e.cy+dy]}
function naturalBounds(){const pts=selected(),xs=pts.map(p=>p[x]),ys=pts.map(p=>p[y]);const g=pair();if(g&&$('regions').checked){for(let delta of [2.30,11.83]){let e=ellipse(g,delta);if(e){const b=ellipseBox(e);xs.push(b[0],b[1]);ys.push(b[2],b[3])}}}
 if(!xs.length)return [0,1,0,1];let lo=Math.min(...xs),hi=Math.max(...xs),bot=Math.min(...ys),top=Math.max(...ys),dx=hi-lo||Math.max(Math.abs(hi)*.05,1e-12),dy=top-bot||Math.max(Math.abs(top)*.05,1e-12);return [lo-.08*dx,hi+.08*dx,bot-.08*dy,top+.08*dy]}
function updateAxisFields(b){for(let [id,v] of [['xmin',b[0]],['xmax',b[1]],['ymin',b[2]],['ymax',b[3]]])$(id).value=Number(v).toPrecision(9)}
function draw(){const b=manual||naturalBounds();updateAxisFields(b);let [lo,hi,bot,top]=b;
 const w=Math.max(canvas.clientWidth,420),h=Math.max(canvas.clientHeight,380),dpr=window.devicePixelRatio||1;canvas.width=Math.round(w*dpr);canvas.height=Math.round(h*dpr);ctx.setTransform(dpr,0,0,dpr,0,0);
 const L=90,R=23,T=25,B=66,pw=w-L-R,ph=h-T-B,X=v=>L+(v-lo)/(hi-lo)*pw,Y=v=>T+ph-(v-bot)/(top-bot)*ph;
 ctx.clearRect(0,0,w,h);ctx.fillStyle='white';ctx.fillRect(L,T,pw,ph);
 ctx.strokeStyle='#dfe5ee';ctx.fillStyle='#566174';ctx.font='11px system-ui';for(let i=0;i<=5;i++){let xx=lo+(hi-lo)*i/5,yy=bot+(top-bot)*i/5;ctx.beginPath();ctx.moveTo(X(xx),T);ctx.lineTo(X(xx),T+ph);ctx.moveTo(L,Y(yy));ctx.lineTo(L+pw,Y(yy));ctx.stroke();ctx.textAlign='center';ctx.fillText(xx.toPrecision(4),X(xx),T+ph+19);ctx.textAlign='right';ctx.fillText(yy.toPrecision(4),L-6,Y(yy)+4)}
 ctx.save();ctx.beginPath();ctx.rect(L,T,pw,ph);ctx.clip();const g=pair();if(g&&$('regions').checked){for(let [delta,color] of [[11.83,'#d3eafa'],[2.30,'#9dd1f5']]){let e=ellipse(g,delta);if(!e)continue;ctx.beginPath();for(let k=0;k<=128;k++){let angle=2*Math.PI*k/128,c=Math.cos(angle),s=Math.sin(angle);let px=e.cx+e.rx*c*Math.cos(e.angle)-e.ry*s*Math.sin(e.angle),py=e.cy+e.rx*c*Math.sin(e.angle)+e.ry*s*Math.cos(e.angle);k?ctx.lineTo(X(px),Y(py)):ctx.moveTo(X(px),Y(py))}ctx.closePath();ctx.fillStyle=color;ctx.fill()}}
 const groups=new Map();for(let p of selected()){let key=p.model+'|'+p.direction;if(!groups.has(key))groups.set(key,[]);groups.get(key).push(p)}
 for(let [k,pts] of groups){pts.sort((a,b)=>a.gamma-b.gamma);let d=k.split('|').pop();ctx.beginPath();pts.forEach((p,i)=>i?ctx.lineTo(X(p[x]),Y(p[y])):ctx.moveTo(X(p[x]),Y(p[y])));ctx.strokeStyle=COLORS[d];ctx.lineWidth=1.5;ctx.stroke();for(let p of pts){ctx.beginPath();ctx.arc(X(p[x]),Y(p[y]),Math.abs(p.gamma-1)<1e-8?3.5:2,0,2*Math.PI);ctx.fillStyle=COLORS[d];ctx.fill()}}
 ctx.restore();ctx.strokeStyle='#556375';ctx.strokeRect(L,T,pw,ph);ctx.fillStyle='#17243b';ctx.font='13px system-ui';ctx.textAlign='center';ctx.fillText(LABELS[x],L+pw/2,h-16);ctx.save();ctx.translate(19,T+ph/2);ctx.rotate(-Math.PI/2);ctx.fillText(LABELS[y],0,0);ctx.restore();
 $('gammaStatus').textContent=selected().length+' visible points (saved range '+gFull[0].toPrecision(3)+'–'+gFull[1].toPrecision(3)+').';
 $('caption').textContent=groups.size+' model–direction trajectories, '+selected().length+' displayed saved points. Blue regions are 2D Gaussian 1σ/3σ contours from supplied covariance (if available). Filters do not delete data.';
 $('status').textContent='Source: '+DATA.source+' | '+(DATA.target.note||'No experimental target metadata');}
function axisZoom(axis,z){let b=(manual||naturalBounds()).slice(),i=axis==='x'?0:2,m=(b[i]+b[i+1])/2;b[i]=m+(b[i]-m)*z;b[i+1]=m+(b[i+1]-m)*z;manual=b;draw()}
$('xIn').onclick=()=>axisZoom('x',.65);$('xOut').onclick=()=>axisZoom('x',1/.65);
$('yIn').onclick=()=>axisZoom('y',.65);$('yOut').onclick=()=>axisZoom('y',1/.65);
$('reset').onclick=$('clearAxes').onclick=()=>{manual=null;draw()};
$('applyAxes').onclick=()=>{let b=['xmin','xmax','ymin','ymax'].map(id=>$(id).value.trim()===''?NaN:Number($(id).value));if(!b.every(Number.isFinite)||b[0]>=b[1]||b[2]>=b[3]){$('axisError').textContent='Require finite min < max for both axes.';return}$('axisError').textContent='';manual=b;draw()};
$('regions').onchange=()=>{if(!manual)draw();else draw()};
canvas.addEventListener('wheel',e=>{e.preventDefault();let b=(manual||naturalBounds()).slice(),z=e.deltaY>0?1.2:0.8;for(let i of [0,2]){let m=(b[i]+b[i+1])/2;b[i]=m+(b[i]-m)*z;b[i+1]=m+(b[i+1]-m)*z}manual=b;draw()},{passive:false});
window.addEventListener('resize',draw);draw();
</script></body></html>'''


def _fallback_target_from_config(config_path: Path):
    """Diagonal Gaussian approximation; no off-diagonal NuFIT information inferred."""
    doc = read_json(config_path)
    centers, errors = doc.get('central_values', {}), doc.get('one_sigma_errors', {})
    names = [k for k in OBS if k in centers and k in errors]
    if not names:
        raise ValueError(f'No supported NuFIT parameters in {config_path}')
    means = [float(centers[k]) for k in names]
    sigmas = [float(errors[k]) for k in names]
    if any(not math.isfinite(s) or s <= 0 for s in sigmas):
        raise ValueError('NuFIT marginal errors must be finite and positive')
    covariance = [[sigmas[i] ** 2 if i == j else 0.0 for j in range(len(names))]
                  for i in range(len(names))]
    return {'observable_names': names, 'central_values': means, 'covariance': covariance,
            'source': doc.get('source', str(config_path)),
            'note': 'Diagonal Gaussian approximation from marginal errors; not official joint contours.'}


def build(source: Path, output: Path, target_config: Path | None = None):
    points, target = gather(Path(source))
    if target is None and target_config is not None:
        target = _fallback_target_from_config(Path(target_config))
    if not points:
        raise RuntimeError(f'No saved successful T3 sensitivity points under {source}')
    payload = {'points': points, 'target': target_payload(target), 'source': str(source)}
    body = HTML.replace('__DATA__', json.dumps(payload, ensure_ascii=False, allow_nan=False).replace('</', '<\\/'))
    body = body.replace('__LABELS__', json.dumps(LABELS, ensure_ascii=False)).replace('__OBS__', json.dumps(OBS))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(body, encoding='utf-8')
    print(f'Wrote {output}: {len(points)} saved successful points')


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path,
                        default=PROJECT_ROOT / 'Reports/output/full/fitted_sensitivity')
    parser.add_argument('--output', type=Path,
                        default=PROJECT_ROOT / 'Reports/output/display/interactive_correlations.html')
    parser.add_argument('--target-config', type=Path,
                        default=PROJECT_ROOT / 'configs/nufit_6_1_ic24_no.json')
    args = parser.parse_args(argv)
    build(args.source, args.output, args.target_config if args.target_config.is_file() else None)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
