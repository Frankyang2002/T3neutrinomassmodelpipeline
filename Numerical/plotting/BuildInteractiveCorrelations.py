"""Presentation-only interactive observable correlations from SAVED T3 scan points.

No RGEs, matching, fits or sensitivity evaluations are performed. Intended to
run after the existing numerical pipeline; missing scans remain missing.

This version upgrades the experimental overlay from marginal rectangles to
Gaussian 2D confidence ellipses whenever the saved target covariance is
available.  The ellipses are derived from the relevant 2x2 covariance submatrix
and use the standard 2-parameter confidence levels

    Delta-chi2 = 2.30  (approximately 68.27%, often called 1 sigma)
    Delta-chi2 = 11.83 (approximately 99.73%, often called 3 sigma)

If no usable covariance is stored, the plot falls back to no ellipse overlay.

Usage (project root):
 python -m Numerical.plotting.BuildInteractiveCorrelations
 python -m Numerical.plotting.BuildInteractiveCorrelations --source Reports/output/full/fitted_sensitivity
 python -m Numerical.plotting.BuildInteractiveCorrelations --source Reports/output/full/common_yukawa_correlation_scans
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
OBS = ('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23')
LABELS = {
    'delta_m21_sq_ev2':'Δm²₂₁ [eV²]',
    'delta_m3l_sq_ev2':'Δm²₃ℓ [eV²]',
    'sin2_theta12':'sin²θ₁₂',
    'sin2_theta13':'sin²θ₁₃',
    'sin2_theta23':'sin²θ₂₃',
}
DIRS = ('T1','T2','T3')
ELLIPSE_DELTAS = {'1sigma': 2.30, '3sigma': 11.83}
ELLIPSE_LEVEL_ORDER = ('3sigma', '1sigma')
ELLIPSE_LEVEL_STYLES = {
    '3sigma': {'fill': '#d3eafa', 'label': 'Gaussian 3σ contour'},
    '1sigma': {'fill': '#9dd1f5', 'label': 'Gaussian 1σ contour'},
}


def finite(v):
    try:
        z = float(v)
        return z if math.isfinite(z) else None
    except (TypeError, ValueError, OverflowError):
        return None


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def gather(source: Path):
    """Read either saved fitted-sensitivity JSON or checkpointed common-Yukawa CSV."""
    items = []
    target = None
    for path in sorted(source.rglob('sensitivity_diagnostics.json')):
        try:
            doc = read_json(path)
            if not isinstance(doc.get('directions'), dict):
                continue
            if target is None and isinstance(doc.get('target'), dict):
                target = doc['target']
            model = path.parent.name
            for d in DIRS:
                points = doc['directions'].get(d, {}).get('points', [])
                for p in points:
                    pred = p.get('prediction', {})
                    gamma = finite(p.get('gamma'))
                    values = {k: finite(pred.get(k)) for k in OBS}
                    if gamma is not None and all(v is not None for v in values.values()):
                        items.append(dict(model=model, direction=d, gamma=gamma, **values))
        except (OSError, ValueError, TypeError, KeyError):
            continue
    for path in sorted(source.rglob('scan_points.csv')):
        with path.open(newline='', encoding='utf-8-sig') as stream:
            for r in csv.DictReader(stream):
                if r.get('status') != 'Success' or r.get('direction') not in DIRS:
                    continue
                gamma = finite(r.get('gamma'))
                values = {k: finite(r.get(k)) for k in OBS}
                if gamma is not None and all(v is not None for v in values.values()):
                    items.append(dict(model=r.get('model') or path.parent.name,
                                      direction=r['direction'], gamma=gamma, **values))
    unique = {}
    for p in items:
        unique.setdefault((p['model'], p['direction'], round(p['gamma'], 12)), p)
    ordered = sorted(unique.values(), key=lambda p: (p['model'], p['direction'], p['gamma']))
    return ordered, target


def target_payload(target):
    """Persist only the verified Gaussian-target metadata needed for plotting."""
    if not isinstance(target, dict):
        return {}
    names = list(target.get('observable_names') or [])
    means = target.get('central_values') or []
    covariance = target.get('covariance') or []
    filtered_names = []
    filtered_means = []
    filtered_covariance = []
    indices = []
    for i, name in enumerate(names):
        if name not in OBS:
            continue
        mean = finite(means[i] if i < len(means) else None)
        if mean is None:
            continue
        indices.append(i)
        filtered_names.append(name)
        filtered_means.append(mean)
    if not filtered_names:
        return {}
    for i in indices:
        row = []
        for j in indices:
            try:
                value = finite(covariance[i][j])
            except (IndexError, TypeError):
                value = None
            row.append(value)
        filtered_covariance.append(row)
    return {
        'observable_names': filtered_names,
        'central_values': filtered_means,
        'covariance': filtered_covariance,
        'source': target.get('source', ''),
        'note': target.get('note', ''),
        'ellipse_delta_chi2': ELLIPSE_DELTAS,
    }


HTML = r'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>T3 interactive neutrino-observable correlations</title><style>
:root{--ink:#17243b;--muted:#556375;--line:#d7dee8;--surface:#f7f9fc;--accent:#225ba8}
*{box-sizing:border-box}body{font-family:system-ui,-apple-system,Segoe UI,Arial,sans-serif;margin:0;background:var(--surface);color:var(--ink)}
header{padding:20px 28px;background:white;border-bottom:1px solid var(--line)}header h2{font-size:20px;margin:0 0 5px}
.small{font-size:12px;color:var(--muted);line-height:1.5}.layout{display:grid;grid-template-columns:275px minmax(0,1fr)}
aside{background:white;border-right:1px solid var(--line);padding:20px 17px}main{padding:20px 22px;min-width:0}
select{width:100%;padding:9px;border:1px solid #bac7d8;border-radius:7px;background:white;margin:5px 0 14px}
button{padding:8px 12px;border:1px solid #ccd5e2;border-radius:7px;background:white;cursor:pointer}
.panel{background:white;border:1px solid var(--line);border-radius:10px;padding:16px}
label{font-size:12px}label.model{display:block;padding:4px 0;overflow-wrap:anywhere}
canvas{width:100%;height:min(70vh,720px);min-height:420px;display:block}
.legend{display:flex;gap:15px;flex-wrap:wrap;font-size:12px;color:var(--muted);margin:12px 0}
.swatch{display:inline-block;width:20px;height:3px;vertical-align:middle;margin-right:5px}
.swatch-box{display:inline-block;width:12px;height:12px;vertical-align:middle;margin-right:5px;border:1px solid #8aaed6}
#caption{background:#f5f8fc;padding:10px 12px;border-left:3px solid #b6cce6;font-size:12px;line-height:1.5}
@media(max-width:830px){.layout{grid-template-columns:1fr}aside{border-right:0;border-bottom:1px solid var(--line)}main{padding:12px}}
</style></head><body><header><h2>T3 observable correlation comparison</h2>
<div class="small">Stored T1/T2/T3 sensitivity trajectories. Presentation only — no fitting, matching or RGE recalculation.</div></header>
<div class="layout"><aside><label for="x">X observable</label><select id="x"></select>
<label for="y">Y observable</label><select id="y"></select><div><button id="all">Show all</button> <button id="none">Hide all</button></div>
<p class="small">Sensitivity directions</p><div id="dirs"></div><p class="small">Models</p><div id="models"></div></aside>
<main><div class="panel"><div><button id="reset">Reset axes</button> <label><input type="checkbox" id="regions" checked> Gaussian 1σ / 3σ regions</label></div>
<canvas id="plot" aria-label="Observable correlation plot"></canvas><div class="legend" id="legend"></div>
<div id="caption"></div><p class="small" id="status"></p></div></main></div>
<script>'use strict';
const DATA=__DATA__;const LABELS=__LABELS__;const OBS=__OBS__;
const COLORS={T1:'#2376c9',T2:'#ed713d',T3:'#269a57'};
const REGION_STYLES={three:{fill:'#d3eafa',label:'Gaussian 3σ contour'},one:{fill:'#9dd1f5',label:'Gaussian 1σ contour'}};
const $=id=>document.getElementById(id);const canvas=$('plot'),ctx=canvas.getContext('2d');
const models=[...new Set(DATA.points.map(p=>p.model))].sort();let enabled=Object.fromEntries(models.map(m=>[m,true]));
let dirs={T1:true,T2:true,T3:true},manual=null;let x=OBS[3],y=OBS[1];
for(let id of ['x','y']){let select=$(id);for(let k of OBS){let o=document.createElement('option');o.value=k;o.textContent=LABELS[k];select.append(o)}select.value=id==='x'?x:y;
select.onchange=()=>{x=$('x').value;y=$('y').value;manual=null;draw()}}
for(let d of ['T1','T2','T3']){let l=document.createElement('label');l.className='model';let i=document.createElement('input');i.type='checkbox';i.checked=true;i.onchange=()=>{dirs[d]=i.checked;draw()};l.append(i,document.createTextNode(' '+d));$('dirs').append(l)}
for(let m of models){let l=document.createElement('label');l.className='model';let i=document.createElement('input');i.type='checkbox';i.checked=true;i.onchange=()=>{enabled[m]=i.checked;draw()};l.append(i,document.createTextNode(' '+m));$('models').append(l)}
$('all').onclick=()=>{models.forEach(m=>enabled[m]=true);$('models').querySelectorAll('input').forEach(i=>i.checked=true);draw()};
$('none').onclick=()=>{models.forEach(m=>enabled[m]=false);$('models').querySelectorAll('input').forEach(i=>i.checked=false);draw()};
$('reset').onclick=()=>{manual=null;draw()};$('regions').onchange=draw;
function gaussianPair(xKey,yKey){
  const target=DATA.target||{};const names=target.observable_names||[];const means=target.central_values||[];const cov=target.covariance||[];
  const ix=names.indexOf(xKey),iy=names.indexOf(yKey);if(ix<0||iy<0)return null;
  const mx=Number(means[ix]),my=Number(means[iy]);
  const sxx=Number(cov[ix]?.[ix]), syy=Number(cov[iy]?.[iy]), sxy=Number(cov[ix]?.[iy]);
  if(!Number.isFinite(mx)||!Number.isFinite(my)||!Number.isFinite(sxx)||!Number.isFinite(syy)||!Number.isFinite(sxy))return null;
  if(sxx<=0||syy<=0)return null;
  const det=sxx*syy-sxy*sxy; if(!(det>0)) return null;
  return {mx,my,sxx,syy,sxy};
}
function ellipseGeometry(pair,delta){
  const a=pair.sxx,b=pair.sxy,d=pair.syy;
  const tr=a+d, disc=Math.sqrt(Math.max(0,(a-d)*(a-d)+4*b*b));
  const l1=(tr+disc)/2, l2=(tr-disc)/2; if(!(l1>0&&l2>0)) return null;
  const angle=0.5*Math.atan2(2*b,a-d);
  return {cx:pair.mx,cy:pair.my,rx:Math.sqrt(delta*l1),ry:Math.sqrt(delta*l2),angle};
}
function boundsFromEllipse(g){
  const c=Math.cos(g.angle), s=Math.sin(g.angle);
  const dx=Math.sqrt(g.rx*g.rx*c*c + g.ry*g.ry*s*s);
  const dy=Math.sqrt(g.rx*g.rx*s*s + g.ry*g.ry*c*c);
  return [g.cx-dx,g.cx+dx,g.cy-dy,g.cy+dy];
}
function drawEllipseRegion(g,fill,X,Y){
  ctx.save(); ctx.translate(X(g.cx),Y(g.cy)); ctx.rotate(-g.angle);
  ctx.beginPath(); ctx.ellipse(0,0,(X(g.cx+g.rx)-X(g.cx)),Math.abs(Y(g.cy+g.ry)-Y(g.cy)),0,0,2*Math.PI);
  ctx.fillStyle=fill; ctx.fill(); ctx.restore();
}
function autoBounds(){
  const points=DATA.points.filter(p=>enabled[p.model]&&dirs[p.direction]); const xs=points.map(p=>p[x]), ys=points.map(p=>p[y]);
  const pair=gaussianPair(x,y);
  if(pair && $('regions').checked){
    for(let key of ['1sigma','3sigma']){
      const g=ellipseGeometry(pair,Number(DATA.target.ellipse_delta_chi2?.[key]));
      if(g){const box=boundsFromEllipse(g); xs.push(box[0],box[1]); ys.push(box[2],box[3]);}
    }
  }
  if(!xs.length) return [0,1,0,1];
  let xmin=Math.min(...xs), xmax=Math.max(...xs), ymin=Math.min(...ys), ymax=Math.max(...ys);
  let dx=xmax-xmin||Math.abs(xmax)||1, dy=ymax-ymin||Math.abs(ymax)||1;
  return [xmin-.08*dx,xmax+.08*dx,ymin-.08*dy,ymax+.08*dy];
}
function draw(){
  const w=Math.max(420,canvas.clientWidth),h=Math.max(420,canvas.clientHeight),dpr=window.devicePixelRatio||1;
  canvas.width=Math.round(w*dpr);canvas.height=Math.round(h*dpr);ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,w,h);
  const groups=new Map();for(let p of DATA.points){if(!enabled[p.model]||!dirs[p.direction])continue;let key=p.model+'|'+p.direction;if(!groups.has(key))groups.set(key,[]);groups.get(key).push(p)}
  let [xmin,xmax,ymin,ymax]=manual||autoBounds();
  const L=92,R=25,T=30,B=75,pw=w-L-R,ph=h-T-B; const X=v=>L+(v-xmin)/(xmax-xmin)*pw, Y=v=>T+ph-(v-ymin)/(ymax-ymin)*ph;
  ctx.fillStyle='#fff';ctx.fillRect(L,T,pw,ph);
  const pair=gaussianPair(x,y);
  if(pair && $('regions').checked){
    const g3=ellipseGeometry(pair,Number(DATA.target.ellipse_delta_chi2?.['3sigma']));
    const g1=ellipseGeometry(pair,Number(DATA.target.ellipse_delta_chi2?.['1sigma']));
    if(g3) drawEllipseRegion(g3,REGION_STYLES.three.fill,X,Y);
    if(g1) drawEllipseRegion(g1,REGION_STYLES.one.fill,X,Y);
  }
  ctx.strokeStyle='#e2e8f0';ctx.fillStyle='#556375';ctx.font='11px system-ui';
  for(let i=0;i<=5;i++){let v=xmin+(xmax-xmin)*i/5,p=X(v);ctx.beginPath();ctx.moveTo(p,T);ctx.lineTo(p,T+ph);ctx.stroke();ctx.textAlign='center';ctx.fillText(v.toPrecision(4),p,T+ph+19)}
  for(let i=0;i<=5;i++){let v=ymin+(ymax-ymin)*i/5,p=Y(v);ctx.beginPath();ctx.moveTo(L,p);ctx.lineTo(L+pw,p);ctx.stroke();ctx.textAlign='right';ctx.fillText(v.toPrecision(4),L-9,p+4)}
  ctx.strokeStyle='#475569';ctx.strokeRect(L,T,pw,ph);
  ctx.save();ctx.beginPath();ctx.rect(L,T,pw,ph);ctx.clip();
  for(let [key,ps] of groups){ps.sort((a,b)=>a.gamma-b.gamma);let d=key.split('|').pop();ctx.strokeStyle=COLORS[d];ctx.lineWidth=1.75;ctx.beginPath();ps.forEach((p,i)=>i?ctx.lineTo(X(p[x]),Y(p[y])):ctx.moveTo(X(p[x]),Y(p[y])));ctx.stroke();
    for(let p of ps){let a=X(p[x]),bb=Y(p[y]);ctx.beginPath();ctx.arc(a,bb,Math.abs(p.gamma-1)<1e-8?3.5:1.8,0,2*Math.PI);ctx.fillStyle=COLORS[d];ctx.fill()}}
  ctx.restore();
  ctx.fillStyle='#17243b';ctx.font='13px system-ui';ctx.textAlign='center';ctx.fillText(LABELS[x],L+pw/2,h-20);
  ctx.save();ctx.translate(21,T+ph/2);ctx.rotate(-Math.PI/2);ctx.fillText(LABELS[y],0,0);ctx.restore();
  $('legend').innerHTML=[
    ...Object.entries(COLORS).map(([k,c])=>'<span><i class="swatch" style="background:'+c+'"></i>'+k+'</span>'),
    '<span><i class="swatch-box" style="background:'+REGION_STYLES.three.fill+'"></i>'+REGION_STYLES.three.label+'</span>',
    '<span><i class="swatch-box" style="background:'+REGION_STYLES.one.fill+'"></i>'+REGION_STYLES.one.label+'</span>'
  ].join('');
  const visible=[...groups.values()].flat();
  $('caption').textContent=groups.size+' model–direction trajectories; '+visible.length+' stored successful points. '+
    (pair ? 'Blue filled ellipses are Gaussian 2D confidence regions derived from the stored covariance submatrix of the selected observable pair. ' : 'No usable experimental covariance is stored for this observable pair, so no ellipse is drawn. ')+
    'Curves connect saved sensitivity points in ascending γ order. γ = 1 points are larger. No extrapolation or recomputation.';
  $('status').textContent='Saved source: '+DATA.source+' | Target source: '+(DATA.target?.source||'')+' | '+(DATA.target?.note||'Gaussian ellipses, not official joint-fit contours.');
}
canvas.addEventListener('wheel',e=>{e.preventDefault();const k=e.deltaY>0?1.2:0.8;let [a,b,c,d]=manual||autoBounds();let mx=(a+b)/2,my=(c+d)/2;manual=[mx+(a-mx)*k,mx+(b-mx)*k,my+(c-my)*k,my+(d-my)*k];draw()},{passive:false});
window.addEventListener('resize',draw);draw();</script></body></html>'''


def _fallback_target_from_config(config_path: Path):
    """Convert the project's diagonal Gaussian NuFIT approximation to covariance form.

    The source has NO cross-observable correlations; off-diagonals are zero by
    assumption, not measured as zero. The resulting 2D ellipses are axis-aligned.
    """
    doc = read_json(config_path)
    centers = doc.get('central_values', {})
    errors = doc.get('one_sigma_errors', {})
    names = [name for name in OBS if name in centers and name in errors]
    if not names:
        raise ValueError(f'No supported NuFIT parameters in {config_path}')
    means = [float(centers[name]) for name in names]
    sigmas = [float(errors[name]) for name in names]
    if any(not math.isfinite(s) or s <= 0 for s in sigmas):
        raise ValueError('NuFIT errors must be finite and positive')
    covariance = [[sigmas[i] ** 2 if i == j else 0.0
                   for j in range(len(names))] for i in range(len(names))]
    return dict(observable_names=names, central_values=means,
                covariance=covariance,
                source=doc.get('source', str(config_path)),
                note=('Diagonal Gaussian approximation inferred from the NuFIT '
                      'one_sigma_errors; no off-diagonal correlations supplied. '
                      'Not an official NuFIT joint confidence surface.'))


def build(source: Path, output: Path, target_config: Path | None = None):
    points, target = gather(source)
    if target is None and target_config is not None:
        target = _fallback_target_from_config(target_config)
    if not points:
        raise RuntimeError(
            'No successful saved sensitivity points under ' + str(source) +
            '. Run an existing scan separately or supply its saved output; full simulation is NOT triggered.'
        )
    payload = {
        'points': points,
        'target': target_payload(target),
        'source': str(source),
    }
    body = HTML.replace('__DATA__', json.dumps(payload, ensure_ascii=False, allow_nan=False).replace('</', '<\\/'))
    body = body.replace('__LABELS__', json.dumps(LABELS, ensure_ascii=False)).replace('__OBS__', json.dumps(OBS))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(body, encoding='utf-8')
    print(f'Wrote {output}: {len(points)} saved successful points, {len(set(p["model"] for p in points))} models')
    if not payload['target']:
        print('Note: no usable target covariance metadata; Gaussian ellipse overlay omitted.')
    else:
        print('Using Gaussian covariance ellipses from saved target metadata.')


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, default=PROJECT_ROOT / 'Reports/output/full/fitted_sensitivity')
    p.add_argument('--output', type=Path, default=PROJECT_ROOT / 'Reports/output/display/interactive_correlations.html')
    p.add_argument('--target-config', type=Path, default=PROJECT_ROOT / 'configs/nufit_6_1_ic24_no.json', help='Gaussian fit-target configuration for saved CSV scans')
    args = p.parse_args(argv)
    build(args.source, args.output, args.target_config if args.target_config.is_file() else None)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
