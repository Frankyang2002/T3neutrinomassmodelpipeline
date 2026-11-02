"""Real UV/EFT1 RGE trajectories and offline HTML for shared-A(-2) benchmark.

No standalone RGE equations or interpolated physics are implemented here.
For the 14 nonmixed T3 models, call the SAME canonical UVRunner and
IntermediateScalarRunner as the general-C5 fitter, with explicit solver save
scales at the shared-Yukawa gamma=1 point. The mixed B_-1/C_-1 models use a
different conditional integrator; omit their dense curves rather than invent
continuations or silently switch to an incompatible RG system.

The produced interactive_RGE.html shows genuinely integrated renormalisable
couplings from 1e7 -> 1e5 GeV (UV) and 1e5 -> 1e3 GeV (EFT1).
No SMEFT Wilson-coefficient running below 1e3 GeV is claimed or shown.
"""
from __future__ import annotations

import hashlib
import html
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from scripts import RunT3GeneralAll16FromReference as campaign

FORMAT = 't3-shared-a2-general-uv-eft1-rge-trajectories-v1'
SERIES = {
    'gY': ('gY', ''), 'g2': ('g2', ''), 'g3': ('g3', ''),
    'lambdaH': ('Higgs quartic λH', ''),
    'lambdaT3_abs': ('|λT3|', ''),
    'lambdaT3_real': ('Re(λT3)', ''),
    'mS1_gev': ('√mS1² [GeV]', 'GeV'),
    'mS2_gev': ('√mS2² [GeV]', 'GeV'),
    'mS1_sq_gev2': ('mS1² [GeV²]', 'GeV²'),
    'mS2_sq_gev2': ('mS2² [GeV²]', 'GeV²'),
    'y1_norm': ('||Y1||F', ''),
    'y2_norm': ('||Y2||F', ''),
    'y1Conj_norm': ('||Y1conj||F', ''),
    'y2Conj_norm': ('||Y2conj||F', ''),
    'muS1_real_gev2': ('Re(μS1²) [GeV²]', 'GeV²'),
    'muS1_imag_gev2': ('Im(μS1²) [GeV²]', 'GeV²'),
    'muS2_real_gev2': ('Re(μS2²) [GeV²]', 'GeV²'),
    'muS2_imag_gev2': ('Im(μS2²) [GeV²]', 'GeV²'),
}


def _finite(value: Any) -> float:
    n = float(value)
    if not math.isfinite(n):
        raise ValueError('Nonfinite sampled RGE value')
    return n


def _norm(value: Any) -> float:
    arr = np.asarray(value, dtype=complex)
    if not (np.isfinite(arr.real).all() and np.isfinite(arr.imag).all()):
        raise ValueError('Nonfinite running Yukawa')
    return _finite(np.linalg.norm(arr))


def state_values(state: Any, *, uv: bool) -> dict[str, float]:
    """Pick only couplings stored in the actual integrated general RG state."""
    result = {name: _finite(getattr(state.sm, name))
              for name in ('gY', 'g2', 'g3', 'lambdaH')}
    lam = complex(state.lambdaT3)
    result.update(lambdaT3_abs=_finite(abs(lam)),
                  lambdaT3_real=_finite(lam.real))
    for label in ('mS1', 'mS2'):
        sq = _finite(getattr(state, label + 'Sq'))
        result[label + '_sq_gev2'] = sq
        if sq >= 0:
            result[label + '_gev'] = _finite(math.sqrt(sq))
    extras = getattr(state, 'general_couplings', {})
    for species, name in [('S1', 'muT3GeneralS1Sq'),
                          ('S2', 'muT3GeneralS2Sq')]:
        if name in extras:
            mu = complex(extras[name])
            result['mu' + species + '_real_gev2'] = _finite(mu.real)
            result['mu' + species + '_imag_gev2'] = _finite(mu.imag)
    if uv:
        for name in ('y1', 'y2'):
            result[name + '_norm'] = _norm(getattr(state, name))
        for name in ('y1Conj', 'y2Conj'):
            if name in extras:
                result[name + '_norm'] = _norm(extras[name])
    return result


def sampled_stage(result: Any, *, uv: bool) -> dict:
    """Use the solver's actual state_at_index; no curve interpolation."""
    scales = [_finite(x) for x in result.mu_gev]
    if len(scales) < 3 or any(scales[i] <= scales[i+1]
                            for i in range(len(scales)-1)):
        raise ValueError('Expected descending, distinct saved solver scales')
    samples = [state_values(result.state_at_index(i), uv=uv)
               for i in range(len(scales))]
    keys = set(samples[0])
    if any(set(row) != keys for row in samples):
        raise ValueError('Inconsistent coupling inventory along RGE trajectory')
    return {'mu_gev': scales,
            'series': {name: [row[name] for row in samples]
                       for name in sorted(keys)},
            'solver_success': bool(result.solver_success),
            'nfev': int(result.nfev)}


def _decode_real(raw: Any) -> float:
    if isinstance(raw, dict):
        if abs(float(raw.get('imag', 0))) > 1e-9 * max(1., abs(float(raw['real']))):
            raise ValueError('Unexpected complex endpoint for real comparison')
        raw = raw['real']
    return _finite(raw)


def validate_endpoints(trajectory: dict, centre: dict) -> None:
    """Check independently re-integrated gamma=1 trajectory reproduces stored endpoints."""
    detail = centre.get('numerical_details') or {}
    for stage, saved in [('uv', 'uv_at_fermion_threshold'),
                         ('eft1', 'eft1_at_scalar_threshold')]:
        expected = detail.get(saved)
        if not isinstance(expected, dict):
            raise ValueError('Saved centre lacks ' + saved)
        series = trajectory[stage]['series']
        for key in ('gY', 'g2', 'g3', 'lambdaH', 'mS1Sq', 'mS2Sq'):
            track = {'mS1Sq': 'mS1_sq_gev2', 'mS2Sq': 'mS2_sq_gev2'}.get(key, key)
            if key not in expected or track not in series:
                raise ValueError('Missing centre endpoint for ' + stage + '/' + key)
            actual = series[track][-1]
            reference = _decode_real(expected[key])
            if not math.isclose(actual, reference, rel_tol=2e-6,
                                abs_tol=2e-9 * max(1., abs(reference))):
                raise ValueError(f'{stage}/{key}: RGE rerun endpoint {actual} differs '
                                 f'from saved gamma=1 endpoint {reference}')


def compute_one(model: str, params: dict, settings: dict,
                shared_reference: dict, conj_source: str, centre: dict,
                n_points: int) -> dict:
    if model in campaign.MIXED:
        raise ValueError('B/C conditional RGE solver has no exported dense trajectory; '
                         'not substituting a different RGE basis')
    if n_points < 3 or n_points > 501:
        raise ValueError('--rge-points must be between 3 and 501')

    # Local import avoids a cycle with the shared-A2 orchestration module.
    from scripts.RunT3SharedA2General16Reports import fixed_for_model
    from Numerical.orchestration.GeneralRunningBenchmark import initial_state
    from Numerical.running.UVRunner import run_uv_segment
    from Numerical.running.FermionThresholdBoundary import build_intermediate_scalar_boundary
    from Numerical.running.IntermediateScalarRunner import run_intermediate_scalar_segment

    reference = campaign.read(campaign.ref_path(model))
    uv_betas, eft_betas, provenance, _kernel, matching_source, _scheme = (
        campaign.nonmixed_inputs(model, str(settings['canonical_run_id']),
                                 generate_general_c5=False))
    fixed = fixed_for_model(model, reference, uv_betas, shared_reference,
                            settings, lambda_scale=1., conj_source=conj_source)
    _, start = initial_state(
        uv_betas, mu_uv=1e7, mu_f=1e5, mu_s=1e3,
        overrides=campaign.override_yukawas(fixed, params))
    uv_run = run_uv_segment(
        start, uv_betas, 1e5,
        save_scales_gev=np.geomspace(1e7, 1e5, n_points),
        rtol=1e-7, atol=1e-10)
    boundary = build_intermediate_scalar_boundary(uv_run.final_state)
    eft_run = run_intermediate_scalar_segment(
        boundary, eft_betas, 1e3,
        save_scales_gev=np.geomspace(1e5, 1e3, n_points),
        rtol=1e-7, atol=1e-10)
    track = {'format': FORMAT, 'model': model, 'status': 'Success',
             'gamma': 1., 'n_save_points_per_stage': n_points,
             'scheme': 'canonical full-general renormalisable UV + EFT1 RGEs',
             'uv': sampled_stage(uv_run, uv=True),
             'eft1': sampled_stage(eft_run, uv=False),
             'matching_source_sha256': hashlib.sha256(matching_source.read_bytes()).hexdigest(),
             'notes': ['Actual saved ODE solution samples; no endpoint interpolation.',
                       'UV: 1e7 to 1e5 GeV; EFT1: 1e5 to 1e3 GeV.',
                       'Yukawas y1,y2 and conjugate Yukawas exist only in UV segment.',
                       'No SMEFT C5 running or matching trajectory is plotted.',
                       'Canonical general matching basis is not physically certified.']}
    validate_endpoints(track, centre)
    return track


# Self-contained interactive RGE dashboard, distinct from the older
# InteractiveModelComparison renderer (which requires a fully saved SMEFT C5
# trajectory and cannot correctly consume current general 16 output).
HTML = r'''<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Shared A(-2) UV and EFT1 interactive RGEs</title>
<style>:root{font-family:system-ui,-apple-system,Segoe UI,Arial,sans-serif;color:#18283b;background:#f5f7fb}
*{box-sizing:border-box}body{margin:0}header{background:white;border-bottom:1px solid #dbe3ee;padding:17px 25px}
h1{font-size:21px;margin:0 0 7px}.note{font-size:12px;color:#536273;line-height:1.55}
main{display:grid;grid-template-columns:250px minmax(0,1fr);gap:16px;padding:16px}
aside,.panel{background:white;border:1px solid #dbe3ee;border-radius:9px;padding:16px}
label{display:block;font-size:13px;margin:9px 0}select{width:100%;padding:9px;border:1px solid #c7d5e4;border-radius:6px}
.model-list{max-height:50vh;overflow:auto}.model{display:block;margin:7px 0;font-size:12px}
#plot{width:100%;height:clamp(380px,67vh,720px);display:block;touch-action:pan-y}
button{padding:8px 12px;background:#fff;border:1px solid #cbd7e5;border-radius:6px;cursor:pointer;margin-right:5px}
.toolbar{display:flex;gap:9px;flex-wrap:wrap;align-items:center}
#tooltip{height:22px;font:12px ui-monospace,monospace;color:#345}
@media(max-width:850px){main{grid-template-columns:1fr}.model-list{max-height:24vh}}
</style></head><body><header><h1>Shared A(-2) Yukawas — interactive UV/EFT1 RGEs</h1>
<div class="note">True solver samples at γ=1. Full-general renormalisable RGBeta running: UV 10⁷→10⁵ GeV; EFT1 10⁵→10³ GeV. No reconstructed C5 evolution. Mixed B(-1)/C(-1) are omitted from continuous RGE curves because their conditional solver does not export comparable dense samples. They remain in the sensitivity oval report.</div></header>
<main><aside><label>Running quantity<select id="quantity"></select></label>
<div class="toolbar"><button id="showAll">All</button><button id="hideAll">None</button></div>
<div id="modelList" class="model-list"></div><p id="coverage" class="note"></p></aside>
<section class="panel"><div class="toolbar">
<label><input id="showUV" type="checkbox" checked> UV</label>
<label><input id="showEFT" type="checkbox" checked> EFT1</label>
<label><input id="logY" type="checkbox"> Log₁₀ |value|</label>
<button id="reset">Reset view</button>
</div><canvas id="plot"></canvas><p id="tooltip">Select a quantity and model.</p>
<p class="note">Dashed vertical line: fermion threshold 10⁵ GeV; the scalar threshold is 10³ GeV. Lines connect actual saved integration samples only (not independent matching calculations). Values are plotted versus log₁₀(μ/GeV). Solid: UV, dashed: EFT1. Gaps mean that a coupling is absent in that EFT. Log-|value| hides zeros; its sign information is lost.</p>
<p class="note"><a href="index.html">Campaign overview</a> · <a href="sensitivity_oval.html">Sensitivity ovals</a> · <a href="rge_status.json">RGE coverage and limitations</a></p></section></main>
<script>'use strict';const D=__DATA__, $=id=>document.getElementById(id), cvs=$('plot'), ctx=cvs.getContext('2d');
const names=Object.keys(D.models), color=['#235aa7','#d76a38','#26936e','#7e5ab6','#bc7b21','#278f9b','#ac496b','#5160a3','#7a943b','#c38c28','#47748a','#975d96','#866d51','#6a82a1','#a94c36','#31917c'];
let enabled=Object.fromEntries(names.map(n=>[n,true])), manual=null, last=[], lastBounds=null;
const qlist=Object.entries(D.series);for(const [k,[lab]] of qlist){if(!names.some(n=>['uv','eft1'].some(st=>D.models[n][st].series[k])))continue;const opt=document.createElement('option');opt.value=k;opt.textContent=lab;$('quantity').append(opt)}
$('quantity').value='lambdaT3_abs';if(!$('quantity').value)$('quantity').selectedIndex=0;
for(const [i,n] of names.entries()){const lab=document.createElement('label');lab.className='model';const c=document.createElement('input');c.type='checkbox';c.checked=true;c.onchange=()=>{enabled[n]=c.checked;manual=null;draw()};lab.append(c);const span=document.createElement('span');span.textContent=' '+n;span.style.color=color[i%color.length];lab.append(span);$('modelList').append(lab)}
$('coverage').textContent=names.length+' genuine sampled model curves; '+D.omitted.length+' excluded. '+D.omitted.map(x=>x.model).join(', ');
for(const k of ['quantity','showUV','showEFT','logY'])$(k).onchange=()=>{manual=null;draw()};
$('showAll').onclick=()=>{for(const n of names)enabled[n]=true;sync();draw()};$('hideAll').onclick=()=>{for(const n of names)enabled[n]=false;sync();draw()};
function sync(){[...$('modelList').querySelectorAll('input')].forEach((c,i)=>c.checked=enabled[names[i]])}
$('reset').onclick=()=>{manual=null;draw()};
function curves(){const q=$('quantity').value, log=$('logY').checked,out=[];for(const [i,n] of names.entries()){if(!enabled[n])continue;for(const st of ['uv','eft1']){if(st==='uv'&&!$('showUV').checked||st==='eft1'&&!$('showEFT').checked)continue;const source=D.models[n][st], a=source.series[q];if(!a)continue;const pts=[];for(let j=0;j<a.length;j++){let y=Number(a[j]);if(log){if(y===0)continue;y=Math.log10(Math.abs(y))}if(Number.isFinite(y))pts.push([Math.log10(source.mu_gev[j]),y,source.mu_gev[j],a[j]])}if(pts.length)out.push({n,st,i,pts})}}return out}
function draw(){const device=window.devicePixelRatio||1,w=cvs.clientWidth,h=cvs.clientHeight;cvs.width=Math.round(w*device);cvs.height=Math.round(h*device);ctx.setTransform(device,0,0,device,0,0);ctx.clearRect(0,0,w,h);
const curves0=curves(), left=85,right=w-22,top=26,bottom=h-57;last=[];
if(!curves0.length){ctx.fillStyle='#667';ctx.fillText('No data for selected quantity/models/stages',left,60);return}
let ymin=Infinity,ymax=-Infinity;for(const c of curves0)for(const p of c.pts){ymin=Math.min(ymin,p[1]);ymax=Math.max(ymax,p[1])}let pad=(ymax-ymin)*.09||Math.abs(ymin)*.12||.15;ymin-=pad;ymax+=pad;
let bounds=manual||[3,7,ymin,ymax];lastBounds=bounds;const [xmin,xmax,y0,y1]=bounds;const x=v=>left+(v-xmin)/(xmax-xmin)*(right-left),y=v=>bottom-(v-y0)/(y1-y0)*(bottom-top);
ctx.font='12px system-ui';ctx.strokeStyle='#e0e6ef';ctx.lineWidth=1;for(let i=0;i<=5;i++){const v=y0+(y1-y0)*i/5;ctx.beginPath();ctx.moveTo(left,y(v));ctx.lineTo(right,y(v));ctx.stroke();ctx.fillStyle='#556';ctx.textAlign='right';ctx.fillText(v.toExponential(2),left-8,y(v)+4)}for(let i=3;i<=7;i++){ctx.beginPath();ctx.moveTo(x(i),top);ctx.lineTo(x(i),bottom);ctx.stroke();ctx.textAlign='center';ctx.fillStyle='#556';ctx.fillText('10^'+i,x(i),bottom+21)}ctx.fillText('μ [GeV] (log scale)',(left+right)/2,h-12);
ctx.strokeStyle='#3c485a';ctx.lineWidth=1;ctx.beginPath();ctx.moveTo(left,top);ctx.lineTo(left,bottom);ctx.lineTo(right,bottom);ctx.stroke();
ctx.setLineDash([3,5]);ctx.strokeStyle='#5e6d7e';ctx.beginPath();ctx.moveTo(x(5),top);ctx.lineTo(x(5),bottom);ctx.stroke();ctx.setLineDash([]);ctx.textAlign='left';ctx.fillStyle='#536273';ctx.fillText('F threshold',Math.min(x(5)+5,right-88),top+12);
for(const c of curves0){ctx.strokeStyle=color[c.i%color.length];ctx.lineWidth=2;ctx.setLineDash(c.st==='uv'?[]:[7,4]);ctx.beginPath();let begun=false;for(const p of c.pts){const px=x(p[0]),py=y(p[1]);if(!begun){ctx.moveTo(px,py);begun=true}else ctx.lineTo(px,py);last.push({px,py,model:c.n,stage:c.st,mu:p[2],value:p[3]})}ctx.stroke()}ctx.setLineDash([]);
}
function nearest(mx,my){let best=null,dist=130;for(const a of last){const d=(mx-a.px)**2+(my-a.py)**2;if(d<dist){best=a;dist=d}}return best}
cvs.onmousemove=e=>{const rect=cvs.getBoundingClientRect(),p=nearest(e.clientX-rect.left,e.clientY-rect.top);$('tooltip').textContent=p?`${p.model} / ${p.stage} / μ=${p.mu.toPrecision(5)} GeV / value=${Number(p.value).toPrecision(8)}`:'Move over a curve to inspect a saved solver point.'};
cvs.addEventListener('wheel',e=>{e.preventDefault();if(!lastBounds)return;const [a,b,c,d]=lastBounds,zoom=e.deltaY>0?1.2:0.8,mx=(a+b)/2,my=(c+d)/2;manual=[Math.max(2.8,mx+(a-mx)*zoom),Math.min(7.2,mx+(b-mx)*zoom),my+(c-my)*zoom,my+(d-my)*zoom];if(manual[1]-manual[0]<.03){manual=[a,b,c,d]}draw()},{passive:false});
cvs.ondblclick=()=>{manual=null;draw()};window.addEventListener('resize',draw);draw();</script></body></html>'''


def write_html(out: Path, statuses: list[dict]) -> bool:
    paths = []
    omitted = []
    for row in statuses:
        model = row['model']
        p = out/campaign.tag(model)/'rge_trajectories.json'
        if p.is_file():
            data = campaign.read(p)
            if data.get('format') == FORMAT and data.get('status') == 'Success' and data.get('model') == model:
                paths.append((model, data))
                continue
        missing = out/campaign.tag(model)/'rge_error.json'
        error = campaign.read(missing).get('error') if missing.is_file() else 'Not generated'
        omitted.append({'model': model, 'reason': error})
    payload = {'format': FORMAT, 'models': {n: d for n,d in paths},
               'series': SERIES, 'omitted': omitted}
    campaign.save(out/'rge_status.json', {
        'format': FORMAT, 'rge_curves_available': len(paths),
        'models_available': [n for n,_ in paths],
        'unavailable': omitted,
        'scope': 'renormalisable UV and EFT1 solver samples only',
        'no_weinberg_trajectory': True})
    if not paths:
        return False
    blob = json.dumps(payload, ensure_ascii=False, allow_nan=False).replace('</','<\\/')
    (out/'interactive_RGE.html').write_text(HTML.replace('__DATA__',blob),encoding='utf-8')
    return True


def generate(out: Path, statuses: list[dict], params: dict, settings: dict,
             reference: dict, conj_source: str, *, n_points: int=41,
             reuse: bool=True) -> bool:
    if not 3 <= n_points <= 501:
        raise ValueError('--rge-points must be in [3,501]')
    for idx,row in enumerate(statuses,1):
        model = row['model']
        model_dir = out/campaign.tag(model)
        model_dir.mkdir(parents=True, exist_ok=True)
        data_path = model_dir/'rge_trajectories.json'
        if model in campaign.MIXED:
            campaign.save(model_dir/'rge_error.json',{'model':model,
                'status':'UnsupportedDenseTrajectory',
                'error':'Conditional B/C mixed-scalar integrator does not expose dense saved RGE states; omitted, not interpolated.'})
            print(f'  RGE [{idx}/{len(statuses)}] {model}: conditional B/C curves unavailable',flush=True)
            continue
        if reuse and data_path.is_file():
            saved = campaign.read(data_path)
            matching_source = campaign.general_coefficient_path(model)
            same_source = (matching_source.is_file() and
                           saved.get('matching_source_sha256') ==
                           hashlib.sha256(matching_source.read_bytes()).hexdigest())
            if (saved.get('format')==FORMAT and saved.get('status')=='Success'
                and saved.get('model')==model and
                saved.get('n_save_points_per_stage')==n_points and same_source):
                print(f'  RGE [{idx}/{len(statuses)}] {model}: reuse dense curve',flush=True)
                continue
        centre_path=model_dir/'centre.json'
        if not centre_path.is_file():
            campaign.save(model_dir/'rge_error.json',{'model':model,
                'error':'No saved gamma=1 centre; cannot cross-check RGE endpoints'})
            continue
        centre=campaign.read(centre_path)
        if centre.get('numerical_details') is None:
            campaign.save(model_dir/'rge_error.json',{'model':model,
                'error':'Gamma=1 centre has no numerical details; cannot cross-check RGE endpoints'})
            continue
        try:
            curve=compute_one(model, params, settings, reference, conj_source,
                              centre, n_points)
            campaign.save(data_path,curve)
            (model_dir/'rge_error.json').unlink(missing_ok=True)
            print(f'  RGE [{idx}/{len(statuses)}] {model}: {n_points} real points per stage',flush=True)
        except Exception as exc:
            campaign.save(model_dir/'rge_error.json',{'model':model,
                'error':f'{type(exc).__name__}: {exc}'})
            print(f'  RGE [{idx}/{len(statuses)}] {model}: unavailable: {exc}',flush=True)
    return write_html(out,statuses)
