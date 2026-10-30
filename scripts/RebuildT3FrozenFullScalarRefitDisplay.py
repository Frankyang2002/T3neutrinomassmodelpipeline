"""Re-evaluate frozen-full-scalar low fitted Yukawas and draw genuine sensitivity scans.

This is a FIXED-SCALAR-BACKGROUND conditional re-evaluation; the original
frozen_general64_v1 pages are never modified. The candidate C5 matching is
not a certified sequential EFT result. Reuses the actual refit Evaluator.
"""
from __future__ import annotations
import argparse
import csv
import html
import json
import math
from pathlib import Path
import hashlib

import numpy as np

from scripts import RefitT3FrozenFullScalarLow as fit
from scripts import GenerateT3FrozenSensitivityScans as directions
from scripts import RunT3FullScalarNeutrino as neutrino
from scripts import T3Frozen64Sources as policy
from Numerical.orchestration import RunSmallLambdaRefit as historical
from Numerical.fitting.OscillationFit import OscillationFitTarget

ROOT = fit.ROOT
BASE = fit.OUTPUT
DEST = ROOT/'Reports/output/display_full_scalar_refitted'
KEYS = fit.NAMES
GAMMA_DEFAULT = '0.95,0.975,1,1.025,1.05'


def parse_grid(text):
    gamma = directions.parse_gamma_grid(text)
    if len(gamma) < 3:
        raise ValueError('need at least three gamma values')
    return gamma


def loaded(refit_dir, source_root, model, target_hash):
    path=refit_dir/(fit.full.model_tag(model)+'_fit.json')
    if not path.is_file():
        raise FileNotFoundError('missing refit: '+str(path))
    data=json.loads(path.read_text(encoding='utf-8-sig'))
    if data.get('format')!=fit.FORMAT or data.get('model')!=model:
        raise ValueError('incompatible refit format/model: '+str(path))
    if data.get('status')!='ConditionalFitAccepted':
        raise ValueError('refit not accepted: '+model+' / '+str(data.get('status')))
    if data.get('target_sha256')!=target_hash:
        raise ValueError('target changed since refit: '+model)
    params=data.get('parameters',{})
    if set(params)!=set(fit.PARAMS):
        raise ValueError('fitted Yukawa inventory mismatch: '+model)
    uv,eft,folder,hashes=fit.load_case(source_root,model)
    if data.get('source_sha256')!=hashes:
        raise ValueError('source EFT1 trajectory/sidecar changed: '+model)
    vector=np.asarray([float(params[k]) for k in fit.PARAMS],dtype=float)
    if not np.isfinite(vector).all():
        raise ValueError('nonfinite Yukawa vector: '+model)
    return data,vector,(uv,eft,folder),path


def assess(model, x, sources, target):
    uv,eft,folder=sources
    ev=fit.Evaluator(model,uv,eft,folder,target)
    residual=ev.evaluate(x)
    if ev.best is None or ev.failures:
        raise RuntimeError('refit evaluator failed: '+model+' '+ev.last_error)
    return float(np.dot(residual,residual)),ev.best['prediction']


def display_plots(rows, model, folder, reference):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    data=[r for r in rows if r['model']==model]
    if not data: return
    pdf=folder/(fit.full.model_tag(model)+'_sensitivity.pdf')
    xref=reference
    with PdfPages(pdf) as pages:
        for kind, ylabel in [('chi2',r'$\chi^2$'),('delta_m3l_sq_ev2',r'$\Delta m_{3\ell}^2\ (\mathrm{eV}^2)$'),('sin2_theta12',r'$\sin^2\theta_{12}$'),('sin2_theta13',r'$\sin^2\theta_{13}$'),('sin2_theta23',r'$\sin^2\theta_{23}$')]:
            fig,ax=plt.subplots(figsize=(8,5))
            for d in ('T1','T2','T3'):
                v=sorted([r for r in data if r['direction']==d], key=lambda r:r['gamma'])
                ax.plot([r['gamma'] for r in v],[r[kind] for r in v],marker='o',label=d)
            if kind in xref:
                ax.axhline(float(xref[kind]),color='0.4',linestyle='--',linewidth=1,label='fit centre')
            ax.set_xlabel(r'Yukawa row scale $\gamma$');ax.set_ylabel(ylabel)
            ax.set_title(model+' | Refitted conditional frozen scalar background')
            ax.grid(alpha=.2);ax.legend();fig.tight_layout();pages.savefig(fig);plt.close(fig)
        fig, ax=plt.subplots(figsize=(7,6))
        for d in ('T1','T2','T3'):
            v=sorted([r for r in data if r['direction']==d],key=lambda r:r['gamma'])
            ax.plot([r['sin2_theta12'] for r in v],[r['sin2_theta23'] for r in v],marker='o',label=d)
        ax.scatter([xref['sin2_theta12']],[xref['sin2_theta23']],marker='x',s=75,color='black',label='fit centre')
        ax.set_xlabel(r'$\sin^2\theta_{12}$');ax.set_ylabel(r'$\sin^2\theta_{23}$')
        ax.set_title(model+' | PMNS sensitivity trajectories (not a covariance ellipse)')
        ax.legend();ax.grid(alpha=.2);fig.tight_layout();pages.savefig(fig);plt.close(fig)
    return pdf


def write_index(out, entries, run_id, gamma):
    cards=[]
    for e in entries:
        if e.get('pdf'):
            link=html.escape(e['pdf'])
            cards.append(f'<li><a href="{link}">{html.escape(e["model"])}</a> &mdash; fitted χ² = {e["chi2"]:.5g}, centred χ² = {e["recomputed_chi2"]:.5g}</li>')
        else:
            cards.append(f'<li>{html.escape(e["model"])} &mdash; {html.escape(e["error"])}</li>')
    markup='''<!doctype html><html lang="en"><meta charset="utf-8"><title>T3 refitted frozen scalar sensitivity</title>
<style>body{font-family:system-ui,sans-serif;max-width:960px;margin:35px auto;padding:0 18px;line-height:1.5}li{margin:8px 0}small{color:#555}</style>
<h1>T3 full-scalar refitted low | conditional sensitivity</h1>
<p>Run: <code>'''+html.escape(run_id)+'''</code>. Gamma: <code>'''+html.escape(','.join(map(str,gamma)))+'''</code>.</p>
<p><strong>Scope:</strong> 18-real-Yukawa refits on the saved frozen scalar EFT1 background. This dashboard recomputes the existing conditional C5 matching and SM+Weinberg evolution for each gamma. It does not rerun Yukawa-dependent scalar RGEs, certify sequential matching, or supply a C12 trajectory. The apparent ovals are sensitivity trajectories, not statistical confidence regions.</p>
<p><a href="sensitivity_points.csv">Raw scan CSV</a> · <a href="summary.json">Source metadata and consistency checks</a></p>
<h2>Models</h2><ul>'''+''.join(cards)+'''</ul></html>'''
    (out/'index.html').write_text(markup,encoding='utf-8')


def run(args):
    models=policy.MODELS if args.models=='all' else tuple(x.strip() for x in args.models.split(','))
    if not models or len(set(models))!=len(models) or any(x not in policy.MODELS for x in models):
        raise ValueError('invalid models')
    gamma=parse_grid(args.gamma_grid)
    refit_dir=BASE/args.refit_run_id
    source_root=(fit.INPUT/'optimised_low').resolve()
    out=DEST/args.display_id
    if out.exists() and not args.plan: raise FileExistsError('output exists; choose new --display-id: '+str(out))
    target_path=historical.TARGET
    target=OscillationFitTarget.from_json(target_path)
    thash=fit.hashfile(target_path)
    inputs=[];errors=[]
    for m in models:
        try:
            data,x,sources,path=loaded(refit_dir,source_root,m,thash)
            inputs.append((m,data,x,sources,path))
            print(m,'READY','fitted chi2=',data['chi2'],flush=True)
        except Exception as exc:
            errors.append({'model':m,'error':type(exc).__name__+': '+str(exc)})
            print(m,'BLOCKED',str(exc),flush=True)
    if args.plan:
        return 0 if not errors else 2
    if errors and not args.allow_incomplete:
        raise RuntimeError('preflight failed for '+str(len(errors))+' models; use --plan to inspect, or --allow-incomplete')
    out.mkdir(parents=True)
    rows=[];summary=[]
    columns=['model','direction','gamma','chi2',*KEYS]
    for m,data,x,sources,path in inputs:
        try:
            chi,pred=assess(m,x,sources,target)
            if not math.isclose(chi,float(data['chi2']),rel_tol=1e-3,abs_tol=1e-7):
                raise ValueError(f'central chi2 disagreement {chi} vs {data["chi2"]}')
            local=[]
            for d in ('T1','T2','T3'):
                for g in gamma:
                    xx=directions.apply_direction(x,d,g)
                    c,v=assess(m,xx,sources,target)
                    local.append({'model':m,'direction':d,'gamma':g,'chi2':c,**{k:float(v[k]) for k in KEYS}})
            rows.extend(local)
            pdf=display_plots(local,m,out,{'chi2':chi,**pred})
            summary.append({'model':m,'chi2':float(data['chi2']),'recomputed_chi2':chi,
                'pdf':pdf.name,'source_fit':str(path.relative_to(ROOT)),'source_fit_sha256':fit.hashfile(path),
                'matching_kind':data.get('matching_kind'), 'matching_certified':False})
            print(m,'DONE',len(local),'scan points',flush=True)
        except Exception as exc:
            summary.append({'model':m,'error':type(exc).__name__+': '+str(exc)})
            print(m,'BLOCKED',str(exc),flush=True)
    with (out/'sensitivity_points.csv').open('w',newline='',encoding='utf-8') as fh:
        writer=csv.DictWriter(fh,fieldnames=columns);writer.writeheader();writer.writerows(rows)
    (out/'summary.json').write_text(json.dumps({'format':'refitted-full-scalar-frozen-low-conditional-display-v1',
        'source_refit_run_id':args.refit_run_id,'gamma_grid':gamma,'models':summary,
        'scalar_background_recomputed_per_yukawa':False,'sequential_matching_certified':False},indent=2)+'\n',encoding='utf-8')
    write_index(out,summary,args.refit_run_id,gamma)
    print('Dashboard:',out/'index.html',flush=True)
    return 0 if all('pdf' in d for d in summary) and not errors else 2


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--refit-run-id',default='frozen_full_scalar_refit_low_pathfix_v2')
    p.add_argument('--display-id',default='frozen_low_refitted_sensitivity_v1')
    p.add_argument('--models',default='all')
    p.add_argument('--gamma-grid',default=GAMMA_DEFAULT)
    p.add_argument('--plan',action='store_true')
    p.add_argument('--allow-incomplete',action='store_true')
    return run(p.parse_args(argv))

if __name__=='__main__':
    raise SystemExit(main())
