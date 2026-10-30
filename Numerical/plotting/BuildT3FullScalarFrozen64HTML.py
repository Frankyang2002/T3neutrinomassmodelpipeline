"""Single HTML entry point for the frozen full-general-scalar 16x4 study.

Uses the original Lambda RGE and sensitivity HTML/canvas sources. Only saved
matched-scenario full-scalar results and complete gamma scans are displayed.
No matching, ODE, fit or gamma evaluations occur in this presentation module.
Historical reduced, bulk48 and six_benchmarks results are NEVER auto-loaded.
"""
from __future__ import annotations

import argparse
import csv
import html
import json
import math
from pathlib import Path

from scripts import T3Frozen64Sources as policy

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'output/studies/general_scalar_extended/full_scalar_frozen_64'
OUTPUT = ROOT / 'Reports/output/display_general_scalar_frozen_64'


def page(title, body):
    return f'''<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(title)}</title>
<style>body{{font:14px system-ui,Arial;max-width:1400px;margin:auto;padding:20px;color:#17243b;background:#f7f9fc}}
a{{color:#225ba8}}table{{border-collapse:collapse;width:100%;font-size:12px;background:white}}
th,td{{padding:7px 10px;border-bottom:1px solid #d7dee8;text-align:left}}
nav{{display:flex;gap:14px;flex-wrap:wrap;margin:15px 0}}
.panel{{background:white;border:1px solid #d7dee8;border-radius:10px;padding:14px;margin:18px 0}}
.small{{color:#556375;font-size:12px}}.scroll{{overflow:auto}}</style>
<body><h1>{html.escape(title)}</h1>{body}</body></html>'''


def esc(s):
    return html.escape(str(s))


def rel(root, path):
    return policy.relative(root, path)


def load_validated(cell, run_id):
    b,m = cell['benchmark'],cell['model']
    r = {}
    if cell['status'] in {'SCALAR_COMPLETE','NEUTRINO_COMPLETE'}:
        expected = (f'output/studies/general_scalar_extended/full_scalar_frozen_64/'
                    f'{run_id}/{b}/{policy.tag(m)}/continuous_trajectory.json')
        policy.require(cell.get('trajectory') == expected, f'{b}/{m}: incorrect source path')
        path=policy.validate_source(ROOT,cell['trajectory'],cell.get('trajectory_sha256'))
        obj=policy.read_json(path)
        policy.require(obj.get('model')==m and obj.get('boundary_mode')=='FROZEN_UV' and
                       obj.get('settings',{}).get('benchmark')==b,
                       f'{b}/{m}: not the requested frozen-general-scalar trajectory')
        policy.require(obj.get('new_extra_quartics_evolved') is not None and
                       obj.get('extra_quadratic_terms_evolved') is not None and
                       obj.get('status')=='CandidateContinuousUVFermionThresholdEFT1Integrated',
                       f'{b}/{m}: incomplete full-scalar candidate')
        policy.require(obj.get('uv_ode_evaluations')==0, f'{b}/{m}: UV is not frozen')
        sources = obj.get('input_sha256')
        policy.require(isinstance(sources,dict) and sources == cell.get('source_sha256'),
                       f'{b}/{m}: source provenance differs from campaign manifest')
        for source_name, source_hash in sources.items():
            policy.validate_source(ROOT, source_name, source_hash)
        for key in ('uv_trajectory','eft1_trajectory'):
            rows=obj.get(key)
            policy.require(isinstance(rows,list) and len(rows)>=2,f'{b}/{m}: incomplete {key}')
            mu=[float(x['mu_gev']) for x in rows]
            policy.require(all(math.isfinite(x) and x>0 for x in mu) and
                           all(a>b for a,b in zip(mu,mu[1:])),f'{b}/{m}: invalid scales')
        r['trajectory']=obj
    if cell['status']=='NEUTRINO_COMPLETE':
        expected=(f'output/studies/general_scalar_extended/full_scalar_frozen_64/'
                  f'{run_id}/{b}/neutrino_v1/{policy.tag(m)}_neutrino.json')
        policy.require(cell.get('neutrino')==expected, f'{b}/{m}: wrong neutrino source')
        path=policy.validate_source(ROOT,cell['neutrino'],cell.get('neutrino_sha256'))
        obj=policy.read_json(path)
        policy.require(obj.get('model')==m and
                       obj.get('source_continuous_sha256')==cell['trajectory_sha256'] and
                       obj.get('full_general_scalar_physical_c5_certified') is False,
                       f'{b}/{m}: neutrino provenance/certification inconsistent')
        r['neutrino']=obj
    return r


def rge_html(benchmark, rows):
    models={}
    for cell,files in rows:
        traj=files.get('trajectory')
        if not traj:
            continue
        stages={}
        for name,key in (('UV','uv_trajectory'),('EFT1','eft1_trajectory')):
            part=[]
            for r in traj[key]:
                p=dict(r)
                # The old canvas expects a flat quartic registry. Preserve a
                # prefix to distinguish original couplings from extra quartics.
                qs=dict(p['extra_quartics'])
                qs.update({'ordinary: '+k:(v if isinstance(v,dict) else
                        {'real':float(v),'imag':0.0})
                           for k,v in p.get('ordinary_quartics',{}).items()})
                p['extra_quartics']=qs
                part.append(p)
            stages[name]=part
        models[cell['model']]={'stages':stages}
    if not models:
        return page('Frozen full-scalar RGE — '+benchmark,
                    '<p class="panel">No verified saved scalar trajectories for this benchmark. '
                    'No RGE data have been manufactured.</p>')
    from Numerical.plotting import BuildGeneralScalarLambdaReports as upstream
    old_models=upstream.MODELS
    try:
        upstream.MODELS=tuple(m for m in policy.MODELS if m in models)
        template=upstream.extract_lambda_template(
            ROOT/'Numerical/plotting/InteractiveModelComparison.py','template','write_dashboard')
        settings={'uv_start_gev':1e7,'fermion_threshold_gev':1e5,'scalar_endpoint_gev':1e3}
        result=upstream.build_rge(models,settings,template)
    finally:
        upstream.MODELS=old_models
    result=result.replace('General-scalar RGE running — continuous_v1',
                          'Frozen general-scalar EFT1 RGE — '+benchmark)
    disclosure=('<div style="padding:12px;background:#fff7e6;border:1px solid #ddc999;margin:12px">'
                '<b>Frozen-UV boundary, not UV RG:</b> the displayed UV interval 10⁷→10⁵ GeV is a '
                'flat boundary guide. Only the EFT1 interval 10⁵→10³ GeV is numerically integrated. '
                'Ordinary and added quartics are distinguished by their selector prefixes. '
                'No time-dependent C5 or PMNS is inferred. All scalar physics is candidate-basis and '
                'full sequential matching is not certified.</div>')
    result=result.replace('<div class="layout">',disclosure+'<div class="layout">',1)
    return result


def scans_for(cell, files, root, scan_id):
    if 'neutrino' not in files:
        return None
    path=root/scan_id/cell['benchmark']/(policy.tag(cell['model'])+'.csv')
    meta=path.with_suffix('.json')
    if not path.is_file() and not meta.is_file():
        return None
    policy.require(path.is_file() and meta.is_file(),'partial gamma scan: '+str(path))
    doc=policy.read_json(meta)
    policy.require(doc['csv_sha256']==policy.digest(path),'gamma scan hash mismatch')
    source=doc.get('source',{})
    policy.require(source.get('model')==cell['model'] and
                   source.get('benchmark')==cell['benchmark'] and
                   source.get('trajectory_sha256')==cell['trajectory_sha256'] and
                   source.get('neutrino_sha256')==cell['neutrino_sha256'] and
                   source.get('generator_sha256') == policy.digest(
                       ROOT / 'scripts/GenerateT3FrozenFullScalarGamma.py'),
                   'gamma scan is from an incompatible numerical source')
    with path.open(newline='',encoding='utf-8-sig') as f:
        lines=list(csv.DictReader(f))
    required={(d,float(g)) for d in policy.DIRECTIONS for g in source['gamma']}
    policy.require(len(lines)==len(required),'incomplete gamma scan')
    indexed={}
    for r in lines:
        k=r['direction'],float(r['gamma'])
        policy.require(k in required and k not in indexed,'duplicate/unknown gamma point')
        policy.require(r['model']==cell['model'] and r['benchmark']==cell['benchmark'] and
                       r['status']=='NeutrinoEvaluated','invalid gamma result')
        p={'gamma':k[1]}
        for field in policy.FIELDS:
            v=float(r[field]);policy.require(math.isfinite(v),'nonfinite gamma response')
            p[field]=v
        indexed[k]=p
    policy.require(set(indexed)==required,'missing gamma point')
    bydir={d:sorted((indexed[d,g] for dd,g in required if dd==d),key=lambda x:x['gamma'])
           for d in policy.DIRECTIONS}
    centre=[next(p for p in bydir[d] if abs(p['gamma']-1)<1e-12) for d in policy.DIRECTIONS]
    for a in centre[1:]:
        for field in policy.FIELDS:
            policy.require(math.isclose(a[field],centre[0][field],rel_tol=1e-7,abs_tol=1e-10),
                           'gamma=1 centres disagree')
    saved=files['neutrino']['oscillation_observables_normal_ordering']
    for field in policy.FIELDS:
        policy.require(math.isclose(centre[0][field],float(saved[field]),rel_tol=1e-6,abs_tol=1e-10),
                       'gamma centre does not reproduce central neutrino')
    return {'model':cell['model'],'label':cell['model'],'status':files['neutrino']['status'],
            'fit_stage':'historical_seed_on_frozen_general_scalar',
            'centre':centre[0], 'directions':bydir}


def sensitivity_html(benchmark,rows,base,scan_id):
    plotted=[]
    missing=[]
    for cell,files in rows:
        item=scans_for(cell,files,base,scan_id)
        if item:
            plotted.append(item)
        else:
            missing.append({'model':cell['model'],'status':'No verified complete full-scalar gamma scan'})
    from Numerical.plotting import BuildFittedModelEnvelopes as original
    ref=policy.read_json(ROOT/'docs/FROZEN_LAMBDA_TARGET_REFERENCE.json')
    target=ref['target']
    for m in plotted:
        # Pure display chi-square: historical diagonal Gaussian NuFIT covariance.
        c=m['centre']; keys=target['observable_names']
        m['chi2']=sum((c[k]-target['central_values'][i])**2/target['covariance'][i][i]
                      for i,k in enumerate(keys))
    payload={'models':plotted,'skipped':missing,
             'target':target,'labels':{k:original.LABELS[k] for k in policy.FIELDS},
             'observables':list(policy.FIELDS),
             'source':'verified frozen full-scalar gamma CSVs: '+scan_id+'/'+benchmark,
             'future_sigma':ref['future_sigma'],'future_source':ref['future_source'],
             'envelope_definition':'2x empirical sample covariance of actual gamma responses; descriptive, not confidence regions'}
    text=original.HTML.replace('__DATA__',json.dumps(payload,allow_nan=False).replace('<','\\u003c'))
    policy.require(text!=original.HTML,'missing sensitivity renderer placeholder')
    text=text.replace('Individually fitted T3 model sensitivity','Frozen general-scalar sensitivity — '+benchmark)
    text=text.replace('Each model uses its own fitted λT3, Yukawa matrices and other fitted parameters. Data from saved scans only — no RGE recomputation.',
                'Genuine T1/T2/T3 perturbations at 100 TeV. Original historical Yukawa seeds are NOT reoptimised for the general scalar potential.')
    text=text.replace('fitted models plotted','frozen-general-scalar candidate models plotted')
    text=text.replace('about each fitted centre','about each saved gamma=1 centre')
    explanation=('<div style="margin:12px;padding:12px;background:#fff7e6;border:1px solid #ddc999">'
                 '<b>Physics status:</b> General-scalar candidate RGE/matching only. '
                 'Six A/D/E self-bilinear C5 calculations are restricted unmixed controls; '
                 'the B/C one-scale mixing result is conditional. No complete sequential EFT matching. '
                 'The blue NuFIT and orange JUNO ovals are historical Gaussian approximations, '
                 'not official likelihood contours. Missing gamma data are not replaced with old scans. '
                 f'Full-gamma models: {len(plotted)}/16.</div>')
    return text.replace('<div class="layout">',explanation+'<div class="layout">',1),len(plotted)


def write_report(run_id: str, scan_id: str, destination: Path):
    base=DATA/run_id
    doc=policy.validate_manifest(policy.read_json(base/'manifest.json'))
    cells={}
    for c in doc['cells']:
        cells[(c['benchmark'],c['model'])]=load_validated(c,run_id)
    destination.mkdir(parents=True,exist_ok=True)
    result=[]
    for b in policy.BENCHMARKS:
        rows=[(c,cells[b,c['model']]) for c in doc['cells'] if c['benchmark']==b]
        folder=destination/b;folder.mkdir(parents=True,exist_ok=True)
        (folder/'interactive_rge.html').write_text(rge_html(b,rows),encoding='utf-8')
        sensitive,covered=sensitivity_html(b,rows,base,scan_id)
        (folder/'sensitivity_ovals.html').write_text(sensitive,encoding='utf-8')
        for c,f in rows:
            result.append({'benchmark':b,'model':c['model'],'status':c['status'],
                           'scalar_rge':bool(f.get('trajectory')),
                           'candidate_neutrino':bool(f.get('neutrino')),
                           'matching_status':f.get('neutrino',{}).get('status') or '',
                           'complete_gamma_scan':bool(scans_for(c,f,base,scan_id))})
        print(f'{b}: scalar trajectories {sum(bool(f.get("trajectory")) for _,f in rows)}/16; '
              f'genuine gamma scans {covered}/16')
    with (destination/'coverage.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(result[0]));writer.writeheader();writer.writerows(result)
    source_manifest=base/'manifest.json'
    (destination/'provenance.json').write_text(json.dumps({
        'run_id':run_id,'scan_id':scan_id,'source_manifest':policy.relative(ROOT,source_manifest),
        'manifest_sha256':policy.digest(source_manifest),
        'physics_recomputed_by_html':False,
        'full_scalar_sequential_matching_certified':False,
        'coverage':{'slots':64,'scalar':sum(x['scalar_rge'] for x in result),
                    'neutrino':sum(x['candidate_neutrino'] for x in result),
                    'sensitivity':sum(x['complete_gamma_scan'] for x in result)}},
        indent=2),encoding='utf-8')
    nav='<nav>' + ''.join(f'<a href="{b}/interactive_rge.html">{esc(b)} RGE</a> <a href="{b}/sensitivity_ovals.html">Ovals</a>'
                          for b in policy.BENCHMARKS)+'</nav>'
    head='<tr>'+''.join('<th>'+esc(k)+'</th>' for k in result[0])+'</tr>'
    body=''.join('<tr>'+''.join('<td>'+esc(row[k])+'</td>' for k in result[0])+'</tr>' for row in result)
    content=('<div class="panel">64 planned cells, all 16 models × 4 frozen-UV cases. '
             'These are candidate general-scalar calculations: historical Yukawa seeds, '
             'all symmetry-allowed additional quartics and mass terms where present, '
             'and EFT1 running. No full sequential C5 matching or new full-scalar Yukawa optimisation '
             'is claimed. Empty slots are deliberately visible.</div>'+nav+
             '<p><a href="coverage.csv">Coverage CSV</a> · '
             '<a href="provenance.json">Provenance hashes</a></p>'+
             '<div class="panel scroll"><table>'+head+body+'</table></div>')
    (destination/'index.html').write_text(page('T3 frozen full-general-scalar — 64 report slots',content),encoding='utf-8')
    return result


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id',default='frozen_general64_v1')
    p.add_argument('--scan-id',default='gamma_scans_v1')
    p.add_argument('--output',type=Path)
    args=p.parse_args(argv)
    policy.require(bool(policy.RUN_ID.fullmatch(args.run_id)), 'unsafe run id')
    policy.require(bool(policy.RUN_ID.fullmatch(args.scan_id)), 'unsafe scan id')
    out=args.output or OUTPUT/args.run_id
    write_report(args.run_id,args.scan_id,out)
    print('HTML report:',out/'index.html')


if __name__=='__main__':
    main()
