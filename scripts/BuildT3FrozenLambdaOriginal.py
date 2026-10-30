"""Build the original Lambda sensitivity dashboard from genuine frozen scans.

Original plotting generators are imported read-only; this does not create scalar
RGE charts or re-label historical trajectories as frozen RGEs.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from Numerical.plotting import BuildFittedModelEnvelopes as lambda_envelopes
from Numerical.plotting import BuildFrozenGeneralScalarLambdaReports as frozen

OBS = tuple(lambda_envelopes.OBS)
DIRECTIONS = ('T1', 'T2', 'T3')


def collect_scans(rows, saved_fits):
    lookup = {(r['model'], r['case']): r for r in saved_fits}
    groups = {}
    for row in rows:
        key = (row['model'], row['case'])
        if key not in lookup or row['direction'] not in DIRECTIONS:
            continue
        if not all(row.get(k) is not None for k in OBS):
            continue
        groups.setdefault(key, {}).setdefault(row['direction'], []).append(row)
    models, skipped = [], []
    for (model, case), directions in sorted(groups.items()):
        if any(len(directions.get(d, [])) < 3 for d in DIRECTIONS):
            skipped.append({'model': model + ' / ' + case, 'status': 'IncompleteGammaScan'})
            continue
        centre = []
        for d in DIRECTIONS:
            points = sorted(directions[d], key=lambda p:p['gamma'])
            points_at_one = [p for p in points if abs(p['gamma']-1.0)<1e-10]
            if len(points_at_one) != 1:
                break
            centre.append(points_at_one[0])
        if len(centre) != 3 or any(any(abs(p[k]-centre[0][k]) > 1e-6*max(abs(p[k]),abs(centre[0][k]),1e-12) for k in OBS) for p in centre[1:]):
            skipped.append({'model':model+' / '+case,'status':'GammaCentreMismatch'})
            continue
        tag = 'small' if case.startswith('small') else 'large'
        label = f'{model} [{tag}]'
        models.append({'model':label,'label':label,'status':'Success',
                       'chi2':lookup[(model,case)]['chi2'],'fit_stage':'frozen_general_scalar',
                       'centre':{k:centre[0][k] for k in OBS},
                       'directions':{d:[{'gamma':p['gamma'], **{k:p[k] for k in OBS}}
                                        for p in sorted(directions[d],key=lambda p:p['gamma'])] for d in DIRECTIONS}})
    return models, skipped


def build(output_dir, run_ids, sensitivity_csv, reference, historical_envelopes=None, historical_rge=None):
    root = frozen.ROOT
    rows, sources = frozen.summary_rows(run_ids)
    scan_rows, scan_state = frozen.load_sensitivity(sensitivity_csv, sources)
    ref_path = root / reference
    ref = json.loads(ref_path.read_text(encoding='utf-8'))
    models, skipped = collect_scans(scan_rows,rows)
    known = {(r['model'],r['case']) for r in rows}
    scanned = {(x['model'].split(' [')[0], 'small_uv_frozen_5e-9' if '[small]' in x['model'] else 'large_uv_frozen_0p5') for x in models}
    coverage=[]
    for case in frozen.CASES:
        for model in frozen.MODELS:
            if (model,case) in scanned:
                status='PlottedFrozenScan'
            elif (model,case) in known:
                status='AcceptedFitNoCompleteScan'
            elif model in {'B_-1','C_-1'}:
                status='SpecialMixedScalarNotScanned'
            elif model in {'A_-4','B_-3','C_-3','E_-4'}:
                status='LegacyC5ParityBlocked'
            else:
                status='NoDiscoveredAcceptedFrozenFit'
            coverage.append({'model':model,'case':case,'status':status})
    payload={'models':models, 'skipped':skipped,
             'target':ref['target'], 'labels':lambda_envelopes.LABELS,
             'observables':list(OBS), 'source':str(sensitivity_csv),
             'future_sigma':ref['future_sigma'],
             'future_source':ref.get('future_source','Historical Lambda reference'),
             'envelope_definition':'2x empirical sample standard-deviation covariance ellipse of genuine gamma points; descriptive only'}
    page=lambda_envelopes.HTML.replace('__DATA__',json.dumps(payload,ensure_ascii=False,allow_nan=False).replace('</','<\\/'))
    if '__DATA__' in page:
        raise RuntimeError('Lambda HTML data marker replacement failed')
    # Only title / provenance text changes; the original JS renderer remains untouched.
    page=page.replace('Individually fitted T3 model sensitivity','Frozen-UV T3 sensitivity — original Lambda renderer')
    page=page.replace('Each model uses its own fitted λT3, Yukawa matrices and other fitted parameters. Data from saved scans only — no RGE recomputation.',
                      'Frozen-UV fitted Yukawas; genuine saved perturbation scans. Experimental contours are the original Lambda-study reference.')
    out=root/'Reports/output/display_frozen_general_scalar'/output_dir
    out.mkdir(parents=True,exist_ok=True)
    (out/'sensitivity_report.html').write_text(page,encoding='utf-8')
    with (out/'coverage.csv').open('w',newline='',encoding='utf-8') as stream:
        w=csv.DictWriter(stream,fieldnames=['model','case','status']);w.writeheader();w.writerows(coverage)
    with (out/'all_results.csv').open('w',newline='',encoding='utf-8') as stream:
        columns=['model','case','chi2','status','run_id','source_path',*frozen.OBS]
        w=csv.DictWriter(stream,fieldnames=columns);w.writeheader();w.writerows(rows)
    links=[('sensitivity_report.html','Frozen sensitivity ovals — original Lambda renderer'),('all_results.csv','Saved fit results'),('coverage.csv','Model/benchmark coverage and exclusions')]
    for src,label,dest in [(historical_envelopes,'Original historical Lambda sensitivity','historical_envelopes.html'),(historical_rge,'Original historical Lambda neutrino RGE','historical_interactive_rge.html')]:
        if src:
            sp=root/src if not Path(src).is_absolute() else Path(src)
            if sp.is_file():
                (out/dest).write_bytes(sp.read_bytes());links.append((dest,label))
    nav=''.join(f'<li><a href="{href}">{label}</a></li>' for href,label in links)
    explanation=(f'Frozen sensitivity series: {len(models)}. Accepted frozen configurations discovered: {len(rows)}. '
                 'Unplotted cases are listed in coverage.csv. The original neutrino RGE requires true running/threshold histories; '
                 'this adapter does not replace them with scalar masses or interpolated endpoint values.')
    (out/'index.html').write_text(f'<!doctype html><html><meta charset="utf-8"><title>T3 Lambda reports</title><body style="font:15px system-ui;max-width:940px;margin:40px auto"><h1>T3 neutrino reports</h1><p>{explanation}</p><ul>{nav}</ul></body></html>',encoding='utf-8')
    (out/'provenance.json').write_text(json.dumps({'source_files':sources,'run_ids':run_ids,'scan_state':scan_state,'plotted':len(models),'accepted_fits':len(rows),'notes':['Original Lambda plot code reused unchanged','No scalar-mass RGE substitution','Missing modes explicitly listed in coverage.csv']},indent=2),encoding='utf-8')
    print(f'Wrote {out}; accepted fits={len(rows)}; plotted scans={len(models)}; missing coverage={sum(x["status"]!="PlottedFrozenScan" for x in coverage)}')
    return out


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-ids',default='frozen_massnorm_a0_v1,frozen_massnorm_remaining_v1')
    p.add_argument('--sensitivity-csv',default='output/studies/general_scalar_extended/frozen_display/sensitivity_v1/frozen_gamma_points.csv')
    p.add_argument('--output-dir',default='frozen_original_lambda_v1')
    p.add_argument('--reference',default='docs/FROZEN_LAMBDA_TARGET_REFERENCE.json')
    p.add_argument('--historical-envelopes',default='')
    p.add_argument('--historical-rge',default='')
    a=p.parse_args(argv)
    build(a.output_dir,[x.strip() for x in a.run_ids.split(',') if x.strip()],a.sensitivity_csv,a.reference,a.historical_envelopes,a.historical_rge)

if __name__=='__main__': main()
