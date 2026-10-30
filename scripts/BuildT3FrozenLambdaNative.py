"""Generate frozen reports using the existing Lambda HTML canvas templates."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import csv
from Numerical.plotting import BuildGeneralScalarLambdaReports as original
from Numerical.plotting import BuildFrozenGeneralScalarLambdaReports as frozen


def render(run_ids, trajectory_dir, sensitivity_csv, output_dir, *, reference=None):
    rows, hashes = frozen.summary_rows(run_ids)
    traj, status = frozen.load_trajectories(trajectory_dir, hashes)
    sens, sstatus = frozen.load_sensitivity(sensitivity_csv, hashes)
    run_models={}
    for entry in traj:
        name = entry['model']+' ['+('small' if entry['case'].startswith('small') else 'large')+']'
        run_models[name]={'stages':{'UV':entry['uv_trajectory'],'EFT1':entry['eft1_trajectory']}}
    if not run_models:
        raise ValueError('No computed frozen RGE trajectories found; generate them first.')
    settings={'uv_start_gev':1e7,'fermion_threshold_gev':1e5,'scalar_endpoint_gev':1e3}
    rge_template=original.extract_lambda_template(original.ROOT/'Numerical/plotting/InteractiveModelComparison.py','template','write_dashboard')
    sens_template=original.extract_lambda_template(original.ROOT/'Numerical/plotting/BuildFittedModelEnvelopes.py','HTML')
    old_models = original.MODELS
    try:
        original.MODELS = tuple(run_models)
        rge = original.build_rge(run_models, settings, rge_template)
    finally:
        original.MODELS = old_models
    # Enforce clearly identified model/case and full five-observable points.
    grouped={}
    for row in sens:
        if not all(row.get(k) is not None for k in original.OBS[3:]):
            continue
        key=row['model']+' ['+('small' if row['case'].startswith('small') else 'large')+']'
        grouped.setdefault(key,{}).setdefault(row['direction'],[]).append(row)
    scans=[]
    for name,dirs in grouped.items():
        if set(dirs)!={'T1','T2','T3'}:continue
        centres=[]
        for d in ('T1','T2','T3'):
            pts=sorted(dirs[d],key=lambda p:p['gamma'])
            c=[p for p in pts if abs(p['gamma']-1)<1e-10]
            if len(c)!=1 or len(pts)<3:break
            centres.append(c[0])
        if len(centres)!=3:continue
        scans.append({'model':name,'obs':set(original.OBS[3:]),'centre':centres[0],
                      'directions':dirs})
    # Reuse the original Lambda HTML renderer, without disabling NuFIT/JUNO controls.
    # Target and projection are taken from the user's original Lambda report.
    ref_path=Path(reference) if reference else frozen.ROOT/'docs/FROZEN_LAMBDA_TARGET_REFERENCE.json'
    if not ref_path.is_absolute():
        ref_path=frozen.ROOT/ref_path
    ref=json.loads(ref_path.read_text(encoding='utf-8'))
    # Avoid falsely overlaying a different experimental fit definition.
    target=ref['target']
    observable_names=target.get('observable_names',[])
    if not all(k in observable_names for k in original.OBS[3:]):
        raise ValueError('The reference NuFIT target lacks the five expected observables.')
    payload={'models':[], 'skipped':[], 'target':target,
             'labels':{k:original.LABELS[k] for k in original.OBS[3:]},
             'observables':list(original.OBS[3:]),
             'source':str(sensitivity_csv),
             'future_sigma':ref['future_sigma'],
             'future_source':ref.get('future_source','JUNO projection from the original Lambda report'),
             'envelope_definition':'2x empirical sample standard-deviation covariance of saved gamma points; descriptive only'}
    for scan in scans:
        row=next((r for r in rows if r['model'] in scan['model'] and
                  ('small' if r['case'].startswith('small') else 'large') in scan['model']),None)
        payload['models'].append({'model':scan['model'],'label':scan['model'],
            'chi2':row['chi2'] if row else None, 'fit_stage':'frozen_general_scalar',
            'status':'Success','centre':{k:scan['centre'][k] for k in original.OBS[3:]},
            'directions':{d:[{'gamma':pt['gamma'], **{k:pt[k] for k in original.OBS[3:]}}
                for pt in sorted(scan['directions'][d], key=lambda p:p['gamma'])]
                for d in ('T1','T2','T3')}})
    sensitive=sens_template.replace('__DATA__',original.json_script(payload))
    if sensitive==sens_template:
        raise ValueError('Original Lambda sensitivity template no longer contains __DATA__')
    sensitive=sensitive.replace('Individually fitted T3 model sensitivity',
                                'Frozen general-scalar sensitivity — original Lambda renderer')
    sensitive=sensitive.replace('Each model uses its own fitted λT3, Yukawa matrices and other fitted parameters. Data from saved scans only — no RGE recomputation.',
                                'Saved frozen general-scalar Yukawa scans. Original Lambda NuFIT/JUNO overlays restored; covariance ovals are descriptive.')
    sensitive=sensitive.replace('<header><h2>Frozen', '<header><a href="index.html">Report index</a><h2>Frozen',1)
    # Use the original UI but distinguish scalar-only RGE data correctly.
    rge=rge.replace('General-scalar RGE running — continuous_v1',
                    'Frozen general-scalar scalar RGEs — original Lambda renderer')
    rge=rge.replace('16 saved general-scalar trajectories:',
                    'Saved frozen general-scalar scalar trajectories:')
    out=frozen.ROOT/'Reports/output/display_frozen_general_scalar'/output_dir
    out.mkdir(parents=True,exist_ok=True)
    (out/'interactive_rge.html').write_text(rge,encoding='utf-8')
    (out/'sensitivity_report.html').write_text(sensitive,encoding='utf-8')
    frozen.build_index(out,rows,status,sstatus)
    frozen.build_table(out,rows)
    # Report coverage for all 16 model IDs in the two frozen benchmarks.
    covered={(r['model'],r['case']) for r in rows}
    rge_covered={(x['model'],x['case']) for x in traj}
    scan_covered={(x['model'],x['case']) for x in sens}
    audit=[]
    for case in frozen.CASES:
        for model in frozen.MODELS:
            audit.append({'model':model,'case':case,'accepted_fit':(model,case) in covered,
                          'rge_data':(model,case) in rge_covered,
                          'sensitivity_data':(model,case) in scan_covered})
    (out/'coverage.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    with (out/'coverage.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=['model','case','accepted_fit','rge_data','sensitivity_data']); w.writeheader(); w.writerows(audit)
    # Explicit legacy reference links: do not mislabel original models as frozen matches.
    historical=out/'historical_lambda_reference.html'
    historical.write_text('''<!doctype html><html><meta charset="utf-8"><title>Historical Lambda comparison</title><body style="font:15px system-ui;max-width:900px;margin:3rem auto"><h1>Historical third benchmark and shared A(-2)</h1><p>These are separate historical analyses. They are not new frozen-general-scalar evaluations.</p><p>The original historical plots are preserved in the existing Lambda study output. Supply its report folder to compare the original fitted benchmark and shared A(-2) results. The current frozen report does not infer or manufacture missing shared-A(-2) data.</p></body></html>''',encoding='utf-8')
    index=(out/'index.html').read_text(encoding='utf-8')
    index=index.replace('</ul></body>', '<li><a href="historical_lambda_reference.html">Historical fitted-λ benchmark and shared A(-2) status</a></li><li><a href="coverage.csv">Coverage audit (16 models × 2 frozen cases)</a></li></ul></body>')
    (out/'index.html').write_text(index,encoding='utf-8')
    (out/'provenance.json').write_text(json.dumps({'input_files':hashes,'generator':'existing Lambda HTML templates','reference':str(ref_path)},indent=2),encoding='utf-8')
    print(f'Wrote {out}; frozen models={len(run_models)}, sensitivity scans={len(scans)}')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-ids',default='frozen_massnorm_a0_v1,frozen_massnorm_remaining_v1')
    p.add_argument('--trajectory-dir',default='output/studies/general_scalar_extended/frozen_display/rge_v1/trajectories_v1')
    p.add_argument('--sensitivity-csv',default='output/studies/general_scalar_extended/frozen_display/sensitivity_v1/frozen_gamma_points.csv')
    p.add_argument('--output-dir',default='frozen_lambda_native_v2')
    p.add_argument('--reference',default='docs/FROZEN_LAMBDA_TARGET_REFERENCE.json')
    a=p.parse_args(argv)
    render([x.strip() for x in a.run_ids.split(',') if x.strip()],a.trajectory_dir,a.sensitivity_csv,a.output_dir,reference=a.reference)

if __name__=='__main__': main()
