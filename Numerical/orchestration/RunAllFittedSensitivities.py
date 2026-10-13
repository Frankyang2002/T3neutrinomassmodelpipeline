"""Resume-safe per-model fitted T3 sensitivity evaluations (no global refitting).

Dry run by default; --run evaluates missing model scans only. Existing scans
are reused only if the gamma grid and selected fitted config SHA256 match.
Older scans without this runner's fingerprint are left intact and reported
as UnverifiedCache; use --new-output-dir instead of overwriting them.

    python -m Numerical.orchestration.RunAllFittedSensitivities
    python -m Numerical.orchestration.RunAllFittedSensitivities --run
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from Numerical.orchestration.ResumePipelineNumericalResults import PROJECT_ROOT, OUTPUT_ROOT, REQUIRED, discover_cases
from Numerical.orchestration.RunFittedSensitivityComparison import best_candidate
from Numerical.fitting.BenchmarkSensitivity import run_t3_sensitivity_scan


def load(p):
    return json.loads(Path(p).read_text(encoding='utf-8-sig'))


def fingerprint(config, target, paths, gamma_min, gamma_max, points):
    digest=hashlib.sha256()
    for p in (config,target,*paths):
        digest.update(str(Path(p).resolve()).encode())
        with Path(p).open('rb') as f:
            for block in iter(lambda:f.read(1024*1024),b''):
                digest.update(block)
    digest.update(json.dumps([gamma_min,gamma_max,points]).encode())
    return digest.hexdigest()


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scenario',default='smallY_smallL')
    p.add_argument('--run',action='store_true')
    p.add_argument('--gamma-min',type=float,default=0.8)
    p.add_argument('--gamma-max',type=float,default=1.2)
    p.add_argument('--points',type=int,default=21)
    p.add_argument('--max-baseline-chi2',type=float,default=0.01)
    p.add_argument('--output-dir',type=Path,default=None)
    a=p.parse_args(argv)
    if not(0<a.gamma_min<1<a.gamma_max) or a.points<3 or a.points%2!=1:
        p.error('Require gamma-min < 1 < gamma-max and odd points >= 3')
    output=a.output_dir or PROJECT_ROOT/'Reports/output/full/fitted_sensitivity'/a.scenario
    target=PROJECT_ROOT/'configs/nufit_6_1_ic24_no.json'
    cases={}
    for scenario,model,folder in discover_cases(OUTPUT_ROOT):
        if scenario==a.scenario:
            if model in cases:raise RuntimeError(f'Duplicate model {model}')
            cases[model]=folder
    rows=[]
    for model,folder in sorted(cases.items()):
        chosen=best_candidate(a.scenario,model)
        row={'model':model,'representation':'','status':'','fit_stage':'','fit_chi2':'',
             'verified_chi2':'','max_abs_prediction_difference':'','scan_json':'','message':''}
        if chosen is None:
            row.update(status='MissingFit',message='No saved successful fit')
        else:
            chi2,stage,config,summary=chosen
            row.update(fit_stage=stage,fit_chi2=chi2)
            paths=[folder/'data'/name for name in REQUIRED]
            missing=[str(x) for x in (config,summary,target,*paths) if not x.is_file()]
            if missing:
                row.update(status='MissingArtifacts',message='; '.join(missing))
            elif chi2>a.max_baseline_chi2:
                row.update(status='FitWarning',message='Saved fit exceeds chi2 threshold')
            else:
                dest=output/model/'sensitivity_diagnostics.json'
                row['scan_json']=str(dest)
                fp=fingerprint(config,target,paths,a.gamma_min,a.gamma_max,a.points)
                scan=None
                if dest.exists():
                    try:
                        scan=load(dest)
                        same=(scan.get('input_fingerprint')==fp)
                        settings=(scan.get('gamma_min')==a.gamma_min and scan.get('gamma_max')==a.gamma_max and scan.get('n_points')==a.points)
                        if not same or not settings or scan.get('status')!='Success':
                            row.update(status='UnverifiedCache',message='Saved scan differs or predates fingerprints; preserved. Select a new --output-dir for a fresh scan.')
                            scan=None
                    except Exception as exc:
                        row.update(status='UnverifiedCache',message=f'Unreadable old scan: {exc}')
                if scan is None and not row['status']:
                    if not a.run:
                        row.update(status='Ready')
                    else:
                        try:
                            payload=load(config)
                            payload['sensitivity']={'enabled':True,'gamma_min':a.gamma_min,'gamma_max':a.gamma_max,'n_points':a.points}
                            scan=run_t3_sensitivity_scan(base_payload=payload,uv_rgbeta_path=paths[0],
                                intermediate_rgbeta_path=paths[1],final_weinberg_path=paths[2],oscillation_target_path=target)
                            scan['input_fingerprint']=fp
                            scan['selection']={'source_config':str(config),'fit_stage':stage,'reported_fit_chi2':chi2}
                            dest.parent.mkdir(parents=True,exist_ok=True)
                            if dest.exists():raise FileExistsError('Refusing to overwrite '+str(dest))
                            with dest.open('x',encoding='utf-8') as f:json.dump(scan,f,indent=2)
                        except Exception as exc:
                            row.update(status='Failed',message=f'{type(exc).__name__}: {exc}')
                if scan is not None and not row['status']:
                    try:
                        centres=[next(q for q in scan['directions'][direction]['points'] if abs(float(q['gamma'])-1)<1e-9) for direction in ('T1','T2','T3')]
                        verified=float(centres[0]['chi2'])
                        reference=load(summary)['best']['prediction']
                        keys=('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23')
                        deviation=max(abs(float(centres[0]['prediction'][k])-float(reference[k])) for k in keys)
                        row.update(verified_chi2=verified,max_abs_prediction_difference=deviation)
                        if not math.isfinite(verified):raise ValueError('Nonfinite baseline chi2')
                        if any(abs(float(c['chi2'])-verified)>max(1e-8,abs(verified)*1e-6) for c in centres):
                            raise ValueError('Inconsistent gamma=1 centres')
                        row['status']='Success' if verified<=a.max_baseline_chi2 else 'VerificationWarning'
                    except Exception as exc:
                        row.update(status='Failed',message=f'Validation: {exc}')
        rows.append(row)
        print(f"{row['status']:19s} {model}",flush=True)
    print('Summary:',{st:sum(r['status']==st for r in rows) for st in sorted(set(r['status'] for r in rows))})
    if a.run:
        output.mkdir(parents=True,exist_ok=True)
        (output/'comparison_summary.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')
        with (output/'comparison_summary.csv').open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    return 1 if any(r['status']=='Failed' for r in rows) else 0

if __name__=='__main__':raise SystemExit(main())
