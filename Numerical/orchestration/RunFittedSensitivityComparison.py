"""Select best feasibility candidates, independently verify gamma=1, scan T1/T2/T3.

Usage: python -m Numerical.orchestration.RunFittedSensitivityComparison [--run] [--limit N]
Never modifies original comparison configurations, fitted candidates or old diagnostics.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path

from Numerical.orchestration.ResumePipelineNumericalResults import PROJECT_ROOT, OUTPUT_ROOT, REQUIRED, discover_cases
from Numerical.fitting.BenchmarkSensitivity import run_t3_sensitivity_scan

SCENARIO = 'smallY_smallL'
STAGES = ('fit_feasibility', 'fit_feasibility_retry', 'fit_feasibility_focused')
OBS = ('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23')

def read_json(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))

def best_candidate(scenario, model):
    candidates=[]
    for stage in STAGES:
        folder=PROJECT_ROOT/'Reports'/'output'/'full'/stage/scenario/model
        summary=folder/'feasibility_summary.json'
        config=folder/'fitted_candidate_config.json'
        if not (summary.is_file() and config.is_file()):
            continue
        try:
            chi2=float(read_json(summary)['best']['chi2'])
            if math.isfinite(chi2):
                candidates.append((chi2,stage,config,summary))
        except (ValueError,KeyError,TypeError):
            continue
    return min(candidates,key=lambda x:x[0]) if candidates else None

def model_class(model):
    # Group by SU(2) representation dimensions; do not infer a multiplicity factor.
    bits=model.split('_')
    return f'({bits[2]},{bits[4]},{bits[6]})'

def plot_comparison(rows, output):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib unavailable; JSON and CSV outputs are still complete.')
        return
    from itertools import combinations
    valid=[r for r in rows if r['status']=='Success']
    if not valid: return
    plotdir=output/'plots'
    plotdir.mkdir(parents=True,exist_ok=True)
    labels={'delta_m21_sq_ev2':r'$\Delta m^2_{21}\ (\mathrm{eV}^2)$',
            'delta_m3l_sq_ev2':r'$\Delta m^2_{3\ell}\ (\mathrm{eV}^2)$',
            'sin2_theta12':r'$\sin^2\theta_{12}$',
            'sin2_theta13':r'$\sin^2\theta_{13}$',
            'sin2_theta23':r'$\sin^2\theta_{23}$'}
    palette=plt.get_cmap('tab20')
    for x,y in combinations(OBS,2):
        fig,ax=plt.subplots(figsize=(8,6))
        for i,row in enumerate(valid):
            scan=read_json(Path(row['scan_json']))
            for direction,style in [('T1','-'),('T2','--'),('T3',':')]:
                points=scan['directions'][direction]['points']
                ax.plot([p['prediction'][x] for p in points],
                        [p['prediction'][y] for p in points],
                        linestyle=style,color=palette(i%20),alpha=.75,lw=1.4,
                        label=row['model'] if direction=='T1' else None)
            p=scan['directions']['T1']['points'][len(scan['directions']['T1']['points'])//2]['prediction']
            ax.scatter([p[x]],[p[y]],color=palette(i%20),s=12)
        target=read_json(Path(valid[0]['scan_json']))['target']
        central=dict(zip(target['observable_names'],target['central_values']))
        if x in central and y in central:
            ax.scatter([central[x]],[central[y]],marker='x',s=70,color='black',zorder=10,label='NuFIT centre')
        ax.set(xlabel=labels[x],ylabel=labels[y],title='Fitted T3 sensitivities: '+SCENARIO)
        ax.grid(alpha=.2)
        ax.legend(fontsize=5,ncol=2,loc='best')
        fig.tight_layout()
        fig.savefig(plotdir/f'{x}_vs_{y}.png',dpi=160)
        plt.close(fig)
    print('Saved 10 observable-pair plots to',plotdir)

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scenario',default=SCENARIO)
    parser.add_argument('--run',action='store_true')
    parser.add_argument('--limit',type=int,default=0)
    parser.add_argument('--force',action='store_true')
    parser.add_argument('--gamma-min',type=float,default=.95)
    parser.add_argument('--gamma-max',type=float,default=1.05)
    parser.add_argument('--points',type=int,default=11)
    parser.add_argument('--max-baseline-chi2',type=float,default=.01,
                        help='Verification threshold, not a fit constraint')
    args=parser.parse_args(argv)
    if args.points<3 or args.points%2!=1 or not 0<args.gamma_min<1<args.gamma_max:
        parser.error('Require odd --points >=3 and gamma-min < 1 < gamma-max')
    if args.limit<0: parser.error('--limit must be >=0')
    cases={}
    for scenario,key,case in discover_cases(OUTPUT_ROOT):
        if scenario==args.scenario:
            if key in cases: parser.error('Multiple artifact directories for '+key)
            cases[key]=case
    output=PROJECT_ROOT/'Reports'/'output'/'full'/'fitted_sensitivity'/args.scenario
    target=PROJECT_ROOT/'configs'/'nufit_6_1_ic24_no.json'
    rows=[]
    for model in sorted(cases)[:args.limit or None]:
        case=cases[model]
        selected=best_candidate(args.scenario,model)
        row={'model':model,'representation':model_class(model),'status':'Ready',
             'fit_stage':'','fit_chi2':'','verified_chi2':'','max_abs_prediction_difference':'',
             'scan_json':'','message':''}
        if selected is None:
            row.update(status='MissingFit',message='No fitted candidate and summary')
        else:
            chi2,stage,config,summary=selected
            row.update(fit_stage=stage,fit_chi2=chi2)
            missing=[str(case/'data'/name) for name in REQUIRED if not (case/'data'/name).is_file()]
            if missing: row.update(status='MissingArtifacts',message='; '.join(missing))
            elif not target.is_file(): row.update(status='MissingTarget',message=str(target))
            else:
                destination=output/model/'sensitivity_diagnostics.json'
                row['scan_json']=str(destination)
                if args.run:
                    try:
                        if destination.is_file() and not args.force:
                            scan=read_json(destination)
                            if scan.get('selection',{}).get('source_config')!=str(config) or scan.get('n_points')!=args.points or scan.get('gamma_min')!=args.gamma_min or scan.get('gamma_max')!=args.gamma_max:
                                raise RuntimeError('Existing scan has different settings; use --force')
                        else:
                            payload=read_json(config)
                            payload['sensitivity']={'enabled':True,'gamma_min':args.gamma_min,
                                                    'gamma_max':args.gamma_max,'n_points':args.points}
                            paths=[case/'data'/name for name in REQUIRED]
                            scan=run_t3_sensitivity_scan(base_payload=payload,
                                uv_rgbeta_path=paths[0],intermediate_rgbeta_path=paths[1],
                                final_weinberg_path=paths[2],oscillation_target_path=target)
                            scan['selection']={'source_config':str(config),'source_summary':str(summary),
                                               'fit_stage':stage,'reported_fit_chi2':chi2}
                            destination.parent.mkdir(parents=True,exist_ok=True)
                            destination.write_text(json.dumps(scan,indent=2)+'\n',encoding='utf-8')
                        centers=[scan['directions'][direction]['points'][args.points//2]
                                 for direction in ('T1','T2','T3')]
                        if any(abs(p['gamma']-1)>1e-12 for p in centers):
                            raise RuntimeError('gamma=1 not present at scan centre')
                        if any(not math.isfinite(p['chi2']) for p in centers):
                            raise RuntimeError('Non-finite gamma=1 chi2')
                        reference=read_json(summary)['best']['prediction']
                        diff=max(abs(centers[0]['prediction'][k]-reference[k]) for k in OBS)
                        verified=centers[0]['chi2']
                        row.update(verified_chi2=verified,max_abs_prediction_difference=diff)
                        if any(abs(p['chi2']-verified)>max(1e-8,abs(verified)*1e-6) for p in centers):
                            raise RuntimeError('gamma=1 predictions differ between directions')
                        row['status']='Success' if verified<=args.max_baseline_chi2 else 'VerificationWarning'
                        if row['status']!='Success': row['message']='gamma=1 chi2 exceeds threshold'
                    except Exception as exc:
                        row.update(status='Failed',message=f'{type(exc).__name__}: {exc}')
        rows.append(row)
        print(f"{row['status']:20} {model} fit_chi2={row['fit_chi2']} verified={row['verified_chi2']}",flush=True)
    if not args.run:
        print('Dry run only; use --run to independently evaluate all gamma points.')
        return 0
    output.mkdir(parents=True,exist_ok=True)
    fields=list(rows[0]) if rows else []
    if fields:
        with (output/'comparison_summary.csv').open('w',newline='',encoding='utf-8') as f:
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
    (output/'comparison_summary.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')
    plot_comparison(rows,output)
    print('Summary:',output/'comparison_summary.csv')
    return 1 if any(r['status']!='Success' for r in rows) else 0

if __name__=='__main__':
    raise SystemExit(main())
