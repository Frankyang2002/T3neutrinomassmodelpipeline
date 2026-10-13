"""Representation-aware 16x16 common-input cross-evaluation using frozen fitted references.

Run: python -m Numerical.orchestration.RunCommonParameterComparison --run
Does not modify source benchmarks, fitting results, earlier scans, or v1 comparison.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path

from Numerical.orchestration.ResumePipelineNumericalResults import PROJECT_ROOT, OUTPUT_ROOT, REQUIRED, discover_cases
from Numerical.core.RGBetaEvaluator import load_rgbeta_payload
from Numerical.running.FinalC5TrajectoryAdapter import FinalC5TrajectoryEvaluator
from Numerical.fitting.OscillationFit import OscillationFitTarget
from Numerical.fitting.ParameterScan import make_t3_point_evaluator
from Numerical.fitting.ScanCLI import build_uv_state_from_config
from types import SimpleNamespace

SCENARIO = 'smallY_smallL'
REF = PROJECT_ROOT / 'Numerical' / 'benchmarks' / 'fitted_reference'
OUT = PROJECT_ROOT / 'Reports' / 'output' / 'full' / 'common_parameter_comparison_v2' / SCENARIO
OBS = ('delta_m21_sq_ev2', 'delta_m3l_sq_ev2', 'sin2_theta12', 'sin2_theta13', 'sin2_theta23')
SHARED = ('y1', 'y2', 'lambdaT3', 'MF', 'mS1Sq', 'mS2Sq')

def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))

def models_and_cases():
    found = {}
    for scenario, model, case in discover_cases(OUTPUT_ROOT):
        if scenario != SCENARIO: continue
        if model in found: raise RuntimeError('Duplicate artifact case: '+model)
        found[model] = case
    configs = {p.stem: p for p in REF.glob('T3_dS1_*.json')}
    if len(configs)!=16 or len(found)!=16 or set(configs)!=set(found):
        raise RuntimeError(f'Expected matching 16 models: configs={len(configs)}, cases={len(found)}, missing_artifacts={sorted(set(configs)-set(found))}, missing_config={sorted(set(found)-set(configs))}')
    return configs, found

def candidate_for_target(source, target):
    """Replace common physical input fields, retaining target representation and SM state."""
    payload = deepcopy(target)
    a = source['base_state']['ordinary']
    b = payload['base_state']['ordinary']
    for key in SHARED:
        if key not in a or key not in b:
            raise ValueError(f'Missing common parameter {key}')
        if key in ('y1', 'y2', 'MF'):
            if len(a[key]['real']) != len(b[key]['real']) or any(len(x) != len(y) for x,y in zip(a[key]['real'],b[key]['real'])):
                raise ValueError(f'Incompatible matrix shape for {key}')
        b[key] = deepcopy(a[key])
    if source['scales'] != target['scales'] or source['base_state']['mu_gev'] != target['base_state']['mu_gev']:
        raise ValueError('Scale mismatch between common candidates')
    if source['base_state']['sm'] != target['base_state']['sm']:
        raise ValueError('SM state mismatch between common candidates')
    payload['benchmark_search'] = {'enabled':False}
    payload['sensitivity'] = {'enabled':False}
    return payload

def evaluate(payload, case, target):
    paths=[case/'data'/name for name in REQUIRED]
    if any(not p.is_file() for p in paths):
        raise FileNotFoundError('Missing artifact(s): '+', '.join(str(p) for p in paths if not p.is_file()))
    uv=load_rgbeta_payload(paths[0]); eft=load_rgbeta_payload(paths[1]); c5=FinalC5TrajectoryEvaluator(paths[2])
    scales=payload['scales']; numerical=payload.get('numerical',{})
    view=SimpleNamespace(raw={**payload,'scan':{'bindings':{}}})
    evaluator=make_t3_point_evaluator(
        uv_state_builder=lambda _params: build_uv_state_from_config(view,{}),
        c5_builder=c5, uv_rgbeta_payload=uv, eft1_rgbeta_payload=eft,
        target=target, mu_fermion_threshold_gev=float(scales['mu_fermion_threshold_gev']),
        mu_scalar_threshold_gev=float(scales['mu_scalar_threshold_gev']),
        mu_low_gev=float(scales['mu_low_gev']), vev_gev=float(numerical.get('vev_gev',246.22)),
        rtol=float(numerical.get('rtol',1e-8)), atol=float(numerical.get('atol',1e-11)))
    fit=evaluator({})
    return float(fit.chi2), dict(zip(fit.observable_names, map(float,fit.prediction)))

def write_new(path, contents):
    with path.open('x',encoding='utf-8',newline='') as f: f.write(contents)

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',action='store_true',help='Evaluate all 256 candidate/model combinations')
    parser.add_argument('--limit-candidates',type=int,default=0,help='Pilot with first N source candidates')
    args=parser.parse_args(argv)
    configs,cases=models_and_cases()
    names=sorted(configs)
    if args.limit_candidates<0: parser.error('--limit-candidates must be nonnegative')
    names_source=names[:args.limit_candidates] if args.limit_candidates else names
    print('Reference files:',REF)
    print('Target models:',len(names),'candidate sources:',len(names_source),'evaluations:',len(names)*len(names_source))
    print('Output:',OUT)
    if not args.run:
        print('Dry run; no outputs written. Add --run to evaluate.')
        return 0
    if OUT.exists():
        raise FileExistsError('Refusing to overwrite an existing common comparison: '+str(OUT))
    target=OscillationFitTarget.from_json(PROJECT_ROOT/'configs'/'nufit_6_1_ic24_no.json')
    sources={name:load(configs[name]) for name in names_source}
    targets={name:load(configs[name]) for name in names}
    rows=[]; ranking=[]
    for source_name,source in sources.items():
        candidate_rows=[]
        for target_name in names:
            try:
                payload=candidate_for_target(source,targets[target_name])
                chi2,pred=evaluate(payload,cases[target_name],target)
                if not math.isfinite(chi2): raise ValueError('Non-finite chi2')
                row={'source_candidate':source_name,'target_model':target_name,'status':'Success','chi2':chi2,'error':''}
                row.update(pred)
            except Exception as exc:
                row={'source_candidate':source_name,'target_model':target_name,'status':'Failed','chi2':'','error':f'{type(exc).__name__}: {exc}'}
                row.update({k:'' for k in OBS})
            candidate_rows.append(row)
            print(f'{source_name} -> {target_name}: {row["status"]} chi2={row["chi2"]}',flush=True)
        rows.extend(candidate_rows)
        successes=[r for r in candidate_rows if r['status']=='Success']
        rank={'source_candidate':source_name,'successful_models':len(successes),'failed_models':16-len(successes),
              'average_chi2':sum(float(r['chi2']) for r in successes)/16 if len(successes)==16 else '',
              'worst_chi2':max(float(r['chi2']) for r in successes) if len(successes)==16 else ''}
        ranking.append(rank)
    ranking.sort(key=lambda r:(r['failed_models'],float(r['average_chi2']) if r['average_chi2']!='' else float('inf')))
    # A unique directory with exclusive file creation prevents accidental replacement of earlier results.
    OUT.mkdir(parents=True,exist_ok=False)
    def csv_output(name, items):
        import io
        s=io.StringIO(); w=csv.DictWriter(s,fieldnames=list(items[0]));w.writeheader();w.writerows(items)
        write_new(OUT/name,s.getvalue())
    csv_output('cross_evaluation.csv',rows)
    csv_output('candidate_ranking.csv',ranking)
    best=ranking[0] if ranking and ranking[0]['failed_models']==0 else None
    summary={'created_utc':datetime.now(timezone.utc).isoformat(),'scenario':SCENARIO,
             'evaluated_candidates':len(sources),'target_models':len(names),
             'objective':'mean chi2 over 16 models; candidate with any failure is ineligible',
             'best_common_candidate':best,'ranking':ranking,
             'note':'Transfers only y1, y2, lambdaT3 and masses. All representation-dependent scalar couplings retained from each target. Frozen fitted references and v1 outputs never modified.'}
    write_new(OUT/'summary.json',json.dumps(summary,indent=2)+'\n')
    if best:
        # A separate new file: do not touch the frozen source configuration.
        write_new(OUT/'best_common_input.json',json.dumps({'source_model':best['source_candidate'],
            'shared_ordinary':{k:sources[best['source_candidate']]['base_state']['ordinary'][k] for k in SHARED},
            'scales':sources[best['source_candidate']]['scales']},indent=2)+'\n')
    if best:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            import numpy as np
            chosen=[r for r in rows if r['source_candidate']==best['source_candidate']]
            center=dict(zip(target.observable_names,map(float,target.central_values)))
            sigma=dict(zip(target.observable_names,np.sqrt(np.diag(target.covariance))))
            for observable in OBS:
                fig,ax=plt.subplots(figsize=(11,5))
                xs=np.arange(len(chosen))
                ys=[float(r[observable]) for r in chosen]
                ax.axhspan(center[observable]-sigma[observable],center[observable]+sigma[observable],alpha=.18,label='NuFIT Gaussian 1 sigma')
                ax.axhline(center[observable],color='black',ls='--',lw=1,label='NuFIT central')
                ax.scatter(xs,ys,s=32,label='Common-input model prediction')
                ax.set_xticks(xs, [r['target_model'] for r in chosen],rotation=90,fontsize=7)
                ax.set_ylabel(observable)
                ax.set_title('Common input from '+best['source_candidate'])
                ax.grid(axis='y',alpha=.2);ax.legend(fontsize=8)
                fig.tight_layout()
                fig.savefig(OUT/(observable+'_common_predictions.png'),dpi=160)
                plt.close(fig)
        except ImportError:
            print('Matplotlib not installed; CSV and JSON outputs complete.')
    print('Best candidate:',best)
    print('Saved:',OUT)
    return 0 if best else 1

if __name__=='__main__':
    raise SystemExit(main())
