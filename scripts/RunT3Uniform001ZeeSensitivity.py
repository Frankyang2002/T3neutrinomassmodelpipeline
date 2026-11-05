"""Recalculate T1/T2/T3 compensated-coupling sensitivity from accepted uniform-0.01 references.

This is a NEW input campaign, but it reuses the existing general-T3 physics:
  - scripts.RunT3GeneralZeeSensitivity.deform/gamma_grid/prediction;
  - scripts.RunT3UniformScalarCampaign.validated_reference;
  - scripts.RunT3UniformScalarBenchmark.ModelFactory (cached beta/C5 kernel per model).

Unlike the old merged-source Zee scan, the only permissible baselines here are
16 accepted *refitted* uniform-0.01 UV configurations. All independent BSM
quartics remain +0.01 and conjugate Yukawas start at UV zero. The nominal
invariance T1/T2/T3 is checked on each deformed UV config. No Yukawa
reoptimisation is performed along gamma: these are actual sensitivity paths.

Output is the SAME points.csv/summary.json input contract consumed by the
UNCHANGED scripts.PlotT3GeneralZeeOvals.py figure/report generator.

Matching remains conditional one-scale C5, NOT certified fully sequential EFT.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REF_ROOT = ROOT / 'output/studies/thesis_qualitative/uniform_scalar_001/campaigns'
OUT_ROOT = ROOT / 'output/studies/thesis_qualitative/uniform_scalar_001/zee_sensitivity'
SCHEMA = 't3-general-zee-compensated-sensitivity-uniform001-v1'
COLUMNS = ('model','direction','gamma','status','chi2',
           'delta_m21_sq_ev2','delta_m3l_sq_ev2',
           'sin2_theta12','sin2_theta13','sin2_theta23','error')
DIRECTIONS = ('T1','T2','T3')
QUARTIC = 0.01


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding='utf-8-sig'))


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def atomic_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    tmp.replace(path)


def atomic_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp')
    with tmp.open('w',newline='',encoding='utf-8') as file:
        writer=csv.DictWriter(file,fieldnames=COLUMNS,extrasaction='ignore')
        writer.writeheader()
        writer.writerows([{k:row.get(k,'') for k in COLUMNS} for row in rows])
    tmp.replace(path)


def row_key(r: dict) -> tuple[str,str,float]:
    return str(r['model']),str(r['direction']),round(float(r['gamma']),12)


def selected_models(value: str, all_models: tuple[str,...]) -> list[str]:
    choice=list(all_models) if value=='all' else [x.strip() for x in value.split(',') if x.strip()]
    if not choice or len(set(choice))!=len(choice) or set(choice)-set(all_models):
        raise ValueError('Unknown/empty/duplicate model selection')
    return [m for m in all_models if m in choice]


def check_fitted_references(folder: Path, models: list[str], *,
                             cutoff: float, canonical_run: str):
    """Read-only acceptance and physical-boundary verification; no matching/RGE run."""
    from scripts import RunT3GeneralZeeSensitivity as zee
    from scripts import RunT3UniformScalarCampaign as camp
    from scripts import RunT3QualitativeScalarRefits as refit
    selected=zee.fit_selections(models,'merged')
    info={}
    for model in models:
        args=argparse.Namespace(fit_run_id='merged',canonical_run_id=canonical_run)
        original,fit,paths=zee.sources(model,args,selected[model])
        saved,cfg=camp.validated_reference(
            folder/'reference'/refit.tag(model),model=model,
            source_paths=paths,manifest_path=zee.MERGED_FILE,
            original_config=original,chi2_cutoff=cutoff)
        if saved['status'] not in camp.REFERENCE_ACCEPTED:
            raise ValueError(f'{model}: reference is not an accepted fit')
        info[model]=(original,fit,paths,saved,cfg)
    return info


def validate_resume(report: dict, *, models: list[str], gamma: list[float],
                    reference_run_id: str, reference_sha: str,
                    canonical_run: str):
    expected={
        'schema':SCHEMA,
        'uniform_quartic':QUARTIC,
        'conjugate_yukawas_uv_zero':True,
        'reference_run_id':reference_run_id,
        'reference_summary_sha256':reference_sha,
        'canonical_run_id':canonical_run,
        'model_order':models,
        'gamma_grid':gamma,
    }
    for k,v in expected.items():
        if report.get(k)!=v:
            raise ValueError(f'Resume provenance mismatch at {k}; do not mix studies')


def stored_points(path: Path, models: list[str], gamma: list[float]) -> dict:
    if not path.is_file():
        return {}
    desired={(m,d,round(float(g),12)) for m in models for d in DIRECTIONS for g in gamma}
    result={}
    with path.open(newline='',encoding='utf-8-sig') as file:
        for row in csv.DictReader(file):
            k=row_key(row)
            if k not in desired or k in result:
                raise ValueError(f'Foreign or duplicate saved scan point {k}')
            if row.get('status') not in ('Success','Failed'):
                raise ValueError(f'Invalid prior point status: {k}')
            if row['status']=='Success':
                for name in ('chi2',*COLUMNS[5:10]):
                    if not math.isfinite(float(row[name])):
                        raise ValueError(f'Invalid saved numerical result {k} {name}')
            result[k]=row
    return result


def summarize(models: list[str], gamma: list[float], store: dict,
              ref_info: dict) -> list[dict]:
    rows=[]
    for model in models:
        expected={(model,d,round(float(g),12)) for d in DIRECTIONS for g in gamma}
        observed=[store[k] for k in expected if k in store]
        successful=sum(r['status']=='Success' for r in observed)
        failures=sum(r['status']!='Success' for r in observed)
        if successful==len(expected):status='Complete'
        elif failures:status='PartialSuccess'
        else:status='Pending'
        rows.append({'model':model,'status':status,'failed_points':failures,
                     'computed_points':len(observed),'successful_points':successful,
                     'matching':'conditional general one-scale C5',
                     'baseline_chi2_saved':float(ref_info[model][3]['chi2']),
                     'new_reference_status':ref_info[model][3]['status']})
    return rows


def report_payload(models: list[str], gamma: list[float], ref_id: str,
                   summary_sha: str, canonical_run: str, model_summaries: list[dict]) -> dict:
    return {
        'schema':SCHEMA,'models':model_summaries,'model_order':models,
        'gamma_grid':gamma,'uniform_quartic':QUARTIC,
        'conjugate_yukawas_uv_zero':True,'reference_run_id':ref_id,
        'reference_summary_sha256':summary_sha,'canonical_run_id':canonical_run,
        'nominal_invariance':'lambdaT3*(Y1 @ Y2.T + Y2 @ Y1.T) in UV',
        'independent_scalar_quartics_fixed':True,
        'zero_conjugate_yukawas_at_UV':True,
        'gamma_directions':{
            'T1':'lambda *= gamma, Y1 *= gamma**-0.5, Y2 *= gamma**-0.5',
            'T2':'lambda *= gamma, Y1 fixed, Y2 *= gamma**-1',
            'T3':'lambda fixed, Y1 *= gamma, Y2 *= gamma**-1'},
        'full_sequential_fermion_then_scalar_matching_performed':False,
        'scope':'conditional general UV+EFT1 RG, one-scale C5',
        'no_Yukawa_refit_during_gamma_scan':True,
    }


def main(argv=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--models',default='all')
    p.add_argument('--reference-run-id',default='uniform001_all16_reference_v1')
    p.add_argument('--canonical-run-id',default='full16_general_v1')
    p.add_argument('--run-id',default='uniform001_zee_all16_v1')
    p.add_argument('--gamma-min',type=float,default=0.99)
    p.add_argument('--gamma-max',type=float,default=1.01)
    p.add_argument('--points',type=int,default=11)
    p.add_argument('--chi2-cutoff',type=float,default=5.86)
    p.add_argument('--preflight',action='store_true',help='Validate all new reference fits, no RGE calls')
    p.add_argument('--resume',action='store_true',help='Continue a checkpointed run; retry failed points')
    args=p.parse_args(argv)
    from scripts import RunT3GeneralZeeSensitivity as zee
    from scripts import RunT3UniformScalarBenchmark as uniform
    try:
        models=selected_models(args.models,zee.MODEL_ORDER)
        if not all(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',v)
                   for v in (args.reference_run_id,args.canonical_run_id,args.run_id)):
            raise ValueError('Invalid run identifier')
        if not math.isfinite(args.chi2_cutoff) or args.chi2_cutoff<=0:
            raise ValueError('Invalid chi2 cutoff')
        gamma=zee.gamma_grid(args.gamma_min,args.gamma_max,args.points)
    except ValueError as exc:
        p.error(str(exc))
    reference_folder=REF_ROOT/args.reference_run_id
    ref_sum=reference_folder/'reference_summary.json'
    if not ref_sum.is_file():
        p.error(f'Accepted uniform reference summary not found: {ref_sum}')
    summary=_read(ref_sum)
    if (summary.get('schema')!='t3-thesis-uniform-001-reference-and-scan-v1'
            or summary.get('stage')!='reference'
            or summary.get('uniform_scalar_quartic_uv')!=QUARTIC
            or summary.get('conjugate_yukawas_uv_zero') is not True):
        p.error('Saved campaign summary is not uniform-0.01 reference stage')
    try:
        info=check_fitted_references(reference_folder,models,
            cutoff=args.chi2_cutoff,canonical_run=args.canonical_run_id)
    except Exception as exc:
        p.error(f'New reference validation failed: {type(exc).__name__}: {exc}')
    for m in models:
        print(f'READY {m}: accepted new 0.01 reference chi2={info[m][3]["chi2"]:.6g}',flush=True)
    if args.preflight:
        print(f'Ready {len(models)}/{len(models)} new uniform-0.01 reference fits; no physics evaluated.')
        return 0
    dest=OUT_ROOT/args.run_id
    if dest.exists() and not args.resume:
        p.error(f'Refusing to overwrite existing sensitivity: {dest}; use --resume or new --run-id')
    if args.resume and not (dest/'summary.json').is_file():
        p.error('Cannot resume without saved summary.json')
    source_hash=_sha(ref_sum)
    if args.resume:
        report=_read(dest/'summary.json')
        try:
            validate_resume(report,models=models,gamma=gamma,
                reference_run_id=args.reference_run_id,
                reference_sha=source_hash,canonical_run=args.canonical_run_id)
            store=stored_points(dest/'points.csv',models,gamma)
        except Exception as exc:
            p.error(f'Resume incompatible with saved inputs: {exc}')
    else:
        dest.mkdir(parents=True,exist_ok=False)
        store={}
    target=zee.target_object()

    def checkpoint():
        rows=[store[k] for k in sorted(store,key=lambda k:(zee.MODEL_ORDER.index(k[0]),
                            DIRECTIONS.index(k[1]),k[2]))]
        atomic_csv(dest/'points.csv',rows)
        models_summary=summarize(models,gamma,store,info)
        atomic_json(dest/'summary.json',report_payload(models,gamma,args.reference_run_id,
                            source_hash,args.canonical_run_id,models_summary))
        return models_summary

    checkpoint()
    for model in models:
        _, original_fit, paths, reference_row, cfg=info[model]
        source_args=info[model][0]
        # Factory reads the old accepted source ONLY for matching, RGBeta,
        # registry and dimensionful bilinear metadata. Its evaluated UV state
        # is REPLACED with the newly refitted uniform-0.01 config below.
        factory=uniform.ModelFactory(model,source_args,original_fit,paths,
                                     target,args.canonical_run_id)
        factory.config=cfg
        baseline_cache=None
        for direction in DIRECTIONS:
            for g in gamma:
                key=(model,direction,round(float(g),12))
                if key in store and store[key].get('status')=='Success':
                    if float(g)==1.0 and baseline_cache is None:
                        baseline_cache={k:float(store[key][k]) for k in ('chi2',*zee.OBS)}
                    continue
                row={'model':model,'direction':direction,'gamma':float(g),
                     'status':'Failed','chi2':'','error':''}
                try:
                    if g==1.0 and baseline_cache is not None:
                        row.update(baseline_cache)
                    else:
                        altered=zee.deform(cfg,direction,g)
                        zee.verify_nominal_invariance(cfg,altered)
                        lam=float(altered['base_state']['ordinary']['lambdaT3']['real'])
                        evaluator,_,_,boundary=factory.build(lam)
                        if float(boundary.get('extra_quartic_value'))!=QUARTIC:
                            raise ValueError('Independent scalar boundary no longer uniform 0.01')
                        fit=evaluator(zee.params_from(altered))
                        chi=float(fit.chi2)
                        if not math.isfinite(chi):
                            raise ValueError('Nonfinite chi2')
                        row.update(zee.prediction(fit),chi2=chi)
                        if g==1.0:
                            baseline_cache={k:row[k] for k in ('chi2',*zee.OBS)}
                    row['status']='Success'
                except Exception as exc:
                    row['error']=f'{type(exc).__name__}: {exc}'[:750]
                store[key]=row
                checkpoint()
                print(f'{model} {direction} gamma={g:g}: {row["status"]}',flush=True)
        centres=[store.get((model,d,1.0),{}) for d in DIRECTIONS]
        if all(r.get('status')=='Success' for r in centres):
            try:
                zee.baseline_consistency(centres,reference_row)
            except Exception as exc:
                # A mismatch is a genuine source/physics blocker: mark centres
                # failed so plotting cannot silently accept mixed baselines.
                for d in DIRECTIONS:
                    key=(model,d,1.0)
                    store[key]={**store[key],'status':'Failed',
                                'error':f'Reference baseline mismatch: {exc}'[:750]}
                checkpoint()
                print(f'{model}: reference baseline mismatch: {exc}',flush=True)
        st=next(s['status'] for s in checkpoint() if s['model']==model)
        print(f'{model}: {st} ({len(gamma)*3} points expected)',flush=True)
    stats=checkpoint()
    n=sum(s['status']=='Complete' for s in stats)
    print(f'{n}/{len(models)} complete uniform-0.01 sensitivity models: {dest}',flush=True)
    return 0 if n==len(models) else 2


if __name__=='__main__':
    raise SystemExit(main())
