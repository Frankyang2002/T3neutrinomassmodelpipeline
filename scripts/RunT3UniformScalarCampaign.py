"""Fresh 16-representation T3 Yukawa fits for uniform BSM quartics = +0.01.

Stage 1 (reference): refit all 18 *real* Y1/Y2 entries at each model's
accepted reference lambdaT3 (0.1 for A-D; 0.05 for E), but replace ALL
permitted BSM scalar quartics other than lambdaT3 by +0.01 at the UV scale.
All symmetry-permitted conjugate Yukawa matrices have UV boundary ZERO.
The original accepted fits are seeds, NOT accepted new-benchmark fits.

Stage 2 (scan): only after an accepted new reference, independently refit
Y1/Y2 at lambdaT3 = (0.0005, 0.005, 0.5) using the NEW reference Yukawas.
Absolute +0.01 scalar quartics stay fixed at every scan point. The original
SM Higgs quartic, SM inputs, masses, thresholds and mixed scalar bilinear
boundary prescription remain unchanged.

Delegates RGE, conditional one-scale C5 matching, oscillation observables and
all 18-param least squares to the existing local physics scripts. This does
NOT certify fully sequential fermion -> scalar -> SMEFT matching.

No GitHub edits; no existing data overwritten. Run from project root:
  python -X utf8 -m scripts.RunT3UniformScalarCampaign --stage preflight --models pilot
  python -X utf8 -m scripts.RunT3UniformScalarCampaign --stage reference --models pilot --run-id new001_pilot_v1 --max-nfev 5
  python -X utf8 -m scripts.RunT3UniformScalarCampaign --stage scan --models pilot --run-id new001_pilot_v1 --max-nfev 5
"""
from __future__ import annotations

import argparse
import csv
from copy import deepcopy
import json
import math
from pathlib import Path
import re
import sys
from typing import Any

import numpy as np

from scripts import RunT3QualitativeScalarRefits as base
from scripts import RunT3UniformScalarBenchmark as uniform

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output/studies/thesis_qualitative/uniform_scalar_001/campaigns'
SCHEMA = 't3-thesis-uniform-001-reference-and-scan-v1'
REFERENCE_ACCEPTED = frozenset(('AcceptedApproxFit', 'AcceptedRescaledSeed'))
CSV_FIELDS = ('model','stage','lambdaT3','status','chi2',*base.OBS,
              'm1_ev','m2_ev','m3_ev',
              'y1_frobenius','y2_frobenius','combined_frobenius',
              'norm_ratio_to_reference','reason')


def select_models(value: str) -> tuple[str, ...]:
    if value == 'pilot':
        selected = ('A_-4','B_-1')
    elif value == 'all':
        selected = base.MODELS
    else:
        selected = tuple(s.strip() for s in value.split(',') if s.strip())
    if not selected or len(set(selected)) != len(selected) or set(selected)-set(base.MODELS):
        raise ValueError('Invalid, empty, repeated or unknown models')
    return selected


def source_fingerprints(paths: dict[str, Path]) -> dict[str, str]:
    return {k:base.sha(v) for k,v in sorted(paths.items())}


def provenance(paths: dict[str, Path], manifest_path: Path) -> dict:
    return {
        'source_sha256':source_fingerprints(paths),
        'source_paths':{k:v.relative_to(ROOT).as_posix() for k,v in sorted(paths.items())},
        'accepted_manifest_sha256':base.sha(manifest_path),
    }


def check_uniform_reference(config: dict, original: dict) -> None:
    """Ensure new reference changes only allowed scalar quartics and Y1/Y2.

    Does not silently accept files with old scalar boundaries or nonzero
    conjugate Yukawa entries. Also checks the UV gauge/SM/mass provenance.
    """
    reference_lam = base.real_lambda(original)
    if not math.isclose(base.real_lambda(config), reference_lam, rel_tol=0,abs_tol=1e-14):
        raise ValueError('Reference lambdaT3 changed')
    expected,_ = uniform.uniform_config(original,reference_lam)
    for key in ('representation','scales'):
        if config.get(key) != expected.get(key):
            raise ValueError(f'Reference {key} changed')
    if config['base_state']['mu_gev'] != expected['base_state']['mu_gev']:
        raise ValueError('Reference UV matching scale changed')
    if config['base_state']['sm'] != expected['base_state']['sm']:
        raise ValueError('Reference SM inputs changed')
    ordinary = config['base_state']['ordinary']
    intended = expected['base_state']['ordinary']
    for key, value in intended.items():
        if key in ('y1','y2'):
            continue
        if key not in ordinary or ordinary[key] != value:
            raise ValueError(f'Non-Yukawa reference boundary changed or missing: {key}')
    if set(ordinary) != set(intended):
        raise ValueError('Reference UV coupling inventory changed')
    # Refitted Yukawas must retain the 18-real-parameter convention.
    x=base.fitted_vector(config)
    if not np.isfinite(x).all():
        raise ValueError('Nonfinite refitted reference Yukawas')


def validated_reference(folder: Path, *, model: str, source_paths: dict[str,Path],
                        manifest_path: Path, original_config: dict,
                        chi2_cutoff: float) -> tuple[dict, dict]:
    row = base.json_read(folder/'result.json')
    cfg = base.json_read(folder/'fitted_config.json')
    if row.get('schema')!=SCHEMA or row.get('stage')!='reference' or row.get('model')!=model:
        raise ValueError('Wrong reference study, stage or model')
    if row.get('status') not in REFERENCE_ACCEPTED:
        raise ValueError(f'New reference fit is not accepted: {row.get("status")}')
    if not (0 <= float(row['chi2']) <= chi2_cutoff and math.isfinite(float(row['chi2']))):
        raise ValueError('New reference fit fails requested chi2 cutoff')
    if row.get('provenance') != provenance(source_paths,manifest_path):
        raise ValueError('Original accepted-source fingerprints changed')
    if cfg.get('thesis_uniform_scalar_benchmark',{}).get('schema') != uniform.SCHEMA:
        raise ValueError('Saved refit configuration lacks uniform-0.01 provenance')
    check_uniform_reference(cfg,original_config)
    got = base.fitted_vector(cfg)
    saved = row.get('fitted_yukawas')
    if not isinstance(saved,dict) or any(cfg['base_state']['ordinary'][k]!=saved.get(k)
                                             for k in ('y1','y2')):
        raise ValueError('Saved Yukawas disagree between reference result and fitted config')
    if abs(float(row['lambdaT3']) - base.real_lambda(cfg)) > 1e-14:
        raise ValueError('Reference fit lambdaT3 mismatch')
    if row.get('seed_skipped_optimizer') and row['status'] != 'AcceptedRescaledSeed':
        raise ValueError('Skipped-optimiser provenance inconsistent')
    return row,cfg


def record_point(folder: Path, row: dict, cfg: dict|None, fixed: dict|None) -> None:
    """Write result first and conditionally save complete, fitted UV input."""
    base.json_save(folder/'result.json',row)
    if row['status'] not in REFERENCE_ACCEPTED and row['status']!='PoorFit':
        return
    if cfg is None:
        raise ValueError('Successful numerical fit has no configuration')
    saved=deepcopy(cfg)
    saved['base_state']['ordinary'].update(row['fitted_yukawas'])
    saved['thesis_uniform_scalar_benchmark'] = {
        'schema':uniform.SCHEMA,'quartic':uniform.QUARTIC,
        'lambdaT3':row['lambdaT3'],'conjugate_yukawas_zero':True,
        'matching_scope':'conditional-one-scale','campaign_schema':SCHEMA,
        'stage':row['stage']}
    base.json_save(folder/'fitted_config.json',saved)
    if fixed is not None:
        from scripts import RunT3GeneralAll16FromReference as general
        full = general.override_yukawas(fixed,base.parameter_dict(base.fitted_vector(saved)))
        base.json_save(folder/'general_uv_overrides.json',full)


def csv_row(r: dict) -> dict:
    flat={k:r.get(k,'') for k in ('model','stage','lambdaT3','status','chi2','reason','norm_ratio_to_reference')}
    flat.update(r.get('prediction',{}))
    flat.update(r.get('yukawa_norms',{}))
    for i,m in enumerate(r.get('neutrino_masses_ev',[]),1):
        flat[f'm{i}_ev']=m
    return {name:flat.get(name,'') for name in CSV_FIELDS}


def save_report(run_folder: Path, stage: str, rows: list[dict], models: list[dict],
                requested: tuple[str,...], args) -> None:
    target = run_folder / f'{stage}_summary.json'
    base.json_save(target, {
        'schema':SCHEMA,'stage':stage,'uniform_scalar_quartic_uv':uniform.QUARTIC,
        'conjugate_yukawas_uv_zero':True,'sm_higgs_quartic_unchanged':True,
        'sequential_eft_matching_certified':False,
        'model_order':list(requested),'model_summaries':models,
        'accepted_points':sum(r.get('status') in REFERENCE_ACCEPTED for r in rows),
        'computed_points':len(rows), 'max_nfev':int(args.max_nfev),
        'accept_seed_enabled':bool(args.accept_seed),
        'chi2_cutoff':float(args.chi2_cutoff),
        'scan_values':list(base.BENCHMARKS) if stage=='scan' else None})
    with (run_folder/f'{stage}_comparison.csv').open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(csv_row(r) for r in rows)


def process_model(model: str, stage: str, folder: Path, args, selection: dict,
                  target: Any, zee: Any) -> tuple[list[dict],dict]:
    """Read accepted source, construct one cached evaluator factory per model."""
    source_args=argparse.Namespace(fit_run_id='merged',canonical_run_id=args.canonical_run_id)
    original_config, original_fit, source_paths = zee.sources(model,source_args,selection[model])
    origin = provenance(source_paths,zee.MERGED_FILE)
    original_lambda = base.real_lambda(original_config)
    if stage=='reference':
        config=original_config
        vector=base.fitted_vector(original_config)
        values=(original_lambda,)
    else:
        row,cfg=validated_reference(folder/'reference'/base.tag(model),model=model,
            source_paths=source_paths,manifest_path=zee.MERGED_FILE,
            original_config=original_config,chi2_cutoff=args.chi2_cutoff)
        config=cfg
        vector=base.fitted_vector(cfg)
        values=base.BENCHMARKS
    factory=uniform.ModelFactory(model,original_config,original_fit,source_paths,
                                  target,args.canonical_run_id)
    factory.config=config
    output=[]
    for lam in values:
        print(f'  {stage}: lambdaT3={lam:g}, BSM quartics={uniform.QUARTIC:g}, Yconj=0',flush=True)
        result,changed,fixed = uniform.run_point(factory,lam,args,vector)
        result.update({'schema':SCHEMA,'stage':stage,'provenance':origin,
                       'original_source_lambdaT3':original_lambda,
                       'original_source_chi2':float(original_fit['chi2_final']),
                       'new_reference_chi2':float(row['chi2']) if stage=='scan' else None,
                       'seed_skipped_optimizer':result['status']=='AcceptedRescaledSeed',
                       'benchmark_scalar_policy':'all-nonSM-quartics-absolute-0.01'})
        if stage=='reference':
            dest=folder/'reference'/base.tag(model)
        else:
            dest=folder/'scan'/base.tag(model)/base.lambda_tag(lam)
        record_point(dest,result,changed,fixed)
        output.append(result)
        chi=f', chi2={result["chi2"]:.6g}' if 'chi2' in result else f', {result.get("reason","")[:160]}'
        print('    '+result['status']+chi,flush=True)
    if stage=='reference':
        summary={'model':model,'stage':stage,
                 'status':'Accepted' if output[0]['status'] in REFERENCE_ACCEPTED else 'Rejected',
                 'reference_lambdaT3':original_lambda,
                 'reference_chi2':output[0].get('chi2'),
                 'original_source_chi2':float(original_fit['chi2_final'])}
    else:
        summary=uniform.model_summary(output)
        summary['stage']='scan'
        summary['new_reference_chi2']=float(row['chi2'])
    return output,summary


def main(argv=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage',choices=('preflight','reference','scan','both'),default='reference')
    p.add_argument('--models',default='pilot',help='pilot, all, or comma-separated 16 model names')
    p.add_argument('--run-id',default='uniform001_pilot_refit_v1')
    p.add_argument('--canonical-run-id',default='full16_general_v1')
    p.add_argument('--max-nfev',type=int,default=5)
    p.add_argument('--bound',type=float,default=0.1)
    p.add_argument('--regularization',type=float,default=1e-3)
    p.add_argument('--chi2-cutoff',type=float,default=5.86)
    p.add_argument('--accept-seed',action='store_true',
                   help='Optional FAST mode: skip optimisation on an already accepted heuristic seed; mark explicitly as unoptimised')
    args=p.parse_args(argv)
    try:
        wanted=select_models(args.models)
    except ValueError as exc:
        p.error(str(exc))
    if (not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.run_id) or
        not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args.canonical_run_id)):
        p.error('Invalid run identifier')
    if (not (math.isfinite(args.bound) and args.bound>0) or
        not (math.isfinite(args.regularization) and 0<=args.regularization<=0.1) or
        not (math.isfinite(args.chi2_cutoff) and args.chi2_cutoff>0) or
        args.max_nfev<1):
        p.error('Invalid optimisation limits')
    from scripts import RunT3GeneralZeeSensitivity as zee
    selection=zee.fit_selections(list(wanted),'merged')
    folder=OUT/args.run_id
    if args.stage=='preflight':
        for model in wanted:
            original,fit,paths=zee.sources(model,
                argparse.Namespace(fit_run_id='merged',canonical_run_id=args.canonical_run_id),
                selection[model])
            reference_lam=base.real_lambda(original)
            boundary,names=uniform.uniform_config(original,reference_lam)
            print(f'{model}: original lambdaT3={reference_lam:g}; '
                  f'new quartics={uniform.QUARTIC:g}; ordinary {names}; '
                  f'Yconj=0; saved chi2={fit["chi2_final"]:.7g}',flush=True)
            if model in base.MIXED:
                from scripts import RunT3BCGeneralFromReference as br
                symbols=uniform.fixed_bc_extra(br.build_reference_source(model))
                print('  independent B/C extra quartics: '+','.join(sorted(symbols)),flush=True)
            else:
                fixed,extra,conj=uniform.uniform_overrides(base.json_read(paths['overrides']))
                print('  UV-layout quartics: '+','.join(extra),flush=True)
                print('  zeroed conjugates: '+(','.join(conj) or '(none)'),flush=True)
        print('Read-only preflight; no physics evaluations.',flush=True)
        return 0
    if args.stage in ('reference','both'):
        if (folder/'reference').exists() or (folder/'reference_summary.json').exists():
            p.error(f'Reference outputs already exist; use new --run-id: {folder}')
    if args.stage=='scan':
        if not (folder/'reference').is_dir():
            p.error(f'No reference fits in {folder}; run --stage reference first')
        if (folder/'scan').exists() or (folder/'scan_summary.json').exists():
            p.error('Scan outputs already exist; choose a fresh campaign --run-id')
    if args.stage=='both' and ((folder/'scan').exists() or (folder/'scan_summary.json').exists()):
        p.error('Scan outputs already exist; choose fresh --run-id')
    folder.mkdir(parents=True,exist_ok=True)
    target=zee.target_object()
    stages=('reference','scan') if args.stage=='both' else (args.stage,)
    overall_ok=True
    for stage in stages:
        (folder/stage).mkdir(parents=True,exist_ok=False)
        rows,model_summaries=[],[]
        for model in wanted:
            print(f'{model}: {stage} under uniform 0.01 UV scalar benchmark',flush=True)
            try:
                chunk,summary=process_model(model,stage,folder,args,selection,target,zee)
                rows.extend(chunk)
                model_summaries.append(summary)
                good=all(r['status'] in REFERENCE_ACCEPTED for r in chunk)
                overall_ok &=good
            except Exception as exc:
                problem={'model':model,'stage':stage,'status':'Blocked',
                         'reason':f'{type(exc).__name__}: {exc}'[:1500]}
                model_summaries.append(problem)
                overall_ok=False
                print('  BLOCKED: '+problem['reason'],flush=True)
            # Flush progressive result summaries even if later models fail.
            save_report(folder,stage,rows,model_summaries,wanted,args)
        print(f'{stage}: accepted {sum(r["status"] in REFERENCE_ACCEPTED for r in rows)}/'
              f'{len(rows)}; outputs: {folder}',flush=True)
    return 0 if overall_ok else 2


if __name__=='__main__':
    raise SystemExit(main())
