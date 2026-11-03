"""All-16 independent Matchete heavy/scalar common-log vs EFT1 mixing audit.

Reuses fitted endpoint snapshots and existing general_seq_05 *independent* LLSS
contraction, but never uses the subtracted hard-matching bookkeeping residual.
For 14 nonmixed models, evaluates original and log-spurion Matchete expressions
with the existing Wolfram flavour-channel evaluator. For B_-1/C_-1, regenerates
only their fitted U=E=1 endpoints and uses their existing compiled one-scale
expression. Optional MF->(1,2,4)*MF scaling tests the hierarchy limit.

No fitted files or existing ablation outputs are modified. Run from repo root:
  python -X utf8 -m scripts.RunT3GeneralMatcheteLogCoefficientAudit --check-only
  python -X utf8 -m scripts.RunT3GeneralMatcheteLogCoefficientAudit --run-id log_check_01
  python -X utf8 -m scripts.RunT3GeneralMatcheteLogCoefficientAudit --run-id log_BC_01 --models B_-1,C_-1
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import re
import shutil
import traceback

from Numerical.running.GeneralT3MatcheteLogCoefficientAudit import (
    MODELS,MIXED,FORMAT,require,read_json,sha256_file,
    saved_case_input,replace_common_mass_logs,analyse_one_model,
)

ROOT=Path(__file__).resolve().parents[1]
FIT_ROOT=ROOT/'output/studies/general_scalar_extended/thesis_general_all16'
CASE_ROOT=ROOT/'output/studies/general_scalar_extended/general_sequential_rge_ablation'
OUT_ROOT=ROOT/'output/studies/general_scalar_extended/general_log_coefficient_audit'
MIXED_SOURCE=ROOT/'output/general_t3_mass_eigen_matchete/one_loop_probe'


def tag(model: str) -> str:
    return model.replace('+','p').replace('-','m')


def coefficient_path(model: str) -> Path:
    if model in MIXED:
        return MIXED_SOURCE / ({'B_-1':'B_m1_final',
                    'C_-1':'C_m1_final'}[model]+
                    '_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')
    from scripts import RunT3GeneralAll16FromReference as general
    return general.general_coefficient_path(model)


def validate_filenames(run: str, what: str) -> None:
    require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',run) is not None,
            f'{what} must contain only safe filename characters')


def inputs(model: str, fit_run: str, case_run: str, *,
           wolframscript: str) -> dict:
    saved=FIT_ROOT/fit_run/tag(model)
    case=CASE_ROOT/case_run/tag(model)/'cases.json'
    coeff=coefficient_path(model)
    paths={'fit_result':saved/'result.json',
           'fitted_config':saved/'fitted_config.json',
           'saved_ablation_cases':case,'one_scale_Matchete_source':coeff}
    errors=[]
    missing=[key for key,p in paths.items() if not p.is_file()]
    if missing:
        errors.append('missing file(s): '+', '.join(f'{k}: {paths[k]}' for k in missing))
    if not missing:
        try:
            saved_case_input(case,model)
            fit=read_json(paths['fit_result'])
            require(fit.get('model')==model,'fit result model mismatch')
            transformed,loginfo=replace_common_mass_logs(coeff.read_text(encoding='utf-8-sig'))
            require(transformed!=coeff.read_text(encoding='utf-8-sig'),
                    'log shift did not change expression')
            if model not in MIXED:
                detail=fit.get('numerical_details',{})
                require(all(k in detail for k in
                        ('uv_at_fermion_threshold','eft1_at_scalar_threshold')),
                        'nonmixed fit missing saved U/E running endpoints')
        except Exception as exc:
            errors.append(f'input validation: {type(exc).__name__}: {exc}')
    if model not in MIXED and not (shutil.which(wolframscript) or Path(wolframscript).is_file()):
        errors.append(f'Wolfram executable not available: {wolframscript}')
    return {'model':model,'ready':not errors,'errors':errors,
            'input_paths':{k:str(v) for k,v in paths.items()},
            'one_scale_source_sha256':sha256_file(coeff) if coeff.is_file() else None}


def row(result: dict) -> dict:
    model=result['model']
    if result.get('status')!='Computed':
        return {'model':model,'status':result.get('status','Blocked'),
                'error':result.get('error','')}
    scan=result['hierarchy_scan']
    head=scan[0]['comparison']
    tail=scan[-1]['comparison']
    return {'model':model,'status':'Computed',
            'matchete_log_occurrences':result['source_provenance']['log_inventory']['log_occurrences'],
            'one_scale_reproduction_relative_residual':result['one_scale_reproduction_relative_residual'],
            'ratio_Matchete_over_EFT_at_MF':head['matchete_over_eft_norm'],
            'best_fit_factor_real_at_MF':head['best_fit_complex_factor']['real'],
            'best_fit_factor_imag_at_MF':head['best_fit_complex_factor']['imag'],
            'relative_matrix_diff_factor1_at_MF':head['relative_matrix_difference_at_assumed_factor_1'],
            'relative_matrix_residual_best_factor_at_MF':head['relative_matrix_residual_after_best_fit_factor'],
            'largest_MF_scale_tested':scan[-1]['heavy_mass_scale_factor'],
            'best_fit_factor_real_at_largest_MF':tail['best_fit_complex_factor']['real'],
            'relative_matrix_diff_factor1_at_largest_MF':tail['relative_matrix_difference_at_assumed_factor_1'],
            'error':''}


def save_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def parser() -> argparse.ArgumentParser:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-run',default='thesis16_merged_v1',
                   help='Unified fitted campaign in thesis_general_all16')
    p.add_argument('--case-run',default='general_seq_05',
                   help='Existing all-16 sequential ablation source')
    p.add_argument('--run-id',default='matchete_log_audit_01',
                   help='New independent audit output basename')
    p.add_argument('--models',default=None,
                   help='Comma-separated T3 model names; default all 16')
    p.add_argument('--hierarchy-factors',default='1,2,4',
                   help='Rescale ONLY heavy MF masses to test large-MF limit')
    p.add_argument('--wolframscript',default='wolframscript')
    p.add_argument('--timeout-sec',type=int,default=300)
    p.add_argument('--check-only',action='store_true')
    return p


def main(argv=None) -> int:
    args=parser().parse_args(argv)
    for name,run in [('source-run',args.source_run),('case-run',args.case_run),
                     ('run-id',args.run_id)]:
        validate_filenames(run,name)
    models=tuple(t.strip() for t in args.models.split(',') if t.strip()) if args.models else MODELS
    require(bool(models) and len(set(models))==len(models) and set(models)<=set(MODELS),
            'invalid models selection')
    try:
        factors=tuple(int(x.strip()) for x in args.hierarchy_factors.split(','))
    except ValueError as e:
        raise ValueError('Hierarchy factors must be positive integers') from e
    require(bool(factors) and factors[0]==1 and len(set(factors))==len(factors)
            and all(x>=1 for x in factors),'factors must start at 1 and not repeat')
    require(args.timeout_sec>0,'invalid Matchete channel timeout')
    require((FIT_ROOT/args.source_run).is_dir(),'fitted campaign source missing')
    require((CASE_ROOT/args.case_run).is_dir(),'ablation source missing')
    case_summary=read_json(CASE_ROOT/args.case_run/'summary.json')
    require(case_summary.get('source_run')==args.source_run,
            'ablation source-run and requested fitted campaign differ')
    assumed_norm=float(case_summary['spinor_factor_assumed'])
    require(math.isfinite(assumed_norm) and assumed_norm!=0,
            'ablation spinor normalization must be nonzero and finite')
    checks=[inputs(m,args.source_run,args.case_run,
                   wolframscript=args.wolframscript) for m in models]
    for chk in checks:
        if args.check_only or not chk['ready']:
            print(chk['model'], 'Ready' if chk['ready'] else 'Missing',
                  '; '.join(chk['errors']),flush=True)
    if args.check_only:
        return 0 if all(chk['ready'] for chk in checks) else 2
    dest=OUT_ROOT/args.run_id
    require(not dest.exists(),f'new audit output already exists: {dest}')
    dest.mkdir(parents=True)
    results=[]
    for check in checks:
        model=check['model']
        if not check['ready']:
            entry={'model':model,'status':'Blocked','error':'; '.join(check['errors'])}
        else:
            try:
                fit=FIT_ROOT/args.source_run/tag(model)
                case=CASE_ROOT/args.case_run/tag(model)/'cases.json'
                coeff=coefficient_path(model)
                entry=analyse_one_model(model,fit,case,coeff,
                           wolframscript=args.wolframscript,
                           timeout_sec=args.timeout_sec,
                           hierarchy_factors=factors)
                save_json(dest/tag(model)/'log_audit.json',entry)
            except Exception as exc:
                entry={'model':model,'status':'Blocked',
                       'error':f'{type(exc).__name__}: {exc}',
                       'traceback':traceback.format_exc(limit=8)}
                save_json(dest/tag(model)/'error.json',entry)
        print(model,entry['status'],entry.get('error','')[:180],flush=True)
        results.append(entry)
        compact=[row(x) for x in results]
        summary={'format':FORMAT,'source_run':args.source_run,
                 'case_run':args.case_run,'requested_models':list(models),
                 'models_processed':len(results),
                 'computed':sum(x['status']=='Computed' for x in results),
                 'hierarchy_factors':list(factors),
                 'normalization_assumed_in_ablation':assumed_norm,
                 'complete_independent_hard_matching':False,
                 'independent_mass_log_coefficient_extracted':True,
                 'model_rows':compact}
        save_json(dest/'summary.json',summary)
    fields=['model','status','matchete_log_occurrences',
            'one_scale_reproduction_relative_residual',
            'ratio_Matchete_over_EFT_at_MF',
            'best_fit_factor_real_at_MF','best_fit_factor_imag_at_MF',
            'relative_matrix_diff_factor1_at_MF',
            'relative_matrix_residual_best_factor_at_MF',
            'largest_MF_scale_tested','best_fit_factor_real_at_largest_MF',
            'relative_matrix_diff_factor1_at_largest_MF','error']
    with (dest/'summary.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore')
        w.writeheader();w.writerows([row(x) for x in results])
    print('Audit output:',dest/'summary.json')
    return 0 if len(results)==len(models) and all(x['status']=='Computed' for x in results) else 2

if __name__=='__main__':
    raise SystemExit(main())
