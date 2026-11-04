"""B(-1)/C(-1) three-quartic Matchete log vs signed EFT1 decomposition.

Read-only use of the EXISTING thesis16_merged_v1 fitted endpoints, signed B/C
ablation cases, and independent common-log audit. No refit or RGE change.

  python -X utf8 -m scripts.RunT3BCQuarticChannelLogAudit --check-only
  python -X utf8 -m scripts.RunT3BCQuarticChannelLogAudit --run-id bc_channel_log_01
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re
import traceback

from Numerical.running.GeneralT3BCQuarticChannelLogAudit import (
    MODELS, CHANNELS, FORMAT, require, run_channel_comparison,
    validate_source_inventory,
)
from Numerical.running.GeneralT3MatcheteLogCoefficientAudit import (
    read_json, sha256_file, decode_matrix, replace_common_mass_logs,
    _mixed_running_endpoint,
)

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'output/studies/general_scalar_extended'
FIT_ROOT = BASE/'thesis_general_all16'
CASE_ROOT = BASE/'general_sequential_rge_ablation'
LOG_ROOT = BASE/'general_log_coefficient_audit'
OUTPUT_ROOT = BASE/'general_bc_quartic_channel_log_audit'
COEFF_ROOT = ROOT/'output/general_t3_mass_eigen_matchete/one_loop_probe'


def tag(model: str) -> str:
    return model.replace('-','m')


def inputs(model: str, args: argparse.Namespace) -> dict:
    require(model in MODELS,'unknown B/C model')
    fit = FIT_ROOT/args.source_run/tag(model)
    case = CASE_ROOT/args.case_run/tag(model)/'cases.json'
    log = LOG_ROOT/args.reference_audit/tag(model)/'log_audit.json'
    coeff = COEFF_ROOT/(tag(model)+'_final_one_loop_probe_weinberg_ordered_coefficient.inputform.txt')
    paths = {'fit_result':fit/'result.json', 'fitted_config':fit/'fitted_config.json',
             'signed_cases':case,'reference_log_audit':log,'one_scale_matchete':coeff}
    problems = [f'{key}: {p}' for key,p in paths.items() if not p.is_file()]
    if not problems:
        try:
            validate_source_inventory(coeff.read_text(encoding='utf-8-sig'))
            data = read_json(case)
            require(data.get('model') == model and data.get('status') == 'Computed'
                    and data.get('cases',{}).get('1111',{}).get('case') == '1111',
                    'saved signed U=E=W=S=1 ablation missing')
            p = data['cases']['1111']['quartic_provenance']
            require(all(p.get('bc_direct_quartic_channel_signs',{}).get(
                        'general:'+q) == -1 for q in
                        ('lambdaHHdagS2S2','lambdaHHdagS1barS1bar')),
                    'signed EFT1 quartic map missing from case')
            old = read_json(log)
            require(old.get('model') == model and old.get('status') == 'Computed',
                    'reference independent Matchete log absent')
            require(old.get('saved_case_record_sha256') == sha256_file(case) and
                    old.get('saved_fit_result_sha256') == sha256_file(fit/'result.json') and
                    old.get('saved_fitted_config_sha256') == sha256_file(fit/'fitted_config.json'),
                    'reference log audit was computed from DIFFERENT fitted/signed files')
        except Exception as exc:
            problems.append(f'failed source validation: {exc}')
    return {'model':model,'ready':not problems,'errors':problems,'paths':paths}


def _save(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def _row(model: str, factor: int, ch: dict) -> dict:
    c = ch['comparison']
    return {'model':model,'mf_scale':factor,
            'quartic':ch['quartic'],'matchete_symbol':ch['matchete_symbol'],
            'matchete_norm':c['matchete_spurion_frobenius_gev_inverse'],
            'eft1_norm':c['eft_spurion_frobenius_gev_inverse'],
            'matchete_over_eft_norm':c['matchete_over_eft_norm'],
            'best_fit_factor_real':c['best_fit_complex_factor']['real'],
            'best_fit_factor_imag':c['best_fit_complex_factor']['imag'],
            'matrix_residual_after_factor':c['relative_matrix_residual_after_best_fit_factor'],
            'relative_difference_factor1':c['relative_matrix_difference_at_assumed_factor_1']}


def parse_args(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-run',default='thesis16_merged_v1')
    p.add_argument('--case-run',default='general_seq_bc_signed_01')
    p.add_argument('--reference-audit',default='log_bc_signed_01')
    p.add_argument('--run-id',default='bc_channel_log_01')
    p.add_argument('--models',default='B_-1,C_-1')
    p.add_argument('--hierarchy-factors',default='1',
                   help='Optional MF scaling, e.g. 1,2,4; default only 1')
    p.add_argument('--check-only',action='store_true')
    return p.parse_args(argv)


def main(argv=None) -> int:
    args=parse_args(argv)
    for s in (args.source_run,args.case_run,args.reference_audit,args.run_id):
        require(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',s) is not None,
                f'invalid safe run-id: {s}')
    models=tuple(x.strip() for x in args.models.split(',') if x.strip())
    require(models and len(models) == len(set(models)) and set(models)<=set(MODELS),
            'expected B_-1 and/or C_-1 without duplicates')
    try:
        scales=tuple(int(s.strip()) for s in args.hierarchy_factors.split(','))
    except ValueError as exc:
        raise ValueError('hierarchy factors must be integers') from exc
    require(scales and scales[0] == 1 and len(set(scales)) == len(scales) and
            all(k>0 for k in scales), 'hierarchy factors must start at 1')
    checks=[inputs(m,args) for m in models]
    for chk in checks:
        print(chk['model'],'Ready' if chk['ready'] else 'Blocked',
              '; '.join(chk['errors']),flush=True)
    if args.check_only:
        return 0 if all(chk['ready'] for chk in checks) else 2
    dest=OUTPUT_ROOT/args.run_id
    require(not dest.exists(),f'refusing to overwrite output: {dest}')
    dest.mkdir(parents=True)
    results=[]; rows=[]
    for chk in checks:
        model=chk['model']
        if not chk['ready']:
            result={'model':model,'status':'Blocked','error':'; '.join(chk['errors'])}
        else:
            try:
                from scripts import RunT3GeneralScalarBCExploratoryFit as bc
                from scripts.RunGeneralT3NumericalC5 import compile_coefficient
                fit=chk['paths']['fit_result'].parent
                case=read_json(chk['paths']['signed_cases'])['cases']['1111']
                prior=read_json(chk['paths']['reference_log_audit'])
                source=chk['paths']['one_scale_matchete'].read_text(encoding='utf-8-sig')
                shifted,log_inventory=replace_common_mass_logs(source)
                original_compiled=compile_coefficient(source)
                shifted_compiled=compile_coefficient(shifted)
                endpoint,heavy,endpoint_info=_mixed_running_endpoint(model,fit)
                result=run_channel_comparison(model,case=case,endpoint=endpoint,
                    heavy=heavy,coefficient_text=source,
                    original_compiled=original_compiled,shifted_compiled=shifted_compiled,
                    prior_matchete=decode_matrix(
                        prior['hierarchy_scan'][0]['matchete_common_log_coefficient']),
                    hierarchy_factors=scales)
                result['provenance']={
                    'fit_result_sha256':sha256_file(chk['paths']['fit_result']),
                    'fitted_config_sha256':sha256_file(chk['paths']['fitted_config']),
                    'signed_cases_sha256':sha256_file(chk['paths']['signed_cases']),
                    'prior_log_audit_sha256':sha256_file(chk['paths']['reference_log_audit']),
                    'one_scale_source_sha256':sha256_file(chk['paths']['one_scale_matchete']),
                    'endpoint_reconstruction':endpoint_info,
                    'log_shift_inventory':log_inventory,
                    'no_existing_outputs_modified':True,
                }
                _save(dest/tag(model)/'channel_audit.json',result)
                for scan in result['hierarchy_scan']:
                    for channel in scan['channels']:
                        rows.append(_row(model,scan['heavy_mass_scale_factor'],channel))
                for channel in result['hierarchy_scan'][0]['channels']:
                    c=channel['comparison']
                    print(f'{model:4s} {channel["quartic"]:28s} '
                          f'factor={c["best_fit_complex_factor"]["real"]:.6g} '
                          f'residual={100*c["relative_matrix_residual_after_best_fit_factor"]:.4f}%',
                          flush=True)
            except Exception as exc:
                result={'model':model,'status':'Blocked','error':f'{type(exc).__name__}: {exc}',
                        'traceback':traceback.format_exc(limit=9)}
                _save(dest/tag(model)/'error.json',result)
                print(model,'Blocked',result['error'],flush=True)
        results.append(result)
        summary={'format':FORMAT,'source_run':args.source_run,
                 'case_run':args.case_run,'reference_audit':args.reference_audit,
                 'requested_models':list(models),
                 'computed':sum(x['status']=='Computed' for x in results),
                 'model_rows':[{'model':x['model'],'status':x['status'],
                      'error':x.get('error',''),
                      'total_at_MF':(x['hierarchy_scan'][0]['total_comparison']
                          if x['status']=='Computed' else None),
                      'closure_at_MF':(x['hierarchy_scan'][0]['checks']
                          if x['status']=='Computed' else None)} for x in results],
                 'conditional_quartic_sign_map_certified':False,
                 'hard_matching_independently_calculated':False}
        _save(dest/'summary.json',summary)
    if rows:
        cols=list(rows[0]);out=dest/'summary.csv'
        with out.open('w',newline='',encoding='utf-8') as fh:
            writer=csv.DictWriter(fh,fieldnames=cols)
            writer.writeheader();writer.writerows(rows)
    print('Output:',dest/'summary.json')
    return 0 if all(x['status']=='Computed' for x in results) else 2


if __name__=='__main__':
    raise SystemExit(main())
