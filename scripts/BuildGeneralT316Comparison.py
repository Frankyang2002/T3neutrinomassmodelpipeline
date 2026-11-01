"""Consolidate the 14 improved direct-matching and 2 conditional archived-C5 T3 results.

This is a reporting-only utility. No RGEs or matching calculations are run.
Run from the repository root:
    python -X utf8 -m scripts.BuildGeneralT316Comparison
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

MODELS = (
    'A_-4','A_-2','A_+0','B_-3','B_-1','B_+1',
    'C_-3','C_-1','C_+1','D_-2','D_+0','D_+2',
    'E_-4','E_-2','E_+0','E_+2',
)
CONDITIONAL = {'B_-1':'B_m1', 'C_-1':'C_m1'}
FIELDS = ['model','calculation_method','status','physical_prediction_certified',
          'mu_low_gev','m1_ev','m2_ev','m3_ev',
          'delta_m21_sq_ev2','delta_m31_sq_ev2','delta_m32_sq_ev2',
          *[f'Uabs_{i}{j}' for i in (1,2,3) for j in (1,2,3)],
          'takagi_residual','source_json']


def read_json(path: Path) -> dict:
    data = path.read_text(encoding='utf-8-sig')
    try:
        value = json.loads(data)
    except json.JSONDecodeError:
        # Earlier 14-model summary writer had a literal backslash-n suffix.
        if not data.endswith('\\n'):
            raise
        value = json.loads(data[:-2])
    if not isinstance(value, dict):
        raise ValueError(f'{path}: expected JSON object')
    return value


def number(value, name):
    if type(value) not in (int, float) or not math.isfinite(float(value)):
        raise ValueError(f'{name}: expected finite number')
    return float(value)


def extract(model: str, path: Path, conditional: bool) -> dict:
    data = read_json(path)
    if data.get('model') != model:
        raise ValueError(f'{path}: incorrect model {data.get("model")!r}')
    if conditional:
        if data.get('status') != 'ConditionalArchivedC5HybridToMZ':
            raise ValueError(f'{model}: wrong conditional result status')
        for flag, expected in [('physical_neutrino_prediction_certified',False),
                               ('mixed_scalar_matching_performed',False),
                               ('archived_C5_assumed_at_mu_s',True)]:
            if data.get(flag) is not expected:
                raise ValueError(f'{model}: unexpected {flag}')
        scale = number(data['scales_gev']['low'], 'MZ scale')
        method = 'ConditionalArchivedC5'
        status = 'SuccessConditional'
    else:
        if data.get('low_energy_rge_status') != 'Success':
            raise ValueError(f'{model}: low-energy RGE unsuccessful')
        scale = number(data['observables_scale_gev'], 'MZ scale')
        method = 'ImprovedDirectMatchingProvisional'
        status = 'Success'
    if not math.isclose(scale, 91.1876, rel_tol=0, abs_tol=1e-6):
        raise ValueError(f'{model}: expected MZ endpoint, got {scale}')
    masses = data['neutrino_masses_ev']
    pmns = data['pmns_abs']
    if not isinstance(masses, list) or len(masses) != 3:
        raise ValueError(f'{model}: expected 3 masses')
    if not isinstance(pmns, list) or len(pmns) != 3 or any(
            not isinstance(r,list) or len(r)!=3 for r in pmns):
        raise ValueError(f'{model}: expected 3x3 PMNS magnitudes')
    row = {'model':model, 'calculation_method':method, 'status':status,
           'physical_prediction_certified':False, 'mu_low_gev':scale,
           'source_json':str(path)}
    for j,m in enumerate(masses,1):
        row[f'm{j}_ev']=number(m,f'{model}.m{j}')
    for key in ('delta_m21_sq_ev2','delta_m31_sq_ev2','takagi_residual'):
        row[key]=number(data[key],f'{model}.{key}')
    # Some earlier 14-model JSONs omit delta_m32; report as missing, not fabricated.
    row['delta_m32_sq_ev2']=(number(data['delta_m32_sq_ev2'],f'{model}.delta_m32')
                               if 'delta_m32_sq_ev2' in data else '')
    for i in range(3):
        for j in range(3):
            row[f'Uabs_{i+1}{j+1}']=number(pmns[i][j],f'{model}.PMNS[{i},{j}]')
    return row


def fmt(value):
    return f'{value:.8g}' if isinstance(value,float) else str(value)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--improved-dir',type=Path,
        default=Path('output/general_numerical_all14_improved/neutrinos'))
    parser.add_argument('--bc-dir',type=Path,
        default=Path('output/general_numerical_BC_archived_hybrid'))
    parser.add_argument('--output-dir',type=Path,
        default=Path('output/general_numerical_16_model_comparison'))
    args=parser.parse_args(argv)
    first=read_json(args.improved_dir/'summary.json')
    second=read_json(args.bc_dir/'summary.json')
    if second.get('physical_predictions_certified') is not False or \
            second.get('scope') != 'conditional_hybrid_archived_fixed_scale_C5_plus_improved_SM_background':
        raise ValueError('Unexpected B/C conditional summary metadata')
    improved={x['model']:x for x in first.get('models',[])}
    bc={x['model']:x for x in second.get('models',[])}
    rows=[]
    for model in MODELS:
        is_bc=model in CONDITIONAL
        rec=(bc if is_bc else improved).get(model)
        expected_status='SuccessConditional' if is_bc else 'Success'
        if rec is None or rec.get('status')!=expected_status:
            raise ValueError(f'{model}: missing {expected_status} summary record: {rec}')
        path=(args.bc_dir/f'{CONDITIONAL[model]}_conditional_neutrinos.json'
              if is_bc else args.improved_dir/model.replace('-','m').replace('+','p')/'neutrinos.json')
        rows.append(extract(model,path,is_bc))
    args.output_dir.mkdir(parents=True,exist_ok=True)
    csv_path=args.output_dir/'t3_16_model_comparison.csv'
    md_path=args.output_dir/'t3_16_model_comparison.md'
    json_path=args.output_dir/'t3_16_model_comparison.json'
    with csv_path.open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    cols=['model','calculation_method','m1_ev','m2_ev','m3_ev',
          'delta_m21_sq_ev2','delta_m31_sq_ev2']
    lines=['# T3 16-model numerical comparison','',
           '**Exploratory, not certified physical neutrino predictions.**','',
           '14 models: improved provisional direct matching and SM+Weinberg running.',
           'B(-1), C(-1): archived fixed-scale C5 imposed at 1 TeV on a separate improved SM background; mixed-scalar matching and consistency not established.',
           'All values are read from saved outputs; no physics is recalculated.','',
           '| '+' | '.join(cols)+' |',
           '| '+' | '.join('---' for _ in cols)+' |']
    for row in rows:
        lines.append('| '+' | '.join(fmt(row[k]) for k in cols)+' |')
    lines.extend(['','## Absolute PMNS matrices','',
                  'Entries are magnitudes |U_ij| (phases are not reported).',''])
    for row in rows:
        lines.extend([f'### {row["model"]} — {row["calculation_method"]}','',
                     '| i / j | 1 | 2 | 3 |','|---|---:|---:|---:|'])
        for i in (1,2,3):
            lines.append('| '+str(i)+' | '+' | '.join(fmt(row[f'Uabs_{i}{j}']) for j in (1,2,3))+' |')
        lines.append('')
    md_path.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    json_path.write_text(json.dumps({'scope':'16_model_mixed_method_exploratory_comparison',
                                    'physical_predictions_certified':False,
                                    'improved_direct_model_count':14,
                                    'conditional_archived_C5_model_count':2,
                                    'models':rows},indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print('Validated 16/16 numerical results (14 provisional direct + 2 conditional archived C5)')
    for path in (csv_path,md_path,json_path):
        print(path)
    return 0

if __name__=='__main__':
    raise SystemExit(main())
