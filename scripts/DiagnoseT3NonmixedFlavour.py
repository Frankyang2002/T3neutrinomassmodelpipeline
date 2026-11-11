"""Diagnose EFT1-induced flavour changes in saved corrected LLSS ablations.

Project path: scripts/DiagnoseT3NonmixedFlavour.py
Read-only: uses previously computed U/E/W/S cases and source-resolved
LLSS channel decomposition; does not rerun matching or evolution.

At fixed UV=1, SMEFT=1, compare EFT1 on/off for W=1 and W=0.
The channel calculations distinguish the quartic and fermion-tree
*alternative* decompositions: never sum across both categories.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'output/studies/general_scalar_extended/corrected_llss_eft1_freeze'
REPORT = ROOT / 'Reports/output/thesis_qualitative/nonmixed_eft1_flavour'
OBS = ('sin2_theta12', 'sin2_theta13', 'sin2_theta23',
       'delta_m21_sq_ev2', 'delta_m3l_sq_ev2')


def matrix(raw):
    if isinstance(raw, dict) and 'real' in raw and 'imag' in raw:
        z = np.asarray(raw['real'], float) + 1j * np.asarray(raw['imag'], float)
    else:
        z = np.asarray(raw, complex)
    if z.shape != (3, 3) or not np.isfinite(z).all():
        raise ValueError('Expected finite complex 3x3 matrix')
    return z


def compare(a, b):
    """Optimal complex factor a ~= r*b and relative Frobenius residual."""
    na = float(np.linalg.norm(a)); nb = float(np.linalg.norm(b))
    if na < 1e-290 or nb < 1e-290:
        return {'scale_real': None, 'scale_imag': None,
                'nonscaling': None, 'on_norm': na, 'off_norm': nb}
    r = np.vdot(b, a) / np.vdot(b, b)
    return {'scale_real': float(r.real), 'scale_imag': float(r.imag),
            'nonscaling': float(np.linalg.norm(a-r*b)/na),
            'on_norm': na, 'off_norm': nb}


def read_cases(path):
    doc = json.loads(path.read_text(encoding='utf-8-sig'))
    if doc.get('status') != 'Computed':
        raise ValueError(f'Not computed: {path}')
    cases = doc['cases']
    for key in ('1111', '1011', '1101', '1001'):
        if key not in cases:
            raise ValueError(f'Missing {key}: {path}')
    return cases


def channels(case):
    groups = {}
    for d in case['per_general_interaction']:
        category = d['category']; name = d['name']
        if category not in ('quartic_decomposition', 'fermion_tree_decomposition'):
            raise ValueError(f'Unknown channel category {category}')
        key = (category, name)
        if key in groups:
            raise ValueError(f'Duplicate channel {key}')
        groups[key] = d
    return groups


def write_csv(path, rows):
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='', encoding='utf-8') as f:
        w=csv.DictWriter(f, fieldnames=keys)
        w.writeheader(); w.writerows(rows)


def evaluate(model, cases):
    output=[]; detail=[]
    for wilson, on_key, off_key in ((1,'1111','1011'), (0,'1101','1001')):
        on=cases[on_key]; off=cases[off_key]
        row={'model':model, 'W':wilson,'on_case':on_key,'off_case':off_key}
        for label, field in [('total','c5_at_scalar_threshold'),
                             ('bubble','c5_joint_scalar_bubble'),
                             ('log','c5_induced_wilson_log')]:
            row.update({label+'_'+k:v for k,v in
                        compare(matrix(on[field]),matrix(off[field])).items()})
        for x in OBS:
            row['delta_'+x]=float(on['prediction'][x])-float(off['prediction'][x])
        output.append(row)
        oc, fc = channels(on), channels(off)
        if set(oc) != set(fc):
            raise ValueError(f'{model}: channel identity changed across EFT1 on/off: '
                             f'{sorted(set(oc)^set(fc))}')
        for (category,name),a in sorted(oc.items()):
            b=fc[(category,name)]
            d={'model':model,'W':wilson,'category':category,'channel':name}
            # 'total' per individual contribution is log+bubble, not the sum
            # over categories. These matrices are source-resolved from the evaluator.
            for label,field in [('total','c5_total'),('bubble','c5_joint_scalar_bubble'),
                                ('log','c5_leading_log')]:
                d.update({label+'_'+k:v for k,v in
                          compare(matrix(a[field]),matrix(b[field])).items()})
            detail.append(d)
    return output, detail


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id',default='corrected_llss16_v2')
    p.add_argument('--models',default='A_-4',help='Comma-separated models, e.g. A_-4,D_-2')
    p.add_argument('--source',type=Path,default=None)
    p.add_argument('--output',type=Path,default=None)
    args=p.parse_args()
    source=(args.source if args.source else BASE/args.run_id).resolve()
    dest=(args.output if args.output else REPORT/args.run_id).resolve()
    if not source.is_dir():p.error(f'Source run missing: {source}')
    if dest == source or source in dest.parents or dest in source.parents:
        p.error('Diagnostic output must be outside source run')
    models=[x.strip() for x in args.models.split(',') if x.strip()]
    if not models or len(models)!=len(set(models)):
        p.error('Specify distinct models')
    rows=[]; channels_out=[]; missing=[]
    for model in models:
        file=source/model.replace('+','p').replace('-','m')/'cases.json'
        if not file.is_file():
            missing.append(model)
            print(f'{model}: no computed cases at {file}')
            continue
        block,detail=evaluate(model,read_cases(file))
        rows.extend(block); channels_out.extend(detail)
        for r in block:
            print(f"{model} W={r['W']} total_nonscaling={r['total_nonscaling']:.6g} "
                  f"bubble_nonscaling={r['bubble_nonscaling']!s} "
                  f"log_nonscaling={r['log_nonscaling']!s} "
                  f"delta_sin2theta12={r['delta_sin2_theta12']:+.6g}")
    write_csv(dest/'matrix_comparison.csv',rows)
    write_csv(dest/'channel_comparison.csv',channels_out)
    (dest/'summary.json').parent.mkdir(parents=True,exist_ok=True)
    (dest/'summary.json').write_text(json.dumps({
        'source':str(source),'requested_models':models,'missing_models':missing,
        'models_processed':len({r['model'] for r in rows}),
        'matrix_rows':len(rows),'channel_rows':len(channels_out),
        'interpretation':'All numbers are read from saved exploratory LLSS cases. '
        'Channel categories are alternative decompositions, not additive across categories. '
        'A nonzero residual alone does not prove a bug.'},indent=2)+'\n',encoding='utf-8')
    print('Reports:',dest)
    return 2 if missing else 0


if __name__ == '__main__':
    raise SystemExit(main())
