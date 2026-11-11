"""Diagnose flavour changes from saved corrected T3 U/E/W/S cases.

Intended project path: scripts/DiagnoseT3EFT1Mixing.py

Reads saved numerical results only; never reruns physics, modifies reference
outputs, or assumes that EFT1 must preserve mixing angles.

Run from repository root:
 python -X utf8 -m scripts.DiagnoseT3EFT1Mixing --run-id corrected_llss16_v1

For each model, compare EFT1 on/off at fixed UV, W and SMEFT switches.
The core diagnostic is the optimal complex common scale and its non-scaling
Frobenius residual for C5 at the scalar threshold. Mixing-angle shifts are
reported in sin^2(theta), the pipeline's native observable convention.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = ROOT / 'output/studies/general_scalar_extended/corrected_llss_eft1_freeze'
DEFAULT_REPORT = ROOT / 'Reports/output/thesis_qualitative/eft1_flavour_diagnostics'
OBS = ('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23')
MATRIX_KEYS = ('c5_at_scalar_threshold','c5_induced_wilson_log','c5_joint_scalar_bubble')


def decode(v):
    if isinstance(v, dict) and 'real' in v and 'imag' in v:
        a = np.asarray(v['real'], dtype=float) + 1j * np.asarray(v['imag'], dtype=float)
    else:
        a = np.asarray(v, dtype=complex)
    if a.shape != (3, 3) or not np.isfinite(a).all():
        raise ValueError('Nonfinite or non-3x3 complex matrix')
    return a


def proportionality(a, b):
    """Compute optimal a ~= r*b. Return r, relative non-scaling residual."""
    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))
    if na == 0 or nb == 0:
        return (complex(float('nan'), float('nan')), float('nan'))
    r = np.vdot(b, a) / np.vdot(b, b)
    eps = float(np.linalg.norm(a-r*b) / na)
    return r, eps


def symmetric_residual(a):
    return float(np.linalg.norm(a-a.T)/max(np.linalg.norm(a), 1e-300))


def get_cases(path):
    raw=json.loads(path.read_text(encoding='utf-8-sig'))
    if raw.get('status') != 'Computed' or not isinstance(raw.get('cases'), dict):
        raise ValueError(f'No Computed cases in {path}')
    cases=raw['cases']
    required = {f'{u}{e}{w}{s}' for u in '01' for e in '01' for w in '01' for s in '01'}
    if set(cases) != required:
        raise ValueError(f'{path}: expected exactly 16 distinct cases, found {sorted(cases)}')
    return cases


def generate(source, output):
    output.mkdir(parents=True, exist_ok=True)
    results=[]
    failures=[]
    case_files=sorted(source.glob('*/cases.json'))
    for file in case_files:
        try:
            cases=get_cases(file)
            model=file.parent.name
            for uv in '01':
                for w in '01':
                    for sm in '01':
                        on = cases[f'{uv}1{w}{sm}']
                        off = cases[f'{uv}0{w}{sm}']
                        record={'model': model, 'uv':uv, 'wilson':w,'smeft':sm,
                                'on_case':f'{uv}1{w}{sm}', 'off_case':f'{uv}0{w}{sm}'}
                        for name in MATRIX_KEYS:
                            a=decode(on[name]); b=decode(off[name])
                            r,eps=proportionality(a,b)
                            record[f'{name}_scale_real']=float(r.real)
                            record[f'{name}_scale_imag']=float(r.imag)
                            record[f'{name}_nonscaling_residual']=eps
                            record[f'{name}_on_norm']=float(np.linalg.norm(a))
                            record[f'{name}_off_norm']=float(np.linalg.norm(b))
                            record[f'{name}_symmetric_residual_on']=symmetric_residual(a)
                        p1=on['prediction'];p0=off['prediction']
                        for key in OBS:
                            record[f'{key}_on']=float(p1[key])
                            record[f'{key}_off']=float(p0[key])
                            record[f'{key}_difference']=float(p1[key])-float(p0[key])
                        results.append(record)
        except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as exc:
            failures.append({'file':str(file),'error':f'{type(exc).__name__}: {exc}'})
    dest=output/'eft1_flavour_nonscaling.csv'
    if results:
        with dest.open('w', newline='', encoding='utf-8') as f:
            writer=csv.DictWriter(f, fieldnames=list(results[0]))
            writer.writeheader();writer.writerows(results)
    summary={'source':str(source),'models_found':len(case_files),
             'models_processed':len(set(r['model'] for r in results)),
             'comparisons':len(results),'failures':failures,
             'interpretation':('C5 non-scaling residual is measured at the scalar threshold; '
              'the angles are downstream observables and need not be invariant. '
              'No physical matching completeness or expected zero angles is asserted.')}
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n', encoding='utf-8')
    print(f"Processed {summary['models_processed']} models, {len(results)} comparisons, {len(failures)} errors")
    print('Source:',source)
    if results:
        print('EFT1 comparisons at UV=1 W=1 SMEFT=1:')
        for r in results:
            if (r['uv'],r['wilson'],r['smeft']) == ('1','1','1'):
                print(f"  {r['model']:8s} C5 nonscaling={r['c5_at_scalar_threshold_nonscaling_residual']:.4g} "
                      f"Δsin²θ12={r['sin2_theta12_difference']:+.5g} "
                      f"Δsin²θ13={r['sin2_theta13_difference']:+.5g} "
                      f"Δsin²θ23={r['sin2_theta23_difference']:+.5g}")
    print('Report:',dest if results else output/'summary.json')
    return 1 if failures or not results else 0


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-id',default='corrected_llss16_v1')
    p.add_argument('--source',type=Path, default=None,
                   help='Full path to existing corrected-run directory')
    p.add_argument('--output',type=Path,default=None,
                   help='Separate directory for diagnostic CSV/JSON')
    args=p.parse_args()
    if not args.run_id or any(x not in 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-' for x in args.run_id):
        p.error('Invalid run-id')
    source=args.source if args.source else DEFAULT_SOURCE/args.run_id
    output=args.output if args.output else DEFAULT_REPORT/args.run_id
    if not source.is_dir():
        p.error(f'Saved run folder not found: {source}')
    if source.resolve()==output.resolve() or source.resolve() in output.resolve().parents:
        p.error('Output must be separate from existing physics results')
    return generate(source,output)


if __name__=='__main__':
    raise SystemExit(main())
