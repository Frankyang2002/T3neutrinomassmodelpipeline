"""Read-only audit of the 16 selected simultaneous-EFT1 cases.json outputs.
Project path: scripts/ValidateT3SimultaneousSaved.py
No RG running; no modification of original output files.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'output/studies/general_scalar_extended'
A4 = BASE / 'experimental_simultaneous_eft1_seven_term/simultaneous_A_m4_v2'
NON = BASE / 'experimental_simultaneous_eft1_general_nonmixed/simultaneous_nonmixed13_v1'
MIX = BASE / 'experimental_simultaneous_eft1_mixed_BC/simultaneous_mixed_BC_v1'
MODELS = ('A_-4','A_-2','A_+0','B_-3','B_-1','B_+1','C_-3','C_-1','C_+1',
          'D_-2','D_+0','D_+2','E_-4','E_-2','E_+0','E_+2')
REQUIRED = ('1111','0111','1001','1110')
OBS = ('delta_m21_sq_ev2','delta_m3l_sq_ev2','sin2_theta12','sin2_theta13','sin2_theta23')
MATRICES = ('c5_at_scalar_threshold','c5_induced_wilson_log','c5_joint_scalar_bubble')
DIAGNOSTICS = ('combined_ode_fermion_symmetry_residual',
               'combined_ode_scalar_symmetry_residual',
               'combined_ode_weinberg_projection_residual',
               'separate_rgbeta_endpoint_max_block_rel',
               'initial_scalar_quartic_source_parity_rel')


def tag(model):
    return model.replace('+','p').replace('-','m')


def decode_matrix(v):
    if isinstance(v, dict) and 'real' in v and 'imag' in v:
        x = np.asarray(v['real'],float) + 1j*np.asarray(v['imag'],float)
    else:
        x = np.asarray(v,complex)
    if x.shape != (3,3) or not np.isfinite(x).all():
        raise ValueError('nonfinite or non-3x3 matrix')
    return x


def walk_nonfinite(v, path='$'):
    if isinstance(v,dict):
        for k,x in v.items():
            yield from walk_nonfinite(x,f'{path}.{k}')
    elif isinstance(v,list):
        for i,x in enumerate(v):
            yield from walk_nonfinite(x,f'{path}[{i}]')
    elif isinstance(v,float) and not math.isfinite(v):
        yield path


def audit(model, path):
    row={'model':model,'file':str(path),'status':'BLOCKED','saved_status':'',
         'cases':0,'required_cases':'','full_nfev':'','c5_symmetry_rel':'',
         'source_revision':'','issues':'','warnings':'','provenance':''}
    errors=[]; warnings=[]
    if not path.is_file():
        errors.append('missing file')
    else:
        try:
            with path.open(encoding='utf-8-sig') as fp:
                doc=json.load(fp,parse_constant=lambda x: (_ for _ in ()).throw(ValueError(f'nonfinite JSON constant {x}')))
            row['saved_status']=str(doc.get('status',''))
            if doc.get('model')!=model:errors.append('model ID mismatch')
            if doc.get('status')!='Computed':errors.append('status != Computed')
            cases=doc.get('cases')
            if not isinstance(cases,dict):
                errors.append('cases not a dictionary'); cases={}
            row['cases']=len(cases)
            if len(cases)!=16:errors.append(f'expected 16 stage combinations, got {len(cases)}')
            missing=[k for k in REQUIRED if k not in cases]
            row['required_cases']=','.join(k for k in REQUIRED if k in cases)
            if missing:errors.append('missing '+','.join(missing))
            for location in walk_nonfinite(doc):errors.append('nonfinite '+location)
            for code,case in cases.items():
                if not isinstance(case,dict):
                    errors.append(f'{code}: not an object');continue
                if case.get('case')!=code:errors.append(f'{code}: case label mismatch')
                if code in REQUIRED:
                    p=case.get('prediction')
                    if not isinstance(p,dict):
                        errors.append(f'{code}: prediction missing');p={}
                    for key in OBS:
                        try:
                            if not math.isfinite(float(p[key])):raise ValueError('nonfinite')
                        except (KeyError,TypeError,ValueError):errors.append(f'{code}: invalid {key}')
                for name in MATRICES:
                    if name not in case:
                        if code in REQUIRED:errors.append(f'{code}: missing {name}')
                        continue
                    try:
                        mat=decode_matrix(case[name])
                        if name=='c5_at_scalar_threshold':
                            sym=np.linalg.norm(mat-mat.T)/max(np.linalg.norm(mat),1e-300)
                            if code=='1111':row['c5_symmetry_rel']=f'{sym:.3e}'
                            if sym>1e-7:errors.append(f'{code}: C5 threshold symmetry residual {sym:.3e}')
                    except (ValueError,TypeError,OverflowError) as e:
                        errors.append(f'{code}: {name}: {e}')
                if code in REQUIRED:
                    for key in DIAGNOSTICS:
                        if key not in case:
                            warnings.append(f'{code}: no {key}');continue
                        try:
                            if not math.isfinite(float(case[key])):raise ValueError()
                        except (TypeError,ValueError):errors.append(f'{code}: invalid {key}')
                    try:
                        n=int(case['combined_ode_nfev'])
                        if n<=0:errors.append(f'{code}: combined_ode_nfev <= 0')
                        if code=='1111':row['full_nfev']=n
                    except (ValueError,TypeError,KeyError):errors.append(f'{code}: invalid combined_ode_nfev')
                    expected={'uv_on':code[0]=='1','eft1_on':code[1]=='1',
                              'wilson_log_on':code[2]=='1','smeft_on':code[3]=='1'}
                    for key,expected_value in expected.items():
                        if case.get(key)!=expected_value:errors.append(f'{code}: {key} switch mismatch')
                    if not case.get('matching_approximation','').startswith(('simultaneous_','mixed_fixed_cartesian_')):
                        errors.append(f'{code}: unexpected matching_approximation')
            if model in ('B_-1','C_-1') and '1111' in cases:
                if cases['1111'].get('fixed_cartesian_wilson_basis') is not True:
                    errors.append('mixed model missing fixed Cartesian basis confirmation')
            prov=doc.get('provenance',{})
            row['provenance']=json.dumps(prov,ensure_ascii=False,default=str)
            if not prov: warnings.append('missing provenance')
            row['source_revision']=str(doc.get('source_revision',doc.get('source_commit','not recorded at top level')))
            for k in ('source_run','canonical_run_id','scales'):
                if k not in doc: warnings.append(f'{k} not at top level; examine provenance/summary')
        except (OSError,ValueError,json.JSONDecodeError) as e:
            errors.append('JSON read failure: '+str(e))
    row['issues']='; '.join(dict.fromkeys(errors))
    row['warnings']='; '.join(dict.fromkeys(warnings))
    row['status']='PASS' if not errors else 'BLOCKED'
    return row


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,default=BASE/'reports/simultaneous_saved_validation_v1')
    args=ap.parse_args()
    sources={m:(A4 if m=='A_-4' else MIX if m in ('B_-1','C_-1') else NON) for m in MODELS}
    rows=[audit(m,sources[m]/tag(m)/'cases.json') for m in MODELS]
    args.output.mkdir(parents=True,exist_ok=True)
    csvpath=args.output/'validation_16_models.csv'
    with csvpath.open('w',newline='',encoding='utf-8') as fp:
        w=csv.DictWriter(fp,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    print('MODEL    STATUS   CASES REQUIRED       NFEV C5_SYMMETRY  ISSUES')
    for r in rows:
        print(f"{r['model']:7} {r['status']:7} {r['cases']:>5} {r['required_cases']:14} "
              f"{str(r['full_nfev']):>6} {str(r['c5_symmetry_rel']):>12}  {r['issues'][:120]}")
    print(f"PASS {sum(r['status']=='PASS' for r in rows)}/16; report: {csvpath}")
    print('Provenance and nonfatal warnings are included in the CSV.')
    return 0 if all(r['status']=='PASS' for r in rows) else 2

if __name__=='__main__':
    raise SystemExit(main())
