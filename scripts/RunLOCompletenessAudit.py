"""All-16 T3 leading-one-loop ledger and optional conditional hard kernel.

python -m scripts.RunLOCompletenessAudit --all --output-dir output/general_observable/lo_certification_completeness/audit_01
python -m scripts.RunLOCompletenessAudit --all --hard-input-dir input/lo_hard_conditional --output-dir output/general_observable/lo_certification_completeness/audit_02

Only creates files under output/general_observable/lo_certification_completeness.
Never edits another workstream, does not call Wolfram or write to GitHub.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
from Numerical.lo_certification.completeness.inventory import build_ledger,ALL_MODELS
from Numerical.lo_certification.completeness.hard_matching import evaluate_hard_scalar_channel

ROOT=Path(__file__).resolve().parents[1]
OUTPUT_PREFIX=Path('output/general_observable/lo_certification_completeness')


def tag(model: str) -> str: return model.replace('-','m').replace('+','p')


def encode(a: np.ndarray) -> dict:
    z=np.asarray(a,dtype=complex)
    return {'real':z.real.tolist(),'imag':z.imag.tolist()}


def compute_conditional(model: str, directory: Path) -> dict:
    p=directory/(tag(model)+'.json')
    if not p.is_file(): return {'model':model,'status':'MISSING_INPUT','physical_C5_certified':False,
                                'missing_input':str(p)}
    d=json.loads(p.read_text(encoding='utf-8'))
    if d.get('model')!=model or d.get('tensor_convention')!='T3_REAL_LLSS_HHSS_WEINBERG_v1' or d.get('scheme')!='MSbar':
        raise ValueError(f'{model}: explicit model, tensor convention and MSbar are mandatory')
    npz=p.parent/d.get('tensor_npz','')
    if npz.parent!=p.parent or not npz.is_file():
        raise ValueError(f'{model}: NPZ must be adjacent and named explicitly')
    with np.load(npz,allow_pickle=False) as z:
        required=('C_LLSS_by_fermion','fermion_masses_gev','lambda_HHSS','scalar_mass_squared_gev2','C_LLSS_total_muF')
        if not set(required).issubset(z.files):
            raise ValueError(f'{model}: missing required arrays {sorted(set(required)-set(z.files))}')
        expected_n=12 if model.startswith('E_') else 8
        if z['scalar_mass_squared_gev2'].shape!=(expected_n,expected_n):
            raise ValueError(f'{model}: expected {expected_n} real scalar coordinates covering BOTH multiplets')
        r=evaluate_hard_scalar_channel(z['C_LLSS_by_fermion'],z['fermion_masses_gev'],
              z['lambda_HHSS'],z['scalar_mass_squared_gev2'],float(d['mu_f_gev']),
              float(d['mu_s_gev']),C_total_muF=z['C_LLSS_total_muF'],
              lambda_channels=z['lambda_HHSS_channels'] if 'lambda_HHSS_channels' in z.files else None,
              channel_names=tuple(d.get('quartic_channel_names',[])))
    out={k:v for k,v in r.items() if k not in ('matrices','per_pair_terms','quartic_channels')}
    out.update({'model':model,'status':'CONDITIONAL_NUMERICAL_NOT_PHYSICS_CERTIFIED',
                'input_provenance':d.get('provenance'),
                'source_normalization_certified':bool(d.get('source_normalization_certified',False)),
                'matrices':{k:encode(v) for k,v in r['matrices'].items()},
                'individual_kernel_terms':r['per_pair_terms'],
                'quartic_channels':{name:{'matrices':{k:encode(v) for k,v in item['matrices'].items()},
                                        'projection_residuals':item['projection_residuals'],
                                        'physical_C5_certified':False}
                                    for name,item in r['quartic_channels'].items()}})
    # Even supplied assertions cannot establish independently verified matching.
    out['physical_C5_certified']=False
    return out


def ensure_owned_output(path:Path) -> None:
    path=path.resolve()
    owned=(ROOT/OUTPUT_PREFIX).resolve()
    if path==owned or owned not in path.parents:
        raise ValueError(f'Output restricted to descendant of {owned}')
    if path.exists(): raise FileExistsError('Output already exists; use a new run name: '+str(path))
    path.mkdir(parents=True,exist_ok=False)


def main(argv: list[str] | None=None) -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--all',action='store_true',help='Generate all 16 rows per contribution')
    p.add_argument('--model',choices=ALL_MODELS,help='Optional single model')
    p.add_argument('--hard-input-dir',type=Path,help='Optional independently prepared mass-eigenstate tensor inputs')
    p.add_argument('--output-dir',type=Path,required=True,help='Fresh owned output directory')
    p.add_argument('--check-only',action='store_true')
    a=p.parse_args(argv)
    if not a.all and a.model is None: p.error('Specify --all or --model')
    if a.all and a.model is not None: p.error('--all and --model cannot be combined')
    target=a.output_dir if a.output_dir.is_absolute() else ROOT/a.output_dir
    allowed=(ROOT/OUTPUT_PREFIX).resolve()
    if target.resolve()==allowed or allowed not in target.resolve().parents:
        p.error('Output outside exclusive ownership or not under a NEW child directory')
    ledger=build_ledger(ROOT)
    selected=list(ALL_MODELS) if a.all else [a.model]
    ledger['rows']=[row for row in ledger['rows'] if row['model'] in selected]
    ledger['models']=selected
    ledger['n_models']=len(selected)
    ledger['n_rows']=len(ledger['rows'])
    from collections import Counter
    ledger['status_counts']=dict(Counter(r['status'] for r in ledger['rows']))
    hard={}
    if a.hard_input_dir:
        source=a.hard_input_dir if a.hard_input_dir.is_absolute() else ROOT/a.hard_input_dir
        for model in selected:
            hard[model]=compute_conditional(model,source)
        ledger['conditional_hard_numerical_available']=sum(r['status'].startswith('CONDITIONAL') for r in hard.values())
    if a.check_only:
        print(json.dumps({'n_models':len(selected),'n_rows':len(ledger['rows']),
                          'hard_counts':{m:v['status'] for m,v in hard.items()}},indent=2))
        return 0
    ensure_owned_output(target)
    (target/'ledger.json').write_text(json.dumps(ledger,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    with (target/'ledger.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=['model','contribution','status','implementation_state',
                                    'physical_certified','source_path','source_supported_reason','external_blockers'])
        w.writeheader()
        for r in ledger['rows']:
            w.writerow({key:json.dumps(r[key]) if key=='external_blockers' else r[key] for key in w.fieldnames})
    if a.hard_input_dir:
        (target/'conditional_hard').mkdir()
        for model,val in hard.items():
            (target/'conditional_hard'/(tag(model)+'.json')).write_text(json.dumps(val,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    summary={'model_count':len(selected),'row_count':len(ledger['rows']),
             'status_counts':ledger['status_counts'],'physical_full_one_loop_certified':False,
             'conditional_hard_numerical_statuses':{m:v['status'] for m,v in hard.items()},
             'eft_sequence':ledger['eft_sequence']}
    (target/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'output':str(target),'summary':summary},indent=2))
    return 0


if __name__=='__main__':
    raise SystemExit(main())
