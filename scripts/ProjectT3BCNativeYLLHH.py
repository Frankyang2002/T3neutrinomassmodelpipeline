"""Reduce verified native finite Y/LLHH InputForm into formal flavour spurions.

Requires the completed v1 finite split and v2 pole/log consistency outputs.
This is an exactly validated *neutral-component coefficient* of the Matchete
operator, NOT a physical C5 matching, EFT loop subtraction, or numerical m_nu.
"""
from __future__ import annotations
import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

from Numerical.experimental.T3BCNativeYLLHHFlavorKernel import (
    KernelError,reduce_finite_y_llhh,check_exact_chirality_conjugacy,
)
from Numerical.experimental.T3BCPoleLogConsistency import polynomial,terms

ROOT=Path(__file__).resolve().parents[1]
MODELS=('B_-1','C_-1')
FINITE=Path('output/audits/t3_bc_scalar_finite_term_extraction_v1')
POLE=Path('output/audits/t3_bc_pole_log_consistency_v2')
OUTPUT=Path('output/audits/t3_bc_native_y_llhh_flavor_kernel_v1')
FORMAT='t3-bc-native-y-llhh-neutral-component-flavor-spurions-v1'
STATUS='ExactNativeNeutralFlavorSpurionsNoPhysicalWeinbergNormalization'


def read(path):
    if not path.is_file():raise KernelError('Missing input evidence: '+str(path))
    return path.read_text(encoding='utf8').strip()


def load_json(path):
    obj=json.loads(read(path))
    if not isinstance(obj,dict):raise KernelError('Invalid JSON report: '+str(path))
    return obj


def safe_rel(x):
    p=Path(x)
    if p.is_absolute() or '..' in p.parts:
        raise KernelError('Unsafe evidence path')
    return p


def sha(data):return sha256(data.encode('utf8')).hexdigest()


def run(models,*,finite_dir,pole_dir,out_dir,root=ROOT):
    finite=load_json(finite_dir/'summary.json')
    if (finite.get('format')!='t3-bc-native-local-matchete-pole-log-split-v1' or
        finite.get('status')!='SymbolicLocalAmplitudePartitionedNoRenormalizedPhysicalMatching'):
        raise KernelError('Unexpected v1 finite audit provenance')
    pole=load_json(pole_dir/'summary.json')
    if (pole.get('format')!='t3-bc-opaque-factor-pole-log-consistency-v2' or
        pole.get('status')!='FormalPoleLogIdentityCheckedNoPhysicalCountertermOrRGE' or
        pole.get('physical_sequential_matching_certified') is not False):
        raise KernelError('Unexpected v2 pole/log audit provenance')
    frows=[r for r in finite['sector_results'] if r['threshold']=='Y' and r['sector']=='LLHH']
    prows=[r for r in pole['sector_results'] if r['threshold']=='Y' and r['sector']=='LLHH']
    if len({r['model'] for r in frows})!=len(frows) or len({r['model'] for r in prows})!=len(prows):
        raise KernelError('Duplicate model/sector evidence')
    frefs={r['model']:r for r in frows};prefs={r['model']:r for r in prows}
    prepared=[]
    for model in models:
        if model not in MODELS or model not in frefs or model not in prefs:
            raise KernelError('Missing or unsupported model')
        f=frefs[model];p=prefs[model]
        if (f['constant_terms']!=32 or f['pole_terms']!=50 or f['log_terms']!=50 or
            p['constant_terms']!=32 or
            p['source_sha256']!=f['source_sha256'] or
            p.get('pole_equals_log_coefficient_in_opaque_factor_algebra') is not True or
            p.get('physical_matching_certified') is not False):
            raise KernelError('Y LLHH provenance/coefficient identity mismatch')
        rel=safe_rel(f['at_mu_equal_heavy_mass_unrenormalized_file'])
        content=read(finite_dir/rel)
        if sha(content)!=f['scale_equal_mass_sha256']:
            raise KernelError('Y LLHH v1 expression SHA256 mismatch')
        if root is not None:
            raw=read(root/safe_rel(f['source_file']))
            if sha(raw)!=f['source_sha256']:
                raise KernelError('Original native Y LLHH source hash mismatch')
        # Verify every quartic fragment from v2 exactly recombines to v1.
        sub=pole_dir/model.replace('_-','_m')
        exprs=[]
        if set(p.get('quartic_groups',{})) != {
            'lambdaT3','lambdaT3ExtraQ020020I1','lambdaT3ExtraQ020200I1'}:
            raise KernelError('Missing or extra finite Y LLHH quartics')
        for label,entry in p['quartic_groups'].items():
            part=read(sub/f'Y_LLHH_{label}_finite_at_mass.inputform.txt')
            if len(terms(part))!=entry['raw_additive_terms']:
                raise KernelError('Quartic source fragment cardinality mismatch')
            exprs.extend(terms(part))
        table={}
        if polynomial(tuple(exprs),table=table)!=polynomial(terms(content),table=table):
            raise KernelError('Finite quartic fragments do not reproduce complete Y LLHH')
        result=reduce_finite_y_llhh(content,model)
        if not check_exact_chirality_conjugacy(result):
            raise KernelError('Missing exact chirality conjugacy')
        result.update({'format':FORMAT,'status':STATUS,
                       'input_native_source_sha256':f['source_sha256'],
                       'input_v1_at_mass_sha256':f['scale_equal_mass_sha256'],
                       'input_v2_quartic_partition_checked':True,
                       'native_PL_PR_spurion_prefactors_exact_conjugates':True,
                       'physical_finite_Y_matching_certified':False,
                       'numerical_C5_modified':False})
        prepared.append((model,result))
    # Input audit is fail-closed for all requested models before any writes.
    out_dir.mkdir(parents=True,exist_ok=True)
    rows=[]
    for model,result in prepared:
        short=model.replace('_-','_m')
        path=out_dir/f'{short}_native_y_llhh_neutral_flavor_kernel.json'
        path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
        row={'model':model,'native_terms':result['source_additive_terms'],
             'PR_monomials':result['chirality_monomial_counts']['PR'],
             'PL_monomials':result['chirality_monomial_counts']['PL'],
             'exact_spurion_chirality_conjugacy':True,
             'physical_matching_certified':False}
        rows.append(row)
        print(f'{model}: neutral finite Y/LLHH = {row["PR_monomials"]} PR + '
              f'{row["PL_monomials"]} PL monomials; exact prefactor conjugacy PASS; '
              'physical C5=NOT CERTIFIED',flush=True)
    with (out_dir/'summary.csv').open('w',encoding='utf8',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summary={'format':FORMAT,'status':STATUS,'models':list(models),'model_results':rows,
             'physical_sequential_matching_certified':False,
             'Y_finite_LLH_H_into_C5_normalization_certified':False,
             'F_first_EOMEFT_equivalence_proven':False,
             'full_counterterm_and_EFT_loop_subtraction_performed':False,
             'numerical_masses_changed':False}
    (out_dir/'summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf8')
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--all',action='store_true')
    group.add_argument('--model',choices=MODELS)
    parser.add_argument('--finite-dir',type=Path,default=FINITE)
    parser.add_argument('--pole-dir',type=Path,default=POLE)
    parser.add_argument('--output-dir',type=Path,default=OUTPUT)
    parser.add_argument('--skip-original-source-validation',action='store_true',
                        help='Use v1+v2 archived hashes only; raw source not needed')
    a=parser.parse_args()
    out=(ROOT/a.output_dir).resolve()
    if not out.is_relative_to(ROOT):raise KernelError('Output must be inside project')
    result=run(MODELS if a.all else (a.model,),
               finite_dir=(ROOT/a.finite_dir).resolve(),
               pole_dir=(ROOT/a.pole_dir).resolve(),out_dir=out,
               root=None if a.skip_original_source_validation else ROOT)
    print('Native Y/LLHH neutral flavour kernel summary:',out/'summary.csv')
    return 0

if __name__=='__main__':raise SystemExit(main())
