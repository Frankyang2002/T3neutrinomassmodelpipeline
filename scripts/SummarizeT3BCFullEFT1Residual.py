"""Validate complete B/C EFT1 residual manifests; NEVER certify physical C5.

Run after Lagrangian/RunT3BCFullEFT1Residual.experimental.wl on both models.
This checks source/candidate/residual files, coupling partition counts,
no-Match status and certification flags. It does not implement SU(2) tensor
reduction, and MUST NOT be used as an equivalence proof.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

MODELS={'B_-1':'B_m1','C_-1':'C_m1'}
FORM='t3-bc-complete-eft1-reconstruction-residual-v1'
STATUS='CompleteSymbolicResidualExportedEquivalenceUnproven'
SUCCESS='EFT1FullResidualExportedBeforeMatching'

class AuditError(ValueError):
    pass

def read_json(p):
    try:
        r=json.loads(p.read_text(encoding='utf-8'))
    except (OSError, ValueError) as e:
        raise AuditError(f'Invalid/missing JSON: {p}') from e
    if not isinstance(r,dict):raise AuditError(f'Not a JSON object: {p}')
    return r

def file_hash(d,name):
    if not isinstance(name,str) or name!=Path(name).name or not name.endswith('.inputform.txt'):
        raise AuditError(f'Invalid output file name: {name}')
    p=d/name
    try: b=p.read_bytes()
    except OSError as e: raise AuditError(f'Missing output: {p}') from e
    if not b:raise AuditError(f'Empty output: {p}')
    return hashlib.sha256(b).hexdigest()

def inspect(d,model,prefix):
    j=read_json(d/f'{prefix}_full_eft1_residual.json')
    a=read_json(d/f'{prefix}_full_eft1_residual_full_residual_manifest.json')
    if j.get('model')!=model or j.get('status')!=SUCCESS or j.get('matching_performed') is not False or j.get('physical_matching_certified') is not False:
        raise AuditError(f'{model}: invalid no-Match status')
    if (a.get('format')!=FORM or a.get('model')!=model or a.get('status')!=STATUS or
        a.get('matching_scope')!='FullEFT' or
        a.get('coupling_group_partition_conserved') is not True or
        a.get('contains_majorana_residual') is not False or
        a.get('CG_and_index_tensor_equivalence_proven') is not False or
        a.get('physical_matching_certified') is not False or
        a.get('renormalized_finite_C5_computed') is not False):
        raise AuditError(f'{model}: invalid manifest or physical flags')
    files=a.get('files',{})
    if set(files)!={'saved_eomeft','reconstructed_candidate','symbolic_residual'}:
        raise AuditError(f'{model}: incomplete full-expression export')
    checks={name:file_hash(d,basename) for name,basename in files.items()}
    groups=a.get('coupling_groups')
    if not isinstance(groups,list): raise AuditError(f'{model}: missing coupling groups')
    n=0
    for g in groups:
        if not isinstance(g,dict) or not isinstance(g.get('structural_sector'),str):
            raise AuditError(f'{model}: malformed group')
        if type(g.get('additive_terms')) is not int or g['additive_terms']<=0:
            raise AuditError(f'{model}: malformed group count')
        n+=g['additive_terms']
        file_hash(d,g.get('file'))
    if n!=a.get('residual_additive_terms'):
        raise AuditError(f'{model}: groups do not cover full residual')
    return {'model':model,'residual_additive_terms':n,'structural_sector_count':len(groups),
            'full_expression_sha256':checks,
            'original_equality_proven': j.get('original_eomeft_equal_to_independent_candidate') is True,
            'reconstruction_matchete_check_passed':j.get('reconstruction_check_passed') is True,
            'CG_tensor_equivalence_proven':False,'physical_matching_certified':False}

def self_test():
    import tempfile
    with tempfile.TemporaryDirectory() as directory:
        d=Path(directory)
        for model,p in MODELS.items():
            for code in ['saved','candidate','residual','part']:
                (d/f'{p}_{code}.inputform.txt').write_text(code,encoding='utf-8')
            main={'model':model,'status':SUCCESS,'matching_performed':False,
                  'physical_matching_certified':False,
                  'original_eomeft_equal_to_independent_candidate':False,
                  'reconstruction_check_passed':True}
            man={'format':FORM,'model':model,'status':STATUS,'matching_scope':'FullEFT',
                 'coupling_group_partition_conserved':True,'contains_majorana_residual':False,
                 'CG_and_index_tensor_equivalence_proven':False,'physical_matching_certified':False,
                 'renormalized_finite_C5_computed':False,
                 'residual_additive_terms':2,
                 'files':{k:f'{p}_{v}.inputform.txt' for k,v in {
                     'saved_eomeft':'saved','reconstructed_candidate':'candidate',
                     'symbolic_residual':'residual'}.items()},
                 'coupling_groups':[{'structural_sector':'lambdaT3__H2_NewScalarX1_NewScalarY1',
                                     'additive_terms':2,'file':f'{p}_part.inputform.txt'}]}
            jf=d/f'{p}_full_eft1_residual.json'
            mf=d/f'{p}_full_eft1_residual_full_residual_manifest.json'
            jf.write_text(json.dumps(main));mf.write_text(json.dumps(man))
            assert inspect(d,model,p)['residual_additive_terms']==2
            man['physical_matching_certified']=True
            mf.write_text(json.dumps(man))
            try:inspect(d,model,p)
            except AuditError:pass
            else:raise AssertionError('False physical certificate accepted')
    print('PASS: full residual inventory and fail-closed physical flags')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory',type=Path,default=Path('output/audits/t3_bc_full_eft1_equivalence_v1'))
    parser.add_argument('--self-test',action='store_true')
    args=parser.parse_args()
    if args.self_test:
        self_test();return 0
    try:rows=[inspect(args.directory,m,p) for m,p in MODELS.items()]
    except AuditError as e:
        print(f'FAIL (not certified): {e}');return 2
    report={'format':'t3-bc-complete-eft1-residual-audit-v1',
            'status':'FullResidualsVerifiedTensorIdentityOutstanding',
            'models':rows,'physical_matching_certified':False}
    path=args.directory/'summary.json'
    path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    for row in rows:
        print(f"{row['model']}: complete residual terms={row['residual_additive_terms']}; "
              f"structural sectors={row['structural_sector_count']}; physical C5=NOT CERTIFIED")
    print(f'PASS: complete residual files validated; report={path}')
    return 0

if __name__=='__main__':raise SystemExit(main())
