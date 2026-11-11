"""Read-only symbolic Laurent audit of native A_-4 and D_-2 scalar-soft sources.

This is an operator-stripped source audit, NOT renormalised one-loop matching.
No assumptions about the UV hard amplitude, CTs, or GammaCC/Weyl dictionary.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
import zipfile
from pathlib import Path
import sympy as sp

NAMES = ('L','LE','W','WE','A','B','C','Q','eps')
SYMS = {name: sp.Symbol(name) for name in NAMES}
SUBS = {
    'Bar[Coupling[lambdaT3ExtraQ020101I1, {}, 0]]':'LE',
    'Bar[Coupling[T3EFT1WExtra, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]':'WE',
    'Bar[Coupling[T3EFT1W12, {Index[d$$1, Flavor], Index[d$$2, Flavor]}, 0]]':'W',
    'Coupling[lambdaT3, {}, 0]':'L',
    'Coupling[MS1, {}, 0]':'A',
    'Coupling[MS2, {}, 0]':'B',
    'Coupling[T32MinusMass, {}, 0]':'B',
    'Coupling[T32PlusMass, {}, 0]':'C',
    'μbar2':'Q', 'ϵ':'eps', 'hbar':'1', 'Sqrt[3]':'sqrt(3)',
}

def parse_exact_source(source):
    s = source.strip().lstrip('\ufeff')
    if not s or 'T3EFT1W12' not in s:
        raise ValueError('Expected original Matchete scalar-soft W12 expression')
    for before, after in SUBS.items():
        s = s.replace(before,after)
    s = s.replace('Log[','log(').replace('[','(').replace(']',')').replace('^','**')
    if not re.fullmatch(r'[A-Za-z0-9_+*/().\s-]+',s):
        raise ValueError('Unsupported Matchete source tokens remain: '+s[:160])
    return sp.sympify(s,locals={**SYMS,'sqrt':sp.sqrt,'log':sp.log},evaluate=True)

def audit_one(source,model):
    expr=parse_exact_source(source)
    eps=SYMS['eps']
    residue=sp.simplify(sp.limit(eps*expr,eps,0))
    finite=sp.simplify(sp.limit(expr-residue/eps,eps,0))
    L,LE,W,WE,A,B,C,Q=[SYMS[x] for x in NAMES[:-1]]
    expected=(-2/sp.sqrt(3))*(L*W if model=='A_-4' else L*W+LE*WE)
    if sp.simplify(residue-expected)!=0:
        raise ValueError(model+': scalar-soft residue mismatches native expectation')
    result={'model':model,'source_soft_pole_residue_formal_hbar1':sp.sstr(residue),
            'source_soft_finite_formal_hbar1':sp.sstr(finite),
            'soft_pole_verified':True,
            'matched_UV_minus_EFT1':False,
            'counterterms_verified':False,
            'gammaCC_to_Weyl_verified':False,
            'physical_Wilson_available':False}
    if model=='A_-4':
        expected_fin=expected*(1+(A*A*sp.log(Q/(A*A))-B*B*sp.log(Q/(B*B)))/(A*A-B*B))
        if sp.simplify(finite-expected_fin)!=0:raise ValueError('A4 finite formula mismatch')
        result['finite_compact']=sp.sstr(expected_fin)
    else:
        coeff={f'{x}_{y}':sp.factor(sp.diff(sp.diff(finite,SYMS[x]),SYMS[y]))
               for x,y in [('W','L'),('W','LE'),('WE','L'),('WE','LE')]}
        if sp.simplify(coeff['W_L']-coeff['WE_LE'])!=0:raise ValueError('D diagonal mismatch')
        if sp.simplify(coeff['W_LE']-coeff['WE_L'])!=0:raise ValueError('D mixing mismatch')
        reconstructed=(coeff['W_L']*(W*L+WE*LE)+coeff['W_LE']*(W*LE+WE*L))
        if sp.simplify(finite-reconstructed)!=0:raise ValueError('D finite reconstruction mismatch')
        result['D_soft_finite_diagonal_coefficient']=sp.sstr(coeff['W_L'])
        result['D_soft_finite_mixed_coefficient']=sp.sstr(coeff['W_LE'])
        result['D_soft_finite_decomposition']='diag*(W*L+WE*LE)+mixed*(W*LE+WE*L)'
        result['extra_quartic_remains_if_WE_zero']=sp.simplify(coeff['W_LE'])!=0
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-zip',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    with zipfile.ZipFile(args.source_zip) as z:
        names=[n for n in z.namelist() if n.endswith('/A_scalar_soft_only_unrenormalized.txt') or n.endswith('/D_scalar_soft_only_unrenormalized.txt')]
        a=[n for n in names if '/A_-4/' in n]
        d=[n for n in names if '/D_-2/' in n]
        if len(a)!=1 or len(d)!=2:raise ValueError(f'Expected one A4 and two D2 soft sources; found {len(a)}, {len(d)}')
        if z.read(d[0])!=z.read(d[1]):raise ValueError('D2 soft exports differ')
        results=[]
        for model,name in [('A_-4',a[0]),('D_-2',d[0])]:
            payload=z.read(name)
            audit=audit_one(payload.decode('utf-8-sig'),model)
            audit['relative_path']=name
            audit['source_sha256']=hashlib.sha256(payload).hexdigest()
            audit['source_bytes']=len(payload)
            results.append(audit)
    doc={'status':'NativeScalarSoftLaurentSymbolicallyChecked_NotFullMatching',
         'source_archive_sha256':hashlib.sha256(args.source_zip.read_bytes()).hexdigest(),
         'D_duplicate_source_identical':True,'results':results}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(doc,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    for x in results:
        print(x['model'], 'soft residue:', x['source_soft_pole_residue_formal_hbar1'])
        if x['model']=='D_-2': print('D mixed finite nonzero:',x['extra_quartic_remains_if_WE_zero'])
    print('Wrote',args.output.resolve())
if __name__=='__main__':main()
