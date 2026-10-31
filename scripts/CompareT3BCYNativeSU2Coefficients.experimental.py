#!/usr/bin/env python3
"""B/C SU(2) component kernel: source LLSS*HHSS vs raw Y/LLHH poles.

Inputs: actual original EFT1 LLSS, full transition-EOMEFT HHSS, native CG
arrays exported from the SAME canonical Y Matchete kernel, and the saved raw
Y/LLHH pole residue. Evaluate external l_SU2=1, H_SU2=2 (Matchete 1-based).
Use rational/surd symbolic arithmetic for all group tensors and couplings.

IMPORTANT: This is a complex-scalar Wick-contraction kernel, NOT Eq. (4.85)
in a converted, fully symmetric REAL scalar tensor basis. It does not compute
an MSbar counterterm, UV/IR separation, or a physical C5. It reports direct
source-to-hard-pole coefficient RATIOS (if calculable), not an equality assumed
by normalization. Never interpret even a universal ratio as CT certification.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import itertools
import json
from pathlib import Path
import re
import sys
from typing import Any

import sympy as sp

TOK = re.compile(r'\s*([A-Za-z$][A-Za-z0-9$]*|\d+|[\[\]\{\}\(\),+\-*/^])')
class AuditError(ValueError): pass

def require(cond,why):
    if not cond: raise AuditError(why)

class Parser:
    def __init__(self,text):
        matches=list(TOK.finditer(text));self.ts=[m.group(1) for m in matches];self.i=0
        require(''.join(m.group() for m in matches).strip()==text.strip(),
                'InputForm includes an unsupported token')
    def pop(self,expect=None):
        require(self.i<len(self.ts),'Unexpected end of InputForm')
        t=self.ts[self.i];self.i+=1
        require(expect is None or t==expect,f'Expected {expect}, got {t}')
        return t
    def peek(self):return self.ts[self.i] if self.i<len(self.ts) else None
    def parse(self):
        x=self.sum();require(self.peek() is None,'Trailing InputForm');return x
    def sum(self):
        first=self.product(); xs=[first]
        while self.peek() in ('+','-'):
            op=self.pop();v=self.product();xs.append(v if op=='+' else ('neg',v))
        return first if len(xs)==1 else ('add',tuple(xs))
    def product(self):
        first=self.power();xs=[first]
        while self.peek() in ('*','/'):
            op=self.pop();v=self.power();xs.append(v if op=='*' else ('inv',v))
        return first if len(xs)==1 else ('mul',tuple(xs))
    def power(self):
        x=self.unary()
        if self.peek()=='^': self.pop();x=('pow',x,self.power())
        return x
    def unary(self):
        if self.peek()=='-':self.pop();return ('neg',self.unary())
        if self.peek()=='+':self.pop();return self.unary()
        return self.atom()
    def atom(self):
        v=self.pop()
        if v=='(':
            x=self.sum();self.pop(')');return x
        if v=='{':
            a=[]
            if self.peek()!='}':
                a.append(self.sum())
                while self.peek()==',':self.pop();a.append(self.sum())
            self.pop('}');return ('list',tuple(a))
        if v.isdecimal():return ('num',int(v))
        require(re.fullmatch(r'[A-Za-z$][A-Za-z0-9$]*',v) is not None,'Bad atom')
        if self.peek()=='[':
            self.pop();a=[]
            if self.peek()!=']':
                a.append(self.sum())
                while self.peek()==',':self.pop();a.append(self.sum())
            self.pop(']');return ('call',v,tuple(a))
        return ('sym',v)

def iscall(x,name):return x[0]=='call' and x[1]==name

def args(x,name):
    require(iscall(x,name),'Expected '+name)
    return x[2]

def val_name(x):
    require(x[0]=='sym','Expected a symbol, not '+str(x)[:100])
    return x[1]

def numerical(x):
    t=x[0]
    if t=='num':return sp.Integer(x[1])
    if t=='sym' and x[1]=='I':return sp.I
    if t=='neg':return -numerical(x[1])
    if t=='inv':return 1/numerical(x[1])
    if t=='add':return sp.Add(*(numerical(y) for y in x[1]))
    if t=='mul':return sp.Mul(*(numerical(y) for y in x[1]))
    if t=='pow':return numerical(x[1])**numerical(x[2])
    if iscall(x,'Sqrt'):return sp.sqrt(numerical(args(x,'Sqrt')[0]))
    if iscall(x,'Rational'):
        a,b=args(x,'Rational');return numerical(a)/numerical(b)
    if iscall(x,'Complex'):
        a,b=args(x,'Complex');return numerical(a)+sp.I*numerical(b)
    raise AuditError('Not an accepted exact number: '+str(x)[:100])

def read_tensor(node):
    if node[0]=='list':return [read_tensor(q) for q in node[1]]
    return numerical(node)

def product(left,right):
    return [(sp.expand(a*b),x+y) for a,x in left for b,y in right]

def expand(x):
    """Fully distribute known numeric/coupling/CG products, never evaluate arbitrary code."""
    k=x[0]
    if k=='add':return [p for y in x[1] for p in expand(y)]
    if k=='neg':return [(-c,a) for c,a in expand(x[1])]
    if k=='mul':
        acc=[(sp.Integer(1),())]
        for y in x[1]:acc=product(acc,expand(y))
        return acc
    if k=='inv':
        a=x[1]
        if a[0]=='mul':
            acc=[(sp.Integer(1),())]
            for part in a[1]:acc=product(acc,expand(('inv',part)))
            return acc
        if iscall(a,'Coupling'):
            label=val_name(args(a,'Coupling')[0]);require(label=='MF',
                'Unexpected inverse coupling '+label)
            return [(sp.Integer(1),())]  # Remove shared 1/MF at all stages
        return [(1/numerical(a),())]
    if k=='pow':
        base, exponent=x[1],x[2]
        exp=numerical(exponent)
        require(exp.is_Integer and -8<=int(exp)<=8,'Unexpected noninteger power')
        if base[0] in ('num','neg','inv') or iscall(base,'Sqrt'):
            return [(numerical(base)**exp,())]
        require(int(exp)>=0,'Inverse nonnumerical power')
        result=[(sp.Integer(1),())]
        for _ in range(int(exp)):result=product(result,expand(base))
        return result
    if k=='num':return [(sp.Integer(x[1]),())]
    if k=='sym' and x[1]=='I':return [(sp.I,())]
    if k=='sym' and x[1]=='hbar':return [(sp.Integer(1),())]
    if iscall(x,'Sqrt') or iscall(x,'Rational') or iscall(x,'Complex'):
        return [(numerical(x),())]
    if iscall(x,'Coupling') and val_name(args(x,'Coupling')[0])=='MF':
        raise AuditError('Unexpected numerator MF')
    if iscall(x,'Bar'):
        # Symbolically barred couplings and fields remain atomic.
        return [(sp.Integer(1),(x,))]
    if k=='call' and x[1] in ('Coupling','CG','Field','NCM'):
        return [(sp.Integer(1),(x,))]
    raise AuditError('Unrecognized source factor: '+str(x)[:180])

def walk_field(x,bar=False):
    if iscall(x,'Bar'):
        return walk_field(args(x,'Bar')[0],not bar)
    if iscall(x,'Field'):
        lab,kind,indices,_=args(x,'Field')
        return [(val_name(lab),val_name(kind),bar,index(indices[1][0]))]
    if x[0]=='call':
        return [v for q in x[2] for v in walk_field(q,bar)]
    if x[0] in ('list','add','mul'):
        return [v for q in x[1] for v in walk_field(q,bar)]
    return []

def unwrap_bar(x):
    if iscall(x,'Bar'):return args(x,'Bar')[0],True
    return x,False

def index(x):
    v,_=unwrap_bar(x);name,kind=args(v,'Index')
    return (val_name(name),str(kind))

def field_scalar(x):
    v,bar=unwrap_bar(x)
    if not iscall(v,'Field'):return None
    label,kind,indices,_=args(v,'Field')
    if val_name(kind)!='Scalar':return None
    return (val_name(label),bar,index(indices[1][0]))

def coupling_factor(x):
    v,bar=unwrap_bar(x)
    if not iscall(v,'Coupling'):return None
    name=val_name(args(v,'Coupling')[0]);require(name!='MF','MF should have been removed')
    return ('bar_' if bar else '')+name

def chirality(raw):
    a=('Proj[1]' in raw);b=('Proj[-1]' in raw)
    require(a!=b,'Missing or ambiguous GammaCC chirality')
    return 'PR' if a else 'PL'

def slice_component_factors(raw: str, tensors: dict, *, external: str):
    monomials=expand(Parser(raw).parse())
    output=defaultdict(lambda:sp.Integer(0))
    for coef, factors in monomials:
        if coef==0:continue
        cgs=[];fields=[];couplings=[]
        for f in factors:
            c=coupling_factor(f)
            if c is not None:couplings.append(c);continue
            x,bar=unwrap_bar(f)
            if iscall(x,'CG'):
                require(not bar,'Unexpected outer-bar CG')
                n,slots=args(x,'CG')
                v,namebar=unwrap_bar(n)
                name=val_name(v) if v[0]=='sym' else ('eps' if iscall(v,'eps') else None)
                require(name is not None,'Unsupported native CG name')
                cgs.append((name,namebar,slots[1]));continue
            if iscall(x,'Field'):
                fields.extend(walk_field(f));continue
            if iscall(f,'NCM'):
                fields.extend(walk_field(f));continue
            raise AuditError('Unexpected symbolic factor: '+str(f)[:130])
        schiral=chirality(raw) if external in ('wilson','pole') else 'NONE'
        if external=='wilson':
            require(sum(z[0]=='l' for z in fields)==2 and
                    sum(z[0] in ('NewScalarX','NewScalarY') for z in fields)==2,
                    'Incorrect LLSS operator inventory')
        elif external=='quartic':
            require(sum(z[0]=='H' for z in fields)==2 and
                    sum(z[0] in ('NewScalarX','NewScalarY') for z in fields)==2,
                    'Incorrect HHSS operator inventory')
        elif external=='pole':
            require(sum(z[0]=='l' for z in fields)==2 and
                    sum(z[0]=='H' for z in fields)==2,
                    'Incorrect LLHH pole operator inventory')
        else:raise AuditError('Bad external kind')
        index_keys=set()
        for _,_,_,ind in fields:
            if 'SU2L' in ind[1] or 'T3BSMd2' in ind[1]:index_keys.add(ind)
        for name,bar,slots in cgs:
            for q in slots:
                ind=index(q)
                if 'SU2L' in ind[1] or 'T3BSMd2' in ind[1]:index_keys.add(ind)
        index_keys=sorted(index_keys)
        for vals in itertools.product((1,2),repeat=len(index_keys)):
            ass=dict(zip(index_keys,vals))
            if any(ass[ind] != (1 if label=='l' else 2) for label,_,_,ind in fields
                   if label in ('l','H')):continue
            fac=coef
            for name,bar,slots in cgs:
                pos=tuple(ass[index(z)]-1 for z in slots)
                if name=='eps':
                    require(len(pos)==2,'Bad eps CG rank')
                    value=sp.Integer(int(pos==(0,1))-int(pos==(1,0)))
                else:
                    require(name in tensors,'CG tensor absent from native registry: '+name)
                    value=tensors[name]
                    for idx in pos:
                        require(isinstance(value,list) and 0<=idx<len(value),
                                f'CG index/shape mismatch for {name}: {pos}')
                        value=value[idx]
                    require(not isinstance(value,list),'Incomplete tensor indexing')
                    if bar:value=sp.conjugate(value)
                fac=sp.expand(fac*value)
                if fac==0:break
            if fac==0:continue
            scalars=tuple(sorted((label,bar,ass[ind]) for label,kind,bar,ind in fields
                                 if kind=='Scalar' and label!='H'))
            hbar=sum(label=='H' and kind=='Scalar' and bar for label,kind,bar,ind in fields)
            key=(tuple(sorted(couplings)),schiral,hbar,scalars)
            output[key]+=fac
    return {k:sp.simplify(v) for k,v in output.items() if sp.simplify(v)!=0}

def read_registry(path: Path, model: str):
    record=json.loads(path.read_text(encoding='utf-8'))
    require(record.get('format')=='t3bc-native-y-cg-component-export-v1','Wrong native CG format')
    require(record.get('model')==model,'Wrong model native CG file')
    require(record.get('native_Y_matching_performed_in_this_invocation') is False,
            'Export unexpectedly involved matching')
    tensors={}
    for t in record['CGRegistry']:
        name=t['name']
        require(re.fullmatch(r'T3[A-Za-z0-9]+',name) is not None,
                'Invalid registered tensor name')
        require(name not in tensors,'Duplicate CG tensor')
        tensors[name]=read_tensor(Parser(t['tensor_inputform']).parse())
    require(bool(tensors),'Empty CG registry')
    return tensors

def scalar_pair(s):
    require(len(s)==2,'Not a scalar pair')
    a,b=s
    pair='YY' if a[0]==b[0]=='NewScalarY' else ('XY' if {a[0],b[0]}==
                                                  {'NewScalarX','NewScalarY'} else 'XX')
    return pair

def contractions(wscalars,qscalars):
    total=0
    for perm in ((0,1),(1,0)):
        if all(wscalars[i][0]==qscalars[perm[i]][0] and
               wscalars[i][1] != qscalars[perm[i]][1] and
               wscalars[i][2] == qscalars[perm[i]][2] for i in range(2)):
            total+=1
    return total

def read_model(root:Path,stem:str,registry_path:Path):
    model={'B_m1':'B_-1','C_m1':'C_-1'}[stem]
    tensors=read_registry(registry_path,model)
    root2=root/'eft1_fermion_first_tree'
    paths={
      'LLSS':root2/f'{stem}_final_eft1_fermion_first_tree_llss_tree.inputform.txt',
      'EOMEFT':root2/f'{stem}_final_eft1_fermion_first_tree_eft1_transition_eom.inputform.txt',
      'pole':Path('output/audits/t3_bc_y_renormalization_ledger_v1')/
             f'{stem}_y_llhh_raw_pole_residue.inputform.txt',
    }
    from collections import Counter
    sources={}
    for k,p in paths.items():
        require(p.is_file(),f'Missing local source: {p}')
        sources[k]=p.read_text(encoding='utf-8')
    # The audited preselection is textual; parsing and arithmetic below are exact.
    from importlib.util import spec_from_file_location,module_from_spec
    sibling=Path(__file__).with_name('AuditT3BCYMixedScalarPoleSupport.experimental.py')
    require(sibling.is_file(),'Previous v2 audit source is required (install prior ZIP)')
    spec=spec_from_file_location('t3_bc_support_v2',sibling)
    module=module_from_spec(spec);spec.loader.exec_module(module)
    ll=[s for s in module.split_additive(sources['LLSS'])
        if module.field_count(s,'NewScalarX')+module.field_count(s,'NewScalarY')==2
        and s.count('Field[l, Fermion')==2]
    qq=[s for s in module.split_additive(sources['EOMEFT'])
        if module.field_count(s,'H')==2 and
        module.field_count(s,'NewScalarX')+module.field_count(s,'NewScalarY')==2
        and 'Field[l, Fermion' not in s]
    pole=module.split_additive(sources['pole'])
    require(len(qq)==42 and len(pole)==50 and len(ll) in (20,24),
            'Not the expected B/C source inventory')
    W=defaultdict(lambda:sp.Integer(0));Q=defaultdict(lambda:sp.Integer(0));P=defaultdict(lambda:sp.Integer(0))
    for s in ll:
        for key,v in slice_component_factors(s,tensors,external='wilson').items():W[key]+=v
    for s in qq:
        for key,v in slice_component_factors(s,tensors,external='quartic').items():Q[key]+=v
    for s in pole:
        for key,v in slice_component_factors(s,tensors,external='pole').items():P[key]+=v
    result=defaultdict(lambda:sp.Integer(0));counts=defaultdict(int)
    for (wc,ch,_,ws),wv in W.items():
        cl=scalar_pair(ws)
        if cl=='XX':continue
        for (qc,_,hbar,qs),qv in Q.items():
            if scalar_pair(qs)!=cl:continue
            n=contractions(ws,qs)
            if not n:continue
            key=(tuple(sorted((*wc,*qc))),ch,hbar)
            result[key]+=n*wv*qv
            counts[cl]+=n
    pole_map=defaultdict(lambda:sp.Integer(0))
    for (pc,ch,hbar,sc),v in P.items():
        require(sc==(), 'Pole unexpectedly contains scalar BSM field')
        pole_map[(pc,ch,hbar)]+=v
    records=[];allkeys=sorted(set(result)|set(pole_map),key=str)
    ratios=set();nonzero=0;unmatched=0
    for key in allkeys:
        r=sp.simplify(result.get(key,0));p=sp.simplify(pole_map.get(key,0))
        ratio=sp.simplify(p/r) if r!=0 else None
        if r!=0 and p!=0:nonzero+=1;ratios.add(str(ratio))
        if r==0 and p!=0:unmatched+=1
        records.append({'coupling_factors':list(key[0]),'chirality':key[1],
          'H_bar_count':key[2], 'tree_contraction_coefficient':str(r),
          'native_hard_pole_residue_coefficient':str(p),
          'pole_over_tree_contraction':str(ratio) if ratio is not None else None})
    return {
      'model':model,'format':'t3bc-y-complex-su2-source-kernel-comparison-v1',
      'source_sha256':{key:hashlib.sha256(value.encode('utf8')).hexdigest()
                       for key,value in sources.items()},
      'native_cg_json_sha256':hashlib.sha256(registry_path.read_bytes()).hexdigest(),
      'external_components':'l_SU2=1,H_SU2=2, native Matchete 1-based components',
      'source_terms':{'LLSS':len(ll),'HHSS':len(qq),'Y_pole':len(pole)},
      'scalar_wick_pairing_count':dict(counts),
      'nonzero_terms_with_direct_and_pole_support':nonzero,
      'native_pole_terms_without_tree_contraction':unmatched,
      'distinct_pole_over_contraction_ratios':sorted(ratios),
      'all_nonzero_ratios_equal':len(ratios)==1 and unmatched==0,
      'contraction_basis':'complex scalar, direct Wick pairing (not real Weyl EFT beta)',
      'physical_C5':None, 'GammaCC_to_Weyl_normalization_certified':False,
      'EFT1_EFT2_counterterm_certified':False,'UV_IR_pole_assignment_certified':False,
      'entries':records,
    }

def self_test():
    r=Parser('(2*Sqrt[3]+1)/2').parse()
    assert numerical(r)==sp.Rational(1,2)+sp.sqrt(3)
    assert len(expand(Parser('Coupling[y1, {}, 0]*(Coupling[y2, {}, 0]+Coupling[y3, {}, 0])').parse()))==2
    assert contractions((('NewScalarY',False,1),('NewScalarY',False,2)),
                        (('NewScalarY',True,1),('NewScalarY',True,2)))==1
    assert contractions((('NewScalarY',False,1),('NewScalarY',False,1)),
                        (('NewScalarY',True,1),('NewScalarY',True,1)))==2
    assert contractions((('NewScalarX',False,1),('NewScalarY',False,2)),
                        (('NewScalarY',True,2),('NewScalarX',True,1)))==1
    tensor={'t':[[sp.Integer(0),sp.Integer(1)],[sp.Integer(-1),sp.Integer(0)]]}
    expr='CG[t, {Index[d$$1, SU2L[fund]], Index[d$$2, T3BSMd2]}]*' \
         'Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*' \
         'Field[H, Scalar, {Index[d$$1, SU2L[fund]]}, {}]*' \
         'Field[NewScalarY, Scalar, {Index[d$$2, T3BSMd2]}, {}]*' \
         'Bar[Field[NewScalarY, Scalar, {Index[d$$2, T3BSMd2]}, {}]]'
    q=slice_component_factors(expr,tensor,external='quartic')
    assert q  # H=2, Y=1 => epsilon21 = -1
    assert all(v==-1 for v in q.values())
    # Full synthetic tree Wilson*quartic slice (not real benchmark data).
    lepton=lambda i:f'Field[l, Fermion, {{Index[{i}, SU2L[fund]], Index[f, Flavor]}}, {{}}]'
    scalar=lambda i:f'Field[NewScalarY, Scalar, {{Index[{i}, T3BSMd2]}}, {{}}]'
    higgs=lambda i:f'Field[H, Scalar, {{Index[{i}, SU2L[fund]]}}, {{}}]'
    src_w=(f'CG[v, {{Index[i, SU2L[fund]], Index[a, T3BSMd2]}}]*'
       f'CG[v, {{Index[j, SU2L[fund]], Index[b, T3BSMd2]}}]*'
       f'Coupling[y1, {{}}, 0]^2*{scalar("a")}*{scalar("b")}*'
       f'NCM[Bar[{lepton("i")}], DiracProduct[GammaCC, Proj[1]],'
       f' Transp[Bar[{lepton("j")}] ] ]/(2*Coupling[MF, {{}}, 0])')
    src_q=(f'CG[v, {{Index[h, SU2L[fund]], Index[a, T3BSMd2]}}]*'
       f'CG[v, {{Index[k, SU2L[fund]], Index[b, T3BSMd2]}}]*'
       f'Coupling[lam, {{}}, 0]*Bar[{scalar("a")}]*Bar[{scalar("b")}]*'
       f'{higgs("h")}*{higgs("k")}')
    full={'v':[[sp.Integer(1),sp.Integer(1)],
               [sp.Integer(1),sp.Integer(1)]]}
    ww=slice_component_factors(src_w,full,external='wilson')
    qq=slice_component_factors(src_q,full,external='quartic')
    assert len(ww)==len(qq)==3
    assert any(contractions(wk[3],qk[3]) for wk in ww for qk in qq)
    bad=expr.replace('CG[t,','CG[unknown,')
    try:slice_component_factors(bad,tensor,external='quartic')
    except AuditError:pass
    else:raise AssertionError('Unknown native CG was not rejected')
    print('PASS: complex CG contractions, dummy-index component projection and Wick multiplicity')
    print('PASS: unknown CG rejection, no beta/CT or physical C5 certification')

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--self-test',action='store_true')
    ap.add_argument('--input-dir',type=Path,
        default=Path('output/audits/t3_bc_y_native_cg_components_v1'))
    ap.add_argument('--root',type=Path,
        default=Path('output/general_t3_mass_eigen_matchete'))
    ap.add_argument('--output',type=Path)
    opts=ap.parse_args()
    try:
        if opts.self_test:self_test()
        if opts.output:
            report={'format':'t3bc-y-native-component-kernel-cross-check-v1',
                    'physical_C5':None,'models':{}}
            for stem in ('B_m1','C_m1'):
                cg=opts.input_dir/f'{stem}_native_CG.json'
                record=read_model(opts.root,stem,cg)
                report['models'][record['model']]=record
                print(f"{record['model']}: {record['nonzero_terms_with_direct_and_pole_support']} "
                      f"coefficient structures compared; distinct ratios="
                      f"{record['distinct_pole_over_contraction_ratios']}")
            opts.output.parent.mkdir(parents=True,exist_ok=True)
            opts.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
            print('Component-level DIAGNOSTIC (not counterterm): '+str(opts.output))
        if not(opts.self_test or opts.output):ap.error('Choose --self-test and/or --output')
        return 0
    except (AuditError,FileNotFoundError,ValueError,IndexError,TypeError,KeyError,AssertionError) as e:
        print('FAIL (no coefficient comparison or matching certificate): '+str(e),file=sys.stderr)
        return 1

if __name__=='__main__':raise SystemExit(main())
