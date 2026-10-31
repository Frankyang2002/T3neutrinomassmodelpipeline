#!/usr/bin/env python3
"""Independent complex-scalar YY/XY bubble finite-structure diagnostic for B/C.

Read the verified v3 native-CG matrices and original EFT1 tree operators.
Recontract the two vertices by scalar channel, compare the resulting UV
kernel with the SAVED native Y/LLHH residue, and compare the XY kernel with
the SAVED formal finite term at mu=M_Y. The latter is tested, NOT assumed.

Compute the dimensionless finite massive two-scalar bubble F(mX^2/mY^2)
for interpretation only; it is not a counterterm, a replacement for Matchete
matching, or a physical Wilson coefficient. Does not execute Wolfram.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

import sympy as sp

BASE = Path('output/general_t3_mass_eigen_matchete')
AUDIT = Path('output/audits')
STEMS = {'B_m1':'B_-1', 'C_m1':'C_-1'}


def require(ok: bool, description: str) -> None:
    if not ok:
        raise ValueError(description)


def load_sibling(filename: str, name: str):
    path = Path(__file__).resolve().with_name(filename)
    require(path.is_file(), 'Required previous-version file missing: '+str(path))
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, 'Cannot load '+str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def finite_bubble_muy(r: float) -> float:
    """MSbar scalar two-propagator bubble finite scalar at mu=mY.

       B0^fin(mY^2,mX^2;mu=mY) = 1 + r log(r)/(1-r), r=mX^2/mY^2.
       This is the full zero-external-momentum bubble, NOT automatically
       the Wilson matching coefficient in the sequential EFT prescription.
    """
    require(math.isfinite(r) and 0 < r <= 1, 'Require 0<MX^2/MY^2<=1')
    d = 1.0-r
    if d < 1.e-5:
        # Exact expansion f= sum_{n>=1} d^n / [n(n+1)]
        return sum(d**n/(n*(n+1)) for n in range(1,7))
    return 1.0+(r*math.log(r))/(1.0-r)


def exact_zero(x) -> bool:
    return sp.simplify(x) == 0


def add_kernel(dictionary, key, value):
    dictionary[key] += value


def get_pair_kernels(v3, v2, root: Path, stem: str, cg: Path):
    t = v3.read_registry(cg, STEMS[stem])
    eft = (root/'eft1_fermion_first_tree'/
           f'{stem}_final_eft1_fermion_first_tree_eft1_transition_eom.inputform.txt')
    llss = (root/'eft1_fermion_first_tree'/
            f'{stem}_final_eft1_fermion_first_tree_llss_tree.inputform.txt')
    require(eft.is_file() and llss.is_file(), f'{stem}: Missing original EFT1 source')
    lltext,efttext = llss.read_text(encoding='utf-8'),eft.read_text(encoding='utf-8')
    ll = [s for s in v2.split_additive(lltext)
          if (v2.field_count(s,'NewScalarX')+v2.field_count(s,'NewScalarY')==2
              and s.count('Field[l, Fermion')==2)]
    qq = [s for s in v2.split_additive(efttext)
          if (v2.field_count(s,'H')==2 and
              v2.field_count(s,'NewScalarX')+v2.field_count(s,'NewScalarY')==2
              and 'Field[l, Fermion' not in s)]
    require(len(ll)==(20 if stem=='B_m1' else 24) and len(qq)==42,
            stem+': Original source term inventory changed')
    W=defaultdict(lambda:sp.Integer(0))
    Q=defaultdict(lambda:sp.Integer(0))
    for src in ll:
        for k,v in v3.slice_component_factors(src,t,external='wilson').items():
            W[k]+=v
    for src in qq:
        for k,v in v3.slice_component_factors(src,t,external='quartic').items():
            Q[k]+=v
    kernels={channel:defaultdict(lambda:sp.Integer(0)) for channel in ('YY','XY')}
    # Two possible Wick pairings. No arbitrary 2! inserted for the mixed pair.
    for (wc,ch,_,ws),wv in W.items():
        pair=v3.scalar_pair(ws)
        if pair=='XX':continue
        for (qc,_,hbar,qs),qv in Q.items():
            if v3.scalar_pair(qs)!=pair:continue
            n=v3.contractions(ws,qs)
            if n:
                k=(tuple(sorted((*wc,*qc))),ch,hbar)
                kernels[pair][k]+=n*wv*qv
    return ({ch:{k:sp.simplify(v) for k,v in vals.items()
                  if not exact_zero(v)} for ch,vals in kernels.items()},
            {'LLSS':len(ll),'HHSS':len(qq)},
            {'LLSS':hashlib.sha256(lltext.encode('utf8')).hexdigest(),
             'EOMEFT':hashlib.sha256(efttext.encode('utf8')).hexdigest(),
             'CG':hashlib.sha256(cg.read_bytes()).hexdigest()})


def native_sector(v3,v2,source:Path,tensors:dict,expected_count:int):
    require(source.is_file(), 'Missing native Y-threshold ledger: '+str(source))
    raw=source.read_text(encoding='utf-8')
    terms=v2.split_additive(raw)
    require(len(terms)==expected_count,
            str(source)+f': expected {expected_count} native terms, found {len(terms)}')
    out=defaultdict(lambda:sp.Integer(0))
    for part in terms:
        for (cs,ch,hbar,scalars),coefficient in v3.slice_component_factors(
                part,tensors,external='pole').items():
            require(scalars==(), 'Unexpected BSM field in native LLHH')
            out[(cs,ch,hbar)]+=coefficient
    return {k:sp.simplify(v) for k,v in out.items() if not exact_zero(v)},hashlib.sha256(raw.encode('utf8')).hexdigest()


def compare(a,b):
    """Compare independently assembled coefficient dictionaries exactly."""
    nonzero=[]
    for key in set(a)|set(b):
        av=sp.simplify(a.get(key,0));bv=sp.simplify(b.get(key,0))
        if not exact_zero(av-bv):
            nonzero.append({'couplings':list(key[0]),'chirality':key[1],
                'H_bar_count':key[2],'lhs':str(av),'rhs':str(bv),
                'difference':str(sp.simplify(av-bv))})
    return sorted(nonzero,key=lambda x:str(x))


def serialize_map(m):
    return [{'couplings':list(k[0]),'chirality':k[1],'H_bar_count':k[2],
             'coefficient':str(v)} for k,v in sorted(m.items(),key=str)]


def run_one(v3,v2,root:Path,stem:str,native_dir:Path,ledger_dir:Path,previous:dict):
    model=STEMS[stem]
    cg=native_dir/f'{stem}_native_CG.json'
    kernels,counts,source_hashes=get_pair_kernels(v3,v2,root,stem,cg)
    t=v3.read_registry(cg,model)
    pole_path=ledger_dir/f'{stem}_y_llhh_raw_pole_residue.inputform.txt'
    finite_path=ledger_dir/f'{stem}_y_llhh_formal_finite_at_muY.NOT_PHYSICAL_C5.inputform.txt'
    P, polehash=native_sector(v3,v2,pole_path,t,50)
    F, finitehash=native_sector(v3,v2,finite_path,t,32)
    result={k:sp.simplify(kernels['YY'].get(k,0)+kernels['XY'].get(k,0))
            for k in set(kernels['YY'])|set(kernels['XY'])}
    result={k:v for k,v in result.items() if not exact_zero(v)}
    uv_diff=compare(result,P)
    finite_diff=compare(kernels['XY'],F)
    # Check the previous v3 report rather than silently trusting its stdout.
    prev=previous.get('models',{}).get(model,{})
    require(prev.get('format')=='t3bc-y-complex-su2-source-kernel-comparison-v1',
            model+': Missing or wrong previous v3 comparison format')
    for provenance_key in ('LLSS','EOMEFT'):
        require(prev.get('source_sha256',{}).get(provenance_key)==
                source_hashes[provenance_key],
                model+': Original EFT1 source differs from previous v3 audit: '+provenance_key)
    require(prev.get('source_sha256',{}).get('pole')==polehash,
            model+': Native pole source differs from previous v3 audit')
    require(prev.get('native_cg_json_sha256')==source_hashes['CG'],
            model+': Native CG export differs from previous v3 audit')
    previous_map={}
    for e in prev.get('entries',[]):
        k=(tuple(e['coupling_factors']),e['chirality'],e['H_bar_count'])
        require(k not in previous_map,model+': Duplicate v3 coefficient key')
        # str -> sympify limited to known arithmetic symbols from self-produced v3 output.
        # It is not used as a source of independent coefficient values or equalities.
        previous_map[k]=e['tree_contraction_coefficient']
    require(set(previous_map)==set(result)|set(P),
            model+': Previous v3 key inventory differs')
    for k,v in result.items():
        require(previous_map[k]==str(v),
                model+': Previously reported coefficient differs: '+str(k)[:120])
    params=root/f'{stem}_final_parameters.json'
    data=json.loads(params.read_text(encoding='utf-8'))
    require(data.get('model')==model and
            data.get('threshold_plan',{}).get('decreasing_scale_order')==
                  ['NewScalarY','NewScalarX'],
            model+': Incorrect mass-basis threshold metadata')
    m2=data['heavy_mass_squared_gev2']
    require(len(m2)==2 and m2[0]>0 and m2[1]>m2[0],
            model+': Invalid heavy mass hierarchy')
    r=float(m2[0])/float(m2[1]);f=finite_bubble_muy(r)
    return {
        'format':'t3bc-y-scalar-bubble-finitemass-diagnostic-v1',
        'model':model,
        'source_sha256':source_hashes|{'pole':polehash,'formal_finite':finitehash},
        'source_terms':counts|{'Y_pole':50,'Y_formal_finite':32},
        'nonzero_keys':{'YY':len(kernels['YY']),'XY':len(kernels['XY']),
                        'total_pole':len(P),'total_formal_finite':len(F)},
        'tree_UV_kernel_equals_native_pole':not uv_diff,
        'tree_XY_kernel_equals_formal_finite_at_muY':not finite_diff,
        'UV_mismatches':uv_diff,
        'finite_mismatches':finite_diff,
        'mass_scales_gev':{'X':float(data['heavy_masses_gev'][0]),
                           'Y':float(data['heavy_masses_gev'][1])},
        'squared_mass_ratio_r':r,
        'scalar_bubble_msbar_finite_at_muY':{
            'YY':0.,'XY_full_massive':f,'XY_light_X_limit':1.,
            'difference_full_minus_light_limit':f-1.,
            'formula':'F_XY(r)=1 + r*ln(r)/(1-r) at mu=mY; F_YY=0'},
        'complex_scalar_tree_kernel_Y_pieces':{
            'YY':serialize_map(kernels['YY']), 'XY':serialize_map(kernels['XY'])},
        'EFT1_EFT2_counterterms_computed':False,
        'Matchete_finite_mass_resummation_applied':False,
        'full_EFT_matching_mass_correction_computed':False,
        'GammaCC_Weyl_normalization_certified':False,
        'physical_C5_Y':None,
        'physical_matching_certified':False,
    }


def self_test():
    for r in (1.,1.-1e-9,0.5,0.1,1e-8):
        f=finite_bubble_muy(r)
        assert 0<=f<=1.00000001
        if r<1:
            fexact=float(1+sp.Float(r)*sp.log(sp.Float(r))/(1-sp.Float(r)))
            assert abs(f-fexact)<2e-7
    assert abs(finite_bubble_muy(0.5)-(1-math.log(2)))<1e-14
    # Independent midpoint Feynman-parameter quadrature at non-degenerate masses.
    for r in (0.76,0.2,0.01):
        n=20000
        numeric=-sum(math.log(r+(1-r)*(i+0.5)/n) for i in range(n))/n
        assert abs(numeric-finite_bubble_muy(r)) < 2e-6
    assert finite_bubble_muy(1.)==0.
    a=(('lam','y','y'),'PR',2)
    assert not compare({a:sp.Rational(1,2)},{a:sp.Rational(2,4)})
    assert len(compare({a:sp.Rational(1,2)},{a:sp.Integer(1)}))==1
    # Identity of two Wick contractions, and exclusion of XX, is reused from
    # the component contraction script rather than restated as a made-up CT.
    print('PASS: exact finite bubble formula, equal-mass and light-X limits')
    print('PASS: independent coefficient-map equality/mismatch detection')
    print('PASS: no EFT counterterms or physical C5 certified')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--self-test',action='store_true')
    ap.add_argument('--output',type=Path)
    ap.add_argument('--root',type=Path,default=BASE)
    ap.add_argument('--native-dir',type=Path,
                    default=AUDIT/'t3_bc_y_native_cg_components_v1')
    ap.add_argument('--ledger-dir',type=Path,
                    default=AUDIT/'t3_bc_y_renormalization_ledger_v1')
    ap.add_argument('--comparison',type=Path,
                    default=AUDIT/'t3_bc_y_native_cg_components_v1'/'comparison.json')
    opt=ap.parse_args()
    try:
        if opt.self_test:self_test()
        if opt.output:
            v3=load_sibling('CompareT3BCYNativeSU2Coefficients.experimental.py','t3_v3_native')
            v2=load_sibling('AuditT3BCYMixedScalarPoleSupport.experimental.py','t3_v2_support')
            require(opt.comparison.is_file(),
                    'Missing completed v3 comparison: '+str(opt.comparison))
            previous=json.loads(opt.comparison.read_text(encoding='utf-8'))
            require(previous.get('format')=='t3bc-y-native-component-kernel-cross-check-v1',
                    'Wrong v3 comparison report format')
            report={'format':'t3bc-y-bubble-mass-finite-audit-v1','models':{},
                    'physical_C5_Y':None,'MSbar_counterterms_certified':False,
                    'physical_matching_certified':False}
            for stem,model in STEMS.items():
                one=run_one(v3,v2,opt.root,stem,opt.native_dir,opt.ledger_dir,previous)
                report['models'][model]=one
                print(f"{model}: pole matches={one['tree_UV_kernel_equals_native_pole']}; "
                      f"formal finite matches XY={one['tree_XY_kernel_equals_formal_finite_at_muY']}; "
                      f"YY/XY={one['nonzero_keys']['YY']}/{one['nonzero_keys']['XY']}; "
                      f"full-mass bubble F_XY={one['scalar_bubble_msbar_finite_at_muY']['XY_full_massive']:.9f}")
            opt.output.parent.mkdir(parents=True,exist_ok=True)
            opt.output.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
            print('Saved independent finite-structure diagnostic: '+str(opt.output))
            if not all(x['tree_UV_kernel_equals_native_pole'] and
                       x['tree_XY_kernel_equals_formal_finite_at_muY']
                       for x in report['models'].values()):
                print('FAIL: component-level pole/finite structure mismatches',file=sys.stderr)
                return 1
        if not (opt.self_test or opt.output):ap.error('Choose --self-test and/or --output')
        return 0
    except (OSError,ValueError,KeyError,TypeError,IndexError,AssertionError) as e:
        print('FAIL (no finite matching certificate): '+str(e),file=sys.stderr)
        return 1

if __name__=='__main__':raise SystemExit(main())
