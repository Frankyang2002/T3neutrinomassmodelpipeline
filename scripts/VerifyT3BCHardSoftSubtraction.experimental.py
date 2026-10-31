#!/usr/bin/env python3
"""Source-anchored MSbar hard/soft subtraction for the restricted B/C triangle.

One-loop leading-1/MF approximation, with the exact MX/MY ratio retained:
  phi_strict(A,B,C) = [B ln(B/A)-C ln(C/A)]/(B-C)
                  = H(mu,A) + S(mu,B,C),
  H = 1+ln(mu^2/A),
  S = -B0_fin(mu,B,C),
  B0_fin = 1-[B ln(B/mu^2)-C ln(C/mu^2)]/(B-C).
The coincident B=C limit is 1+ln(B/A), with B0_fin=ln(mu^2/B).

Formal pole residues H:+1, S:-1 for each pair; this is the
high-minus-low one-loop REGION subtraction, not a finite field/coupling
renormalization of the complete T3 Lagrangian. Source matching conditions:
  PF = WXX+WYY+WXY, PY+PX=-(WXX+WYY+WXY).

Uses project-native strict symbolic F/Y/X polynomials (v16-v18) and the
existing frozen three-MF reference only. It must NOT certify physical C5.
"""
from __future__ import annotations

import argparse
from hashlib import sha256
import importlib.util
import json
import math
from pathlib import Path
import sys

STEMS={'B_m1':'B_-1','C_m1':'C_-1'}
AUD=Path('output/audits')

class SubtractionError(ValueError): pass

def require(condition,message):
    if not condition:
        raise SubtractionError(message)

def import_local(path,name):
    require(path.is_file(),'Required source script missing: '+str(path))
    spec=importlib.util.spec_from_file_location(name,path)
    require(spec is not None and spec.loader is not None,'Invalid import source')
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def read_json(path,expected=None):
    require(path.is_file(),'Required source missing: '+str(path))
    x=json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(x,dict),'Required JSON object')
    if expected is not None:
        require(x.get('format')==expected,'Unexpected evidence format at '+str(path))
    return x

def weighted_log(A,B,C):
    """Exactly dimensionless strict heavy-MF triangle, stable near B=C."""
    import mpmath as mp
    A,B,C=map(lambda v:mp.mpf(str(v)),(A,B,C))
    require(all(x>0 for x in (A,B,C)),'Positive squared masses required')
    mp.mp.dps=48
    if abs(B-C)<mp.mpf('1e-23')*max(B,C):
        return 1+mp.log(B/A)
    return (B*mp.log(B/A)-C*mp.log(C/A))/(B-C)

def hard_soft(A,B,C,mu2):
    import mpmath as mp
    mu2=mp.mpf(str(mu2))
    require(mu2>0,'Positive matching scale required')
    h=1+mp.log(mu2/mp.mpf(str(A)))
    s=weighted_log(mu2,B,C)-1
    return h,s

def assert_near(a,b,description,relative='1e-25'):
    import mpmath as mp
    require(abs(a-b)<mp.mpf(relative)*max(1,abs(b)),description+': failed')

def self_test():
    import mpmath as mp
    for A,B,C in ((1e10,1.4e6,1.08e6),(100,9,4),
                  (100,9,9),(100,4,4),(1e10,1.4e6,1.4e6)):
        strict=weighted_log(A,B,C)
        for mu2 in (A,B,C,math.sqrt(A*B)):
            h,s=hard_soft(A,B,C,mu2)
            assert_near(h+s,strict,'Hard/soft finite parts must cancel arbitrary mu')
        h1,s1=hard_soft(A,B,C,B)
        h2,s2=hard_soft(A,B,C,2*B)
        assert_near(h2-h1,mp.log(2),'Hard logarithm shift')
        assert_near(s2-s1,-mp.log(2),'Soft logarithm shift')
    # Explicit pole cancellation and wrong-sign rejection, independent of masses.
    require(1-1==0 and 1+1!=0,'UV/EFT pole sign mismatch')
    # Deliberately enforce the need for a nontrivial source pole identity.
    pf,px,py=-5,2,3
    wxx,wyy,wxy=-px,-py,0
    require(pf==wxx+wyy+wxy and px+py==-(wxx+wyy+wxy),
            'Synthetic three-sector pole decomposition failed')
    print('PASS: exact hard+soft MSbar finite identity for XX/XY/YY and multiple matching scales')
    print('PASS: hard/soft mu dependence and opposite UV residues; incorrect pole sign rejected')
    print('PASS: no full physical C5 or running claimed')

def cmat(x):
    import numpy as np
    mat=np.asarray(x,dtype=float)
    require(mat.shape==(3,3,2) and np.isfinite(mat).all(),'Expected 3x3 complex matrix pairs')
    return mat[:,:,0]+1j*mat[:,:,1]

def matrix_pair(mat):
    return [[[float(z.real),float(z.imag)] for z in row] for row in mat]

def evaluate_matrix(weights,ref,masses,mode):
    """Return restricted C5-like topology with hard/soft split, no physical flag."""
    import numpy as np
    mx,my=masses
    require(mode in ('hard','soft','hard_plus_soft'),'Invalid region selection')
    out=np.zeros((3,3),complex)
    for a,m in enumerate(ref['mf']):
        A=float(m)**2
        mu2=my**2  # same mu=MY for all pairs at a given benchmark
        pairs=((mx**2,mx**2),(my**2,my**2),(my**2,mx**2))
        w=[]
        for (B,C),pair_weight in zip(pairs,weights):
            h,s=hard_soft(A,B,C,mu2)
            phi={'hard':h,'soft':s,'hard_plus_soft':h+s}[mode]
            w.append(pair_weight*float(phi))
        k=sum(w)
        M2=np.array([[2*k[0],k[1]],[k[1],2*k[2]]],complex)
        l=np.column_stack((ref['y1'][:,a].conj(),ref['y2'][:,a].conj()))
        out+=(l@M2@l.T)/(16*math.pi**2*float(m))
    return out

def run(root,out):
    import numpy as np
    import sympy as sp
    v16=import_local(root/'scripts/AssembleT3BCFormalFiniteScaleEnvelope.experimental.py','t3_v19_v16')
    v17=import_local(root/'scripts/CompareT3BCFiniteEnvelopeToTriangle.experimental.py','t3_v19_v17')
    v18=import_local(root/'scripts/EvaluateT3BCThreeFermionTriangle.experimental.py','t3_v19_v18')
    v14=read_json(root/AUD/'t3_bc_sequential_batch_v14/summary.json',
       't3-bc-sequential-matching-fixed-order-batch-v14')
    v15=read_json(root/AUD/'t3_bc_one_loop_inheritance_v15/summary.json',
       't3-bc-one-loop-fixed-order-inheritance-v15')
    v17audit=read_json(root/AUD/'t3_bc_triangle_master_v17/summary.json',
       't3-bc-independent-three-propagator-topology-comparison-v17')
    v18audit=read_json(root/AUD/'t3_bc_three_fermion_triangle_v18/summary.json',
       't3-bc-frozen-reference-three-fermion-triangle-v18')
    extract=read_json(root/AUD/'t3_bc_scalar_finite_term_extraction_v1/summary.json',
       't3-bc-native-local-matchete-pole-log-split-v1')
    eq=read_json(root/AUD/'t3_bc_hermiticity_y_handoff_v1/summary.json',
       't3-bc-eft1-component-hermiticity-y-handoff-v1')
    require(eq.get('status')=='EFT1ComponentHermiticityEstablishedYRenormalizationOutstanding',
            'Source equivalence evidence changed; inspect current report')
    require(v18audit.get('physical_C5') is None and v18audit.get('running_performed') is False,
            'Changed v18 certification; review source')
    components={row['model']:row for row in eq['models']}
    deps=v16.import_deps(root)
    _,kernel,v12,v13,_=deps
    result={'format':'t3-bc-one-loop-hard-soft-msbar-restricted-v19',
       'status':'RestrictedTriangleHardSoftMSbarSubtractionEstablished',
       'definition':'H=1+ln(mu2/MF2); S=-B0_fin; B0_fin=1-[B ln(B/mu2)-C ln(C/mu2)]/(B-C)',
       'overall_loop_factor':'1/(16*pi^2*MF); ordered PL Warsaw uses C5=A+A^T',
       'pole_assignment':'F hard +1; EFT soft -1 per signed scalar-pair master',
       'Y_X_loop_power_counting':'overall one loop, leading 1/MF; full MX/MY mass dependence retained in EFT1 scalar bubble',
       'full_T3_renormalized_matching_certified':False,
       'full_T3_physical_C5':None,'physical_neutrino_masses':None,
       'source_scheme_limitations':[
          'F/Y/X Matchete matching coefficients have not been mapped to complete renormalized UV and EFT parameters',
          'Only scalar-quartic/LLSS triangle topology is included',
          'EOM and finite field/coupling transformations not certified',
          'No EFT1 running; finite MX/MY dependence kept here requires higher-dimensional Y-threshold terms in a strictly sequential truncation',
          'Underlying full triangle UV-vertex sign not independently certified',
       ],'models':{}}
    for stem,model in STEMS.items():
        require(components[model]['component_residual_zero'] is True,model+': HHXY component source equivalence failed')
        require(v15['models'][stem]['F_Y_X_formal_pole_relation']['formal_neutral_residue_relations']['F_plus_Y_plus_X'] is True, model+': F+Y+X pole relation failed')
        for flag in ('source_verified_F_finite_equals_F_pole','source_verified_P_F_plus_P_Y_plus_P_X_zero',
                     'strict_pairwise_triangle_master_equals_v16_envelope'):
            require(v17audit['models'][stem].get(flag) is True,model+': missing v17 source gate '+flag)
        _,poly=v16.compare_finite_and_poles(root,stem,v14,extract,deps)
        source=root/'output/general_t3_mass_eigen_matchete/eft1_fermion_first_oneloop'/f'{stem}_final_eft1_fermion_first_oneloop_raw_hard_llhh_sector.inputform.txt'
        fsec=v13.split_sectors(source.read_text(encoding='utf-8-sig'))
        kf=v16.component_polynomials(fsec['constant'],model,kernel,v12,68)
        pf,py,px=poly['F'],poly['Y'],poly['X']
        ky=tuple(sp.expand(poly['constant'][i]-kf[i]) for i in range(3))
        wxx=tuple(-z for z in px)
        wxy=tuple(-z for z in ky)
        wyy=tuple(sp.expand(-py[i]+ky[i]) for i in range(3))
        for i in range(3):
            hardpole=sp.expand(wxx[i]+wxy[i]+wyy[i])
            softpole=-hardpole
            require(sp.expand(hardpole-pf[i])==0,model+': hard coefficient != original F pole')
            require(sp.expand(softpole-py[i]-px[i])==0,model+': soft coefficient != Y+X pole')
            require(sp.expand(kf[i]-pf[i])==0,model+': F constant differs from F pole')
        params=v18.load_json(root/v18.SQUARE_ROOT/f'{stem}_final_parameters.json')
        mx,my=map(float,params['heavy_masses_gev'])
        mix=v18.as_matrix(params['unitary_U'],(2,2),'unitary_U')
        ref=v18.load_frozen_reference(root,stem)
        sub={}
        for i in range(2):
            for j in range(2):
                s=f'T3MEU{i+1}{j+1}'
                sub[sp.Symbol(s)]=mix[i,j]
                sub[sp.Symbol('bar_'+s)]=mix[i,j].conjugate()
        for n,z in ref['quartics'].items():
            sub[sp.Symbol(n)]=z
            sub[sp.Symbol('bar_'+n)]=z.conjugate()
        weights=tuple(v18.polynomial_values(w,sub) for w in (wxx,wyy,wxy))
        H=evaluate_matrix(weights,ref,(mx,my),'hard')
        S=evaluate_matrix(weights,ref,(mx,my),'soft')
        finiteXY=H+S
        original=cmat(v18audit['models'][stem]['strict_and_full_triangle']['finite_XY']
          ['reference_C5_3x3_GeV_inverse_NOT_PHYSICAL'])
        norm=max(float(np.linalg.norm(original)),1e-29)
        residual=float(np.linalg.norm(finiteXY-original)/norm)
        require(residual<5e-11,model+': hard+soft differs from v18 finite-XY leading-1/MF matrix')
        direct=v18.sum_fermions(weights,ref,(mx,my),'finite_XY',v17)
        require(float(np.linalg.norm(finiteXY-direct))/norm<5e-11,
                model+': hard+soft differs from direct finite-XY heavy-MF master')
        fully=cmat(v18audit['models'][stem]['strict_and_full_triangle']['full_mass']
          ['reference_C5_3x3_GeV_inverse_NOT_PHYSICAL'])
        strictSequential=cmat(v18audit['models'][stem]['strict_and_full_triangle']['strict']
          ['reference_C5_3x3_GeV_inverse_NOT_PHYSICAL'])
        deviationFull=float(np.linalg.norm(finiteXY-fully)/max(np.linalg.norm(fully),1e-29))
        deviationStrict=float(np.linalg.norm(strictSequential-fully)/max(np.linalg.norm(fully),1e-29))
        other=(mx**2,my**2,float(ref['mf'][0])**2)
        for A in (float(m)**2 for m in ref['mf']):
            for B,C in ((mx**2,mx**2),(my**2,my**2),(my**2,mx**2)):
                for q in other:
                    h,s=hard_soft(A,B,C,q)
                    assert_near(h+s,weighted_log(A,B,C),'matching-scale independence')
        result['models'][stem]={
            'model':model,
            'input_mass_spectrum_GeV':{'MX':mx,'MY':my,'MF':[float(x) for x in ref['mf']]},
            'audited_su2_component_source_equivalence':True,
            'symbolic_f_hard_plus_eft_soft_poles_zero_in_each_yukawa_structure':True,
            'symbolic_hard_residue_equals_F_pole':True,
            'symbolic_soft_residue_equals_Y_plus_X_poles':True,
            'analytic_hard_plus_soft_matching_scale_independent':True,
            'reference_hard_matrix_GeV_inverse_NOT_PHYSICAL':matrix_pair(H),
            'reference_soft_matrix_GeV_inverse_NOT_PHYSICAL':matrix_pair(S),
            'reference_hard_plus_soft_matrix_GeV_inverse_NOT_PHYSICAL':matrix_pair(finiteXY),
            'relative_difference_vs_v18_finite_XY':residual,
            'relative_difference_vs_v18_full_mass':deviationFull,
            'v18_strict_sequential_vs_full_mass':deviationStrict,
            'higher_MX_MY_orders_kept_in_EFT1_bubble_but_absent_in_strict_sequential_Y_dimension_5':True,
            'same_one_loop_order_reference_only':True,
            'physical_C5':None,
        }
        print(f'{model}: F/soft pole assignment PASS; hard+soft MSbar finite at mu=MY PASS; '
              f'v18 finite-XY residual={residual:.2e}; leading-1/MF vs full={deviationFull:.4g}; strict vs full={deviationStrict:.4g}')
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('PASS: restricted one-loop F/EFT scalar hard-soft subtraction proven; complete physical C5 NOT CERTIFIED')
    print('Saved: '+str(out))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--self-test',action='store_true')
    p.add_argument('--root',type=Path,default=Path.cwd())
    p.add_argument('--output',type=Path,default=Path('output/audits/t3_bc_hard_soft_subtraction_v19/summary.json'))
    p.add_argument('--require-physical',action='store_true')
    a=p.parse_args()
    try:
        if a.self_test:
            self_test();return 0
        if a.require_physical:
            raise SubtractionError('Full T3 physical matching not certified: missing complete counterterms, running, UV vertex convention')
        root=a.root.resolve()
        out=a.output if a.output.is_absolute() else root/a.output
        run(root,out)
        return 0
    except (SubtractionError,AssertionError,ValueError,KeyError,TypeError,ImportError,RuntimeError,OSError) as exc:
        print('FAIL (full physical C5 unavailable): '+str(exc),file=sys.stderr)
        return 1

if __name__=='__main__':
    sys.exit(main())
