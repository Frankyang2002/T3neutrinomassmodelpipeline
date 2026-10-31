#!/usr/bin/env python3
"""Independent signed one-loop triangle master test of the v16 *formal* F/Y/X envelope.

This is a source-validated topology comparison, NOT full renormalized T3 matching.
The UV three-propagator zero-momentum loop is evaluated independently of
Matchete in Euclidean signature. Its overall negative sign is selected so its
strict hard/soft expansion can be compared with the signed Matchete source.

Let P_F+P_Y+P_X=0 (proved in v15) and K_F,K_Y,K_X be native PL ordered
constant sectors with common hbar/MF removed. With K_X=0 and *if verified*
K_F=P_F, the strict source envelope decomposes into scalar-pair sectors:
    W_XX = -P_X, W_XY = -K_Y, W_YY = -P_Y+K_Y
and
    K_strict = W_XX * [1+log(MX^2/MF^2)]
             + W_YY * [1+log(MY^2/MF^2)]
             + W_XY * [log(MY^2/MF^2)].
These identities are checked in the actual B/C polynomial algebra; no
external vertex sign or physical Wilson-coefficient identification is implied.

The *full massive one-loop triangle* is a distinct reference. The leading
MF^-1 expansion with finite MX gives an additional XY correction
    deltaK = K_Y * r log(r)/(1-r),  r=MX^2/MY^2.
That term is formally subleading only in the MX/MY expansion, which is poor
for the current B/C masses. Do NOT blindly append this correction to the
sequential dim-5 matching: doing so would mix EFT truncation orders.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys

ROOT=Path.cwd()
MODELS={'B_m1':'B_-1','C_m1':'C_-1'}
AUD=Path('output/audits')
BASE=Path('output/general_t3_mass_eigen_matchete')

class TriangleError(ValueError):pass

def require(ok,msg):
    if not ok:raise TriangleError(msg)

def read_json(path,fmt=None):
    require(path.is_file(),'Missing local source or audit '+str(path))
    d=json.loads(path.read_text(encoding='utf-8-sig'))
    require(isinstance(d,dict),'Expected JSON object in '+str(path))
    if fmt:require(d.get('format')==fmt,'Changed source evidence format '+str(path))
    return d

def script_import(path,name):
    require(path.is_file(),'Missing v16 script '+str(path))
    spec=importlib.util.spec_from_file_location(name,path)
    require(spec is not None and spec.loader is not None,'Script import failed')
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def signed_triangle_master(A,B,C):
    """-A times positive Euclidean radial three-propagator integral.

    Exact analytic expression, dimensionless. B and C are scalar mass^2,
    A is heavy fermion mass^2. Uses mpmath to avoid degeneracy loss.
    """
    import mpmath as mp
    A,B,C=mp.mpf(str(A)),mp.mpf(str(B)),mp.mpf(str(C))
    require(all(v>0 for v in (A,B,C)),'Master requires positive squared masses')
    mp.mp.dps=48
    def g(x):
        z=(A-x)/A
        if abs(z)<mp.mpf('1e-22'):
            # x*ln(A/x)/(A-x) -> 1 as x -> A
            return (1-z)*(-mp.log1p(-z))/z if z else mp.mpf(1)
        return x*mp.log(A/x)/(A-x)
    if abs(B-C)<=mp.mpf('1e-22')*max(B,C):
        if abs(A-B)<=mp.mpf('1e-22')*A:
            return -mp.mpf('0.5')
        return -A*(A*mp.log(A/B)-(A-B))/(A-B)**2
    return -A*(g(B)-g(C))/(B-C)

def triangle_quad(A,B,C):
    import mpmath as mp
    A,B,C=map(lambda v:mp.mpf(str(v)),(A,B,C));mp.mp.dps=42
    return -A*mp.quad(lambda t:t/((t+A)*(t+B)*(t+C)),[0,mp.inf])

def psi_strict(A,B,C):
    """Leading 1/MF master at finite MX, not a full-threshold coefficient."""
    import mpmath as mp
    A,B,C=map(lambda v:mp.mpf(str(v)),(A,B,C))
    if abs(B-C)<mp.mpf('1e-22')*max(B,C):
        return 1+mp.log(B/A)
    return (B*mp.log(B/A)-C*mp.log(C/A))/(B-C)

def soft_factor(r):
    r=float(r)
    require(0<r<1 and math.isfinite(r),'Require 0<MX^2/MY^2<1')
    return r*math.log(r)/(1-r)

def check_three_mass_integral():
    import mpmath as mp
    for A,B,C in ((100.,9.,4.),(100.,9.,9.),(100.,9.,1e-5),
                  (4.,4.,4.),(1e8,3.,2.)):
        a=signed_triangle_master(A,B,C); b=triangle_quad(A,B,C)
        require(abs(a-b)<mp.mpf('1e-25')*max(1,abs(b)),
                f'Analytic and integral master mismatch at {A,B,C}: {a}, {b}')
    # The strict leading result must converge for all three pair types.
    for B,C in ((9.,9.),(9.,4.),(4.,4.)):
        exact=signed_triangle_master(1e11,B,C)
        strict=psi_strict(1e11,B,C)
        require(abs(exact-strict)<mp.mpf('1e-7'),
                'Heavy-fermion master asymptotic mismatch')
    # Independent algebraic pair-decomposition test, including tamper rejection.
    import sympy as sp
    px,py,ky,LX,LY=sp.symbols('px py ky LX LY')
    pf=-px-py
    wxx=-px; wxy=-ky; wyy=-py+ky
    source=pf+ky+pf*LX+py*LY
    master=wxx*(1+LX)+wyy*(1+LX-LY)+wxy*(LX-LY)
    require(sp.expand(source-master)==0,'Synthetic signed-pair decomposition mismatch')
    require(sp.expand((source+1)-master)!=0,'Wrong finite constant not rejected')
    r=0.759189881
    correction=soft_factor(r)
    require(abs((1+correction)-0.131434492)<2e-6,
            'Independent XY finite-bubble benchmark mismatch')
    print('PASS: full three-mass triangle analytic = independent convergent integral')
    print('PASS: strict MF^-1 and finite MX/MY asymptotics; v7 XY benchmark reproduced')
    print('PASS: independent XX/YY/XY signed master identity; wrong finite constant rejected')

def source_run(root,out,numeric_input):
    import sympy as sp
    import mpmath as mp
    v16=script_import(root/'scripts/AssembleT3BCFormalFiniteScaleEnvelope.experimental.py','t3bc_v17_v16')
    mods=v16.import_deps(root)
    pole,kernel,v12,v13,v14=mods
    v14report=v16.read_json(root/AUD/'t3_bc_sequential_batch_v14/summary.json',
                    't3-bc-sequential-matching-fixed-order-batch-v14')
    v15report=v16.read_json(root/AUD/'t3_bc_one_loop_inheritance_v15/summary.json',
                    't3-bc-one-loop-fixed-order-inheritance-v15')
    v16report=v16.read_json(root/AUD/'t3_bc_formal_finite_envelope_v16/summary.json',
                    't3-bc-formal-finite-common-scale-envelope-v16')
    extract=v16.read_json(root/AUD/'t3_bc_scalar_finite_term_extraction_v1/summary.json',
                    't3-bc-native-local-matchete-pole-log-split-v1')
    require(v16report.get('physical_matching_certified') is False and
            v16report.get('mass_log_scale_cancellation_exact') is True,
            'Unexpected v16 certification or missing scale identity')
    report={'format':'t3-bc-independent-three-propagator-topology-comparison-v17',
        'scope':'strict signed one-loop scalar-quartic/LLSS triangle topology',
        'status':'IndependentTriangleComparison_PendingModelEvaluation',
        'source_renormalized_matching_certified':False,
        'external_fermion_vertex_sign_independently_certified':False,
        'mass_power_counting':'strict dimension-five/MF first, then optional MX/MY expansion',
        'physical_C5':None,'physical_neutrino_mass_matrix':None,'models':{}}
    LF,LY=sp.symbols('LF LY') # LF=ln(MX^2/MF^2), LY=ln(MX^2/MY^2)
    polynomials={}
    for stem,model in MODELS.items():
        require(v15report['models'][stem]['F_Y_X_formal_pole_relation']
                ['formal_neutral_residue_relations']['F_plus_Y_plus_X'] is True,
                model+': missing v15 pole cancellation')
        stages,poly=v16.compare_finite_and_poles(root,stem,v14report,extract,mods)
        fpath=root/BASE/'eft1_fermion_first_oneloop'/f'{stem}_final_eft1_fermion_first_oneloop_raw_hard_llhh_sector.inputform.txt'
        fsec=v13.split_sectors(fpath.read_text(encoding='utf-8-sig'))
        KF=v16.component_polynomials(fsec['constant'],model,kernel,v12,68)
        PF,PY,PX=poly['F'],poly['Y'],poly['X']
        KT=poly['constant']
        KY=tuple(sp.expand(KT[i]-KF[i]) for i in range(3))  # K_X=0 strictly verified in v16
        fminus=tuple(sp.expand(KF[i]-PF[i]) for i in range(3))
        bad=[i for i,v in enumerate(fminus) if v!=0]
        # This is an additional, nontrivial source-level equality, NOT a priori assumed.
        require(not bad,model+': native F finite != native F pole; need inspect actual subtraction')
        WXX=tuple(-x for x in PX)
        WXY=tuple(-y for y in KY)
        WYY=tuple(sp.expand(-PY[i]+KY[i]) for i in range(3))
        for i in range(3):
            src=sp.expand(KT[i]+PF[i]*LF+PY[i]*LY)
            tri=sp.expand(WXX[i]*(1+LF)+
                    WYY[i]*(1+LF-LY)+WXY[i]*(LF-LY))
            require(sp.expand(src-tri)==0,
                model+f': independent pair-master decomposition failed at Yukawa {i}')
        params=read_json(root/BASE/f'{stem}_final_parameters.json')
        masses=params['heavy_masses_gev'];require(len(masses)==2,'Missing scalar mass eigenvalues')
        MX,MY=map(float,masses)
        r=(MX/MY)**2; delta=soft_factor(r)
        # Signed exact-leading finite-MX triangle differs from r->0 expansion
        # by delta*K_Y, using the source-specific XY finite sign.
        polynomials[stem]=(WXX,WYY,WXY,KY,poly,params)
        report['models'][stem]={
            'model':model,'source_verified_F_finite_equals_F_pole':True,
            'source_verified_K_X_constant_zero':True,
            'source_verified_P_F_plus_P_Y_plus_P_X_zero':True,
            'strict_pairwise_triangle_master_equals_v16_envelope':True,
            'ordered_Yukawa_basis':['y1*y1','y1*y2','y2*y2'],
            'signed_W_XX_3spurions':[str(p) for p in WXX],
            'signed_W_XY_3spurions':[str(p) for p in WXY],
            'signed_W_YY_3spurions':[str(p) for p in WYY],
            'finite_XY_source_K_Y_3spurions':[str(p) for p in KY],
            'r_MX2_over_MY2':r,
            'XY_delta_full_mass_minus_strict_limit_multiplying_K_Y':delta,
            'XY_full_mass_finite_bubble_factor':1+delta,
            'soft_corrected_but_NOT_matching':'K_strict + K_Y * r*log(r)/(1-r)',
            'full_renormalized_subtraction_done':False,
            'Y_X_soft_overlap_included_at_same_EFT_truncation':False,
            'physical_C5':None}
        print(f'{model}: F constant=pole PASS; signed XX/YY/XY triangle decomposition PASS; '
              f'XY omitted soft remainder={delta:.9f} times K_Y')
    if numeric_input is not None:
        # Optional full-mass independent reference; does NOT replace matching.
        import numpy as np
        data=read_json(numeric_input)
        if set(data)==set(MODELS.values()):
            records={stem:data[model] for stem,model in MODELS.items()}
        else:
            require(data.get('model') in MODELS.values(),'Unrecognized numeric model')
            stem=next(k for k,v in MODELS.items() if v==data['model'])
            records={stem:data}
        for stem,rec in records.items():
            WXX,WYY,WXY,KY,poly,params=polynomials[stem]
            formal=v16.numeric_evaluate(root,stem,rec,poly,params,v12)
            U=np.array(params['unitary_U']['real'],float)+1j*np.array(params['unitary_U']['imag'],float)
            quartics={k:v12.parse_complex(v,k) for k,v in rec['quartics'].items()}
            sub={}
            for a in range(2):
                for b in range(2):
                    n=f'T3MEU{a+1}{b+1}'
                    sub[sp.Symbol(n)]=U[a,b]
                    sub[sp.Symbol('bar_'+n)]=U[a,b].conjugate()
            for n,v in quartics.items():
                sub[sp.Symbol(n)]=v
                sub[sp.Symbol('bar_'+n)]=v.conjugate()
            def val(vec):
                for e in vec:
                    require(not(e.free_symbols-set(sub)),
                            'Unexpected benchmark spurion in numerical master')
                return np.array([complex(sp.N(e.subs(sub),16)) for e in vec])
            mf=float(rec['MF_GeV']); MX,MY=map(float,params['heavy_masses_gev'])
            A,B,C=mf**2,MY**2,MX**2
            Pxx=1+math.log(C/A);Pyy=1+math.log(B/A);Pxy=math.log(B/A)
            dxy=soft_factor(C/B)
            weights=(val(WXX),val(WYY),val(WXY))
            s=weights[0]*Pxx+weights[1]*Pyy+weights[2]*Pxy
            sc=weights[0]*Pxx+weights[1]*Pyy+weights[2]*(Pxy-dxy)
            full=weights[0]*float(signed_triangle_master(A,C,C))+ \
                 weights[1]*float(signed_triangle_master(A,B,B))+ \
                 weights[2]*float(signed_triangle_master(A,B,C))
            ys=[np.array([v12.parse_complex(z,n) for z in rec[n]],complex)
                for n in ('y1','y2')]
            L=np.column_stack([ys[0].conjugate(),ys[1].conjugate()])
            def to3(k):
                B2=np.array([[2*k[0],k[1]],[k[1],2*k[2]]],complex)
                return v16.to_complex_matrix((1/(16*math.pi**2*mf))*L@B2@L.T)
            recout={'formal_v16':formal,
               'strict_triangle_3x3_GeVinv_NOT_PHYSICAL':to3(s),
               'soft_XY_corrected_3x3_GeVinv_NOT_PHYSICAL':to3(sc),
               'full_mass_triangle_topology_3x3_GeVinv_NOT_PHYSICAL':to3(full),
               'full_mass_requires_independent_UV_Lagrangian_vertex_check':True,
               'full_matching_scheme_and_EFT_order_not_certified':True,
               'physical_C5':None}
            rep=report['models'][stem]
            rep['optional_numeric_triangle_comparison']=recout
            require(np.max(np.abs(s-(val(poly['constant'])+
                2*math.log(MX/mf)*val(poly['F'])+
                2*math.log(MX/MY)*val(poly['Y']))))<1e-9*
                max(1,float(np.max(np.abs(s)))), 'Numeric strict master mismatch')
            print(f'{MODELS[stem]}: optional strict/soft/full-mass 3x3 topology reference written (NOT C5)')
    report['status']='ExactStrictPairMasterMatchedToNativeSources_FullTriangleReferenceAvailable'
    report['numeric_input_sha256']=hashlib.sha256(numeric_input.read_bytes()).hexdigest() if numeric_input else None
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print('PASS: independent triangle topology checked; physical C5 NOT CERTIFIED')
    print('Saved: '+str(out))

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--self-test',action='store_true')
    ap.add_argument('--root',type=Path,default=Path.cwd())
    ap.add_argument('--output',type=Path,
        default=Path('output/audits/t3_bc_triangle_master_v17/summary.json'))
    ap.add_argument('--numeric-input',type=Path)
    ap.add_argument('--require-physical',action='store_true')
    args=ap.parse_args()
    try:
        if args.self_test:check_three_mass_integral();return 0
        if args.require_physical:
            raise TriangleError('Not a physical C5: EFT full finite subtraction and scheme missing')
        root=args.root.resolve();out=args.output if args.output.is_absolute() else root/args.output
        ni=args.numeric_input
        if ni is not None and not ni.is_absolute():ni=root/ni
        source_run(root,out,ni);return 0
    except (TriangleError,ValueError,TypeError,KeyError,OSError,ImportError,IndexError) as exc:
        print('FAIL (physical C5 unchanged): '+str(exc),file=sys.stderr)
        return 1

if __name__=='__main__':sys.exit(main())
