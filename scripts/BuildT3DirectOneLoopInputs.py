"""Build PROVISIONAL all-16 UV Weyl-component inputs from existing Matchete CG audits.

Reads frozen benchmark JSONs and repository-generated group-component audits;
never invents UV masses or Yukawa entries. Source tensor rotation and charge
selection are exact, but 4-component -> two-component spinor phases and the
source mass/potential normalization still require independent certification.
Therefore every generated input is marked UNVERIFIED and cannot pass the
strict RunT3DirectOneLoop gate without an explicit convention audit.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
from itertools import product
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Numerical.direct_one_loop.core import encoded # noqa: E402
from RunT3DirectOneLoop import models,CLASSES,atomic_json # noqa: E402

FLAGS={'T3Y1CG':(True,True,False),
       'T3Y2CG':(True,False,True),
       'T3MixCG':(False,False,False,True)}


def _spmat(raw):
    return sp.Matrix([[sp.sympify(z,locals={'I':sp.I}) for z in row] for row in raw])


def _cnumber(z):
    return complex(sp.N(z,16))


def _number(v):
    return complex(v['real'],v['imag']) if isinstance(v,dict) else complex(v)


def _matrix(v):
    if isinstance(v,dict):
        return np.asarray(v['real'],float)+1j*np.asarray(v['imag'],float)
    return np.asarray(v,complex)


def _rotations(analysis):
    if analysis.get('status') not in ('generator_algebra_verified_CG_basis_unverified','generator_algebra_verified','all_representations_verified',None):
        # Individual representation rows are the authoritative check.
        raise ValueError('Unknown generator audit status')
    us={1:sp.eye(1)}
    for e in analysis['representations']:
        if e.get('status')!='generator_algebra_verified':
            raise ValueError('Unverified SU(2) generator basis: '+str(e.get('dimension')))
        d=int(e['dimension'])
        U=_spmat(e['U_columns_are_weight_eigenvectors'])
        if sp.simplify(U.H*U-sp.eye(d))!=sp.zeros(d):
            raise ValueError('Nonunitary SU(2) generator diagonalization')
        us[d]=U
    if set(us)!={1,2,3}:
        raise ValueError('Expected SU(2) 1,2,3 generator rotations')
    return us


def _tensor(vertex, flags, us):
    dims=vertex['dimensions']
    if len(flags)!=len(dims):raise ValueError('Wrong CG tensor dimensions')
    maps=[us[d].conjugate() if c else us[d] for d,c in zip(dims,flags)]
    data=np.zeros(tuple(dims),complex)
    for component in vertex['components']:
        old=tuple(i-1 for i in component['indices'])
        val=sp.sympify(component['cg'],locals={'I':sp.I})
        for new in product(*(range(d) for d in dims)):
            coeff=val
            for axis in range(len(dims)):
                coeff*=maps[axis][old[axis],new[axis]]
            if coeff != 0:
                data[new]+=_cnumber(coeff)
    return data


def _charge(d,y,index):
    return (d-1)/2-index+y


def _source_hash(path:Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def construct(case,benchmark,us,model,source_benchmark:Path,source_cg_audit:Path,
              generator_audit:Path, covariant:bool) -> dict:
    """Source-auditable provisional basis translation for the restricted fit.

    This does NOT assert equivalence of spinor mass/phase conventions to the
    live Matchete UV Lagrangian. It is deliberately not a certified result.
    """
    if not covariant:raise ValueError('Source CG tensors not SU2 covariant')
    m=case['model'];ds1,ds2,df,a=(int(m[k]) for k in ('dS1','dS2','dF','alpha'))
    info=benchmark['representation']
    if (ds1,ds2,df,a)!=(info['d_s1'],info['d_s2'],info['d_f'],info['alpha']):
        raise ValueError('Source CG model and source benchmark disagree')
    cg={name:_tensor(case['vertices'][name],f,us) for name,f in FLAGS.items()}
    o=benchmark['base_state']['ordinary']
    y1=_matrix(o['y1']);y2=_matrix(o['y2']);mf=_matrix(o['MF'])
    if y1.shape!=(3,3) or y2.shape!=(3,3) or mf.shape!=(3,3):
        raise ValueError('Invalid source flavour shape')
    v=float(benchmark['numerical']['vev_gev'])
    if v<=0:raise ValueError('Invalid source vev')
    maj=(df%2==1 and a==-1)
    nc=3*df;nf=nc if maj else 2*nc;ns=ds1+ds2
    mm=np.zeros((nf,nf),complex)
    metric=np.asarray(us[df].T*us[df],complex)
    # Matchete's self-conjugate integer-isospin fermions have a nontrivial
    # charge-pairing metric when expressed in a T3 weight basis.
    if maj:
        for r in range(3):
            for t in range(3):
                mm[r*df:(r+1)*df,t*df:(t+1)*df]=mf[r,t]*metric
    else:
        for r in range(3):
            for t in range(3):
                for w in range(df):
                    mm[r*df+w,nc+t*df+w]=mf[r,t]
                    mm[nc+t*df+w,r*df+w]=mf[r,t]
    if np.linalg.norm(mm-mm.T)>1e-9*max(1.,np.linalg.norm(mm)):
        raise ValueError('Input fermion mass not symmetric in constructed Weyl basis')
    m1=float(o['mS1Sq']);m2=float(o['mS2Sq'])
    # Source Lagrangian uses +lambda_Hi (H^dag H)(Si^dag Si), with V=-L.
    # This exact sign is provisional until tested with the UV quadratic action.
    A=np.diag([m1-float(o['lambdaH1'])*v*v/2]*ds1+
              [m2-float(o['lambdaH2'])*v*v/2]*ds2).astype(complex)
    lam=_number(o['lambdaT3'])
    for s in range(ds1):
        for t in range(ds2):
            if abs(cg['T3MixCG'][1,1,s,t])<1e-13:continue
            if abs(_charge(ds1,a/2,s)-_charge(ds2,(a+2)/2,t))>1e-12:
                raise ValueError('VEV mixing connects scalars of different electric charge')
            # L includes +lambdaT3 H^0 H^0 S1 S2^dag +h.c.
            Vcoeff=-lam*v*v*cg['T3MixCG'][1,1,s,t]/2
            A[ds1+t,s]+=Vcoeff
            A[s,ds1+t]+=Vcoeff.conjugate()
    B=np.zeros((ns,ns),complex)
    # Source frozen benchmarks contain no extra quadratic inputs. The
    # generated point is the restricted benchmark with those set to zero.
    for extra in ('lambdaHHdagS2S2','lambdaHHdagS1barS1bar',
                  'lambdaHHdagS1barS2barCross'):
        if extra in o and abs(_number(o[extra]))>0:
            raise ValueError('Nonzero additional Higgs portal requires a source CG component map: '+extra)
    for extra in ('lambdaH1Adj','lambdaH2Adj','lambdaS1Adj','lambdaS2Adj'):
        if extra in o and abs(_number(o[extra]))>0:
            raise ValueError('Nonzero adjoint portal requires component normalization: '+extra)
    # The HC of Bar[L] F^c S1 and Bar[L] F S2^dag is represented in an
    # all-left-handed Weyl basis. Complex-conjugation and chiral phases need
    # a further independent Feynman-rule check against native Matchete.
    y1dag=np.zeros((3,nf,ns),complex)
    y2plain=np.zeros((3,nf,ns),complex)
    for r in range(3):
        for f in range(df):
            q=_charge(df,(a+1)/2,f)
            for s in range(ds1):
                factor=cg['T3Y1CG'][0,f,s].conjugate()
                if abs(factor)<1e-14:continue
                if abs(q-_charge(ds1,a/2,s))>1e-11:
                    raise ValueError('Charge mismatch at Yukawa1 vertex')
                index=r*df+f if maj else nc+r*df+f
                y1dag[:,index,s]+=np.conjugate(y1[:,r])*factor
            for s in range(ds2):
                factor=cg['T3Y2CG'][0,f,s].conjugate()
                if abs(factor)<1e-14:continue
                if abs(q-_charge(ds2,(a+2)/2,s))>1e-11:
                    raise ValueError('Charge mismatch at Yukawa2 vertex')
                if maj:
                    for j in range(df):
                        y2plain[:,r*df+j,ds1+s]+=np.conjugate(y2[:,r])*factor*metric[f,j]
                else:
                    y2plain[:,r*df+f,ds1+s]+=np.conjugate(y2[:,r])*factor
    provenance={
        'benchmark_path':str(source_benchmark),
        'benchmark_sha256':_source_hash(source_benchmark),
        'cg_support_path':str(source_cg_audit),
        'cg_support_sha256':_source_hash(source_cg_audit),
        'generator_analysis_path':str(generator_audit),
        'generator_analysis_sha256':_source_hash(generator_audit),
        'cg_generator_covariance_checked':True,
        'validated_uv_component_conventions':False,
        'unverified_items':['source lambdaT3 mixing normalization/sign',
                            'Matchete four-spinor to Weyl phase convention',
                            'self-conjugate fermion charge-pairing phase',
                            'extra symmetry-allowed Yukawas/quadratics not supplied by the restricted benchmark'],
        'input_assumptions':['frozen fitted_reference benchmark evaluated at base_state.mu_gev',
                             'all unspecified UV general-coupling extras set to zero',
                             'source Yukawas at the same UV scale; no RGE running']}
    return {'model':model,'renormalization_scale_gev':float(benchmark['base_state']['mu_gev']),
            'provenance':provenance,'fermion_mass_gev':encoded(mm),
            'scalar_A_gev2':encoded(A),'scalar_B_gev2':encoded(B),
            'vertices':[{'name':'y1_HC','conjugate':encoded(y1dag)},
                        {'name':'y2_HC','field':encoded(y2plain)}],
            'suppression_mass_ratios':{'max_scalar_mass_over_min_fermion_mass':
                float(np.sqrt(max(np.linalg.eigvalsh(A)))/min(np.linalg.svd(mf)[1])),
                'mS1Sq_over_MFminSq':m1/min(np.linalg.svd(mf)[1])**2,
                'mS2Sq_over_MFminSq':m2/min(np.linalg.svd(mf)[1])**2}}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--project-root',type=Path,default=ROOT)
    ap.add_argument('--audit-dir',type=Path,default=None)
    ap.add_argument('--output',type=Path,default=ROOT/'input'/'lo_direct_one_loop')
    args=ap.parse_args()
    proj=args.project_root
    audit=args.audit_dir or (proj/'Reports'/'output'/'full'/'multiplicity')
    cp=audit/'component_support_audit.json'
    gp=audit/'matchete_generator_analysis.json'
    cv=audit/'cg_covariance_audit.json'
    for p in (cp,gp,cv):
        if not p.is_file():ap.error('Missing source audit file: '+str(p))
    cases=json.loads(cp.read_text(encoding='utf-8'))
    analysis=json.loads(gp.read_text(encoding='utf-8'))
    covariance=json.loads(cv.read_text(encoding='utf-8'))
    if cases.get('errors') or len(cases.get('cases',[]))!=16:
        ap.error('CG source does not have 16 error-free models')
    if covariance.get('status')!='cg_generator_covariance_verified':
        ap.error('CG covariance was not verified')
    us=_rotations(analysis)
    lookup={tuple(case['model'][k] for k in ('dS1','dS2','dF','alpha')):case
            for case in cases['cases']}
    valid={tuple((m['model'][k] for k in ('dS1','dS2','dF','alpha'))):m['all_covariant']
           for m in covariance['models']}
    problems=[]
    for short,key in models():
        dims=CLASSES[short[0]][:3]+(int(short[2:]),)
        bench=proj/'Numerical'/'benchmarks'/'fitted_reference'/(key+'.json')
        try:
            case=lookup[dims]
            base=json.loads(bench.read_text(encoding='utf-8-sig'))
            obj=construct(case,base,us,short,bench,cp,gp,bool(valid.get(dims)))
            atomic_json(args.output/(short+'.json'),obj)
            print(short,'provisional component input built; NOT source-certified')
        except (ValueError,KeyError,IndexError,OSError,TypeError) as e:
            problems.append({'model':short,'error':str(e)})
            print(short,'BLOCKED:',e)
    atomic_json(args.output/'build_status.json',{'count':16,'failed':problems,
                'all_inputs_certified':False})
    return 2 if problems else 0

if __name__=='__main__':raise SystemExit(main())
