"""Independent gauge-invariant SM L/eC/H stress test for the seven-term master beta.

This intentionally FAILS CLOSED (exit 2) while the beta violates Ward identity.
Use --diagnostic for a zero-exit numerical report. No physical A(-4) inputs.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from RGE.running.intermediate.GeneralEFT1SevenTermFlow import Couplings,master_beta,total_beta,TERMS
TOL=2e-11

def real_generator(t):
    out=np.zeros((4,4),complex)
    for p in range(2):
        for q in range(2):
            z=t[p,q]
            out[2*p,2*q] += 1j*z.imag
            out[2*p,2*q+1] += 1j*z.real
            out[2*p+1,2*q] -= 1j*z.real
            out[2*p+1,2*q+1] += 1j*z.imag
    return out

def example():
    """All inputs are gauge invariant: L Y=-1/2, eC Y=+1, H Y=+1/2."""
    sigma=[np.array([[0,1],[1,0]],complex)/2,
           np.array([[0,-1j],[1j,0]],complex)/2,
           np.diag([.5,-.5]).astype(complex)]
    gauge=((.65,tuple((np.pad(t,((0,1),(0,1))),real_generator(t)) for t in sigma)),
           (.36,((np.diag([-.5,-.5,1]).astype(complex),real_generator(.5*np.eye(2))),)))
    y=np.zeros((3,3,4),complex)
    for n in range(2):
        y[2,n,2*n]=y[n,2,2*n]=(.30+.10j)/np.sqrt(2)
        y[2,n,2*n+1]=y[n,2,2*n+1]=-1j*(.30+.10j)/np.sqrt(2)
    lam=np.zeros((4,)*4)
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    lam[a,b,c,d]=.4*((a==b)*(c==d)+(a==c)*(b==d)+(a==d)*(b==c))
    v=(np.array([0,0,1,1j])/np.sqrt(2),-np.array([1,1j,0,0])/np.sqrt(2))
    C=np.zeros((3,3,4,4),complex)
    for i in range(2):
        for j in range(2):
            C[i,j]=.25*(1+.13j)*(np.outer(v[i],v[j])+np.outer(v[j],v[i]))
    return C,Couplings(y,lam,gauge)

def ward_y(y,t,s):
    return (np.einsum('ki,kja->ija',t,y)+np.einsum('kj,ika->ija',t,y)
            +np.einsum('ba,ijb->ija',s,y))

def ward_l(l,s):
    return (np.einsum('ea,ebcd->abcd',s,l)+np.einsum('eb,aecd->abcd',s,l)
            +np.einsum('ec,abed->abcd',s,l)+np.einsum('ed,abce->abcd',s,l))

def ward_c(c,t,s):
    return (np.einsum('ki,kjab->ijab',t,c)+np.einsum('kj,ikab->ijab',t,c)
            +np.einsum('ca,ijcb->ijab',s,c)+np.einsum('cb,ijac->ijab',s,c))

def diagnose():
    C,x=example(); generators=[(t,s) for g,ps in x.gauge for t,s in ps]
    maximum=lambda arrays:max(float(np.max(np.abs(a))) for a in arrays)
    ins={'yukawa':maximum(ward_y(x.yukawa,t,s) for t,s in generators),
         'quartic':maximum(ward_l(x.quartic,s) for t,s in generators),
         'wilson':maximum(ward_c(C,t,s) for t,s in generators)}
    if max(ins.values())>TOL: raise RuntimeError('Toy inputs are not Ward invariant')
    variants={}
    for name,flag in [('repository_two_simultaneous_exchanges',False),
                      ('batch34_four_independent_exchanges',True)]:
        terms=master_beta(C,x,corrected_crossed=flag); beta=total_beta(terms)
        numer=maximum(ward_c(beta,t,s) for t,s in generators)
        den=float(np.max(np.abs(beta)))
        variants[name]={
          'full_beta_relative_ward_residual':numer/max(den,1e-300),
          'full_beta_ward_absolute':numer,
          'per_term_absolute_ward_residual':{
              key:maximum(ward_c(terms[key],t,s) for t,s in generators) for key in TERMS},
          'per_term_frobenius_norm':{key:float(np.linalg.norm(terms[key])) for key in TERMS}}
    failures=[k for k,v in variants.items() if v['full_beta_relative_ward_residual']>TOL]
    return {'format':'t3-batch34-ward-gate-v1',
            'status':'BLOCKED_GAUGE_WARD_VIOLATION' if failures else 'PASSED',
            'source_revision':'12128de602ac6bf109dbad50dbac8446ed629434',
            'toy':'gauge-invariant SM L,eC,H Yukawa/quartic/Weinberg C, not physical A(-4)',
            'input_max_absolute_ward_residuals':ins,
            'variants':variants,'tolerance':TOL,'failed_variants':failures,
            'physical_Aminus4_certified':False,
            'explanation':'Crossed-Yukawa term/gauge index conventions remain unresolved; do not infer a correction merely from this test.'}

def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path)
    ap.add_argument('--diagnostic',action='store_true',help='Return zero after writing a diagnostic failure')
    a=ap.parse_args(argv); result=diagnose()
    if a.output:
        if a.output.exists(): raise FileExistsError(a.output)
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf8')
    print(result['status'])
    for name,v in result['variants'].items():
        print(f"{name}: relative Ward residual {v['full_beta_relative_ward_residual']:.9g}")
    return 0 if a.diagnostic or not result['failed_variants'] else 2
if __name__=='__main__': raise SystemExit(main())
