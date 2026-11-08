"""Standalone 45-Weyl, 12-real-scalar SYNTHETIC scaling check.

Only dimensions match three-generation SM+A(-4). Artificial couplings have
no SU(3)xSU(2)xU(1) selection-rule guarantee. NO physical SMEFT C5.
"""
from __future__ import annotations

from itertools import permutations
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))
from RGE.running.intermediate.GeneralEFT1SevenTermFlow import (
    Couplings,TensorBasis,evolve,save_results)


def symmetric_set(c,i,j,a,b,v):
    for f,g in ((i,j),(j,i)):
        for x,y in ((a,b),(b,a)):
            c[f,g,x,y]=v


def example():
    nf,ns=45,12
    basis=TensorBasis(nf,ns,((0,1),(2,3),(4,5)),(0,1,2,3))
    c=np.zeros((nf,nf,ns,ns),complex)
    symmetric_set(c,0,0,4,6,2e-11+0.7e-11j)
    symmetric_set(c,0,2,4,6,1e-11-0.4e-11j)
    symmetric_set(c,2,4,5,7,1.4e-11+0.3e-11j)
    y=np.zeros((nf,nf,ns),complex)
    for i,j,a,v in ((0,0,4,.13+.07j),(44,44,4,.11+.01j),
                    (0,1,2,.04+.02j),(2,2,5,.06+.03j),
                    (2,4,7,.03-.02j)):
        y[i,j,a]=y[j,i,a]=v
    lam=np.zeros((ns,)*4)
    for key,v in (((2,2,4,6),.07),((3,3,4,6),-.07),
                  ((2,2,5,7),-.03),((0,0,4,6),.07),
                  ((4,4,6,6),.11)):
        for q in set(permutations(key)):
            lam[q]=v
    low=Couplings(y*.82,lam*.90,())
    high=Couplings(y,lam,())
    return basis,c,[2000.,100000.],[low,high]


def main():
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-root',type=Path,required=True)
    p.add_argument('--beta-mode',choices=('source_v1_unstarred','ward_covariant_candidate','fonseca_eq22_six_term'),
                   default='source_v1_unstarred')
    a=p.parse_args()
    b,c,t,s=example()
    out=evolve(b,c,100000.,2000.,t,s,points=3,atol=1e-20,rtol=1e-7,beta_mode=a.beta_mode)
    out['authority']='synthetic'
    out['source_ids']=['SYN_SM_45_large_closure']
    out['fermion_field_blocks']=[
        {'name':f'SYN_L{q}','indices':[2*q-1,2*q]} for q in (1,2,3)] + [
        {'name':'SYN_39_unlabelled_SM_spectators','indices':list(range(7,46))}]
    out['scalar_field_blocks']=[{'name':'H','role':'Higgs','indices':[1,2,3,4]},
         {'name':'S1','role':'heavy','indices':[5,6]},
         {'name':'S2','role':'heavy','indices':list(range(7,13))}]
    root=a.output_root/'A_m4_full45_synthetic'
    save_results(out,str(root/'seven_term_trajectory.npz'),str(root/'seven_term_summary.json'))
    print('beta mode:',a.beta_mode,'(unverified); nf=45,ns=12,independent=',b.independent_size,'generated=',len(out['generated_indices']),
          'nfev=',out['solver']['nfev'])

if __name__=='__main__':
    main()
