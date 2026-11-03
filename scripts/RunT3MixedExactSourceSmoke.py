"""Read-only source-exact B(-1)/C(-1) Majorana F-tree CG/flavour smoke check.

Checks the archived reference native CG reconstruction, then evaluates every
PL source with independent test 3-family y1/y2 and a nontrivial unitary X/Y
scalar mixing matrix. No RGE integration, Matchete rerun or saved output edits.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

from Numerical.running.GeneralT3MixedExactFermionWilson import (
    MIXED, read_native, archived_reference_reproduction, project_source,
)

ROOT=Path(__file__).resolve().parents[1]


def generic_yukawas():
    y1=1e-4*np.array([[1,2,3],[-1,2,-2],[2,1,-1]],complex)
    y2=1e-4*np.array([[2,1j,1],[3,1,-2j],[1,3,1]],complex)
    return {'y1':y1,'y2':y2}


def run(models):
    theta=0.24
    U=np.array([[math.cos(theta),math.sin(theta)],
                [-math.sin(theta),math.cos(theta)]],complex)
    mf=np.array([1.01e5,1.13e5,1.27e5])
    for model in models:
        source,pl,cg,ordered=read_native(model,ROOT)
        archived_res=archived_reference_reproduction(model,pl,cg,ordered)
        by,indices,details,audit=project_source(model,source,pl,cg,mf,generic_yukawas(),U)
        assert by.shape==(len(pl),3,3,2,2,8,8) and np.isfinite(by).all()
        assert len(set(indices))==len(pl)
        classes={tuple(v['yukawa_pair']) for v in details}
        assert {('y1','y1'),('y1','y2'),('y2','y2')} <= classes
        bridge={str(k):float(np.linalg.norm(v)) for k,v in zip(indices,by)
                if model=='C_-1' and k in (23,24)}
        print(f'{model} Verified native PL={len(indices)} CG reference residual='
              f'{archived_res:.3e} dynamic LLSS norm={np.linalg.norm(by.sum(axis=0)):.7g}'
              f' bridge_norms={bridge}',flush=True)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--models',default='B_-1,C_-1')
    args=p.parse_args(argv)
    models=tuple(x.strip() for x in args.models.split(',') if x.strip())
    if not models or len(set(models))!=len(models) or not set(models)<=set(MIXED):
        p.error('models must be distinct B_-1 and/or C_-1')
    run(models)
    return 0

if __name__=='__main__':
    raise SystemExit(main())
