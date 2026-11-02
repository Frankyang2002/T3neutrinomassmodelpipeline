"""Synthetic, NON-PHYSICAL tensor fixtures for all 16 model dimensions.

Their amplitudes are deliberately engineered to project onto the Weinberg
tensor, so these tests validate contractions, mixing, degeneracy, phases,
finite scale dependence and serialization, NOT gauge-invariant T3 couplings.
"""
from __future__ import annotations
import math
import numpy as np
from .core import b0bar,weinberg_tensor
from .adapter import InputSet,normalize_model,scalar_groups,MIXED,SELF_SPLIT,ALL_MODELS


def make_fixture(model: str) -> InputSet:
    model=normalize_model(model)
    groups=scalar_groups(model)
    n=len(groups); n1=groups.count("S1")
    mu_s=1000.0
    masses = np.array([925.+31*i for i in range(n1)]+[1175.+27*i for i in range(n-n1)])
    D=np.diag(masses*masses)
    # Engineer six independent *real* Higgs-quartic tensor structures.
    # For an internal pair (a,b), C carries an SU2 weak tensor while q carries
    # a symmetric real Higgs tensor; their product reconstructs O5.
    T=weinberg_tensor()
    flavor=np.array([[1.+.37j, .2+.08j, -.11+.04j],
                     [.2+.08j, .8-.2j, .13-.1j],
                     [-.11+.04j,.13-.1j,.6+.24j]],complex)*1.0e-10 * (1+0.03*ALL_MODELS.index(model)) * np.exp(1j*.04*ALL_MODELS.index(model))
    C=np.zeros((3,3,2,2,n,n),complex)
    Q=np.zeros((4,4,n,n))
    Qchannels=np.zeros((6,4,4,n,n))
    channels=[]
    for base in (0,n1):
        channels.extend([(base,base),(base+1,base+1),(base,base+1)])
    for k,(a,b) in enumerate(channels):
        weak=(0,0) if k < 2 else (0,1) if k < 4 else (1,1)
        coeff = (T[weak].real if k%2==0 else T[weak].imag)
        weak_matrix=np.zeros((2,2),complex)
        w=1.0 if k%2==0 else 1j
        weak_matrix[weak]=w
        weak_matrix[weak[::-1]]=w
        kernel=b0bar(D[a,a],D[b,b],mu_s)
        if abs(kernel)<1e-8:raise RuntimeError("Fixture bubble unexpectedly zero")
        # Tuned solely to make non-physical toy data a pure Weinberg tensor.
        C[...,a,b]=np.einsum("fg,ij->fgij",flavor,weak_matrix)
        C[...,b,a]=C[...,a,b]
        q=coeff / ((2 if a!=b else 1)*kernel)
        Q[...,a,b]=q
        Q[...,b,a]=q
        Qchannels[k,...,a,b]=q
        Qchannels[k,...,b,a]=q
    # Real orthogonal mass mixing of S1 and S2 (B/C), including all tensors.
    O=np.eye(n)
    if model in MIXED:
        angle=0.23
        i,j=0,n1
        O[i,i]=O[j,j]=math.cos(angle)
        O[i,j]=math.sin(angle);O[j,i]=-math.sin(angle)
        i,j=1,n1+1
        O[i,i]=O[j,j]=math.cos(angle)
        O[i,j]=-math.sin(angle);O[j,i]=math.sin(angle)
    mass2=O@D@O.T
    C=np.einsum("fgijab,ka,lb->fgijkl",C,O,O,optimize=True)
    Q=np.einsum("ABab,ka,lb->ABkl",Q,O,O,optimize=True)
    Qchannels=np.einsum("tABab,ka,lb->tABkl",Qchannels,O,O,optimize=True)
    C5_direct=1.5e-14*flavor/1e-10
    return InputSet(model,5000.,mu_s,1100.,mass2,C,Q,C5_direct,
                    {"source":"synthetic algebra-only fixture (not Chat 1, not physical T3)",
                     "rank_checks":"engineered Weinberg external tensor"},True,
                    "synthetic_algebra_only",Qchannels,
                    tuple(f"synthetic_higgs_tensor_channel_{i}" for i in range(6)))
