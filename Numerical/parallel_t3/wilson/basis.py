"""Canonical left Weyl + interleaved real scalar RGBeta reference coordinates.

Reproduces label/generator conventions inspected in Aminus4EFT1Couplings.py,
SHA-1 GitHub blob e1797b35446e7d98cbf18ce2b47b4bc86c28603d.
The representation is a *reference*, not a certified Matchete change of basis.
"""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np

DIMENSIONS = {"A": (1, 3, 2), "B": (2, 2, 1), "C": (2, 2, 3),
              "D": (3, 1, 2), "E": (3, 3, 2)}
ALPHAS = {"A": (-4, -2, 0), "B": (-3, -1, 1),
          "C": (-3, -1, 1), "D": (-2, 0, 2), "E": (-4, -2, 0, 2)}
ALL_MODELS = tuple(f"{k}_{a:+d}" for k, aa in ALPHAS.items() for a in aa)


def spin_generators(n: int) -> tuple[np.ndarray, ...]:
    s = (n - 1) / 2
    m = [s - i for i in range(n)]
    raising = np.zeros((n, n), complex)
    for j, mm in enumerate(m):
        if mm+1 in m:
            raising[m.index(mm+1), j] = math.sqrt((s-mm)*(s+mm+1))
    lowering = raising.conj().T
    return ((raising+lowering)/2, (raising-lowering)/(2j),
            np.diag(m).astype(complex))


def su3_generators() -> tuple[np.ndarray, ...]:
    s = math.sqrt(3)
    return tuple(np.asarray(g, complex)/2 for g in (
        [[0,1,0],[1,0,0],[0,0,0]], [[0,-1j,0],[1j,0,0],[0,0,0]],
        [[1,0,0],[0,-1,0],[0,0,0]], [[0,0,1],[0,0,0],[1,0,0]],
        [[0,0,-1j],[0,0,0],[1j,0,0]], [[0,0,0],[0,0,1],[0,1,0]],
        [[0,0,0],[0,0,-1j],[0,1j,0]], [[1/s,0,0],[0,1/s,0],[0,0,-2/s]]))


def real_scalar_generator(t: np.ndarray) -> np.ndarray:
    """Interleaved (R,I) generator on real scalars; Hermitian, imaginary/antisymmetric."""
    n = t.shape[0]
    out = np.zeros((2*n, 2*n), complex)
    for i in range(n):
        for j in range(n):
            a, b = t[i,j].real, t[i,j].imag
            out[2*i:2*i+2, 2*j:2*j+2] = 1j*np.array([[b,a],[-a,b]])
    return out


def _fermions() -> tuple[str, ...]:
    labels = []
    for sp in ('Q', 'uc', 'dc', 'L', 'ec'):
        for p in (1,2,3):
            if sp == 'Q':
                labels.extend(f'Q[{p},{c},{w}]' for c in (1,2,3) for w in (1,2))
            elif sp in ('uc','dc'):
                labels.extend(f'{sp}[{p},{c}]' for c in (1,2,3))
            elif sp == 'L':
                labels.extend(f'L[{p},{w}]' for w in (1,2))
            else:
                labels.append(f'ec[{p}]')
    return tuple(labels)

F_LABELS = _fermions()
assert len(F_LABELS) == 45

@dataclass(frozen=True)
class ModelBasis:
    model: str
    dims: tuple[int, int, int]
    alpha: int
    hypercharges: tuple[float, float, float]
    fermion_labels: tuple[str, ...]
    scalar_labels: tuple[str, ...]

    @property
    def nf(self): return len(self.fermion_labels)
    @property
    def ns(self): return len(self.scalar_labels)
    @property
    def leptons(self):
        lookup = {x:i for i,x in enumerate(self.fermion_labels)}
        return tuple((lookup[f'L[{p},1]'],lookup[f'L[{p},2]']) for p in (1,2,3))
    @property
    def higgs(self): return (0,1,2,3)

    def to_json(self) -> dict:
        return {"basis_id": "T3_RGBeta_reference_Weyl45_real_interleaved_v1",
                "model": self.model, "dims_S1_S2_F": list(self.dims),
                "alpha":self.alpha, "hypercharges_S1_S2_F": list(self.hypercharges),
                "fermion_labels": list(self.fermion_labels),
                "scalar_labels": list(self.scalar_labels),
                "fermion_axis": ["i","j"],"scalar_axis": ["a","b"],
                "wilson_axes": ["weyl_i","weyl_j","real_scalar_a","real_scalar_b"],
                "all_fermions_are_left_handed_Weyl":True,
                "fermion_reps": {"Q": [3,2,1/6],"uc":["anti3",1,-2/3],
                                  "dc":["anti3",1,1/3],"L":[1,2,-.5],"ec":[1,1,1]},
                "heavy_scalars_Z2":-1, "SM_fields_Z2":1,
                "complex_scalar_definition":"z=(R+iI)/sqrt(2)",
                "H_index_weak1":"charged", "H_index_weak2":"neutral",
                "all_heavy_scalars_present_in_EFT1":True,
                "Matchete_RGBeta_change_of_basis_certified":False}

    def generators(self) -> dict[str, tuple[tuple[np.ndarray,np.ndarray], ...]]:
        fi = {name:i for i,name in enumerate(self.fermion_labels)}
        ff = {'SU2': np.zeros((3,self.nf,self.nf),complex),
              'U1':np.zeros((1,self.nf,self.nf),complex),
              'SU3':np.zeros((8,self.nf,self.nf),complex)}
        q = {'Q':1/6,'uc':-2/3,'dc':1/3,'L':-.5,'ec':1}
        weak2, color3 = spin_generators(2), su3_generators()
        for sp in ('Q','uc','dc','L','ec'):
            for p in (1,2,3):
                if sp in ('Q','L'):
                    for color in ((1,2,3) if sp=='Q' else (None,)):
                        idx = [fi[f'Q[{p},{color},{w}]' if sp=='Q' else f'L[{p},{w}]'] for w in (1,2)]
                        for k in range(3): ff['SU2'][k][np.ix_(idx,idx)] = weak2[k]
                if sp in ('Q','uc','dc'):
                    for w in ((1,2) if sp=='Q' else (None,)):
                        idx = [fi[f'Q[{p},{col},{w}]' if sp=='Q' else f'{sp}[{p},{col}]'] for col in (1,2,3)]
                        for k in range(8):
                            ff['SU3'][k][np.ix_(idx,idx)] = color3[k] if sp=='Q' else -color3[k].conj()
        for i, lab in enumerate(self.fermion_labels):
            ff['U1'][0,i,i] = q[lab.split('[')[0]]
        ss = {'SU2':np.zeros((3,self.ns,self.ns),complex),
              'U1':np.zeros((1,self.ns,self.ns),complex),
              'SU3':np.zeros((8,self.ns,self.ns),complex)}
        off = 0
        for d,y in ((2,.5),(self.dims[0],self.hypercharges[0]),
                    (self.dims[1],self.hypercharges[1])):
            ids = list(range(off,off+2*d)); off += 2*d
            for k,t in enumerate(spin_generators(d)):
                ss['SU2'][k][np.ix_(ids,ids)] = real_scalar_generator(t)
            ss['U1'][0][np.ix_(ids,ids)] = real_scalar_generator(y*np.eye(d))
        return {group:tuple((ff[group][k],ss[group][k]) for k in range(len(ff[group])))
                for group in ('SU2','U1','SU3')}


def build_basis(model: str) -> ModelBasis:
    try:
        k, a = model.split('_', 1); alpha = int(a)
    except (ValueError,AttributeError) as e:
        raise ValueError('Expected model A_-4, B_-1, etc') from e
    if k not in ALPHAS or alpha not in ALPHAS[k]:
        raise ValueError(f'Unsupported ordinary T3 model {model}')
    dims = DIMENSIONS[k]
    species = (('H',2),('S1',dims[0]),('S2',dims[1]))
    scalars = tuple(f'{name}[{m},{ri}]' for name,d in species
                    for m in range(1,d+1) for ri in ('R','I'))
    return ModelBasis(model,dims,alpha,(alpha/2,(alpha+2)/2,(alpha+1)/2),
                      F_LABELS,scalars)
