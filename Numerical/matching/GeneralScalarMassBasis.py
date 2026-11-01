"""Physical real-scalar basis of a general complex quadratic potential.

V2 = phi^dagger A phi + (phi^T B phi + h.c.)/2,
A=A^dagger, B=B^T.  phi=(x+i y)/sqrt(2).

B is a *quadratic-form coefficient* and NOT automatically an RGBeta
muT3General... coupling.  Translate the project's explicit SU(2) CG and
PlusHc normalization before constructing A and B.

The result defines phi = (T_x + i T_y) R / sqrt(2), where R contains
canonically normalized real mass eigenfields, and V2 = R^T diag(m2) R / 2.
This handles same-field complex bilinears and cross-field mixing together.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from typing import Any
import numpy as np


def _read_complex(value: Any) -> complex:
    if isinstance(value, dict):
        if set(value) != {'real','imag'}:
            raise ValueError('Complex objects must have only real, imag')
        return complex(float(value['real']),float(value['imag']))
    return complex(value)


def _matrix(raw: Any, name: str) -> np.ndarray:
    if not isinstance(raw,list) or not raw or any(not isinstance(row,list) for row in raw):
        raise ValueError(name+' must be a square nested list')
    n=len(raw)
    if any(len(row)!=n for row in raw):
        raise ValueError(name+' must be square')
    result=np.array([[_read_complex(v) for v in row] for row in raw],dtype=complex)
    if not np.all(np.isfinite(result.real)) or not np.all(np.isfinite(result.imag)):
        raise ValueError(name+' contains nonfinite entries')
    return result


def real_quadratic_hessian(A: np.ndarray,B: np.ndarray) -> np.ndarray:
    """Return real Hessian in ordering (Re phi_1..n, Im phi_1..n)."""
    A=np.asarray(A,complex)
    B=np.asarray(B,complex)
    if A.ndim!=2 or A.shape[0]!=A.shape[1] or B.shape!=A.shape:
        raise ValueError('A and B must be square with equal dimensions')
    if not np.allclose(A,A.conj().T,rtol=1e-11,atol=1e-11):
        raise ValueError('A must be Hermitian')
    if not np.allclose(B,B.T,rtol=1e-11,atol=1e-11):
        raise ValueError('B must be symmetric')
    # xTx(A_R+B_R)x/2 + yTy(A_R-B_R)y/2
    #  - xT(A_I+B_I)y, where A_I is antisymmetric.
    top=np.column_stack((A.real+B.real, -A.imag-B.imag))
    bottom=np.column_stack((A.imag-B.imag, A.real-B.real))
    H=np.concatenate((top,bottom),axis=0)
    return (H+H.T)/2


def mass_basis(A: np.ndarray,B: np.ndarray, *,tachyon_tolerance:float=1e-10) -> dict[str,Any]:
    H=real_quadratic_hessian(A,B)
    eigenvalues,O=np.linalg.eigh(H)
    scale=max(np.linalg.norm(H,ord=2),1.)
    if np.min(eigenvalues)<-tachyon_tolerance*scale:
        raise ValueError('Tachyonic quadratic eigenvalue: '+str(np.min(eigenvalues)))
    if np.any(eigenvalues < 0):
        eigenvalues=np.maximum(eigenvalues,0.)
    n=A.shape[0]
    U=(O[:n,:]+1j*O[n:,:])/np.sqrt(2.)
    check=np.linalg.norm(O.T@H@O-np.diag(eigenvalues))/max(np.linalg.norm(H),1.)
    return {'mass_squared_gev2':eigenvalues.tolist(),
            'masses_gev':np.sqrt(eigenvalues).tolist(),
            'real_rotation':O.tolist(),
            'complex_field_from_real_eigenstate':_encode(U),
            'residual':float(check),
            'field_convention':'phi = U @ real_eigenstates; V2 = 0.5 * sum(m2_i * real_i**2)',
            'matching_status':'PhysicalQuadraticBasisOnly_NoMatcheteMatching'}


def _encode(obj):
    if isinstance(obj,np.ndarray): return _encode(obj.tolist())
    if isinstance(obj,list): return [_encode(x) for x in obj]
    if isinstance(obj,complex):return {'real':float(obj.real),'imag':float(obj.imag)}
    return obj


def run(args=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--quadratic-json',type=Path,required=True,
                   help='JSON with Hermitian A and symmetric B, in GeV^2, V=phi^dag A phi + (phi^T B phi+h.c.)/2')
    p.add_argument('--output-json',type=Path,required=True)
    a=p.parse_args(args)
    data=json.loads(a.quadratic_json.read_text(encoding='utf-8-sig'))
    A=_matrix(data['A'],'A'); B=_matrix(data['B'],'B')
    if A.shape!=B.shape: raise ValueError('A and B dimensions differ')
    result=mass_basis(A,B)
    result['input_fields']=data.get('fields',[])
    result['quadratic_input']=str(a.quadratic_json)
    a.output_json.parent.mkdir(parents=True,exist_ok=True)
    a.output_json.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Diagonalized',2*len(A),'real scalar modes ->',a.output_json)
    return 0

if __name__=='__main__':raise SystemExit(run())
