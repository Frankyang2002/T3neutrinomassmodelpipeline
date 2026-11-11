"""Batch 16: source-convention candidate real scalar mass Hessians at muS.

Project path: Numerical/experimental/Batch16ScalarMassHessians.py

Uses the Batch15 saved scalar-threshold audit. The mass eigenvalue spectra
follow the explicitly declared quadratic potential with RGBeta-native
mass-parameter normalisation; cross-package scalar CG and mass-basis maps
are NOT certified. Neither the Wilson boundary nor a physical C5 is computed.

For complex scalar components S=(R+i I)/sqrt(2), with real coordinates
q=(R_1,I_1,...), construct Hessian H of the quadratic potential V=1/2 q^T H q.

Potential assumptions:
  V_diag = m1^2 S1†S1 + m2^2 S2†S2
  V_self = sign [mu (S_k*)^T (S_k*) + h.c.]
  V_mixed= sign [mu (S1† epsilon S2*) + h.c.]
  epsilon = [[0,1],[-1,0]]
sign=+1 for the literal RGBeta V2 convention in the export metadata;
sign=-1 for the source Matchete L_int=+mu*invariant+h.c. convention.

The sign choice changes the eigenvectors/field map even when the spectrum
is unchanged. The triplet Cartesian identity metric is a *candidate*
component interpretation; never identify its eigenvectors with Matchete
source tensors until the SU(2) invariant map is certified.
"""
from __future__ import annotations

import json
import math
from typing import Any, Mapping
import numpy as np

DIMS = {"A": (1,3), "B": (2,2), "C": (2,2), "D": (3,1), "E": (3,3)}
SPECIAL = {
    "A_+0": "self_s1", "A_-2": "self_s2",
    "B_-1": "mixed", "C_-1": "mixed",
    "D_+0": "self_s1", "D_-2": "self_s2",
    "E_+0": "self_s1", "E_-2": "self_s2",
}
MASS_KEY = {
    "self_s1": "muT3GeneralS1Sq",
    "self_s2": "muT3GeneralS2Sq",
    "mixed": "muT3GeneralS1S2",
}
CONVENTIONS = {
    "rgbeta_potential_plus": +1,
    "matchete_lint_plus": -1,
}
EPS = np.array([[0.,1.],[-1.,0.]])


def _complex(raw: Any, name: str) -> complex:
    if not isinstance(raw, Mapping) or set(raw) != {"real", "imag"}:
        raise ValueError(f"{name}: expected explicit real/imag object")
    z = complex(float(raw["real"]), float(raw["imag"]))
    if not (np.isfinite(z.real) and np.isfinite(z.imag)):
        raise ValueError(f"{name}: nonfinite complex parameter")
    return z


def _model_info(report: Mapping[str, Any]) -> tuple[str,int,int,str,float,float,complex]:
    if report.get("status") != "CandidateScalarThresholdHandoffAudited":
        raise ValueError("Requires successful Batch15 scalar-threshold audit")
    if report.get("physical_neutrino_mass_prediction") is not False:
        raise ValueError("Input unexpectedly claims a physical neutrino-mass prediction")
    model = report.get("model")
    if not isinstance(model,str) or model not in {
        f"{fam}_{a}" for fam, alphas in {
            "A": ("+0","-2","-4"), "B": ("+1","-1","-3"),
            "C": ("+1","-1","-3"), "D": ("+0","+2","-2"),
            "E": ("+0","+2","-2","-4"),
        }.items() for a in alphas
    }:
        raise ValueError(f"Unknown T3 model {model!r}")
    d1,d2 = DIMS[model[0]]
    audit = report["scalar_threshold"]
    masses = audit["diagonal_mass_squared_gev2"]
    m1,m2 = float(masses["mS1Sq"]),float(masses["mS2Sq"])
    if not all(math.isfinite(x) and x>0 for x in (m1,m2)):
        raise ValueError("Diagonal mass coefficients must be finite and positive")
    kind = SPECIAL.get(model,"none")
    symbols = audit["additional_quadratic_symbols"]
    vals = audit["additional_quadratic_values"]
    expected = [] if kind=="none" else [MASS_KEY[kind]]
    if sorted(symbols)!=expected or set(vals)!=set(expected):
        raise ValueError(f"{model}: quadratic inventory mismatch: expected {expected}, got {symbols}")
    mu = 0j if kind=="none" else _complex(vals[expected[0]],expected[0])
    return model,d1,d2,kind,m1,m2,mu


def _complex_components(q: np.ndarray, d1:int,d2:int) -> tuple[np.ndarray,np.ndarray]:
    if q.shape != (2*(d1+d2),):
        raise ValueError("Incorrect real scalar component vector shape")
    z = (q[0::2]+1j*q[1::2])/math.sqrt(2.)
    return z[:d1],z[d1:]


def _quadratic_potential(q:np.ndarray, d1:int,d2:int,kind:str,
                         m1:float,m2:float,mu:complex,sign:int) -> float:
    s1,s2 = _complex_components(q,d1,d2)
    result = m1*np.vdot(s1,s1).real + m2*np.vdot(s2,s2).real
    if kind == "self_s1":
        result += sign*2*np.real(mu*np.dot(s1.conj(),s1.conj()))
    elif kind == "self_s2":
        result += sign*2*np.real(mu*np.dot(s2.conj(),s2.conj()))
    elif kind == "mixed":
        result += sign*2*np.real(mu*np.vdot(s1,EPS@s2.conj()))
    elif kind != "none":
        raise ValueError(f"Unknown quadratic kind {kind}")
    return float(result)


def real_mass_hessian(d1:int,d2:int,kind:str,m1:float,m2:float,
                      mu:complex,*,convention:str) -> np.ndarray:
    """Exact polarization of the declared quadratic potential; no finite differences."""
    if convention not in CONVENTIONS:
        raise ValueError(f"Unknown quadratic convention {convention!r}")
    sign=CONVENTIONS[convention]
    n=2*(d1+d2)
    basis=np.eye(n)
    pot=lambda q: _quadratic_potential(q,d1,d2,kind,m1,m2,mu,sign)
    vdiag=np.array([2*pot(e) for e in basis])
    result=np.diag(vdiag)
    for i in range(n):
        for j in range(i+1,n):
            hij=pot(basis[i]+basis[j])-pot(basis[i])-pot(basis[j])
            result[i,j]=result[j,i]=hij
    return result


def analytic_mass_spectrum(d1:int,d2:int,kind:str,m1:float,m2:float,mu:complex) -> np.ndarray:
    if kind=="none":
        return np.sort([m1]*(2*d1)+[m2]*(2*d2))
    if kind in ("self_s1","self_s2"):
        d = d1 if kind=="self_s1" else d2
        m = m1 if kind=="self_s1" else m2
        other_d = d2 if kind=="self_s1" else d1
        other_m = m2 if kind=="self_s1" else m1
        return np.sort([m-2*abs(mu)]*d+[m+2*abs(mu)]*d+[other_m]*(2*other_d))
    if kind=="mixed":
        avg=0.5*(m1+m2)
        gap=math.hypot(0.5*(m1-m2),abs(mu))
        return np.sort([avg-gap]*4+[avg+gap]*4)
    raise ValueError(f"Unknown quadratic kind: {kind}")


def audit_mass_spectrum(report:Mapping[str,Any]) -> dict[str,Any]:
    model,d1,d2,kind,m1,m2,mu = _model_info(report)
    analytic=analytic_mass_spectrum(d1,d2,kind,m1,m2,mu)
    result={}
    for conv in CONVENTIONS:
        h=real_mass_hessian(d1,d2,kind,m1,m2,mu,convention=conv)
        vals,vecs=np.linalg.eigh(h)
        norm=max(abs(m1),abs(m2),abs(mu),1.)
        if not np.allclose(vals,analytic,atol=1e-11*norm,rtol=1e-11):
            raise RuntimeError(f"{model}: Hessian and analytic spectrum disagree for {conv}")
        residual=np.max(np.abs(h@vecs-vecs*vals[None,:]))/norm
        if residual>1e-10 or not np.allclose(vecs.T@vecs,np.eye(h.shape[0]),atol=1e-10):
            raise RuntimeError(f"{model}: mass basis decomposition failed")
        if not np.isclose(np.trace(h),2*d1*m1+2*d2*m2,rtol=1e-11):
            raise RuntimeError(f"{model}: mass-squared trace mismatch")
        mu_s=float(report["mu_s_gev"])
        mu_f=float(report["mu_f_gev"])
        if not 0<mu_s<mu_f:
            raise ValueError("Invalid simultaneous threshold scales")
        positive=bool(np.min(vals)>0)
        hierarchy=bool(positive and np.sqrt(vals[0])>mu_s
                       and np.sqrt(vals[-1])<mu_f)
        result[conv]={
            "mass_squared_eigenvalues_gev2":vals.tolist(),
            "masses_gev":[float(np.sqrt(x)) if x>0 else None for x in vals],
            "tachyon_free_candidate":positive,
            "eigenmass_threshold_hierarchy_candidate":hierarchy,
            "diagonalization_max_relative_residual":float(residual),
            "real_mass_hessian_gev2":h.tolist(),
            # Real eigenvectors are convention-dependent and cannot be silently
            # transferred to Matchete/LLSS fields.
            "real_eigenvector_columns_candidate":vecs.tolist(),
        }
    max_shift=float(np.max(np.abs(analytic-
        np.sort([m1]*(2*d1)+[m2]*(2*d2)))))
    return {
        "format":"t3-batch16-scalar-mass-convention-audit-v1",
        "model":model,
        "status":"CandidateNativeQuadraticHessianComputed",
        "quadratic_kind":kind,
        "quadratic_parameter":MASS_KEY.get(kind),
        "quadratic_value":{"real":mu.real,"imag":mu.imag} if kind!="none" else None,
        "scalar_component_count_real":2*(d1+d2),
        "spectra_independent_of_bilinear_sign":True,
        "maximum_shift_from_diagonal_proxy_gev2":max_shift,
        "conventions":result,
        "scalar_CG_basis_mapping_certified":False,
        "matchete_rgbeta_mass_parameter_mapping_certified":False,
        "mass_eigenvectors_authoritative_for_wilson_matching":False,
        "simultaneous_scalar_wilson_matching_performed":False,
        "physical_neutrino_mass_prediction":False,
        "note":"Declared native scalar quadratic convention only; no certified Matchete-to-RGBeta mass/CG map.",
    }
