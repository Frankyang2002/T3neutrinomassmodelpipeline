"""Source-conditioned RGBeta quartic real-scalar action audit, NOT native certification.

The 6*Matrix prefactor and AvgPermutations convention are reported by RGBeta's
FieldsAndCouplings.m.  For four complex legs the component projection is
z=(x+iy)/sqrt(2).  Construct a *real-component* Lam tensor in a declared
Cartesian frame, symmetrize it exactly as an ordered quartic tensor, and check
(1/4!) Lam[x,x,x,x] against c*I(z)+conjugate.  The check is an independent
expansion of the same source invariant; it is NOT native RGBeta execution.
"""
from __future__ import annotations

import ast
from hashlib import sha256
import itertools
import json
from pathlib import Path
import re
from typing import Any
import numpy as np

from .basis import DIMENSIONS
from .matchete_rgbeta_dictionary import (
    DictionaryEvidenceError, MAT_FORMAT, RGB_CONDITIONAL_FORMAT,
    ALL16, calibrated_triplet_cg, _target_tensor, _mat_tensor,
    _species, _axis_permutation,
)

STATUS = "Conditional_RGBetaRealActionAlgebraPassed_NativeScalarBasisUnverified"


def _permutations(raw: str) -> list[tuple[int, ...]]:
    if not isinstance(raw, str) or len(raw) > 3000:
        raise DictionaryEvidenceError("Missing/oversized native permutations")
    try:
        x=ast.literal_eval(raw.replace("{", "[").replace("}", "]"))
    except (ValueError, SyntaxError) as exc:
        raise DictionaryEvidenceError("Cannot parse native permutations") from exc
    if (not isinstance(x,list) or not 1<=len(x)<=24 or
            any(type(p) is not list or len(p)!=4 or sorted(p)!=[1,2,3,4] for p in x) or
            len({tuple(p) for p in x})!=len(x)):
        raise DictionaryEvidenceError("Invalid/duplicate RGBeta UniqueArrangements")
    return [tuple(v-1 for v in p) for p in x]


def _projectors(species: list[str], dims: list[int]) -> tuple[list[np.ndarray],int]:
    """Four matrices carrying non-normalized (x+/-iy) real projections."""
    mapping={};offset=0
    for s,d in zip(species,dims):
        key=s[4:-1] if s.startswith("Bar[") else s
        if key in mapping:
            if mapping[key][1]!=d: raise DictionaryEvidenceError("Inconsistent species dimension")
        else:
            mapping[key]=(offset,d);offset+=2*d
    projectors=[]
    for s,d in zip(species,dims):
        key=s[4:-1] if s.startswith("Bar[") else s
        start,_=mapping[key]; P=np.zeros((offset,d),dtype=np.complex128)
        for k in range(d):
            P[start+2*k,k]=1
            P[start+2*k+1,k]=(-1j if s.startswith("Bar[") else 1j)
        projectors.append(P)
    return projectors,offset


def complex_field_from_real(species: list[str], dims: list[int],x:np.ndarray) -> dict[str,np.ndarray]:
    result={}; off=0
    for s,d in zip(species,dims):
        k=s[4:-1] if s.startswith("Bar[") else s
        if k not in result:
            result[k]=(x[off:off+2*d:2]+1j*x[off+1:off+2*d:2])/np.sqrt(2.)
            off+=2*d
    if off!=len(x): raise DictionaryEvidenceError("Invalid real-field dimension")
    return result


def monomial(tensor:np.ndarray,species:list[str],fields:dict[str,np.ndarray]) -> complex:
    legs=[]
    for s in species:
        bar=s.startswith("Bar["); key=s[4:-1] if bar else s
        v=fields[key]; legs.append(v.conj() if bar else v)
    return complex(np.einsum("abcd,a,b,c,d->",tensor,*legs,optimize=True))


def real_lam(tensor:np.ndarray,species:list[str],c:complex,perms:list[tuple[int,...]],
             *,self_conjugate:bool=False) -> np.ndarray:
    """Symmetric real 4-tensor, with RGBeta 6 and explicit complex H.c."""
    projectors,n=_projectors(species,list(tensor.shape))
    target=np.einsum("abcd,ia,jb,kc,ld->ijkl",tensor,*projectors,optimize=True)
    total=np.zeros((n,n,n,n),dtype=np.complex128)
    # The RGBeta AvgPermutations factor is 1/len(unique arrangements).
    # Every permutation represents the same commuting scalar polynomial.
    for p in perms:
        total+=np.transpose(target,p)/len(perms)
    raw=6*(c*total if self_conjugate else c*total+(c*total).conj())
    # The full symmetric quartic Lam action is unaltered by symmetrization;
    # explicitly symmetrize to check the commuting-variable convention.
    symmetric=sum(np.transpose(raw,p) for p in itertools.permutations(range(4)))/24
    if not np.allclose(symmetric.imag,0,atol=1e-12):
        raise DictionaryEvidenceError("RGBeta real-scalar Lam tensor not real")
    return symmetric.real


def baseline_invariants(d1:int,d2:int)->list[tuple[str,list[str],np.ndarray]]:
    """Six explicitly source-defined norm quartics; no higher portal CGs."""
    if any(type(d)!=int or d not in (1,2,3) for d in (d1,d2)):
        raise DictionaryEvidenceError('Unsupported scalar dimension')
    dims={'H':2,'S1':d1,'S2':d2}
    out=[]
    for key,first,second,factor in (
        ('lambdaH','H','H',.5),('lambdaS1','S1','S1',.5),
        ('lambdaS2','S2','S2',.5),('lambdaH1','H','S1',1.),
        ('lambdaH2','H','S2',1.),('lambda12','S1','S2',1.),
    ):
        a,b=dims[first],dims[second]
        tensor=factor*np.einsum('ij,kl->ijkl',np.eye(a),np.eye(b))
        species=[f'Bar[{first}]',first,f'Bar[{second}]',second]
        out.append((key,species,tensor))
    return out


def audit_baseline_all16(*,seed:int=1019,samples:int=3)->dict:
    rng=np.random.default_rng(seed)
    models={};total=0;worst=0.
    for tag in sorted(ALL16):
        d1,d2,_=DIMENSIONS[tag[0]]
        statuses={}
        for symbol,species,tensor in baseline_invariants(d1,d2):
            residual=0.
            for _ in range(samples):
                _,nr=_projectors(species,list(tensor.shape))
                x=rng.normal(size=nr); coeff=float(rng.uniform(-1.5,1.5))
                # SelfConjugate=True: RGBeta does not add a second H.c.
                lam=real_lam(tensor,species,coeff,[tuple(range(4))],self_conjugate=True)
                lhs=float(np.einsum('ijkl,i,j,k,l->',lam,x,x,x,x,optimize=True)/24.)
                fields=complex_field_from_real(species,list(tensor.shape),x)
                rhs=float((coeff*monomial(tensor,species,fields)).real)
                residual=max(residual,abs(lhs-rhs)/max(1.,abs(lhs),abs(rhs)))
            if residual>1e-10:raise DictionaryEvidenceError(f'{tag} {symbol}: baseline mismatch')
            statuses[symbol]={'relative_residual':residual,'self_conjugate':True,
                               'matched_explicit_source_norm_monomial':True,
                               'all_independent_portal_CGs_included':False}
            worst=max(worst,residual);total+=1
        models[tag]=statuses
    return {'model_count':16,'baseline_norm_coupling_slots':total,
            'max_baseline_norm_action_residual':worst,'models':models,
            'covered_couplings':['lambdaH','lambdaS1','lambdaS2','lambdaH1','lambdaH2','lambda12'],
            'NOT_covered':['lambdaH1Inv2','lambdaH2Inv2','lambda12Cross',
                           'other RG-closed independent quartics','A_-4 anti-Hermitian Higgs triplet portal']}


def audit(mdoc:dict,rdoc:dict,*,sha_m:str,sha_r:str,seed:int=3241,samples:int=3,tolerance:float=1e-9)->dict:
    if (mdoc.get("format")!=MAT_FORMAT or rdoc.get("format")!=RGB_CONDITIONAL_FORMAT or
            set(mdoc.get("models",{}))!=ALL16 or set(rdoc.get("models",{}))!=ALL16):
        raise DictionaryEvidenceError("Missing complete 16-model source evidence")
    if type(seed)!=int or not 1<=samples<=20 or not 0<tolerance<=1e-4:
        raise DictionaryEvidenceError("Invalid audit settings")
    cg,_=calibrated_triplet_cg(mdoc)
    rng=np.random.default_rng(seed)
    overall=0.;count=0;models={}
    for tag in sorted(ALL16):
        mb=mdoc['models'][tag]['tensors']; rb=rdoc['models'][tag]['records']
        if len(mb)!=len(rb):raise DictionaryEvidenceError(f"{tag}: inconsistent inventories")
        worst=0.;audit_records=[]
        for rec in rb:
            name=rec['symbol'];species=rec['species']; dims=rec['dimensions']
            tensor=_target_tensor(rec,cg)
            if tuple(tensor.shape)!=tuple(dims):raise DictionaryEvidenceError("Invalid tensor shape")
            perms=_permutations(rec['native_unique_arrangements'])
            if not re.fullmatch(r"6\s*\*\s*Matrix\["+re.escape(name)+r"\]\[##1\]\s*&\s*",
                                rec['native_coupling_function']):
                raise DictionaryEvidenceError("Native RGBeta complex quartic coefficient !=6")
            if rec['native_self_conjugate'] not in (False,'False'):
                raise DictionaryEvidenceError("Self-conjugate quartic not in extra inventory")
            residual=0.
            for _ in range(samples):
                _,nr=_projectors(species,dims)
                x=rng.normal(size=nr)
                c=complex(*rng.normal(size=2))
                lam=real_lam(tensor,species,c,perms)
                left=np.einsum('ijkl,i,j,k,l->',lam,x,x,x,x,optimize=True)/24
                fields=complex_field_from_real(species,dims,x)
                right=2*(c*monomial(tensor,species,fields)).real
                rel=float(abs(left-right)/max(1.,abs(left),abs(right)))
                residual=max(residual,rel)
            worst=max(worst,residual)
            audit_records.append({'symbol':name,'unique_arrangements':len(perms),
                                  'max_real_action_relative_residual':residual,
                                  'source_quartic_coefficient':6,
                                  'source_action_denominator':24,
                                  'complex_leg_projector':'(x_re +/- I*x_im)/sqrt(2)',
                                  'physical_normalization_certified':False})
            count+=1
        models[tag]={'n_extra_quartics':len(rb),'max_real_action_relative_residual':worst,
                     'records':audit_records}
        overall=max(overall,worst)
    if count!=50 or overall>=tolerance:
        raise DictionaryEvidenceError(f"Not all 50 pass real quartic-action check: count={count}, residual={overall}")
    baseline=audit_baseline_all16(seed=seed+1,samples=samples)
    return {'format':'t3-parallel-a-conditional-real-quartic-action-audit-v4',
            'status':STATUS,'models':models,'model_count':16,'extra_quartic_slots':50,
            'max_real_action_relative_residual':overall,'baseline_norm_audit':baseline,
            'source_evidence_sha256':{'matchete':sha_m,'rgbeta':sha_r},
            'source_RGBeta_quartic_prefactor':6,
            'source_RGBeta_quartic_action_convention':'Lam_abcd phi_a phi_b phi_c phi_d / 4!',
            'scope':{'numerical_check':'independent real-component action contraction in chosen Matchete-aligned frame',
                     'permutation_averaging_verified_algebraically':True,
                     'native_RGBeta_real_scalar_projectors_executed':False,
                     'native_RGBeta_scalar_coordinate_basis_certified':False,
                     'six_baseline_norm_monomials_audited':True,
                     'baseline_quartic_tensor_dictionary_certified':False,
                     'A_minus4_anti_Hermitian_portal_phase_certified_by_this_run':False,
                     'physical_scalar_dictionary_certified':False,
                     'physical_wilson_normalization_certified':False}}


def audit_from_paths(mat:Path,rgb:Path,**options)->dict:
    mb=mat.read_bytes();rb=rgb.read_bytes()
    return audit(json.loads(mb.decode('utf-8-sig')),json.loads(rb.decode('utf-8-sig')),
                 sha_m=sha256(mb).hexdigest(),sha_r=sha256(rb).hexdigest(),**options)
