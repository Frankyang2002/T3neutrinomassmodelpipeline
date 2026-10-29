"""T3 A/D/E self-bilinears: tensor-checked anomalous one-loop spectral term.

At one quartic insertion, the two-Higgs extra operator is the conjugate of
an independent Stage22 Hdag Hdag S S vertex. In the native Matchete triplet
basis the self-scalar invariant metric is delta_ab (not the declared RGBeta
spherical metric). The rotated anomalous scalar propagator contributes

  K_A(p) = (mu2/abs(mu2))**eta * [G(p,m_plus)-G(p,m_minus)]/2,

where eta=-1 for alpha=0 (S1dag S1dag), +1 for alpha=-2 (S2 S2).
The loop is evaluated by spectral replacement in the *existing* original
one-loop C5 kernel, after setting its lambdaT3 to one. The source tensors
must prove proportionality of the new conjugate HH tensor to T3MixCG.

Only a CANDIDATE matching: the Matchete/RGBeta coupling identification is
an explicit sign assumption, and complete sequential Wilson matching is not
claimed. No independent one-loop formula or group factor is invented.
"""
from __future__ import annotations

import re
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np

from Numerical.experimental.T3SelfBilinearMassBasis import parse_matchete_constant

MODELS=('A_+0','A_-2','D_+0','D_-2','E_+0','E_-2')
EXTRA={0:'lambdaT3GeneralQ021010I1',-2:'lambdaT3GeneralQ020101I1'}
MATCHE_TE={0:'lambdaT3ExtraQ021010I1',-2:'lambdaT3ExtraQ020101I1'}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError('T3 self-bilinear anomalous C5: '+message)


def model_alpha(model: str) -> int:
    _require(model in MODELS, 'unsupported model '+str(model))
    return int(model.split('_')[1])


def native_self_metric(dim: int) -> np.ndarray:
    """Matchete real Cartesian singlet/triplet metric; not RGBeta spherical C."""
    _require(dim in (1,3), 'only singlet and triplet real irreps are supported')
    return np.eye(dim,dtype=complex)


def _number_list(raw: str) -> list[int]:
    _require(bool(re.fullmatch(r'\d+(?:\s*,\s*\d+)*',raw.strip())), 'invalid sparse-array index list')
    return [int(x) for x in raw.split(',')]


def parse_matchete_sparse(raw: str) -> np.ndarray:
    """Restricted parser for the exact Stage group-factor SparseArray InputForm.

    The source stores CSR with one-based column-index tuples; the first
    (row) index is encoded by row pointers. A general Mathematica parser is
    deliberately not used here.
    """
    pattern=(r'^SparseArray\[Automatic,\s*\{(?P<dims>[0-9, ]+)\},\s*0,\s*'
             r'\{1,\s*\{\{(?P<ptr>[0-9, ]+)\},\s*\{(?P<cols>(?:\{[0-9, ]+\},?\s*)+)\}\},\s*'
             r'\{(?P<vals>.*?)\}\}\]$')
    match=re.fullmatch(pattern,raw.strip(),flags=re.S)
    _require(match is not None,'T3MixCG is not a supported native sparse tensor')
    dims=_number_list(match.group('dims'))
    ptr=_number_list(match.group('ptr'))
    _require(len(dims) in (2,3,4) and all(0<d<=3 for d in dims), 'invalid CG dimensions')
    _require(len(ptr)==dims[0]+1 and ptr[0]==0 and all(a<=b for a,b in zip(ptr,ptr[1:])),
             'invalid CG CSR row pointers')
    coords=[_number_list(x[1:-1]) for x in re.findall(r'\{[0-9, ]+\}',match.group('cols'))]
    value_strings=[x.strip() for x in match.group('vals').split(',')]
    _require(ptr[-1]==len(coords)==len(value_strings) and
             all(len(c)==len(dims)-1 and all(1<=v<=d for v,d in zip(c,dims[1:])) for c in coords),
             'CG CSR coordinates/values mismatch')
    tensor=np.zeros(dims,dtype=complex)
    for row in range(dims[0]):
        for index in range(ptr[row],ptr[row+1]):
            tup=(row,)+tuple(v-1 for v in coords[index])
            _require(tensor[tup]==0,'duplicate nonzero CG sparse index')
            tensor[tup]=parse_matchete_constant(value_strings[index])
    _require(np.linalg.norm(tensor)>0,'zero T3 mixing tensor')
    return tensor


def topology_tensor(group_record: dict, ds1: int,ds2: int) -> np.ndarray:
    _require(group_record.get('status')=='Success','invalid source SU2 CG registry')
    tensors=[entry for entry in group_record.get('CGRegistry',[])
             if entry.get('Name')=='T3MixCG']
    _require(len(tensors)==1,'missing/ambiguous T3MixCG')
    t=parse_matchete_sparse(tensors[0]['TensorInputForm'])
    expected=tuple(d for d in (2,2,ds1,ds2) if d!=1)
    _require(t.shape==expected,'CG topology has unexpected SU2 representation dimensions')
    return t.reshape((2,2,ds1,ds2))


def extra_hh_tensor(stage22: dict,model: str,ds1: int,ds2: int) -> tuple[np.ndarray,str]:
    _require(stage22.get('format')=='stage22-exact-tensors-v1' and
             stage22.get('source')=='MatcheteInvariantTensors' and
             stage22.get('quartic_slot_count')==50 and stage22.get('model_count')==16,
             'unexpected Stage22 Matchete tensor catalogue')
    row=stage22.get('models',{}).get(model,{})
    alpha=model_alpha(model)
    symbol=MATCHE_TE[alpha]
    match=[x for x in row.get('tensors',[]) if x.get('symbol')==symbol]
    _require(len(match)==1,'missing/ambiguous independent HH extra quartic: '+symbol)
    obj=match[0]
    species=(['Hdag','Hdag','S1','S2'] if alpha==0
             else ['Hdag','Hdag','S1dag','S2dag'])
    _require(obj.get('species')==species and
             obj.get('dimensions')==[2,2,ds1,ds2] and
             obj.get('potential_prefactor')=='1',
             'wrong fields, tensor axes or quartic prefactor')
    components=obj.get('components')
    _require(isinstance(components,list) and len(components)==4*ds1*ds2,
             'unexpected native HH component count')
    t=np.array([parse_matchete_constant(x) for x in components],dtype=complex)
    return t.reshape((2,2,ds1,ds2)),symbol


def verify_native_self_metric(stage22:dict,model:str,dim:int) -> float:
    """Use the independent 4-self-field Matchete tensor to determine delta.

    Its source-contracted invariant is sym(delta_ab delta_cd), not the
    spherical C_ab; this is an explicit basis check, not a RGBeta map.
    """
    alpha=model_alpha(model)
    symbol=('lambdaT3ExtraQ000400I1' if alpha==0
            else 'lambdaT3ExtraQ000004I1')
    obj=[x for x in stage22['models'][model]['tensors'] if x['symbol']==symbol]
    _require(len(obj)==1 and obj[0].get('dimensions')==[dim]*4,
             'missing native self-quartic metric witness')
    vals=np.asarray([parse_matchete_constant(x) for x in obj[0]['components']],dtype=complex)
    source=vals.reshape((dim,dim,dim,dim))
    eye=native_self_metric(dim)
    expected=(np.einsum('ab,cd->abcd',eye,eye)+
              np.einsum('ac,bd->abcd',eye,eye)+
              np.einsum('ad,bc->abcd',eye,eye))/3
    err=float(np.linalg.norm(source-expected)/np.linalg.norm(expected))
    _require(err<2e-12,'Matchete self-tensor does not certify native delta metric')
    return err


def tensor_ratio(stage22:dict,group_record:dict,model:str,ds1:int,ds2:int)->dict:
    """Require the anomalous HH tensor and old topology be proportional.

    Conjugation changes Hdag Hdag to H H and one scalar orientation;
    anomalous propagator supplies native delta_ab. No fitted factor.
    """
    alpha=model_alpha(model)
    _require((ds1,ds2) in ((1,3),(3,1),(3,3)),'unexpected A/D/E dimensions')
    dim=ds1 if alpha==0 else ds2
    metric_res=verify_native_self_metric(stage22,model,dim)
    old=topology_tensor(group_record,ds1,ds2)
    candidate,symbol=extra_hh_tensor(stage22,model,ds1,ds2)
    candidate=candidate.conj()
    # The self scalar index is 2 for alpha=0, 3 for alpha=-2.
    axis=2 if alpha==0 else 3
    candidate=np.moveaxis(np.tensordot(candidate,native_self_metric(dim),
                                      axes=([axis],[0])),-1,axis)
    den=np.vdot(old,old)
    _require(den.real>0,'zero topology normalization')
    ratio=complex(np.vdot(old,candidate)/den)
    res=float(np.linalg.norm(candidate-ratio*old)/np.linalg.norm(candidate))
    _require(res<2e-12,'extra HH tensor not proportional to original topology')
    _require(abs(ratio.imag)<1e-12 and abs(ratio.real-1)<1e-12,
             'unexpected HH anomalous/topology CG ratio')
    return {'ratio':ratio,'relative_tensor_residual':res,
            'native_self_metric_residual':metric_res,
            'matchete_extra_symbol':symbol,'rgbeta_extra_symbol':EXTRA[alpha],
            'metric_basis':'native Matchete Cartesian delta; RGBeta spherical map not certified'}


def anomalous_phase(mu_sq:complex,alpha:int)->complex:
    """The <Sdag Sdag> (alpha=0) or <S S> (alpha=-2) phase."""
    _require(alpha in (0,-2),'invalid holomorphic mass model alpha')
    z=complex(mu_sq)
    _require(np.isfinite([z.real,z.imag]).all(),'nonfinite holomorphic mass')
    if z==0:return 0j
    p=z/abs(z)
    return p.conjugate() if alpha==0 else p


def anomalous_from_unit_kernels(unit_plus:np.ndarray,unit_minus:np.ndarray,
                                 *, mu_sq:complex,alpha:int,
                                 matchete_extra:complex,tensor_factor:complex=1+0j)->np.ndarray:
    """One-anomalous-scalar-line spectral reweighting of the source C5 kernel."""
    a,b=(np.asarray(x,dtype=complex) for x in (unit_plus,unit_minus))
    _require(a.shape==b.shape==(3,3),'expected two 3x3 unit-topology kernels')
    _require(np.isfinite(a).all() and np.isfinite(b).all(),'nonfinite kernels')
    _require(np.allclose(a,a.T,rtol=1e-10,atol=1e-18) and
             np.allclose(b,b.T,rtol=1e-10,atol=1e-18),'non-symmetric kernels')
    phase=anomalous_phase(mu_sq,alpha)
    lam=complex(matchete_extra)
    fac=complex(tensor_factor)
    _require(np.isfinite([lam.real,lam.imag,fac.real,fac.imag]).all(),
             'nonfinite quartic or group factor')
    result=(lam.conjugate()*phase*fac/2)*(a-b)
    _require(np.isfinite(result).all() and np.allclose(result,result.T,atol=1e-18),
             'nonfinite or asymmetric anomalous result')
    return result


def unit_topology_kernels(kernel:Any,uv:Any,eft:Any,basis:Any,normal_saved:np.ndarray):
    """Evaluate source λT3-linear kernel at one, not divide by tiny λT3.

    Enforce exact linearity with the saved original topology C5 result.
    """
    _require(basis.species in ('S1','S2'),'bad self species')
    attr='mS1Sq' if basis.species=='S1' else 'mS2Sq'
    lam=complex(eft.lambdaT3)
    expected=np.asarray(normal_saved,dtype=complex)
    _require(expected.shape==(3,3),'missing saved normal C5')
    unit=[]
    for msq in basis.eigen_masses_squared:
        new=replace(eft,**{attr:float(msq),'lambdaT3':1+0j}).validated()
        val=np.asarray(kernel.evaluate(uv,new),dtype=complex)
        _require(val.shape==(3,3) and np.isfinite(val).all(),
                 'invalid source C5 at unit topology coupling')
        unit.append(val)
    reconstructed=lam*0.5*(unit[0]+unit[1])
    relative=float(np.linalg.norm(reconstructed-expected)/max(np.linalg.norm(expected),1e-35))
    _require(relative<3e-8,'existing kernel is not linear in λT3 at mass-split endpoint')
    return unit[0],unit[1],relative
