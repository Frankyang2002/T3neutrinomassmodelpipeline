#!/usr/bin/env python3
"""Track B: independent broken-phase T3 A/D/E Dirac self-energies.

Strict data ordering: source benchmark + native Matchete CGs + source operator
coefficients -> independent scalar Hessians -> loop integrals -> Takagi ->
(optional, final) provisional comparison. NEVER imports the production loop
engine or any precomputed provisional tensor/matrix in its calculation.

A present fresh `TrackBADENativeAction.wl` witness checks local native actions.
Absent fresh actions: archived *actual native CG exports* may be used only for
explicitly labelled conditional calculations (never certified).
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime
from functools import lru_cache
from hashlib import sha256
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import traceback

import numpy as np
from scipy.integrate import quad
from sympy.parsing.mathematica import parse_mathematica
import mpmath as mp

MODELS = {
    "A_-4": (1,3,2,-4), "A_-2": (1,3,2,-2),
    "D_-2": (3,1,2,-2), "D_+0": (3,1,2,0), "D_+2": (3,1,2,2),
    "E_-4": (3,3,2,-4), "E_-2": (3,3,2,-2), "E_+0": (3,3,2,0),
    "E_+2": (3,3,2,2),
}
SRC = [
    "Lagrangian/model/T3ModelCatalog.wl",
    "Lagrangian/model/T3Fields.wl",
    "Lagrangian/model/SU2Invariants.wl",
    "Lagrangian/model/LagrangianBuilder.wl",
    "Lagrangian/interactions/T3Topology.wl",
    "Lagrangian/interactions/ScalarPotential.wl",
]
CG_BASE = Path("output/general_observable/lo_direct_one_loop/source_cg_all16")
PORTAL_BASE = Path("output/general_observable/lo_direct_one_loop/source_portal_cg_all16")
FROZEN_BASE = Path("output/general_observable/lo_direct_one_loop/provisional")
OUT_BASE = Path("output/general_observable/physical_normalisation")
WITNESS_FORMAT = "t3-track-b-ade-native-action-v1"
MAX_TENSOR_RESIDUAL = 5e-11
MATRIX_TOL = 2e-12
LOOP_TOL = 5e-10


def require(test, why):
    if not test:
        raise ValueError(why)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    require(path.is_file(), f"Missing file: {path}")
    return json.loads(path.read_text(encoding="utf-8-sig"))


def complex_json(x):
    """Source complex number, vector or matrix (real/imag object)."""
    if isinstance(x, dict) and set(x.keys()) == {"real", "imag"}:
        a = np.asarray(x["real"], dtype=float) + 1j*np.asarray(x["imag"], dtype=float)
    elif isinstance(x, (int,float)):
        a = np.asarray(x, dtype=complex)
    else:
        raise ValueError("Unexpected complex datum encoding")
    require(np.isfinite(a).all(), "Nonfinite source parameter")
    return a


def encode(z):
    a = np.asarray(z, dtype=complex)
    return {"real": a.real.tolist(), "imag": a.imag.tolist()}


def numeric_tensor_from_inputform(source):
    """Parse native Matchete InputForm; no eval and no provisional CG builder.

    Supported: dense Lists; Mathematica SparseArray CSR with the first index
    as outer row and the remaining indices as 1-based coordinate tails.
    """
    x = parse_mathematica(source)
    if getattr(x.func, "__name__", "") != "SparseArray":
        def arrify(t):
            return [arrify(v) for v in t] if isinstance(t, (tuple, list)) else complex(t.evalf(45))
        a = np.asarray(arrify(x), dtype=complex)
    else:
        dims = tuple(int(v) for v in x.args[1]); data = x.args[3]
        require(int(data[0]) == 1, "Nonstandard SparseArray version")
        rowptr, coords = data[1]
        vals = data[2]
        require(len(rowptr) == dims[0]+1 and int(rowptr[-1]) == len(vals) == len(coords),
                "Malformed native CSR offsets")
        a = np.zeros(dims, dtype=complex)
        for r in range(dims[0]):
            for n in range(int(rowptr[r]), int(rowptr[r+1])):
                tail = tuple(int(v)-1 for v in coords[n])
                ix = (r,) + tail
                require(len(ix) == len(dims) and all(0<=v<d for v,d in zip(ix,dims)),
                        "Native CSR out-of-range tensor index")
                a[ix] = complex(vals[n].evalf(45))
    require(np.isfinite(a).all(), "Nonfinite native CG")
    return a


def fresh_tensor(obj):
    sh = tuple(obj["shape"])
    real=np.asarray(obj["real_flat"],float)
    imag=np.asarray(obj["imag_flat"],float)
    require(real.size==imag.size==math.prod(sh),"Malformed fresh tensor payload")
    return (real+1j*imag).reshape(sh)


def pad_singlets(a, expected):
    nontrivial=tuple(d for d in expected if d>1)
    require(a.shape==nontrivial or a.shape==expected,
            f"Native dimensions {a.shape} incompatible with source field dimensions {expected}")
    return a.reshape(expected)


def weights(d):
    """Exact eigenvector phase convention for the saved native SU(2) generator.

    U3 columns are (+1,0,-1) eigenvectors of source T3 with explicit
    T3_cart=[[0,0,-i],[0,0,0],[i,0,0]]. This basis choice must be checked
    against a live generator witness before physical absolute certification.
    """
    if d==1:
        return np.eye(1,dtype=complex), [0]
    if d==3:
        u=np.array([[-1j,0,1j],[0,1,0],[1,0,1]],complex)/np.sqrt(2)
        u[1,1]=1.0 # do not divide neutral-axis normalization by sqrt(2)
        t=np.array([[0,0,-1j],[0,0,0],[1j,0,0]],complex)
        require(np.linalg.norm(u.conj().T@u-np.eye(3))<1e-13, "Bad native triplet U")
        require(np.linalg.norm(u.conj().T@t@u-np.diag([1,0,-1]))<1e-13,
                "Native triplet Cartan diagonalization failed")
        return u, [1,0,-1]
    raise ValueError("Only source 1,3 scalar irreps allowed")


def file_key(d1,d2,df,alpha):
    a=("m"+str(-alpha)) if alpha<0 else ("p"+str(alpha))
    return f"T3_dS1_{d1}_dS2_{d2}_dF_{df}_alpha_{a}"


def field_inputs(root, model, dims, native_dir=None, require_fresh=False):
    d1,d2,df,alpha=dims
    stem=file_key(*dims)
    bp=root/"Numerical/benchmarks/fitted_reference"/(stem+".json")
    bench=read(bp)
    require(tuple(bench["representation"][k] for k in ("d_s1","d_s2","d_f","alpha"))==dims,
            f"{model}: benchmark representation mismatch")
    require(bench.get("base_state",{}).get("ordinary") is not None,
            f"{model}: ordinary source benchmark absent")
    o=bench["base_state"]["ordinary"]
    y1,y2,mf=[complex_json(o[z]) for z in ("y1","y2","MF")]
    require(y1.shape==y2.shape==mf.shape==(3,3),f"{model}: flavour shape")
    require(np.linalg.norm(mf-np.diag(np.diag(mf)))<1e-10 and
            np.max(np.abs(np.diag(mf).imag))<1e-10 and np.min(np.diag(mf).real)>0,
            f"{model}: only positive diagonal benchmark MF accepted")
    # These are the actual source benchmark values, NOT chosen benchmark fit CGs.
    lam=complex(complex_json(o["lambdaT3"]).item())
    require(abs(lam)>0,"Benchmark lambdaT3 is zero")
    for k in ("lambdaH1Adj","lambdaH2Adj","lambda12Adj", "lambdaHHdagS2S2",
              "lambdaHHdagS1barS1bar", "lambdaHHdagS1barS2barCross"):
        if k in o:
            require(abs(complex_json(o[k]).item())<1e-14,
                    f"{model}: unsupported extra scalar coupling is nonzero: {k}")

    fresh_path = root/Path(native_dir)/f"{model.replace('+','p').replace('-','m')}.json" if native_dir else None
    fresh = read(fresh_path) if fresh_path is not None and fresh_path.exists() else None
    if fresh is None and require_fresh:
        raise ValueError(f"{model}: FRESH_NATIVE_ACTION_REQUIRED at {fresh_path}")
    if fresh is not None:
        require(fresh.get("format")==WITNESS_FORMAT and fresh.get("model")==model and
                fresh.get("status")=="NativeRestrictedADEActionExtracted" and
                fresh.get("native_matchete_checked_lagrangian") is True,
                f"{model}: invalid fresh native action witness")
        require(tuple(fresh["dimensions_alpha"])==dims, "Fresh native representation mismatch")
        checks=fresh["checks"]
        require(checks and all(v is True for v in checks.values()),
                f"{model}: native action checks failed")
        require(fresh["native_action"]["fermion_mass_coefficient"]=="-1" and
                fresh["native_action"]["flavour_mass_indexed"] is True,
                f"{model}: native Dirac mass/index contract failed")
        for path,h in fresh.get("source_sha256",{}).items():
            require((root/path).is_file() and digest(root/path).lower()==h.lower(),
                    f"{model}: source changed since native witness: {path}")
        t={key:fresh_tensor(obj) for key,obj in fresh["native_tensors"].items()}
        cg_origin="FRESH_VALIDATED_Matchete_action_and_CG"
    else:
        oldcg=read(root/CG_BASE/stem/"uv_cg_registry.json")
        oldp=read(root/PORTAL_BASE/stem/"native_higgs_portal_registry.json")
        require(oldcg["status"]==oldp["status"]=="Success" and
                all(int(oldcg["metadata"][k])==v for k,v in zip(("dS1","dS2","dF","alpha"),dims)),
                f"{model}: bad archived native source metadata")
        require(all(int(oldp["metadata"][k])==v for k,v in zip(("dS1","dS2","dF","alpha"),dims)),
                f"{model}: bad archived portal metadata")
        t={r["Name"]:numeric_tensor_from_inputform(r["TensorInputForm"])
           for r in oldcg["CGRegistry"]+oldp["CGRegistry"]}
        cg_origin="ARCHIVED_ACTUAL_Matchete_CG_ONLY_native_action_not_measured_for_model"
    require(all(x in t for x in ("T3Y1CG","T3Y2CG","T3MixCG")),
            f"{model}: topology CG absent")
    t["T3Y1CG"]=pad_singlets(t["T3Y1CG"],(2,2,d1))
    t["T3Y2CG"]=pad_singlets(t["T3Y2CG"],(2,2,d2))
    t["T3MixCG"]=pad_singlets(t["T3MixCG"],(2,2,d1,d2))
    for k in (1,2):
        d=(d1,d2)[k-1]
        names=[f"T3HiggsPortal{k}CGInv{i}" for i in range(1,3 if d==3 else 2)]
        for n in names:
            require(n in t,f"{model}: missing full portal invariant {n}")
            t[n]=pad_singlets(t[n],(2,2,d,d))
    return bench,y1,y2,np.diag(mf).real,lam,t,cg_origin,{
        "benchmark":str(bp.relative_to(root)).replace('\\','/'),
        "benchmark_sha256":digest(bp),
        "native_witness":str(fresh_path.relative_to(root)) if fresh is not None else None,
        "archived_cg":str((CG_BASE/stem/"uv_cg_registry.json")).replace('\\','/') if fresh is None else None,
        "archived_portal":str((PORTAL_BASE/stem/"native_higgs_portal_registry.json")).replace('\\','/') if fresh is None else None,
        "local_source_sha256":{p:digest(root/p) for p in SRC if (root/p).exists()},
    }


def identity_decomposition(t,which,d):
    names=[f"T3HiggsPortal{which}CGInv{i}" for i in range(1,3 if d==3 else 2)]
    tensors=[t[n] for n in names]
    goal=np.einsum('ij,ab->ijab',np.eye(2),np.eye(d)).astype(complex)
    cols=np.column_stack([z.ravel() for z in tensors])
    arr=np.vstack((cols.real,cols.imag))
    v=np.r_[goal.ravel().real,goal.ravel().imag]
    sol, *_=np.linalg.lstsq(arr,v,rcond=1e-13)
    resid=float(np.linalg.norm(cols@sol-goal.ravel()) / max(np.linalg.norm(goal),1e-50))
    require(resid<MAX_TENSOR_RESIDUAL, f"Portal {which}: full tensor identity not in native invariant span; residual={resid}")
    return names,sol,resid


def build_mass_blocks(dims,t,o,lam,vev,portal):
    d1,d2,df,alpha=dims
    U1, w1=weights(d1);U2,w2=weights(d2)
    n1,c1,r1=identity_decomposition(t,1,d1)
    n2,c2,r2=identity_decomposition(t,2,d2)
    require(portal in ('zero-portal','conditional-benchmark-identity'),"Unknown portal point")
    if portal=='zero-portal':
        src1=np.zeros(len(c1));src2=np.zeros(len(c2))
    else:
        src1=float(o['lambdaH1'])*c1
        src2=float(o['lambdaH2'])*c2
    p1=sum((src1[i]*t[name][1,1,:,:] for i,name in enumerate(n1)),
           np.zeros((d1,d1),complex))
    p2=sum((src2[i]*t[name][1,1,:,:] for i,name in enumerate(n2)),
           np.zeros((d2,d2),complex))
    # Source Lint has +1/2 PlusHc[lambdaH Hbar H Sbar S CG]
    # => physical Hessian A = mS^2 I - v^2/2 sum(lambdaH_i CG_i).
    H1=float(o['mS1Sq'])*np.eye(d1)-vev**2/2*(U1.conj().T@p1@U1)
    H2=float(o['mS2Sq'])*np.eye(d2)-vev**2/2*(U2.conj().T@p2@U2)
    cm=t['T3MixCG'][1,1,:,:]
    mix_weight=U1.T@cm@U2.conj()
    # Source Lint contains +lambdaT3 HH S1 S2bar + h.c. (no 1/2).
    # Lquadratic=-Sdag A S therefore A_S2,S1=-lambda v^2/2*Cmix.
    A21=-lam*vev**2/2*mix_weight.T
    YS1=alpha/2;YS2=(alpha+2)/2
    charges1=[w+YS1 for w in w1]; charges2=[w+YS2 for w in w2]
    max_wrong=max([abs(A21[j,i]) for i in range(d1) for j in range(d2)
                   if abs(charges1[i]-charges2[j])>1e-10] or [0.0])
    require(max_wrong<1e-10*max(1.0,np.linalg.norm(A21)),
            f"Charge-nonconserving source mixing survives native projection: {max_wrong}")
    n=d1+d2
    H=np.zeros((n,n),complex)
    H[:d1,:d1]=H1;H[d1:,d1:]=H2
    H[d1:,:d1]=A21;H[:d1,d1:]=A21.conj().T
    rel=float(np.linalg.norm(H-H.conj().T)/max(np.linalg.norm(H),1e-50))
    require(rel<1e-12,"Source Hessian non-Hermitian")
    allQ=charges1+charges2
    sectors={}
    for Q in sorted(set(allQ)):
        inds=[k for k,q in enumerate(allQ) if abs(q-Q)<1e-8]
        block=H[np.ix_(inds,inds)]
        masses=np.linalg.eigvalsh(block)
        require(np.min(masses)>0, f"Non-positive scalar mass in charge sector Q={Q}: {masses}")
        sectors[str(int(Q)) if float(Q).is_integer() else str(Q)]={
            'Q':Q,'indices':inds,'masses_squared_gev2':masses.tolist(),
            'H_gev2':encode(block)}
    return H, sectors,(U1,U2,w1,w2,charges1,charges2),{
        "identity_full_tensor_relative_residual":{"S1":r1,"S2":r2},
        "native_identity_decomposition":{"S1":dict(zip(n1,c1.tolist())),
                                         "S2":dict(zip(n2,c2.tolist()))},
        "source_portal_coefficients":{"S1":dict(zip(n1,src1.tolist())),
                                      "S2":dict(zip(n2,src2.tolist()))},
        "max_wrong_charge_mixing_gev2":max_wrong,
        "full_Hessian_gev2":encode(H),
        "scalar_mass_positive":True}


@lru_cache(maxsize=8192)
def divide_loop(M,slow,shigh):
    """Two genuinely separate finite one-loop scalar B0 divided differences.

    Finite B0 at p²=0, mu²=M²: 1+s*ln(s/M²)/(M²-s).
    Integral: -∫ dx log[1+(1-x)(shigh-slow)/(x*M²+(1-x)*slow)]/(shigh-slow).
    Analytic side uses 80-digit mpmath; quadrature side scipy QUADPACK.
    """
    require(M>0 and slow>0 and shigh>slow,"Invalid masses for loop")
    M2=M*M;diff=shigh-slow
    def f(x):
        return -math.log1p((1-x)*diff/(x*M2+(1-x)*slow))/diff
    val,err=quad(f,0.0,1.0,epsabs=5e-21,epsrel=4e-13,limit=400)
    with mp.workdps(75):
        mm=mp.mpf(str(M))**2;a=mp.mpf(str(slow));b=mp.mpf(str(shigh))
        def B0(s):
            # mpmath handles even accidental M² ~= s without float cancellation
            z=s/mm
            if abs(z-1)<mp.mpf('1e-40'):
                return mp.mpf(0)
            return 1+s*mp.log(z)/(mm-s)
        ana=float((B0(b)-B0(a))/(b-a))
    relative=abs(val-ana)/max(abs(ana),1e-50)
    require(relative<LOOP_TOL,f"Analytic B0 / Feynman quadrature mismatch={relative}")
    return ana,val,relative,err


def evaluate_matrix(dims,t,bench,y1,y2,masses,lam,portal,perform_independent=True):
    o=bench['base_state']['ordinary'];d1,d2,df,alpha=dims
    v=float(bench['numerical']['vev_gev'])
    H,sectors,weights_data,meta=build_mass_blocks(dims,t,o,lam,v,portal)
    U1,U2,w1,w2,q1,q2=weights_data
    fy=[(alpha+1)/2+0.5,(alpha+1)/2-0.5]
    m=np.zeros((3,3),complex)
    channel_rows=[];pole_max=0.0;quad_max=0.0
    y1g=t['T3Y1CG'];y2g=t['T3Y2CG']
    for fi,Qf in enumerate(fy):
        idx1=[i for i,q in enumerate(q1) if abs(q-Qf)<1e-8]
        idx2=[j for j,q in enumerate(q2) if abs(q-Qf)<1e-8]
        require(len(idx1)<=1 and len(idx2)<=1, "Non-unique scalar weight in a charge sector")
        if not idx1 or not idx2:
            channel_rows.append({'fermion_SU2_index_1based':fi+1,'fermion_Q':int(Qf),
                                 'active':False,'reason':'No matching scalar in one multiplet'})
            continue
        i,j=idx1[0],idx2[0]
        c1=np.dot(y1g[0,fi,:],U1[:,i])
        c2=np.dot(y2g[0,fi,:],U2[:,j].conj())
        a21=H[d1+j,i]
        # Selected physical 2x2 block in charge-Q sector; do not use frozen A/B matrices.
        two=H[np.ix_([i,d1+j],[i,d1+j])]
        # Additional same-charge states (not expected for 1/3 irreps) are prohibited.
        matching=[k for k,q in enumerate(q1+q2) if abs(q-Qf)<1e-8]
        require(len(matching)==2, f"Unexpected degenerate SU(2) charge block Q={Qf}")
        sl,sh=np.linalg.eigvalsh(two)
        # The UV pole multiplies sum_i U_{2i} U*_{1i}=0, unlike individual B0 terms.
        w,uu=np.linalg.eigh(two)
        pole=abs(np.sum(uu[1,:]*uu[0,:].conj()))
        pole_max=max(pole_max,pole)
        require(pole<2e-12, f"UV pole does not cancel: {pole}")
        coefficients=[]
        for r,M in enumerate(masses):
            if abs(a21)<1e-40 or abs(c1*c2)<1e-40:
                kernel=0.0;rel=0.0
            else:
                kernel,_,rel,_=divide_loop(float(M),float(sl),float(sh))
            quad_max=max(quad_max,rel)
            g1=-np.conj(y1[:,r]*c1)
            g2=-np.conj(y2[:,r]*c2)
            # Two distinguishable chiral vertices, one Dirac mass insertion;
            # sum the two external neutrino orientations, no invented 1/2.
            coeff=M*a21*kernel/(16*np.pi**2)
            m+=coeff*(np.outer(g1,g2)+np.outer(g2,g1))
            coefficients.append({'flavour_index':r+1,'M_gev':float(M),
                                 'finite_divided_difference_gevminus2':float(kernel),
                                 'B0_vs_Feynman_relative_residual':float(rel)})
        channel_rows.append({'fermion_SU2_index_1based':fi+1,'fermion_Q':int(Qf),
                             'active':True,'S1_weight':w1[i],'S2_weight':w2[j],
                             'S1_native_selected_CG':encode(c1),'S2bar_native_selected_CG':encode(c2),
                             'mixing_A_S2_S1_gev2':encode(a21),
                             'charge_block_eigenmass_squared_gev2':[float(sl),float(sh)],
                             'uv_pole_cross_weight_absolute':pole,'flavours':coefficients})
    sym=float(np.linalg.norm(m-m.T)/max(np.linalg.norm(m),1e-50))
    require(sym<MATRIX_TOL, f"Majorana matrix asymmetry {sym}")
    m_ev=m*1e9
    takagi=np.linalg.svd(m_ev,compute_uv=False)
    require(np.isfinite(takagi).all() and min(takagi)>=0,"Nonphysical Takagi singular values")
    sv2=np.sqrt(np.maximum(0,np.linalg.eigvalsh(m_ev.conj().T@m_ev)))
    sv2=sv2[::-1]
    svres=np.linalg.norm(takagi-sv2)/max(np.linalg.norm(takagi),1e-40)
    require(svres<2e-7, f"Takagi two-method singular value residual {svres}")
    require(pole_max<2e-12 and quad_max<LOOP_TOL,"UV validation failed")
    return m_ev, {
        'scalar':meta,'charge_sectors':sectors,'channels':channel_rows,
        'checks':{'matrix_relative_asymmetry':sym,'maximum_channel_uv_pole_absolute':pole_max,
                  'maximum_integral_relative_residual':quad_max,
                  'takagi_svd_vs_gram_relative_residual':float(svres)},
        'neutrino_takagi_masses_ev_descending':takagi.tolist(),
        'matrix_ev':encode(m_ev)}


def rel(a,b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b),1e-40))


def phase_tests(dims,t,bench,y1,y2,masses,lam):
    portal='conditional-benchmark-identity'
    a=y1.copy();b=y2.copy();a[0,1]*=np.exp(0.37j);a[1,2]*=np.exp(-0.19j)
    b[2,0]*=np.exp(0.41j);b[0,2]*=np.exp(-0.23j)
    z=lam*np.exp(0.29j)
    orig,_=evaluate_matrix(dims,t,bench,a,b,masses,z,portal)
    phi=np.array([.12,-.28,.39])
    phase=np.exp(-1j*phi[:,None])
    new,_=evaluate_matrix(dims,t,bench,a*phase,b*phase,masses,z,portal)
    P=np.diag(np.exp(1j*phi))
    ext=rel(new,P@orig@P)
    ph1,ph2,phL=.21,-.17,.43
    rot,_=evaluate_matrix(dims,t,bench,a*np.exp(1j*ph1),b*np.exp(1j*ph2),
                          masses,z*np.exp(1j*phL),portal)
    spur=rel(rot,orig*np.exp(1j*(phL-ph1-ph2)))
    wrong,_=evaluate_matrix(dims,t,bench,a.conj(),b,masses,z,portal)
    wrongres=rel(wrong,orig)
    zero_y1,_=evaluate_matrix(dims,t,bench,np.zeros_like(a),b,masses,z,portal)
    zero_lam,_=evaluate_matrix(dims,t,bench,a,b,masses,0j,portal)
    nz=max(np.linalg.norm(orig),1e-40)
    require(ext<2e-11 and spur<2e-11,"Complex rephasing covariance failed")
    require(wrongres>1e-4,"Negative control: wrong Yukawa conjugation not detected")
    require(np.linalg.norm(zero_y1)<1e-12*nz and np.linalg.norm(zero_lam)<1e-12*nz,
            "Negative control: lepton-number violating coupling zero limit failed")
    return {'complex_nonuniform_matrix_ev':encode(orig),
            'complex_nonuniform_takagi_ev_descending':np.linalg.svd(orig,compute_uv=False).tolist(),
            'external_neutrino_congruence_relative_residual':ext,
            'source_spurion_phase_relative_residual':spur,
            'wrong_y1_conjugation_negative_control_residual':wrongres,
            'zero_y1_relative_frobenius':float(np.linalg.norm(zero_y1)/nz),
            'zero_lambdaT3_relative_frobenius':float(np.linalg.norm(zero_lam)/nz),
            'phase_choices':{'y1_entries':{'0,1':.37,'1,2':-.19},
                             'y2_entries':{'2,0':.41,'0,2':-.23},
                             'lambdaT3':.29,'external_flavour':phi.tolist(),
                             'spurion_y1':ph1,'spurion_y2':ph2,'spurion_lambdaT3':phL},
            'absolute_Matchete_GammaCC_phase_verified':False}


def compare_last(root,model,our):
    f=root/FROZEN_BASE/(model+'_direct_one_loop.json')
    if not f.is_file():
        return {'status':'frozen_comparison_pending','reason':'frozen reference file missing'}
    j=read(f)
    x=j.get('total',{}).get('complex_symmetric_matrix_gev')
    require(x is not None, "Frozen comparison JSON lacks matrix")
    frozen=complex_json(x)*1e9
    require(frozen.shape==(3,3),"Frozen comparison wrong matrix shape")
    return {'status':'comparison_only_not_certification','frozen_file':str(f.relative_to(root)).replace('\\','/'),
            'frozen_sha256':digest(f),'relative_matrix_residual':rel(our,frozen),
            'relative_matrix_residual_up_to_global_minus':min(rel(our,frozen),rel(our,-frozen)),
            'frozen_takagi_ev_descending':np.linalg.svd(frozen,compute_uv=False).tolist()}


def model_result(root,model,dims,args):
    bench,y1,y2,masses,lam,t,origin,prov=field_inputs(root,model,dims,args.native_dir,args.native_only)
    zero,mz=evaluate_matrix(dims,t,bench,y1,y2,masses,lam,'zero-portal')
    identity,mi=evaluate_matrix(dims,t,bench,y1,y2,masses,lam,'conditional-benchmark-identity')
    controls=phase_tests(dims,t,bench,y1,y2,masses,lam)
    # Comparison intentionally happens AFTER both matrices and validation.
    comparison=(compare_last(root,model,identity) if not args.no_compare
                else {"status":"not_requested"})
    fresh=origin.startswith('FRESH_')
    cert=False # physical independent spinor sign and portal-coupling dictionary outstanding
    return {'model':model,'status':'independent_chosen_Weyl_numerical_result_NOT_physically_certified',
        'native_source_status':origin,
        'scalar_fermion_normalisation_status':(
            'native_scalar_action_and_indexed_Dirac_mass_checked; absolute spinor dictionary pending'
            if fresh else 'native CG registry only; nine-model scalar action and indexed mass pending'),
        'independent_matrix_status':('source_action_checked_convention_dependent_matrix' if fresh
                                    else 'archived_native_CGs_source_operator_convention_dependent_matrix'),
        'certification_eligible':cert,
        'certification_blockers':(
            ['concrete installed GammaCC relative-spinor convention not independently measured',
             'native portal/RGBeta benchmark dictionary conditional'] +
            ([] if fresh else ['fresh native nine-model scalar/Yukawa action and flavour-indexed MF absent'])),
        'source_field_dimensions_alpha':list(dims),
        'source_field_hypercharges':{'S1':dims[3]/2,'S2':(dims[3]+2)/2,'F':(dims[3]+1)/2},
        'source_mass_eigenvalues_gev':masses.tolist(),
        'source_lambdaT3':encode(lam),
        'source_provenance':prov,
        'zero_portal_source_point':mz,
        'conditional_benchmark_identity_point':mi,
        'complex_and_negative_controls':controls,
        'frozen_comparison_performed_last':comparison,
        'assumed_two_component_convention':(
            "F=(chi,bar(xi)); C=+i gamma^2 gamma^0; "
            "L_Y=-nu(chi*g1*S1bar+xi*g2*S2)+h.c.; "
            "g1=-conj(y1*C1), g2=-conj(y2*C2), "
            "A_S2,S1=-lambdaT3*v^2*Cmix/2; "
            "L_mass=-MF[r] bar(F[r])F[r] (no Dirac factor 1/2). "
            "Overall/relative installed Matchete GammaCC phase not certified."),
        'comparison_does_not_set_loop_normalisation':True,
        'no_EFT_or_running':True}


def git_head(root):
    try:
        r=subprocess.run(['git','-C',str(root),'rev-parse','HEAD'],text=True,
                         capture_output=True,check=False,timeout=8)
        return r.stdout.strip() if r.returncode==0 else 'UNAVAILABLE'
    except Exception:return 'UNAVAILABLE'


def write_json(path,obj):
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def update_roadmap(root, out, results):
    # A runtime copy only. NEVER edit/overwrite the installed handoff roadmap.
    parent=Path('Numerical/lo_certification/physical_normalisation/parallel_B_ADE/ROADMAP.md')
    installed=root/parent
    template=installed.read_text(encoding='utf-8') if installed.is_file() else '# ROADMAP not installed\n'
    appendix='\n## B4 runtime status update (generated this run, not a source edit)\n\n'
    appendix+=f'Run directory: `{out.relative_to(root)}`. B0 source inventory checked; B1 and B2 numerical gates per row below; B3 numerical controls tested for completed rows; B4 consolidated reports written. Physical eligibility remains false until unresolved native source and spinor gates are independently measured.\n\n'
    appendix+='| Model | Native action | Matrix | Takagi descending (eV) | Comparison residual | Eligible |\n|---|---|---|---|---|---|\n'
    for m,row in results.items():
        if row.get('error'):
            appendix+=f'| {m} | pending | blocked: {str(row["error"]).replace("|","/")} | pending | pending | no |\n';continue
        ms=row['conditional_benchmark_identity_point']['neutrino_takagi_masses_ev_descending']
        cmp=row['frozen_comparison_performed_last'].get('relative_matrix_residual')
        appendix+=f'| {m} | {"passed" if row["native_source_status"].startswith("FRESH_") else "pending (archived CG only)"} | evaluated, convention dependent | {", ".join(f"{x:.8g}" for x in ms)} | {cmp:.4g} | no |\n' if cmp is not None else f'| {m} | pending | evaluated | {", ".join(f"{x:.8g}" for x in ms)} | pending | no |\n'
    (out/'ROADMAP_RUNTIME_STATUS.md').write_text(template+appendix,encoding='utf-8')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root',type=Path,default=Path('.'))
    parser.add_argument('--native-dir',type=Path,default=Path('output/general_observable/physical_normalisation/parallel_B_ADE_native'))
    parser.add_argument('--native-only',action='store_true',help='Block all rows lacking fresh verified Matchete native action')
    parser.add_argument('--timestamp',default=None,help='Optional unique filename suffix for reproducibility')
    parser.add_argument('--no-compare',action='store_true',help='Do not inspect provisional comparison reference')
    parser.add_argument('--run-self-tests',action='store_true')
    args=parser.parse_args(argv)
    root=args.repo_root.resolve()
    if args.run_self_tests:
        self_test()
        return 0
    require((root/'Lagrangian/model/T3ModelCatalog.wl').exists(),
            'Repository root missing; install archive to project root')
    require((root/OUT_BASE).is_dir(), 'Existing output/general_observable/physical_normalisation directory not present')
    stamp=args.timestamp or datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    require(re.fullmatch(r'[A-Za-z0-9_]+',stamp) is not None,'Invalid timestamp suffix')
    out=root/OUT_BASE/('parallel_B_ADE_'+stamp)
    out.mkdir(parents=False,exist_ok=False)
    rows={}
    for model,dims in MODELS.items():
        try:
            result=model_result(root,model,dims,args)
            rows[model]=result
        except Exception as exc:
            rows[model]={'model':model,'status':'BLOCKED',
                         'error':str(exc),'exact_failing_gate':str(exc),
                         'native_source_status':'pending_or_failed',
                         'scalar_fermion_normalisation_status':'pending_or_failed',
                         'independent_matrix_status':'blocked',
                         'takagi_masses_ev_descending':None,
                         'comparison_residual':None,'certification_eligible':False}
        write_json(out/(model.replace('+','p').replace('-','m')+'_report.json'),rows[model])
    table=[]
    for model,row in rows.items():
        good=row.get('status')!='BLOCKED'
        masses=row['conditional_benchmark_identity_point']['neutrino_takagi_masses_ev_descending'] if good else None
        c=row['frozen_comparison_performed_last'] if good else {}
        table.append({'model':model,
                      'native_source_status':row['native_source_status'],
                      'scalar_fermion_normalisation_status':row['scalar_fermion_normalisation_status'],
                      'independent_matrix_status':row['independent_matrix_status'],
                      'takagi_masses_ev_descending':masses,
                      'comparison_residual':c.get('relative_matrix_residual'),
                      'certification_eligible':False,
                      'failing_gate':row.get('exact_failing_gate','; '.join(row.get('certification_blockers',[])))})
    summary={'format':'t3-track-b-ade-independent-uv-nine-batch-v1',
             'status':'nine_rows_complete_with_gates' if len(table)==9 else 'INTERNAL_ERROR',
             'local_git_head':git_head(root),'comparison_is_read_only_and_last':True,
             'all_models':table,
             'evaluated_models':sum(r['status']!='BLOCKED' for r in rows.values()),
             'fresh_native_action_models':sum(r.get('native_source_status','').startswith('FRESH_') for r in rows.values()),
             'physical_certification_count':0,
             'no_import_of_provisional_loop_engine':True,
             'no_github_writes_or_EFT_RGE':True,
             'output_relative':str(out.relative_to(root)).replace('\\','/')}
    write_json(out/'batch_summary.json',summary)
    with (out/'nine_model_ledger.csv').open('w',newline='',encoding='utf-8') as f:
        wr=csv.DictWriter(f,fieldnames=list(table[0]));wr.writeheader()
        for r in table:wr.writerow({k:json.dumps(v) if isinstance(v,list) else v for k,v in r.items()})
    lines=['# Track B A/D/E independent nine-model batch','',
           'All models use the same independent scalar-source/Hessian and two-loop-integral-check code. '
           'Archived native CGs alone do not certify assembled Matchete actions. '
           'Neither a close comparison to frozen masses nor a correct Takagi spectrum is a physical certificate.','',
           '| Model | Native source | Fermion/scalar | Independent matrix | Takagi descending (eV) | Frozen residual | Eligible |',
           '|---|---|---|---|---|---|---|']
    for row in table:
        m=row['takagi_masses_ev_descending'];mass=', '.join(f'{v:.8g}' for v in m) if m else 'pending'
        cres=row['comparison_residual'];cstr=f'{cres:.4g}' if cres is not None else 'pending'
        lines.append(f"| {row['model']} | {'fresh action' if row['native_source_status'].startswith('FRESH_') else 'archived CG / action pending'} | {'native mass checked' if row['native_source_status'].startswith('FRESH_') else 'conditional'} | {'evaluated' if m else 'blocked'} | {mass} | {cstr} | no |")
    lines.extend(['','## Interpretation','',
       'Zero-portal and conditional benchmark-identity matrices are in each model JSON. '
       'The benchmark identity is solved from the complete H-H-S†-S native invariant tensors, '
       'not assigned to Inv1/Inv2 by name.','',
       'For every model, the charged internal fermion components are included whenever both '
       'scalar multiplets have the corresponding electric-charge component. In family E, '
       'the two distinct charge sectors are both considered.','',
       'The result is convention-dependent until an independently inspected concrete spinor '
       'charge-conjugation dictionary and native action are accepted. '
       'Do not interpret small frozen comparison residual as certification.',''])
    (out/'BATCH_SUMMARY.md').write_text('\n'.join(lines),encoding='utf-8')
    update_roadmap(root,out,rows)
    print(json.dumps({'output':str(out),'evaluated':summary['evaluated_models'],
                      'physical_certification_count':0,'nine_models':len(table)},indent=2))
    return 0 if summary['evaluated_models']==len(MODELS) else 2


def self_test():
    # Never use fitted/unfitted CG or frozen masses as a reference in this test.
    sparse='SparseArray[Automatic,{2,2,3},0,{1,{{0,2,3},{{1,1},{2,3},{1,2}}},{3^(-1/4),(-I)/3^(1/4),1/Sqrt[2]}}]'
    ar=numeric_tensor_from_inputform(sparse)
    require(ar.shape==(2,2,3) and abs(ar[0,0,0]-3**(-.25))<1e-14 and
            abs(ar[0,1,2]+1j*3**(-.25))<1e-14,'Native SparseArray parser test')
    for M,a,b in ((1e5,1e6,1.3e6),(1.1e5,1e6,1.01e6),(1200.0,1e6,1.4e6)):
        z,q,r,_=divide_loop(M,a,b)
        require(r<LOOP_TOL and z<0 and q<0,"Loop independent test")
    for d in (1,3):
        U,w=weights(d)
        require(np.linalg.norm(U.conj().T@U-np.eye(d))<1e-13,'Weight basis unitary')
    # Synthetic charge-conserving native-shaped tensors: numerical tests of
    # the ENGINE, explicitly NOT a substitute for live model-specific CG data.
    test_y=np.array([[1, .2, -.1],[.2, .8, .3],[-.1, .3, .6]],complex)*1e-4
    test_z=np.array([[.7,-.3,.2],[.4,1,-.2],[.1,.3,.8]],complex)*1e-4
    obs={'mS1Sq':1306122.448979592,'mS2Sq':1e6,'lambdaH1':.02,'lambdaH2':.03}
    synthetic={'base_state':{'ordinary':obs},'numerical':{'vev_gev':246.22}}
    masses=np.array([1e5,1.1e5,1.2e5])
    for _model,(d1,d2,_df,alpha) in MODELS.items():
        fam=_model[0]
        U1,w1=weights(d1);U2,w2=weights(d2)
        tensors={'T3Y1CG':np.zeros((2,2,d1),complex),
                 'T3Y2CG':np.zeros((2,2,d2),complex),
                 'T3MixCG':np.zeros((2,2,d1,d2),complex)}
        for which,d in ((1,d1),(2,d2)):
            if d==1:
                tensors[f'T3HiggsPortal{which}CGInv1']=np.einsum(
                    'ij,ab->ijab',np.eye(2),np.eye(d))
            else:
                tensors[f'T3HiggsPortal{which}CGInv1']=np.zeros((2,2,d,d))
                tensors[f'T3HiggsPortal{which}CGInv2']=np.einsum(
                    'ij,ab->ijab',np.eye(2),np.eye(d))
        targets=([(1.,-1j,1j)] if fam!='E' else
                 [(-1j,1.,-1j),(-.7,-1j,1j)])
        count=0
        for f in range(2):
            charge=(alpha+1)/2 +(.5 if f==0 else -.5)
            r1=[i for i,w in enumerate(w1) if abs(w+alpha/2-charge)<1e-12]
            r2=[j for j,w in enumerate(w2) if abs(w+(alpha+2)/2-charge)<1e-12]
            if not r1 or not r2:continue
            i,j=r1[0],r2[0];c1,c2,cm=targets[count];count+=1
            tensors['T3Y1CG'][0,f,:]=c1*U1[:,i].conj()
            tensors['T3Y2CG'][0,f,:]=c2*U2[:,j]
            tensors['T3MixCG'][1,1,:,:]+=cm*np.outer(U1[:,i].conj(),U2[:,j])
        require(count==(2 if fam=='E' else 1),'Synthetic charge inventory')
        raw,zmeta=evaluate_matrix((d1,d2,2,alpha),tensors,synthetic,
                                  test_y,test_z,masses,.1,'zero-portal')
        require(np.linalg.norm(raw)>0,'Synthetic UV matrix zero')
        require(len([c for c in zmeta['channels'] if c['active']])==count,
                'Charged channel missing from synthetic self test')
        check=phase_tests((d1,d2,2,alpha),tensors,synthetic,test_y,test_z,masses,.1)
        require(check['source_spurion_phase_relative_residual']<2e-11,
                'Complex source test failure')
    print('Track B independent math/parser/loop/NINE synthetic charge spectra/complex self-tests PASS')


if __name__=='__main__':
    try:sys.exit(main())
    except Exception as e:
        print('Track B batch failed before per-model result creation:',e,file=sys.stderr)
        traceback.print_exc()
        sys.exit(3)
