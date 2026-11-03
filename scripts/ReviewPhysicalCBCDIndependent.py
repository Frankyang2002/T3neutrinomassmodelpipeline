#!/usr/bin/env python3
"""Track C: charge-resolved B/C Dirac direct UV loop. No production imports.

Reads existing GitHub-native CG/portal exports, fitted-reference source JSONs, and
optionally NEW installed-Matchete action witnesses from the matching Track C WL
script. Archived registries are not passed off as newly executed native actions.

All four models are one batch. Failures are isolated per model. Never overwrites
an existing output directory. No production, RGE, EFT, matching or Git writes.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import re
import sys
import traceback
from datetime import datetime
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

MODELS = {
    "B_-3": ("B", -3, 1, "m3"),
    "B_+1": ("B", +1, 1, "p1"),
    "C_-3": ("C", -3, 3, "m3"),
    "C_+1": ("C", +1, 3, "p1"),
}
SOURCE_PATHS = (
    "Lagrangian/model/T3ModelCatalog.wl",
    "Lagrangian/model/T3Fields.wl",
    "Lagrangian/model/LagrangianBuilder.wl",
    "Lagrangian/model/SU2Invariants.wl",
    "Lagrangian/interactions/T3Topology.wl",
    "Lagrangian/interactions/ScalarPotential.wl",
)
TENSOR_NAMES = ("T3Y1CG", "T3Y2CG", "T3MixCG",
                "T3HiggsPortal1CGInv1", "T3HiggsPortal1CGInv2",
                "T3HiggsPortal2CGInv1", "T3HiggsPortal2CGInv2")
PORTAL_NAMES = TENSOR_NAMES[3:]
PARSER = re.compile(
    r"^SparseArray\[Automatic, \{([\d, ]+)\}, 0, \{1, \{\{([\d, ]+)\}, \{(.*?)\}\}, \{(.*)\}\}\]$"
)
mp.mp.dps = 55


def require(x, msg):
    if not x:
        raise ValueError(msg)


def sha256(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, data):
    require(not path.exists(), f"Refusing overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def enc_matrix(m):
    a = np.asarray(m)
    return {"real": np.real(a).tolist(), "imag": np.imag(a).tolist()}


def expr(s):
    # Read-only, source-tracked Matchete scalar literals, not arbitrary Python eval.
    val = sp.sympify(s.strip().replace("^", "**").replace("Sqrt[", "sqrt(").replace("]", ")"),
                     locals={"I": sp.I, "sqrt": sp.sqrt}, evaluate=True)
    require(not val.free_symbols, f"Non-numerical Matchete tensor value: {s}")
    return complex(val.evalf(35))


def sparse_registry(raw: str):
    m = PARSER.fullmatch(raw)
    require(m is not None, f"Unexpected archived Matchete SparseArray format: {raw[:70]}")
    shape = tuple(int(x.strip()) for x in m.group(1).split(','))
    ptr = [int(x.strip()) for x in m.group(2).split(',')]
    coords = [tuple(map(int, x.split(','))) for x in re.findall(r"\{([\d, ]+)\}", m.group(3))]
    vals = [expr(x) for x in m.group(4).split(',')]
    require(len(coords) == len(vals) == ptr[-1], "Matchete sparse coordinate/value lengths")
    require(len(ptr) == shape[0] + 1, "Matchete CSR first-axis offsets")
    arr = np.zeros(shape, complex)
    for i in range(shape[0]):
        for k in range(ptr[i], ptr[i+1]):
            idx = (i,) + tuple(z-1 for z in coords[k])
            arr[idx] = vals[k]
    return arr


def native_tensor(d):
    shape = tuple(d["shape"])
    return (np.asarray(d["real_flat"], float) + 1j*np.asarray(d["imag_flat"], float)).reshape(shape)


def load_inputs(root, mid, native_dir, source_only):
    cls, alpha, dF, suf = MODELS[mid]
    stem = f"T3_dS1_2_dS2_2_dF_{dF}_alpha_{suf}"
    bpath = root / "Numerical/benchmarks/fitted_reference" / (stem + ".json")
    cr = root / "output/general_observable/lo_direct_one_loop/source_cg_all16" / stem / "uv_cg_registry.json"
    pr = root / "output/general_observable/lo_direct_one_loop/source_portal_cg_all16" / stem / "native_higgs_portal_registry.json"
    paths = [bpath, cr, pr] + [root / p for p in SOURCE_PATHS]
    require(all(p.is_file() for p in paths), "Missing required source: " + ", ".join(str(p) for p in paths if not p.is_file()))
    bench = json.loads(bpath.read_text(encoding="utf-8"))
    cg = json.loads(cr.read_text(encoding="utf-8"))
    portal = json.loads(pr.read_text(encoding="utf-8"))
    rep = bench["representation"]
    require((rep["d_s1"],rep["d_s2"],rep["d_f"],rep["alpha"]) == (2,2,dF,alpha), "Benchmark model mismatch")
    require(cg["status"] == portal["status"] == "Success", "Archived native registry not successful")
    require((int(cg["metadata"]["dF"]), int(cg["metadata"]["alpha"])) == (dF,alpha), "Native CG metadata mismatch")
    reg = {o["Name"]: sparse_registry(o["TensorInputForm"]) for o in (cg["CGRegistry"] + portal["CGRegistry"])}
    require(set(reg) == set(TENSOR_NAMES), "Missing/extra CG tensor")
    flag = "archived_GitHub_native_tensor_registry_only"
    native_path = None
    if native_dir:
        native_path = native_dir / (mid + "_native_action.json")
        if native_path.is_file():
            act = json.loads(native_path.read_text(encoding="utf-8"))
            require(act.get("model") == mid and act.get("status") == "native_action_and_tensors_extracted_spinor_absolute_pending", "Native action status wrong")
            require(act.get("full_CheckLagrangian_valid") is True, "Native action validation absent")
            require(all(v is True for v in act["native_checks"].values()), "Native source checks failed")
            require(act.get("mass",{}).get("coefficients") == ["-1"], "Charged Dirac mass must have -1 coefficient")
            require(act.get("restricted_flags_effective",{}) == {
                "T3_EXPERIMENTAL_GENERAL_FERMIONS":"0",
                "T3_EXPERIMENTAL_GENERAL_QUADRATICS":"0",
                "T3_EXPERIMENTAL_GENERAL_QUARTICS":"0",
                "T3_EXPERIMENTAL_FLAVOR_MASSES":"1"}, "Wrong native flags")
            for src, h in act.get("source_sha256", {}).items():
                require(sha256(root/src) == h, f"Native action/source hash diverged: {src}")
            for name in TENSOR_NAMES:
                a = native_tensor(act["native_tensors"][name]); require(a.shape == reg[name].shape, name+" shape differs")
                require(np.max(np.abs(a - reg[name])) < 2e-12, name+" changed from archived source registry")
            flag = "fresh_installed_Matchete_action_validated"
        else:
            require(source_only, f"Missing fresh native action: {native_path}")
    else:
        require(source_only, "--native-dir required for certification review; use explicit --source-only for archived-only evaluation")
    return bench,reg,flag,{str(p.relative_to(root)):sha256(p) for p in paths},str(native_path) if native_path else None


def read_flavour(obj):
    return np.asarray(obj["real"],dtype=float) + 1j*np.asarray(obj["imag"],dtype=float)


def all_portal_map(reg):
    # Full four-index least squares. No chosen portal identity coefficients.
    ident = np.einsum("ij,ab->ijab",np.eye(2),np.eye(2)).ravel()
    maps={}; projections={}
    for index in (1,2):
        p = reg[f"T3HiggsPortal{index}CGInv1"]
        q = reg[f"T3HiggsPortal{index}CGInv2"]
        basis = np.stack([p.ravel(),q.ravel()],axis=1)
        co,_,rank,_=np.linalg.lstsq(basis,ident,rcond=None)
        err = np.linalg.norm(basis@co-ident)/np.linalg.norm(ident)
        require(rank == 2 and err < 1e-13, f"S{index} no exact full identity decomposition")
        require(np.max(np.abs(co-np.array([0.5,math.sqrt(3)/2]))) < 1e-12,"Unexpected native portal basis identity")
        maps[f"S{index}"] = {"Inv1":float(co[0].real),"Inv2":float(co[1].real),"full_tensor_residual":float(err)}
        projections[f"S{index}"]={"Inv1":enc_matrix(p[1,1,:,:]),"Inv2":enc_matrix(q[1,1,:,:])}
    return maps,projections


def generator_data(root):
    path=root/"output/general_observable/lo_direct_one_loop/source_audits/matchete_generator_analysis.json"
    require(path.is_file(), "Missing archived Matchete generator basis report")
    data=json.loads(path.read_text(encoding="utf-8"))
    reps={e["dimension"]:e for e in data["representations"]}
    return path,reps


def expr_mat(data):
    return np.array([[expr(str(x)) for x in row] for row in data],complex)


def charge_contract(reg,dF,alpha,reps):
    """Read actual Matchete T3/U, select all external L[1] charged components."""
    v1,v2=reg["T3Y1CG"],reg["T3Y2CG"]
    t3f=(np.zeros((1,1),complex) if dF == 1 else expr_mat(reps[dF]["T3_in_original_matchete_basis"]))
    U=(np.eye(1,dtype=complex) if dF == 1 else expr_mat(reps[dF]["U_columns_are_weight_eigenvectors"]))
    weights=np.array([0] if dF == 1 else [1,0,-1],float)
    require(np.max(np.abs(U.conj().T@U-np.eye(dF)))<1e-13,"Native rep rotation not unitary")
    require(np.max(np.abs(U.conj().T@t3f@U-np.diag(weights)))<1e-13,"Native rep Cartan weights invalid")
    require(np.max(np.abs(U.conj().T@np.eye(dF)@U-np.eye(dF)))<1e-13,"Dirac hermitian mass metric invalid")
    # The real Cartesian triplet is real AS AN SU(2) irrep; charged YF makes F complex.
    # CConj[F_cart]=U* CConj[F_weight], vs F_cart=U F_weight.
    if dF == 1:
        a=v1[0][np.newaxis,:]
        b=v2[0][np.newaxis,:]
    else:
        a=U.conj().T @ v1[0,:,:]
        b=U.T @ v2[0,:,:]
    # Check matrix conversion: sum_cart CG_cart,f U*_{f,w} = U^dag_{w,f} CG_cart,f.
    # y2: sum_cart CG_cart,f U_{f,w} = U^T_{w,f} CG_cart,f.
    yf=(alpha+1)/2
    ys1=alpha/2;ys2=(alpha+2)/2
    qf=weights+yf
    qs1=np.array([ys1+.5,ys1-.5]);qs2=np.array([ys2+.5,ys2-.5])
    nonzero1=[]; nonzero2=[]
    for fi in range(dF):
        for si in range(2):
            if abs(a[fi,si])>1e-12:
                require(abs(-qf[fi]+qs1[si])<1e-12, "y1 CG violates electric charge")
                nonzero1.append({"F_weight_T3":int(weights[fi]),"F_Q":float(qf[fi]),"S1_index":si+1,"S1_Q":float(qs1[si]),"CG":enc_matrix(np.array([[a[fi,si]]]))})
            if abs(b[fi,si])>1e-12:
                require(abs(qf[fi]-qs2[si])<1e-12, "y2 CG violates electric charge")
                nonzero2.append({"F_weight_T3":int(weights[fi]),"F_Q":float(qf[fi]),"S2_index":si+1,"S2_Q":float(qs2[si]),"CG":enc_matrix(np.array([[b[fi,si]]]))})
    cm=reg["T3MixCG"][1,1,:,:]
    pairs=[]
    for i in range(2):
        for j in range(2):
            if abs(cm[i,j])>1e-12:
                require(abs(qs1[i]-qs2[j])<1e-12,"H0 mixing changes scalar electric charge")
                # y1 and y2 must select the SAME chirality-flipping Dirac-mass channel.
                common=[f for f in range(dF) if abs(a[f,i])>1e-12 and abs(b[f,j])>1e-12]
                pairs.append((i,j,common))
    require(len(pairs)==1 and pairs[0][:2]==(0,1) and pairs[0][2]==[0 if dF==1 else 1],
            f"Contributing Dirac charge channels unexpected: {pairs}")
    fidx=pairs[0][2][0]
    # Directly evaluate ONLY the unique charge-selected CG vertex pair.
    c1=a[fidx,0];c2=b[fidx,1]
    require(abs(c1-(1 if dF==1 else 3**(-.25)))<1e-12,"Source y1 native selected CG mismatch")
    require(abs(c2-(1 if dF==1 else -3**(-.25)))<1e-12,"Source y2 native selected CG mismatch")
    return c1,c2,{
        "SU2_T3_cart":enc_matrix(t3f),"SU2_cart_to_weight_U":enc_matrix(U),
        "dirac_mass_metric_cart":enc_matrix(np.eye(dF)),
        "dirac_mass_metric_weight":enc_matrix(U.conj().T@U),
        "symmetric_transpose_weight_metric_is_NOT_Dirac":enc_matrix(U.T@U),
        "all_y1_Lnu_components":nonzero1,"all_y2_Lnu_components":nonzero2,
        "mixing_pairs_S1_S2_zero_based": [[i,j] for i,j,_ in pairs],
        "sole_mass_loop_weight_T3":int(weights[fidx]),
        "sole_mass_loop_electric_charge":float(qf[fidx]),
        "mixing_electric_charge":float(qs1[0]),
        "number_mass_channels":len(pairs[0][2]),
        "source_mixing_CG_H0_S1_1_S2_2":enc_matrix(np.array([[cm[0,1]]]))
    }


def blocks_from_full_tensors(bench,reg,portal_coeff):
    """Full 4-complex-scalar quadratic Hessian, block organised by electric charge."""
    d=bench["base_state"]["ordinary"]
    v=float(bench["numerical"]["vev_gev"])
    alpha=int(bench["representation"]["alpha"])
    lam=complex(d["lambdaT3"]["real"],d["lambdaT3"]["imag"])
    # Native LInt mixing is +lambdaT3 H H S1 Bar[S2] C + H.c.
    # V=-L -> A[S2,S1]= -lambdaT3*v^2/2 * C[H0,H0,S1,S2bar].
    mix=reg["T3MixCG"][1,1,:,:]
    p1=sum(portal_coeff["S1"][k]*reg[f"T3HiggsPortal1CGInv{k}"][1,1,:,:] for k in (1,2))
    p2=sum(portal_coeff["S2"][k]*reg[f"T3HiggsPortal2CGInv{k}"][1,1,:,:] for k in (1,2))
    # Native quadratic source: V_mass = mS1Sq S1dag S1 + mS2Sq S2dag S2.
    a=np.zeros((4,4),complex)
    a[:2,:2] = d["mS1Sq"]*np.eye(2)-.5*v*v*p1
    a[2:,2:] = d["mS2Sq"]*np.eye(2)-.5*v*v*p2
    # The two scalar pairs are (S1[0/1],S2[0/1]).
    a[2:,:2]=-.5*v*v*lam*mix.T
    a[:2,2:]=a[2:,:2].conj().T
    require(np.max(abs(a-a.conj().T))<1e-7,"Scalar Hessian not Hermitian")
    ys1=alpha/2;ys2=(alpha+2)/2
    charges=[ys1+.5,ys1-.5,ys2+.5,ys2-.5]
    for i in range(4):
        for j in range(4):
            if abs(a[i,j])>1e-8 and i!=j:
                require(charges[i]==charges[j], "Scalar mass block violates electric charge")
    blocks={}
    for q in sorted(set(charges)):
        indices=[i for i in range(4) if charges[i]==q]
        vals=np.linalg.eigvalsh(a[np.ix_(indices,indices)])
        require(np.min(vals)>0,"Non-positive scalar squared-mass block")
        blocks[str(q)]={"field_indices_1_based":[k+1 for k in indices],
                        "charge_Q":float(q),"matrix_GeV2":enc_matrix(a[np.ix_(indices,indices)]),
                        "eigenvalues_GeV2":vals.tolist()}
    return a,blocks


def finite_difference(M, s0,s1):
    M2=mp.mpf(str(M))**2
    lo,hi=map(lambda t:mp.mpf(str(float(t))),(s0,s1))
    require(M>0 and lo>0 and hi>lo,"Nonpositive loop mass")
    delta=hi-lo
    integral=-mp.quad(lambda x:mp.log1p((1-x)*delta/(x*M2+(1-x)*lo))/delta,
                      [mp.mpf('0'),mp.mpf('0.5'),mp.mpf('1')])
    def bfin(s):
        return 1+s*mp.log(s/M2)/(M2-s)
    analytic=(bfin(hi)-bfin(lo))/delta
    require(abs(integral-analytic)<mp.mpf('1e-40')*max(abs(integral),mp.mpf('1e-30')),
            "Independent parameter integral/analytic B0 finite difference mismatch")
    return float(integral),float(analytic)


def make_matrix(bench,reg,c1,c2,portal,masses, y1=None,y2=None,lam=None):
    d=bench["base_state"]["ordinary"]
    y1=read_flavour(d["y1"]) if y1 is None else y1
    y2=read_flavour(d["y2"]) if y2 is None else y2
    if lam is None:lam=complex(d["lambdaT3"]["real"],d["lambdaT3"]["imag"])
    local=json.loads(json.dumps(bench))
    local["base_state"]["ordinary"]["lambdaT3"]={"real":lam.real,"imag":lam.imag}
    h,blocks=blocks_from_full_tensors(local,reg,portal)
    # Mixed charge block indexed (S1[1], S2[2]) = (0,3) in full basis.
    x=h[np.ix_([0,3],[0,3])]
    eig,rot=np.linalg.eigh(x)
    require(np.all(eig>0),"Mixed scalar eigenvalues non-positive")
    pole=abs(sum(rot[1,k]*np.conj(rot[0,k]) for k in range(2)))
    require(pole < 1e-13,"B0 UV pole fails cancellation")
    A21=x[1,0]
    require(abs(A21+.5*bench["numerical"]["vev_gev"]**2*lam*reg["T3MixCG"][1,1,0,1])<1e-7,
            "Portal/mixing source matrix orientation inconsistent")
    g1=-np.conj(y1*c1);g2=-np.conj(y2*c2)
    m=np.zeros((3,3),complex)
    each=[]
    for r,M in enumerate(masses):
        D,D2=finite_difference(M,eig[0],eig[1])
        coefficient=M*A21*D/(16*math.pi**2)*1e9  # GeV -> eV
        m+=coefficient*(np.outer(g1[:,r],g2[:,r])+np.outer(g2[:,r],g1[:,r]))
        each.append({"flavour_index":r+1,"M_GeV":float(M),
                     "finite_D_GeV_minus2":D,"analytic_D_GeV_minus2":D2,
                     "coefficient_eV_over_two_Yukawas":enc_matrix(np.array([[coefficient]]))})
    require(np.all(np.isfinite(m)),"Mass matrix contains NaN/Inf")
    return m, {"scalar_full_basis":"S1[1], S1[2], S2[1], S2[2]",
               "scalar_blocks":blocks,"mixed_mass_matrix_GeV2":enc_matrix(x),
               "scalar_mixed_eigenvalues_GeV2":eig.tolist(),
               "A_S2_2_S1_1_GeV2":enc_matrix(np.array([[A21]])),
               "B0_pole_offdiagonal_sum":float(pole),
               "loop_flavour_terms":each}


def takagi(m):
    values=np.linalg.svd(m,compute_uv=False)
    e,V=np.linalg.eigh(m.conj().T@m)
    order=np.argsort(e)[::-1];V=V[:,order]
    z=V.T @ m @ V
    theta=np.exp(-.5j*np.angle(np.diag(z)))
    U=V*theta[np.newaxis,:]
    residual=np.linalg.norm(U.T@m@U-np.diag(values))/max(np.linalg.norm(m),1e-300)
    # Degenerate SVD directions are not relevant for the nondegenerate fitted points.
    require(residual<3e-8,f"Takagi reconstruction residual {residual}")
    return values.tolist(),float(residual)


def negative_controls(bench,reg,c1,c2,portal,masses,reference):
    d=bench["base_state"]["ordinary"]
    y1=read_flavour(d["y1"]);y2=read_flavour(d["y2"])
    lam=complex(d["lambdaT3"]["real"],d["lambdaT3"]["imag"])
    # Genuine complex deformation, not an overall phase of the real matrix.
    y1[0,1]*=np.exp(.37j);y1[1,2]*=np.exp(-.19j)
    y2[2,0]*=np.exp(.41j);y2[0,2]*=np.exp(-.23j)
    lam*=np.exp(.29j)
    m,_=make_matrix(bench,reg,c1,c2,portal,masses,y1,y2,lam)
    u1,u2,up=.21,-.17,.43
    mr,_=make_matrix(bench,reg,c1,c2,portal,masses,y1*np.exp(1j*u1),
                     y2*np.exp(1j*u2),lam*np.exp(1j*up))
    phase=np.exp(1j*(up-u1-u2))
    cov=float(np.linalg.norm(mr-phase*m)/np.linalg.norm(m))
    require(cov<3e-12,"Complex source-coupling phase transport failed")
    wrong,_=make_matrix(bench,reg,c1,c2,portal,masses,np.conj(y1),y2,lam)
    neg=float(np.linalg.norm(wrong-m)/np.linalg.norm(m))
    require(neg>1e-4,"Wrong Yukawa conjugation negative control failed")
    zero,_=make_matrix(bench,reg,c1,c2,portal,masses,y1,y2,0j)
    require(np.linalg.norm(zero)==0,"lambdaT3=0 should give exactly zero")
    # A wrong charged scalar (S1[2],S2[1]) has no source H0 mixing.
    require(abs(reg["T3MixCG"][1,1,1,0])<1e-14,"Wrong-charge scalar mixing did not vanish")
    masses0=np.linalg.svd(m,compute_uv=False)
    masses1=np.linalg.svd(mr,compute_uv=False)
    spin=float(np.max(abs(masses1-masses0))/max(masses0))
    require(spin<3e-12,"Takagi masses not rephasing invariant")
    # wrong-sign relative vertex changes the matrix sign, although singular masses
    # are unchanged; report this as a convention test, NOT physical sign resolution.
    wrong_sign=-reference
    sign_change=float(np.linalg.norm(wrong_sign-reference)/np.linalg.norm(reference))
    return {"genuine_complex_point_matrix_eV":enc_matrix(m),
            "global_source_rephasing_theta_y1_y2_lambda":[u1,u2,up],
            "expected_overall_matrix_phase_rad":up-u1-u2,
            "complex_phase_covariance_relative_residual":cov,
            "Takagi_rephasing_relative_residual":spin,
            "wrong_y1_conjugation_relative_change":neg,
            "lambdaT3_zero_mass_norm_eV":float(np.linalg.norm(zero)),
            "unmixed_charge_channel_CG_H0":float(abs(reg["T3MixCG"][1,1,1,0])),
            "wrong_relative_vertex_sign_matrix_change":sign_change,
            "wrong_relative_vertex_sign_changes_Takagi_masses":False}


def compare_frozen(m, reference, mid):
    """Optional comparison performed *after* the independent numerical matrix.

    Accept the repository's official `github_provisional_matrix_reference.json`
    list-of-records format, plus a model-keyed (model -> mnu_eV) mapping.
    Never fit, rescale, rephase or change coefficients to improve agreement.
    """
    if reference is None:
        return None, "not_requested"
    data = json.loads(reference.read_text(encoding="utf-8"))
    recs = data.get("models", data)
    if isinstance(recs, list):
        matching = [e for e in recs if e.get("model") == mid]
        require(len(matching) == 1, f"Missing/duplicate frozen model {mid}")
        rec = matching[0]
    elif isinstance(recs, dict):
        rec = recs.get(mid)
    else:
        rec = None
    if rec is None:
        return None, "model_missing_in_reference"
    factor = 1.
    if isinstance(rec, dict) and "mnu_3x3_gev" in rec:
        rec = rec["mnu_3x3_gev"]
        factor = 1e9
    else:
        for key in ("mnu_eV", "mnu_eigenbasis_ev", "mnu_matrix_ev", "mnu_ev",
                    "mass_matrix_eV", "mnu_matrix_eV", "matrix_ev"):
            if isinstance(rec, dict) and key in rec:
                rec = rec[key]
                break
    if isinstance(rec, dict) and "real" in rec and "imag" in rec:
        f = read_flavour(rec) * factor
    elif isinstance(rec, list) and np.asarray(rec).shape == (3, 3):
        f = np.asarray(rec, complex) * factor
    else:
        return None, "unknown_frozen_schema_no_coefficient_fit"
    require(f.shape == (3, 3) and np.linalg.norm(f) > 0,
            "Frozen reference lacks a nonzero 3x3 matrix")
    return float(np.linalg.norm(m-f)/np.linalg.norm(f)), "provisional_comparison_after_independent_loop_not_certification"


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo-root",type=Path,default=Path("."))
    ap.add_argument("--native-dir",type=Path,help="Fresh RunPhysicalCBCDNative.wl outputs")
    ap.add_argument("--source-only",action="store_true",help="Allow archived GitHub-native tensors; NOT newly executed Matchete")
    ap.add_argument("--output-dir",required=True,type=Path,help="Must NOT already exist")
    ap.add_argument("--frozen-reference",type=Path,help="Comparison ONLY after independent evaluation")
    a=ap.parse_args()
    root=a.repo_root.resolve();out=a.output_dir.resolve()
    require(not out.exists(),"Output directory already exists: "+str(out))
    require(root.is_dir(),"Repo root not found")
    require(a.native_dir is not None or a.source_only,"Require --native-dir or explicit --source-only")
    if a.frozen_reference:
        require(a.frozen_reference.is_file(),"Frozen comparison explicitly requested but missing")
    rows=[];models={}
    genpath, reps=generator_data(root)
    for mid,(cls,alpha,dF,suf) in MODELS.items():
        rec={"model":mid,"native_evidence":"blocked","scalar_normalisation":"blocked",
             "fermion_normalisation":"blocked","independent_matrix_eV":None,
             "Takagi_masses_eV":None,"comparison_residual":None,
             "comparison_status":"not_computed","remaining_blockers":[]}
        try:
            bench,reg,status,hashes,npath=load_inputs(root,mid,a.native_dir,a.source_only)
            rec["native_evidence"]=status
            require(dF == 1 or reps[dF]["status"]=="generator_algebra_verified", "Generator witness incomplete")
            c1,c2,charges=charge_contract(reg,dF,alpha,reps)
            portal_map,projectors=all_portal_map(reg)
            d=bench["base_state"]["ordinary"]
            MF=read_flavour(d["MF"])
            require(np.linalg.norm(MF-np.diag(np.diag(MF)))<1e-9 and
                    np.max(abs(np.imag(np.diag(MF))))<1e-12 and np.min(np.real(np.diag(MF)))>0,
                    "Benchmark masses must be real positive diagonal")
            masses=np.real(np.diag(MF)).tolist()
            require(len(masses)==3,"Three benchmark Dirac masses required")
            require(read_flavour(d["y1"]).shape==read_flavour(d["y2"]).shape==(3,3),
                    "Yukawa flavour sizes must match native action")
            # Model-specific source couplings; identity maps solved from FULL tensors.
            zero={"S1":{1:0.,2:0.},"S2":{1:0.,2:0.}}
            identity={"S1":{i:d["lambdaH1"]*portal_map["S1"][f"Inv{i}"] for i in (1,2)},
                      "S2":{i:d["lambdaH2"]*portal_map["S2"][f"Inv{i}"] for i in (1,2)}}
            zmat,zextra=make_matrix(bench,reg,c1,c2,zero,masses)
            imat,iextra=make_matrix(bench,reg,c1,c2,identity,masses)
            mcommon,cextra=make_matrix(bench,reg,c1,c2,zero,[masses[0]]*3)
            icommon,icextra=make_matrix(bench,reg,c1,c2,identity,[masses[0]]*3)
            spectra={}
            for key,mass in (("zero_diagonal",zmat),("identity_diagonal",imat),
                             ("zero_native_common",mcommon),("identity_native_common",icommon)):
                values,tkres=takagi(mass)
                sym=float(np.linalg.norm(mass-mass.T)/max(np.linalg.norm(mass),1e-300))
                require(sym<1e-13,"Majorana neutrino matrix not symmetric")
                spectra[key]={"matrix_eV":enc_matrix(mass),"takagi_descending_eV":values,
                              "symmetry_residual":sym,"takagi_reconstruction_residual":tkres}
            controls=negative_controls(bench,reg,c1,c2,identity,masses,imat)
            resid,cstatus=None,"deferred_until_all_four_independent_matrices_computed"
            rec.update({
                "scalar_normalisation":"native_CG_full_Hessian_positive_all_charge_blocks" if status.startswith("fresh") else "source_registry_conditional_full_Hessian_positive",
                "fermion_normalisation":"native_Dirac_minus_one_flavour_diagonal_verified_absolute_GammaCC_pending" if status.startswith("fresh") else "source_Dirac_minus_one_conditioned_on_live_action",
                "independent_matrix_eV":enc_matrix(imat),"Takagi_masses_eV":spectra["identity_diagonal"]["takagi_descending_eV"],
                "comparison_residual":resid,"comparison_status":cstatus,
                "remaining_blockers":([
                    "Fresh local installed Matchete four-model action not run"] if not status.startswith("fresh") else []) + [
                    "Track A owns absolute common GammaCC spinor phase certification",
                    "Benchmark identity portal is a selected native Inv1/Inv2 direction, not a canonical native coupling identification",
                    "Benchmark 100/110/120 TeV assignment is not numerically inserted into Matchete native action"],
            })
            models[mid]={"C0": {"representation":[2,2,dF],"hypercharges":[alpha/2,(alpha+2)/2,(alpha+1)/2],
                                 "bench_mass_GeV":masses,"source_sha256":hashes,"fresh_native_path":npath,
                                 "archived_generator_reference":str(genpath.relative_to(root))},
                         "C1":{"source_tensor_selected_y1":enc_matrix(np.array([[c1]])),
                                "source_tensor_selected_y2":enc_matrix(np.array([[c2]])),
                                "charge_contract":charges,"full_portal_identity_solutions":portal_map,
                                "native_portal_projectors":projectors,
                                "selected_identity_couplings":identity,
                                "source_zero_portal_couplings":zero,
                                "triplet_is_complex_charged_Dirac_not_Majorana":True},
                         "C2":{"points":spectra,"zero_diagonal_scalar":zextra,
                                "identity_diagonal_scalar":iextra,
                                "zero_native_common_scalar":cextra,
                                "identity_native_common_scalar":icextra},
                         "C3":{"controls":controls,"frozen_comparison_relative_Frobenius":resid,
                                "frozen_status":cstatus},"C4":{"row":rec}}
        except Exception as exc:
            rec["remaining_blockers"].append(str(exc))
            models[mid]={"status":"blocked","reason":str(exc),"traceback":traceback.format_exc(limit=4),"C4":{"row":rec}}
        rows.append(rec)
    # FROZEN COMPARISON GATE: no frozen bytes are opened until ALL FOUR independent
    # model calculations/negative controls have finished, even if one was blocked.
    for rec in rows:
        if rec["independent_matrix_eV"] is None:
            continue
        mid=rec["model"]
        matrix=read_flavour(rec["independent_matrix_eV"])
        comparison,status=compare_frozen(matrix,a.frozen_reference,mid)
        rec["comparison_residual"]=comparison
        rec["comparison_status"]=status
        models[mid]["C3"]["frozen_comparison_relative_Frobenius"]=comparison
        models[mid]["C3"]["frozen_status"]=status
    summary={"format":"track-c-bc-independent-source-loop-v1",
             "timestamp_local":datetime.now().isoformat(),
             "origin":"source-only numerical calculation unless native fresh witness supplied",
             "installed_Matchete_executed_in_this_python_session":False,
             "model_order":list(MODELS),"source_only":a.source_only,
             "frozen_ever_read_before_direct_evaluation":False,
             "common_spinor_absolute_convention_certified":False,
             "C4_four_rows":rows,"models":models}
    out.mkdir(parents=True,exist_ok=False)
    write_json(out/"track_c_four_models.json",summary)
    with (out/"track_c_four_rows.csv").open("x",newline="",encoding="utf-8") as f:
        cols=list(rows[0].keys())
        writer=csv.DictWriter(f,fieldnames=cols);writer.writeheader()
        for row in rows:
            writer.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in row.items()})
    succeeded=sum(bool(r["independent_matrix_eV"]) for r in rows)
    print(f"Track C: {succeeded}/4 source-derived independent mass matrices; native evidence status in JSON")
    print(out/"track_c_four_rows.csv")
    return 0 if succeeded==4 else 1


if __name__ == "__main__":
    sys.exit(main())
