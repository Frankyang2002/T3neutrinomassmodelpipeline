"""Workstream 4: exact Matchete/RGBeta basic Higgs-portal embedding.

Compares registered Matchete CGs after the native real-coupling
1/2 PlusHc projection against delta(Hbar,H) delta(Sbar,S). The embedding
may be non-unique; a canonical representative plus nullspace are reported.
Never mistakes a dependent Hermitian operator family for a unique coupling
map, or infers the full physical source-to-Weyl normalisation.
"""
from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import json
import re
from typing import Any

import numpy as np
import sympy as sp
from sympy.parsing.mathematica import parse_mathematica

# Independently reproduce the already source-verified class inventory. Importing
# the production numerical engine here would defeat independent auditing.
CLASSES = {
    "A": (1,3,2,(-4,-2,0)),"B": (2,2,1,(-3,-1,1)),
    "C": (2,2,3,(-3,-1,1)),"D": (3,1,2,(-2,0,2)),
    "E": (3,3,2,(-4,-2,0,2)),
}


def all_models():
    for cls,(d1,d2,df,alphas) in CLASSES.items():
        for a in alphas:
            yield f"{cls}_{a:+d}", (
                f"T3_dS1_{d1}_dS2_{d2}_dF_{df}_alpha_{'m' if a<0 else 'p'}{abs(a)}")


def decode_matrix(obj:Any,name:str)->np.ndarray:
    try:
        if isinstance(obj,dict):
            a=np.asarray(obj['real'],dtype=float)+1j*np.asarray(obj['imag'],dtype=float)
        else:
            a=np.asarray(obj,dtype=complex)
    except (KeyError,TypeError,ValueError) as exc:
        raise PortalDictionaryError(f"Invalid complex matrix {name}") from exc
    require(a.ndim==2 and a.shape[0]==a.shape[1] and np.all(np.isfinite(a)),
            f"Not a finite square matrix: {name}")
    return a

NATIVE_FORMAT = "t3-ws4-native-portal-registered-tensors-five-classes-v1"
FORMAT = "t3-ws4-registered-native-portal-singlet-dictionary-v2"
FIVE = {"A": (1,3,2,0),"B": (2,2,1,-1),"C": (2,2,3,-1),
        "D": (3,1,2,0),"E": (3,3,2,0)}
SOURCE_REL = (
    "Lagrangian/model/T3ModelCatalog.wl",
    "Lagrangian/model/LagrangianBuilder.wl",
    "Lagrangian/model/SU2Invariants.wl",
    "Lagrangian/model/T3Fields.wl",
    "Lagrangian/interactions/T3Topology.wl",
    "Lagrangian/interactions/ScalarPotential.wl",
)


class PortalDictionaryError(ValueError):
    pass


def require(ok: bool, why: str) -> None:
    if not ok:
        raise PortalDictionaryError(why)


def read(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding="utf-8-sig"))
    require(isinstance(obj,dict),f"Not a JSON object: {path}")
    return obj


def hashfile(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def source_contract(root: Path, native: dict) -> dict:
    """Require matching native source bytes and independently inspect RGBeta.

    Guards check exact relevant source operators. They are not a full parser.
    """
    imported = native.get("native_source_sha256")
    require(isinstance(imported,dict) and set(imported)==set(SOURCE_REL),
            "Native JSON lacks complete source SHA256 evidence")
    hashes = {}
    for name in SOURCE_REL:
        p=root/name
        require(p.is_file(),f"Source missing: {p}")
        h=hashfile(p)
        require(re.fullmatch(r"[0-9a-fA-F]{64}",str(imported[name])) is not None and
                str(imported[name]).lower()==h,
                f"Native registered CG JSON is stale against {name}; re-export")
        hashes[name]=h
    rgbpath=root/"RGE/running/rgbeta/T3RGBetaModel.wl"
    buildpath=root/"scripts/BuildT3DirectOneLoopInputs.py"
    for p in (rgbpath,buildpath):
        require(p.is_file(),f"Missing source: {p}")
    rgb=rgbpath.read_text(encoding="utf-8")
    src=(root/SOURCE_REL[-1]).read_text(encoding="utf-8")
    builder=buildpath.read_text(encoding="utf-8")
    squash=lambda s:re.sub(r"\s+","",s)
    rb=squash(rgb);sc=squash(src);bu=squash(builder)
    for needle in (
        "T3RGBetaPairInvariant[1,i_,j_]:=1;",
        "T3RGBetaPairInvariant[2,i_,j_]:=del[SU2L@fund,i,j];",
        "T3RGBetaPairInvariant[3,i_,j_]:=del[SU2L@S2,i,j];",
        "AddQuartic[lambdaH1,{Bar@H,H,Bar@S1,S1},GroupInvariant->(del[SU2L@fund,#1,#2]*T3RGBetaPairInvariant[dS1,#3,#4]&)",
        "AddQuartic[lambdaH2,{Bar@H,H,Bar@S2,S2},GroupInvariant->(del[SU2L@fund,#1,#2]*T3RGBetaPairInvariant[dS2,#3,#4]&)",
    ):
        require(needle in rb,f"RGBeta source operator guard changed: {needle[:70]}")
    require("1/2PlusHc[coupling[]fieldProductIf[cg===None,1,CG[cg,labels]]]" in sc,
            "Matchete source portal 1/2 PlusHc convention changed")
    for needle in ("-float(o['lambdaH1'])*v*v/2","-float(o['lambdaH2'])*v*v/2"):
        require(needle in bu,f"Provisional input portal sign changed: {needle}")
    hashes[str(rgbpath.relative_to(root))]=hashfile(rgbpath)
    hashes[str(buildpath.relative_to(root))]=hashfile(buildpath)
    return {"status":"native_source_bytes_match_and_portal_forms_guarded",
            "source_sha256":hashes,
            "signed_portal_Lagrangian_source":"+lambdaHiInvN (Hdag H)(Si dag Si) times native CG (1/2 PlusHc)",
            "RGBeta_basic_portal_group_invariant":"delta(Hbar,H) delta(Sbar,S)"}


def parse_exact(x: str) -> sp.Expr:
    require(isinstance(x,str) and len(x)<=500,"Invalid native CG component literal")
    try:
        e=sp.simplify(parse_mathematica(x))
    except Exception as exc:
        raise PortalDictionaryError(f"Cannot parse Mathematica CG number {x!r}: {exc}") from exc
    require(not e.free_symbols and not e.atoms(sp.Function) and
            not e.atoms(sp.Float) and
            bool(e.is_number) and bool(e.is_finite),
            f"CG component is not a finite exact number: {x}")
    return e


def portal_tensor(record: dict,d: int,which: int,index: int) -> sp.Matrix:
    expected=[2,2]+([d,d] if d>1 else [])
    require(record.get("dimensions")==expected,
            f"Unexpected registered portal tensor dims for H{which}Inv{index}")
    require(record.get("scalar")==which and record.get("invariant")==index and
            record.get("coupling")==f"lambdaH{which}Inv{index}" and
            record.get("name")==f"HiggsPortal{which}Inv{index}" and
            record.get("cg_symbol")==f"T3HiggsPortal{which}CGInv{index}" and
            record.get("tensor_accessor")=="Matchete`GetCGTensor" and
            record.get("registered_in_assembled_interaction") is True,
            f"Native portal registration evidence invalid for H{which}Inv{index}")
    values=record.get("components_inputform")
    require(isinstance(values,list) and len(values)==4*d*d,
            f"Incorrect number of CG elements for H{which}Inv{index}")
    entries=[parse_exact(v) for v in values]
    # Flat tensor uses Mathematica row-major Flatten[{hBar,h,sBar,s}].
    # Rows are (hBar,h); columns are (sBar,s); works also when d=1.
    out=sp.Matrix(4,d*d,entries)
    require(out!=sp.zeros(4,d*d),"Native registered tensor is identically zero")
    return out


def identity_tensor(d: int) -> sp.Matrix:
    return sp.Matrix(4,d*d,[int(i==j and a==b)
        for i in range(2) for j in range(2)
        for a in range(d) for b in range(d)])


def hermitian_portal_tensor(tensor: sp.Matrix, d: int) -> sp.Matrix:
    """Effective tensor of +1/2 PlusHc[lambda*CG*Hbar H Sbar S].

    Hermitian conjugation swaps barred/unbarred component positions. This
    cannot be ignored when an invariant CG basis is not itself Hermitian.
    """
    return sp.Matrix(4,d*d,[sp.simplify((tensor[2*i+j,d*a+b] +
        sp.conjugate(tensor[2*j+i,d*b+a]))/2)
        for i in range(2) for j in range(2)
        for a in range(d) for b in range(d)])


def solve_singlet_basis(tensors: list[sp.Matrix], d: int, *,
                        label: str = "unlabelled") -> dict:
    """Embed the RGBeta singlet in the REAL Hermitian native portal span.

    Native ScalarPotential.wl defines *real* lambdaHiInvN and includes
    1/2 PlusHc. Consequently independently registered raw complex CGs need
    not remain independent after Hermitian projection. Rank deficiency is
    not evidence of a unique coupling dictionary. In that case, retain the
    exact nullspace, and give a deterministic pivot-column representative
    while explicitly marking its coefficients NON-UNIQUE.

    Never accept a rank-deficient RAW registered basis, a missing RGBeta
    singlet, complex coefficients for real couplings, or an inexact identity.
    """
    require(1 <= len(tensors) <= 2,
            f"{label}: unexpected native portal basis size")
    require(all(t.shape == (4, d*d) for t in tensors),
            f"{label}: inconsistent native portal tensor dimensions")
    n = len(tensors)
    raw = sp.Matrix.hstack(*[sp.Matrix(t).reshape(4*d*d, 1)
                              for t in tensors])
    raw_rank = raw.rank()
    require(raw_rank == n,
            f"{label}: RAW Matchete registered CG tensors already dependent "
            f"(raw_rank={raw_rank}, count={n}); cannot trust native basis")

    projected = [hermitian_portal_tensor(t, d) for t in tensors]
    complex_matrix = sp.Matrix.hstack(*[
        sp.Matrix(t).reshape(4*d*d, 1) for t in projected])
    target = sp.Matrix(identity_tensor(d)).reshape(4*d*d, 1)

    # Couplings are real (DefineDynamicRealCoupling), so the identity must
    # hold over real coefficients, not just over the complex-number field.
    real_matrix = complex_matrix.applyfunc(lambda z: sp.simplify(sp.re(z)))
    imaginary_matrix = complex_matrix.applyfunc(lambda z: sp.simplify(sp.im(z)))
    system = real_matrix.col_join(imaginary_matrix)
    rhs = target.applyfunc(lambda z: sp.simplify(sp.re(z))).col_join(
        target.applyfunc(lambda z: sp.simplify(sp.im(z))))
    rank = system.rank()
    augmented_rank = system.row_join(rhs).rank()
    zero_projected = [k + 1 for k in range(n)
        if all(sp.simplify(v) == 0 for v in system[:, k])]
    require(rank > 0,
            f"{label}: all native real-coupling Hermitian portal operators "
            f"vanish (raw_rank={raw_rank}, projected_rank=0, count={n})")
    require(augmented_rank == rank,
            f"{label}: RGBeta basic singlet is OUTSIDE native real-Hermitian "
            f"portal span (raw_rank={raw_rank}, projected_rank={rank}, "
            f"augmented_rank={augmented_rank}, count={n}, "
            f"zero_projected_indices={zero_projected}). "
            "Do not infer a source coupling dictionary; inspect native "
            "tensor axes and PlusHc conjugation")

    # Pivot columns are an exact independent subfamily. Setting all other
    # coefficients to zero is a reproducible representative of the affine
    # solution space, NOT a unique physical dictionary in the degenerate case.
    pivots = tuple(system.rref()[1])
    require(len(pivots) == rank,
            f"{label}: inconsistent exact pivot rank")
    reduced = system[:, list(pivots)]
    chosen, free = reduced.gauss_jordan_solve(rhs)
    require(len(free) == 0,
            f"{label}: independent pivot solve retained an unexpected freedom")
    coefficients = [sp.Integer(0)] * n
    for col, value in zip(pivots, chosen):
        coefficients[col] = sp.simplify(value)
    require(all(z.is_real is True for z in coefficients),
            f"{label}: real native coupling basis yielded a non-real coefficient")
    residual = complex_matrix * sp.Matrix(coefficients) - target
    require(all(sp.simplify(z) == 0 for z in residual),
            f"{label}: exact native vs RGBeta operator residual is nonzero")

    null_vectors = system.nullspace()
    require(len(null_vectors) == n - rank,
            f"{label}: exact nullspace rank inconsistent")
    for v in null_vectors:
        require(all(sp.simplify(z) == 0 for z in system*v),
                f"{label}: native nullspace witness failed")
    unique = (rank == n)
    return {
        "status": ("exact_full_four_index_operator_identity_unique" if unique
                   else "exact_full_four_index_operator_identity_nonunique"),
        "label": label,
        "dimension": d,
        "raw_registered_CG_rank_over_complex": raw_rank,
        "real_Hermitian_operator_rank": rank,
        "native_invariant_count": n,
        "rank": rank,
        "augmented_real_rank": augmented_rank,
        "unique_coupling_dictionary": unique,
        "free_real_coupling_directions": n-rank,
        "native_invariant_redundancy_after_PlusHc": not unique,
        "identically_zero_Hermitian_projection_indices_one_based": zero_projected,
        "independent_pivot_indices_one_based": [k+1 for k in pivots],
        "nullspace_basis_real_coupling_shifts": [
            [str(sp.simplify(z)) for z in vec] for vec in null_vectors],
        "coefficients_native_lambdaHiInvN_per_RGBeta_lambdaHi": [
            str(z) for z in coefficients],
        "coefficient_approx": [float(sp.N(z, 17)) for z in coefficients],
        "coefficient_interpretation": (
            "unique exact operator identity" if unique else
            "ONE canonical pivot-column representative; infinitely many "
            "equivalent native real-coupling assignments, related by nullspace"),
        "symbolic_identity_residual_zero": True,
        "tensor_relation": (
            "sum_k c[k]*HermitianPart(T_native[k])[Hbar,H,Sbar,S]="
            "delta_H*delta_S; c[k] real"),
        "hermitian_conjugate_included": True,
        "potential_shift_at_H0_v_over_sqrt2":
            "-lambdaHi*v^2/2 * identity (conditional on source convention)",
        "not_certified": "RGBeta internal physical action normalisation, "
            "scalar free-mass convention, full source-to-Weyl dictionary",
    }


def validate_native_models(native:dict) -> dict:
    require(native.get("format")==NATIVE_FORMAT and
            native.get("status")=="five_classes_native_registered_portal_tensors_extracted",
            "Unexpected native portal CG report schema/status")
    env=native.get("restricted_environment",{})
    require(env.get("isolated_to_wolframscript_process") is True and
            all(env.get("effective",{}).get(n)=="0" for n in (
                "T3_EXPERIMENTAL_GENERAL_FERMIONS",
                "T3_EXPERIMENTAL_GENERAL_QUADRATICS",
                "T3_EXPERIMENTAL_GENERAL_QUARTICS",
                "T3_EXPERIMENTAL_FLAVOR_MASSES")),
            "Native Matchete build was not restricted")
    require(native.get("native_rgbeta_portal_map_certified") is False and
            native.get("matching_or_running_performed") is False,
            "Unexpected source report certification or matching status")
    rows=native.get("models")
    require(isinstance(rows,list) and len(rows)==5,"Expected five native classes")
    results={}
    for row in rows:
        cls=row.get("class")
        require(cls in FIVE and cls not in results,f"Missing or repeated native class {cls}")
        d1,d2,df,alpha=FIVE[cls]
        require(all(row.get(name)==v for name,v in (
                ("model",f"{cls}_{alpha:+d}"),("alpha",alpha),("dS1",d1),
                ("dS2",d2),("dF",df),("native_build_status","Success"),
                ("native_full_lagrangian_validation",True))),
                f"Native model metadata mismatch: {cls}")
        portals=row.get("portals")
        require(isinstance(portals,list) and len(portals) in (2,3,4),
                f"Malformed portal list: {cls}")
        results[cls]={}
        for which,d in ((1,d1),(2,d2)):
            selected=[r for r in portals if r.get("scalar")==which]
            require(len(selected)==(1 if d==1 else 2),
                    f"Wrong native portal count: {cls}, S{which}")
            ordered=sorted(selected,key=lambda r:r.get("invariant",-1))
            matrices=[portal_tensor(rec,d,which,k)
                      for k,rec in enumerate(ordered,1)]
            results[cls][f"S{which}"]=solve_singlet_basis(
                matrices,d,label=f"{cls}_{alpha:+d}/S{which}")
        require(sum(len([r for r in portals if r.get("scalar")==i])
                    for i in (1,2))==len(portals),"Unexpected portal record")
    require(set(results)==set(FIVE),"Five class portal export incomplete")
    return results


def _as_real(z:Any,name:str)->float:
    if isinstance(z,dict):
        z=complex(z["real"],z["imag"])
    c=complex(z)
    require(np.isfinite(c.real) and np.isfinite(c.imag) and abs(c.imag)<1e-12,
            f"{name} is not finite real")
    return float(c.real)


def frozen_inputs(root:Path,dictionary:dict) -> dict:
    """Verify all 16 unchanged provisional scalar blocks vs RGBeta singlet.

    Does NOT recompute UV one-loop matrices or assert physical matching.
    """
    rows=[]
    for model,key in all_models():
        cls=model[0]
        d1,d2,df,alphas=CLASSES[cls]
        inp_path=root/"output/general_observable/lo_direct_one_loop/provisional_inputs"/f"{model}.json"
        bench_path=root/"Numerical/benchmarks/fitted_reference"/f"{key}.json"
        require(inp_path.is_file(),f"Frozen provisional input missing: {inp_path}")
        require(bench_path.is_file(),f"Frozen reference missing: {bench_path}")
        inp=read(inp_path)
        bench=read(bench_path)
        require(inp.get("model")==model,f"Input model mismatch: {model}")
        rep=bench["representation"]
        require(all(rep[k]==v for k,v in (
            ("d_s1",d1),("d_s2",d2),("d_f",df),("alpha",int(model[2:])))),
            f"Frozen benchmark representation mismatch: {model}")
        ordinary=bench["base_state"]["ordinary"]
        v=_as_real(bench["numerical"]["vev_gev"],"vev")
        A=decode_matrix(inp["scalar_A_gev2"],"scalar_A")
        require(A.shape==(d1+d2,d1+d2),f"Scalar A dimension mismatch: {model}")
        norm=float(np.linalg.norm(A))
        err=[]
        for which,d,offset in ((1,d1,0),(2,d2,d1)):
            lam=_as_real(ordinary[f"lambdaH{which}"],f"lambdaH{which}")
            mass=_as_real(ordinary[f"mS{which}Sq"],f"mS{which}Sq")
            target=np.eye(d,dtype=complex)*(mass-lam*v*v/2)
            block=A[offset:offset+d,offset:offset+d]
            residual=float(np.linalg.norm(block-target)/max(1.,norm))
            require(residual<1e-11,
                    f"{model}: frozen diagonal scalar A disagrees with source basic portal")
            err.append(residual)
        rows.append({"model":model,"input_sha256":hashfile(inp_path),
            "benchmark_sha256":hashfile(bench_path),"scalar_block_residuals":err,
            "source_basic_portal_blocks_verified":True,
            "native_portal_dictionary_class":cls,
            "physical_source_to_Weyl_certified":False})
    return {"status":"all_16_restricted_basic_portal_blocks_agree_with_source_builder",
            "models_checked":len(rows),"max_relative_block_residual":max(max(r["scalar_block_residuals"]) for r in rows),
            "models":rows}


def audit(root:Path,native_path:Path,check_frozen:bool=True) -> dict:
    native=read(native_path)
    src=source_contract(root,native)
    dictionaries=validate_native_models(native)
    frozen=frozen_inputs(root,dictionaries) if check_frozen else None
    all_unique = all(
        row["unique_coupling_dictionary"]
        for cls in dictionaries.values() for row in cls.values())
    degeneracies = [{"class": cls, "scalar": scalar,
                     "native_invariant_count": row["native_invariant_count"],
                     "real_Hermitian_rank": row["real_Hermitian_operator_rank"],
                     "zero_projected_indices": row[
                         "identically_zero_Hermitian_projection_indices_one_based"]}
                    for cls,entries in dictionaries.items()
                    for scalar,row in entries.items()
                    if not row["unique_coupling_dictionary"]]
    return {"format":FORMAT,
        "status": ("native_portal_operator_identity_unique_conditional_physical" if all_unique
                   else "native_portal_operator_identity_exact_but_dictionary_nonunique"),
        "source_evidence":src,"native_report_sha256":hashfile(native_path),
        "native_portal_representation_classes_checked":5,
        "portal_operator_identity_solved_for_five_classes":True,
        "models_covered_by_class_transfer":16,
        "portal_dictionary_by_class":dictionaries,
        "frozen_input_audit":frozen,
        "scope":"Restricted basic Higgs portals only, not opt-in adjoint or holomorphic scalar operators",
        "basic_RGBeta_to_Matchete_operator_dictionary_verified":True,
        "exact_operator_embedding_verified":True,
        "unique_basic_portal_coupling_dictionary_verified":all_unique,
        "nonunique_operator_subspaces":degeneracies,
        "native_RGBeta_quantum_action_prefactors_independently_certified":False,
        "physical_source_to_Weyl_model_certified":0,
        "still_missing":[
            "native RGBeta internal scalar action prefactor vs Matchete physical fields",
            "Matchete GammaCC and RGBeta heavy fermion epsilon phase with both Yukawas and MF",
            "T3 scalar mixing CG coefficient and Higgs vev sign in common components",
            "mass-free-sector scalar quadratic action and all omitted allowed couplings",
            "hard one-loop matching of ALL heavy BSM fermions at muF",
            "integration of BOTH BSM scalars at muS into SMEFT"],
        "EFT_sequence":"UV -> all heavy fermions -> EFT1 (both BSM scalars and SM fermions) -> both scalars -> SMEFT",
        "matching_or_running_performed":False}
