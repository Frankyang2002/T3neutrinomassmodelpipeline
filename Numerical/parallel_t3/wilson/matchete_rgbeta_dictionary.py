"""Source-conditioned, *NOT physically certified*, T3 extra-quartic dictionary.

Input: the user's actual 50-slot Matchete CG export and the conditional RGBeta
symbolic-invariant reconstruction (v2).  This code does NOT interpret the v2
numbers as native RGBeta components, and never equates the two frames by name.

A reference Cartesian triplet CG is inferred from the explicit A_+0 Matchete
Hdag Hdag S1 S2 tensor.  We substitute it into the **source-exported** RGBeta
invariant InputForm expressions, separately applying a conjugate CG to a pair
of unbarred doublets when an internal triplet is contracted.  The resulting
component polynomials are expressed in a deliberately Matchete-compatible
*conditional* frame; no native RGBeta transformation is asserted.

For a common, explicit complex scalar-coordinate frame, solve
    V_M(x; m) = V_R(x; r),  r = A @ m,
where each non-self-adjoint quartic is included with its Hermitian conjugate.
The numeric A matrices are only for *extra* quartic monomial columns, with
potential prefactors set to unity conditionally, NOT complete native coupling
certificates.  They do not cover baseline quartics, masses, Yukawas, the
A(-4) anti-Hermitian portal, or GammaCC-to-Weyl spinor normalization.

No eval/exec or unrestricted SymPy parsing of untrusted JSON is used.
"""
from __future__ import annotations

import ast
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import product
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping

import numpy as np

MAT_FORMAT = "stage22-exact-tensors-v1"
RGB_CONDITIONAL_FORMAT = "t3-parallel-a-rgbeta-declared-invariants-conditional-frame-v2"
OUT_FORMAT = "t3-parallel-a-extra-quartic-polynomial-dictionary-conditional-v3"
ALL16 = frozenset((
    "A_-4", "A_-2", "A_+0", "B_-3", "B_-1", "B_+1",
    "C_-3", "C_-1", "C_+1", "D_-2", "D_+0", "D_+2",
    "E_-4", "E_-2", "E_+0", "E_+2",
))
BARED = {"Hdag": "Bar[H]", "S1dag": "Bar[S1]", "S2dag": "Bar[S2]"}
PARTIAL_SCOPE = "50 additional quartic invariant slots ONLY (not full scalar potential)"


class DictionaryEvidenceError(ValueError):
    """A required component or independently supplied convention is missing."""


def _ast_number(node: ast.AST) -> complex:
    if isinstance(node, ast.Constant) and type(node.value) in (int, float):
        return complex(node.value)
    if isinstance(node, ast.Name) and node.id == "I":
        return 1j
    if isinstance(node, ast.UnaryOp) and type(node.op) in (ast.UAdd, ast.USub):
        ans = _ast_number(node.operand)
        return -ans if isinstance(node.op, ast.USub) else ans
    if isinstance(node, ast.BinOp) and type(node.op) in (
        ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow,
    ):
        a, b = _ast_number(node.left), _ast_number(node.right)
        if isinstance(node.op, ast.Add): return a + b
        if isinstance(node.op, ast.Sub): return a - b
        if isinstance(node.op, ast.Mult): return a * b
        if isinstance(node.op, ast.Div):
            if b == 0: raise DictionaryEvidenceError("Zero denominator")
            return a / b
        if abs(b.imag) > 1e-13 or abs(b.real) > 16 or abs(a) > 1e8:
            raise DictionaryEvidenceError("Unsupported component exponent")
        return a ** b.real
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "sqrt":
        if node.keywords or len(node.args) != 1:
            raise DictionaryEvidenceError("Malformed numeric Sqrt")
        z = _ast_number(node.args[0])
        if abs(z.imag) > 1e-13 or z.real < 0:
            raise DictionaryEvidenceError("Only real nonnegative Sqrt is allowed")
        return complex(math.sqrt(z.real))
    raise DictionaryEvidenceError("Unsupported Matchete numerical component")


def matchete_number(item: Any) -> complex:
    if not isinstance(item, str) or len(item) > 150:
        raise DictionaryEvidenceError("Matchete component must be an InputForm string")
    source = item.strip().replace("Sqrt[", "sqrt(").replace("]", ")").replace("^", "**")
    try:
        out = _ast_number(ast.parse(source, mode="eval").body)
    except (SyntaxError, TypeError, OverflowError, ZeroDivisionError) as exc:
        raise DictionaryEvidenceError("Could not parse Matchete component") from exc
    if not np.isfinite(out):
        raise DictionaryEvidenceError("Nonfinite Matchete component")
    return out


def _mat_tensor(rec: Mapping[str, Any]) -> np.ndarray:
    if not isinstance(rec, Mapping) or not isinstance(rec.get("symbol"), str):
        raise DictionaryEvidenceError("Malformed Matchete tensor record")
    dims = rec.get("dimensions")
    if (not isinstance(dims, list) or len(dims) != 4 or
            any(type(d) is not int or d < 1 or d > 3 for d in dims)):
        raise DictionaryEvidenceError("Invalid Matchete tensor shape")
    values = rec.get("components")
    if not isinstance(values, list) or len(values) != math.prod(dims):
        raise DictionaryEvidenceError("Incomplete Matchete components")
    if rec.get("potential_prefactor") != "1":
        raise DictionaryEvidenceError("Changed Matchete raw invariant prefactor")
    tensor = np.array([matchete_number(v) for v in values], dtype=np.complex128).reshape(dims)
    if not np.any(abs(tensor) > 1e-12):
        raise DictionaryEvidenceError("Zero Matchete invariant")
    return tensor


def calibrated_triplet_cg(mat_doc: Mapping[str, Any]) -> tuple[np.ndarray, float]:
    """Orthonormal 2x2x3 Cartesian CG obtained from ACTUAL Matchete A_+0.

    This is a *coordinate convention* for comparing polynomial subspaces,
    not an independently measured or native RGBeta triplet basis.
    """
    entries = mat_doc["models"]["A_+0"]["tensors"]
    anchors = [x for x in entries if x.get("symbol") == "lambdaT3ExtraQ021010I1" and
               x.get("species") == ["Hdag", "Hdag", "S1", "S2"]]
    if len(anchors) != 1 or anchors[0]["dimensions"] != [2, 2, 1, 3]:
        raise DictionaryEvidenceError("Missing unique Matchete A_+0 triplet anchor")
    raw = _mat_tensor(anchors[0]).reshape(2, 2, 3)
    col_norms = np.linalg.norm(raw.reshape(4, 3), axis=0)
    scale = float(col_norms[0])
    if not 0.1 < scale < 10 or not np.allclose(col_norms, scale, rtol=1e-11):
        raise DictionaryEvidenceError("Matchete CG columns have inconsistent norms")
    cg = raw / scale
    gram = np.einsum("ija,ijb->ab", cg.conj(), cg, optimize=True)
    if not np.allclose(gram, np.eye(3), atol=1e-12) or not np.allclose(
            cg, cg.swapaxes(0, 1), atol=1e-12):
        raise DictionaryEvidenceError("Matchete 2x2x3 triplet anchor is not orthonormal/symmetric")
    return cg, scale


def _normalise_expr(native: str) -> ast.Expression:
    if not isinstance(native, str) or len(native) > 1600 or not native.strip().endswith("&"):
        raise DictionaryEvidenceError("Need source-exported RGBeta pure invariant InputForm")
    s = native.strip()[:-1].strip()
    s = re.sub(r"\bdel\[", "delta[", s)
    s = s.replace("SU2L[fund]", "fund").replace("SU2L[S2]", "triplet")
    s = re.sub(r"#([1-4])(?![0-9])", r"idx\1", s)
    s = s.replace("[", "(").replace("]", ")").replace("^", "**")
    try:
        return ast.parse(s, mode="eval")
    except SyntaxError as exc:
        raise DictionaryEvidenceError("Unsupported native invariant grammar") from exc


def _native_invariant(expr: str, component: tuple[int, int, int, int],
                      species: tuple[str, ...], cg: np.ndarray) -> complex:
    """Evaluate one coefficient in chosen Cartesian triplet component frame.

    The asymmetric CG treatment is explicit: delS2 acting on an *unbarred*
    doublet pair uses conjugate Cartesian CG.  It is needed to contract the
    B/C Hbar Hbar S2 S2 polynomial with the Matchete SU(2) singlet.  It is a
    CONDITIONAL convention, not an RGBeta-native confirmation.
    """
    node = _normalise_expr(expr)
    env: dict[str, Any] = {f"idx{i+1}": value for i, value in enumerate(component)}
    env.update({"SU2L": "SU2L", "fund": "fund", "triplet": "triplet"})

    def delta(rep, a, b):
        if rep not in ("fund", "triplet") or type(a) is not int or type(b) is not int:
            raise DictionaryEvidenceError("Malformed SU(2) delta")
        n = 2 if rep == "fund" else 3
        if not (1 <= a <= n and 1 <= b <= n): raise DictionaryEvidenceError("Delta out of range")
        return int(a == b)

    def eps(rep, a, b):
        if rep != "fund" or type(a) is not int or type(b) is not int or not (
            1 <= a <= 2 and 1 <= b <= 2
        ): raise DictionaryEvidenceError("Malformed SU(2) epsilon")
        return int((a, b) == (1, 2)) - int((a, b) == (2, 1))

    def cg_triplet(group, a, i, j):
        if group != "SU2L" or any(type(z) is not int for z in (a, i, j)) or not (
            1 <= a <= 3 and 1 <= i <= 2 and 1 <= j <= 2
        ): raise DictionaryEvidenceError("Malformed delS2 indices")
        return complex(cg[i-1, j-1, a-1])

    def triple(a, b, c):
        if any(type(z) is not int or not 1 <= z <= 3 for z in (a, b, c)):
            raise DictionaryEvidenceError("Malformed triplet invariant")
        value = 0j
        for i, j, k, l, mm, n in product((1, 2), repeat=6):
            value += (cg_triplet("SU2L", a, i, j) * cg_triplet("SU2L", b, k, l) *
                      cg_triplet("SU2L", c, mm, n) * eps("fund", i, k) *
                      eps("fund", j, mm) * eps("fund", l, n))
        return value

    primitives = {"delta": delta, "eps": eps, "delS2": cg_triplet,
                  "T3RGBetaS2TripleInvariant": triple}
    def compute(n: ast.AST):
        if isinstance(n, ast.Constant) and type(n.value) in (int, float):
            return n.value
        if isinstance(n, ast.Name) and n.id in env:
            return env[n.id]
        if isinstance(n, ast.UnaryOp) and type(n.op) in (ast.USub, ast.UAdd):
            v = compute(n.operand)
            if not isinstance(v, (int, float, complex, np.number)):
                raise DictionaryEvidenceError("Unusable unary operand")
            return -v if isinstance(n.op, ast.USub) else v
        if isinstance(n, ast.BinOp) and type(n.op) in (ast.Add, ast.Sub, ast.Mult, ast.Div):
            a, b = compute(n.left), compute(n.right)
            if type(a) in (str, list, tuple) or type(b) in (str, list, tuple):
                raise DictionaryEvidenceError("Nonnumeric contraction")
            if isinstance(n.op, ast.Add): return a + b
            if isinstance(n.op, ast.Sub): return a - b
            if isinstance(n.op, ast.Mult): return a * b
            if b == 0: raise DictionaryEvidenceError("Zero denominator")
            return a / b
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id in primitives:
            if n.keywords:
                raise DictionaryEvidenceError("Keyword calls forbidden")
            args = [compute(v) for v in n.args]
            value = primitives[n.func.id](*args)
            if n.func.id == "delS2" and len(n.args) == 4:
                # Arguments #1,#2 are the triplet/group and #3,#4 the two
                # fundamental slots. Complex conjugation is needed for an
                # unbarred pair in the Matchete-aligned Cartesian convention.
                pair = n.args[2:]
                if all(isinstance(v, ast.Name) and
                       v.id in ("idx1", "idx2", "idx3", "idx4") for v in pair):
                    refs = [int(v.id[3:])-1 for v in pair]
                    bar = [species[ref].startswith("Bar[") for ref in refs]
                    if bar == [False, False]: return complex(value).conjugate()
                    if bar not in ([True, True], [False, False]):
                        raise DictionaryEvidenceError("Mixed bar/unbar pair needs native CG evidence")
            return value
        raise DictionaryEvidenceError("Unsupported node in native invariant expression")

    dummy = [z for z in ast.walk(node) if isinstance(z, ast.Name) and z.id == "a"]
    if dummy:
        result = 0j
        for a in range(1, 4):
            env["a"] = a
            result += compute(node.body)
    else:
        result = compute(node.body)
    if not np.isfinite(result):
        raise DictionaryEvidenceError("Nonfinite invariant component")
    return complex(result)


def _target_tensor(rec: Mapping[str, Any], cg: np.ndarray) -> np.ndarray:
    if rec.get("status") != "Constructed_ConditionalSU2Frame_NativeNormalizationUncertified":
        raise DictionaryEvidenceError("RGBeta record has unexpected certification/status")
    if rec.get("prefactor_certified") is not False or rec.get("native_triplet_frame_certified") is not False:
        raise DictionaryEvidenceError("Unexpected prior certification: re-audit")
    dims, species = rec.get("dimensions"), rec.get("species")
    if (not isinstance(dims, list) or len(dims) != 4 or
            any(type(x) is not int or not 1 <= x <= 3 for x in dims) or
            not isinstance(species, list) or len(species) != 4):
        raise DictionaryEvidenceError("Malformed RGBeta component record")
    if not isinstance(rec.get("components"), list) or len(rec["components"]) != math.prod(dims):
        raise DictionaryEvidenceError("Invalid source RGBeta conditional array")
    matched = re.fullmatch(
        r"6\s*\*\s*Matrix\[([A-Za-z_0-9]+)\]\[##1\]\s*&\s*",
        rec.get("native_coupling_function", ""),
    )
    if not matched or matched.group(1) != rec.get("symbol"):
        raise DictionaryEvidenceError("RGBeta native quartic 6 Matrix factor/symbol changed")
    permutations = rec.get("native_unique_arrangements")
    if not isinstance(permutations,str) or len(permutations)>1500:
        raise DictionaryEvidenceError("Missing/oversized RGBeta native arrangement list")
    try:
        permutations = ast.literal_eval(permutations.replace("{","[").replace("}","]"))
    except (SyntaxError, ValueError, TypeError) as exc:
        raise DictionaryEvidenceError("Malformed native quartic permutations") from exc
    if (not isinstance(permutations,list) or not 1<=len(permutations)<=24 or
        any(not isinstance(row,list) or sorted(row)!=[1,2,3,4] for row in permutations)):
        raise DictionaryEvidenceError("Native quartic permutation metadata inconsistent")
    if rec.get("native_self_conjugate") not in ("False", False):
        raise DictionaryEvidenceError("This fit requires complex non-self-conjugate quartics")
    tensor = np.array([
        _native_invariant(rec["native_invariant_inputform"], tuple(indices), tuple(species), cg)
        for indices in product(*(range(1, d+1) for d in dims))
    ], dtype=np.complex128).reshape(dims)
    if not np.any(abs(tensor) > 1e-12):
        raise DictionaryEvidenceError("Zero conditionally reconstructed RGBeta tensor")
    return tensor


def _species(rec: Mapping[str, Any], kind: str) -> tuple[str, ...]:
    name = rec.get("species")
    if not isinstance(name, list) or len(name) != 4:
        raise DictionaryEvidenceError("Invalid field species")
    if kind == "matchete":
        return tuple(BARED.get(v, v) for v in name)
    if kind == "rgbeta":
        return tuple(name)
    raise ValueError("Unknown species kind")


def _axis_permutation(reference: tuple[str, ...], target: tuple[str, ...]) -> tuple[int, ...]:
    """Reorder target tensor axes to the Matchete monomial order."""
    remaining = list(range(4)); order = []
    for name in reference:
        candidates = [i for i in remaining if target[i] == name]
        if not candidates: raise DictionaryEvidenceError("Species multiset does not match")
        k = candidates[0]; order.append(k); remaining.remove(k)
    if remaining:
        raise DictionaryEvidenceError("Unmatched target axes")
    return tuple(order)


def _encode_complex(value: np.ndarray) -> dict[str, Any]:
    return {"real": value.real.tolist(), "imag": value.imag.tolist()}


def _sample_monomial(tensor: np.ndarray, species: tuple[str, ...],
                     fields: Mapping[str, np.ndarray]) -> complex:
    legs = []
    for name in species:
        conjugated = name.startswith("Bar[")
        key = name[4:-1] if conjugated else name
        arr = fields[key]
        legs.append(arr.conj() if conjugated else arr)
    return complex(np.einsum("abcd,a,b,c,d->", tensor, *legs, optimize=True))


def _sample_holdout(mat_cols: np.ndarray, rgb_cols: np.ndarray,
                    species: tuple[str, ...], dims: tuple[int, ...],
                    amap: np.ndarray, seed: int, samples: int = 25) -> float:
    rng = np.random.default_rng(seed)
    species_dims: dict[str, int] = {}
    for name, d in zip(species, dims):
        key = name[4:-1] if name.startswith("Bar[") else name
        if key in species_dims and species_dims[key] != d:
            raise DictionaryEvidenceError("Inconsistent field dimension")
        species_dims[key] = d
    ns, nm = rgb_cols.shape[1], mat_cols.shape[1]
    worst = 0.0
    for _ in range(samples):
        fields = {key: (rng.normal(size=d)+1j*rng.normal(size=d))/math.sqrt(2)
                  for key, d in species_dims.items()}
        cm = rng.normal(size=nm)+1j*rng.normal(size=nm)
        cr = amap@cm
        vm = 0j;vr=0j
        for j in range(nm):
            vm += cm[j]*_sample_monomial(mat_cols[:,j].reshape(dims),species,fields)
        for k in range(ns):
            vr += cr[k]*_sample_monomial(rgb_cols[:,k].reshape(dims),species,fields)
        # Include the explicit Hermitian conjugate of each complex quartic.
        lhs, rhs = 2*vm.real, 2*vr.real
        worst = max(worst, abs(lhs-rhs)/max(1.,abs(lhs),abs(rhs)))
    return worst


def _compare_model(model: str, mraw: list[dict], rraw: list[dict],
                   cg: np.ndarray, *, seed: int, tol: float) -> dict:
    if len(mraw) != len(rraw):
        raise DictionaryEvidenceError(f"{model}: mismatched quartic inventory")
    nm = len(mraw)
    if not nm:
        return {"status": "NoAdditionalQuartics", "n_extra_quartics": 0,
                "matchete_symbols":[], "rgbeta_symbols":[],
                "map_rgbeta_from_matchete":{"real":[],"imag":[]},
                "rank_matchete":0,"rank_rgbeta":0,
                "max_component_fit_residual":0.,"max_holdout_potential_residual":0.,
                "groups":[]}
    mnames = [x["symbol"] for x in mraw]
    rnames = [x["symbol"] for x in rraw]
    if len(set(mnames)) != nm or len(set(rnames)) != nm:
        raise DictionaryEvidenceError(f"{model}: repeated coupling in one-model inventory")
    mtensor = [_mat_tensor(x) for x in mraw]
    rtensor = [_target_tensor(x,cg) for x in rraw]
    mspec = [_species(x,"matchete") for x in mraw]
    rspec = [_species(x,"rgbeta") for x in rraw]
    mk = [tuple(sorted(x)) for x in mspec]
    rk = [tuple(sorted(x)) for x in rspec]
    if Counter(mk)!=Counter(rk):
        raise DictionaryEvidenceError(f"{model}: source/target quartic species groups differ")
    matrix = np.zeros((nm,nm),dtype=np.complex128)
    groups=[]
    maxfit=0.; maxhold=0.;rankm=rankr=0
    for idx,key in enumerate(sorted(set(mk))):
        mi=[i for i,k in enumerate(mk) if k==key]
        ri=[i for i,k in enumerate(rk) if k==key]
        reference=mspec[mi[0]]
        t0=mtensor[mi[0]]; dims=t0.shape
        if any(mspec[i]!=reference for i in mi):
            # Required by the current 50-slot source; this is not an algebraic
            # restriction, but prevents an unnoticed different tensor ordering.
            raise DictionaryEvidenceError("Source monomials in same group have reordered species")
        mt = np.stack([mtensor[i].reshape(-1) for i in mi],axis=1)
        rt = np.stack([rtensor[i].transpose(_axis_permutation(reference,rspec[i])).reshape(-1)
                       for i in ri],axis=1)
        if mt.shape!=rt.shape:
            raise DictionaryEvidenceError("Mismatch after species-axis alignment")
        a,_,rank,sv=np.linalg.lstsq(rt,mt,rcond=1e-12)
        rank_m=int(np.linalg.matrix_rank(mt,tol=1e-11*np.linalg.norm(mt,ord=2)))
        rank_r=int(np.linalg.matrix_rank(rt,tol=1e-11*np.linalg.norm(rt,ord=2)))
        if rank_m != len(mi) or rank_r != len(ri) or rank != len(ri):
            raise DictionaryEvidenceError(f"{model}: insufficient invariant rank in {key}")
        residual=float(np.linalg.norm(rt@a-mt)/max(np.linalg.norm(mt),1e-300))
        for row,j in enumerate(ri):
            for col,i in enumerate(mi):
                matrix[j,i]=a[row,col]
        hold=_sample_holdout(mt,rt,reference,dims,a,seed+idx)
        maxfit=max(maxfit,residual);maxhold=max(maxhold,hold)
        groups.append({"species_multiset":list(key), "source_indices":[i+1 for i in mi],
                       "target_indices":[i+1 for i in ri],"rank_matchete":rank_m,
                       "rank_rgbeta":rank_r,"component_fit_relative_residual":residual,
                       "independent_complex_potential_holdout_relative_residual":hold,
                       "map_rgbeta_from_matchete":_encode_complex(a),
                       "singular_values_rgbeta":[float(v) for v in sv],
                       "condition_number_rgbeta":float(sv[0]/sv[-1])})
        rankm+=rank_m;rankr+=rank_r
    complete=(maxfit<tol and maxhold<tol and rankm==rankr==nm)
    return {"status":("Candidate_RawInvariantSpacesIdentical_ConditionalFrame" if complete
                       else "Blocked_ConditionalInvariantSpaceMismatch"),
            "n_extra_quartics":nm,"matchete_symbols":mnames,"rgbeta_symbols":rnames,
            "map_rgbeta_from_matchete":_encode_complex(matrix),
            "rank_matchete":rankm,"rank_rgbeta":rankr,
            "max_component_fit_residual":maxfit,
            "max_holdout_potential_residual":maxhold,
            "groups":groups}


def solve_from_uploaded_evidence(mdata: Mapping[str, Any], rdata: Mapping[str, Any], *,
                                 msha: str = "unknown", rsha: str = "unknown",
                                 tol: float = 1e-9, seed: int = 1129) -> dict:
    """Compute *conditional raw-invariant* maps, never native physical ones."""
    if (mdata.get("format")!=MAT_FORMAT or mdata.get("source")!="MatcheteInvariantTensors" or
            rdata.get("format")!=RGB_CONDITIONAL_FORMAT or
            rdata.get("status")!="ConditionalFrameConstructed_NotRGBetaNativeOrMatcheteCertified"):
        raise DictionaryEvidenceError("Wrong raw Matchete/conditional RGBeta evidence formats")
    if not math.isfinite(tol) or not 0<tol<=1e-3 or type(seed) is not int:
        raise DictionaryEvidenceError("Invalid fitting controls")
    models_m=mdata.get("models");models_r=rdata.get("models")
    if (not isinstance(models_m,dict) or not isinstance(models_r,dict) or
            set(models_m)!=ALL16 or set(models_r)!=ALL16 or
            mdata.get("quartic_slot_count")!=50 or rdata.get("quartic_slot_count")!=50):
        raise DictionaryEvidenceError("Require the 16 model / 50 slot inventory")
    if rdata.get("interpretation",{}).get("physical_matchete_to_RGBeta_dictionary_certified") is not False:
        raise DictionaryEvidenceError("Unexpected physical-certification field")
    cg, cg_scale = calibrated_triplet_cg(mdata)
    outcomes={};total=0
    for i,tag in enumerate(sorted(ALL16)):
        mblock=models_m[tag];rblock=models_r[tag]
        mrecords=mblock.get("tensors");rrecords=rblock.get("records")
        if (not isinstance(mrecords,list) or not isinstance(rrecords,list) or
                rblock.get("n_records") != len(rrecords)):
            raise DictionaryEvidenceError(f"{tag}: invalid model record list")
        total+=len(mrecords)
        outcomes[tag]=_compare_model(tag,mrecords,rrecords,cg,seed=seed+100*i,tol=tol)
    if total!=50 or sum(x["n_extra_quartics"] for x in outcomes.values())!=50:
        raise DictionaryEvidenceError("Quartic slot count not 50")
    failed=[tag for tag,v in outcomes.items() if v["status"].startswith("Blocked_")]
    maxfit=max((x["max_component_fit_residual"] for x in outcomes.values()),default=0.)
    maxhold=max((x["max_holdout_potential_residual"] for x in outcomes.values()),default=0.)
    return {
        "format":OUT_FORMAT,
        "status":("All50RawExtraQuarticInvariantsMatched_ConditionalFrameOnly" if not failed else
                  "Blocked_ConditionalExtraQuarticDictionary"),
        "model_count":len(outcomes),"additional_quartic_slots":total,
        "n_failed_models":len(failed),"failed_models":failed,
        "triplet_frame":{
            "origin":"A_+0 Matchete Hdag Hdag S1 S2 CG; normalised for conditional RGBeta contraction",
            "anchor_scale":cg_scale,
            "calibrated_cg_2x2x3":_encode_complex(cg),
            "cg_orthonormality_residual":float(np.max(abs(np.einsum(
                "ija,ijb->ab",cg.conj(),cg)-np.eye(3)))),
            "unbarred_doublet_pair_uses_conjugate_cg":True,
            "triplet_frame_is_native_RGBeta_certified":False,
            "intertwiner_phase_is_independently_verified":False,
        },
        "source_evidence_sha256":{"matchete":msha,"conditional_rgbeta":rsha},
        "max_component_fit_residual":maxfit,
        "max_independent_potential_holdout_residual":maxhold,
        "models":outcomes,
        "scope":{
            "operator_inventory":PARTIAL_SCOPE,
            "coefficients_complex":True,
            "matchete_extra_quartic_plus_hc":True,
            "rgbeta_quartic_plus_hc":True,
            "unit_raw_invariant_potential_prefactors_assumed":True,
            "native_RGBeta_Coupling_and_UniqueArrangements_normalization_certified":False,
            "native_scalar_basis_verified":False,
            "native_mass_and_yukawa_dictionary_verified":False,
            "A_minus4_anti_Hermitian_Higgs_triplet_portal_preserved_but_not_in_this_export":True,
            "signed_GammaCC_to_Weyl_verified":False,
            "correct_complete_one_loop_ADM_certified":False,
            "physical_scalar_dictionary_certified":False,
        },
    }


def solve_from_paths(mat_path: Path, rgb_path: Path, *, tol: float=1e-9,
                     seed: int=1129) -> dict:
    for p in (mat_path,rgb_path):
        if not p.is_file(): raise DictionaryEvidenceError(f"Missing input file: {p}")
    mbytes=mat_path.read_bytes();rbytes=rgb_path.read_bytes()
    md=json.loads(mbytes.decode("utf-8-sig"));rd=json.loads(rbytes.decode("utf-8-sig"))
    return solve_from_uploaded_evidence(md,rd,msha=sha256(mbytes).hexdigest(),
                                        rsha=sha256(rbytes).hexdigest(),tol=tol,seed=seed)
