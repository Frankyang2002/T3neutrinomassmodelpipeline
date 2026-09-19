from __future__ import annotations

"""Scalar-quartic basis mapping and invariant-completeness validation."""

import sys
import argparse, json, re, sys
from collections import Counter
from dataclasses import dataclass
from math import factorial
from pathlib import Path
from typing import Iterable
import sympy as sp
from RGE.general.GaugeGenerators import su2_complex_generators
from RGE.group_factors.core.MixingQuarticTensorAlgebra import basis_tensor, tensor_inner_product
from RGE.running.EFT1TensorAdapters import SparseQuarticTensor, load_and_build_eft1_quartic_tensor
import argparse
import json

# ---------------------------------------------------------------------------
# Scalar-quartic basis mapping
# Former source: ScalarQuarticBasisMap.py
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


@dataclass(frozen=True)
class FieldBlock:
    name: str
    dimension: int
    first_real_index: int
    def real_index(self, c: int) -> int: return self.first_real_index + 2*c
    def imag_index(self, c: int) -> int: return self.first_real_index + 2*c + 1

def _text(x): return sp.sstr(sp.factor(sp.simplify(x)))

def _blocks(d1, d2):
    return {
        "H": FieldBlock("H", 2, 1),
        "S1": FieldBlock("S1", d1, 5),
        "S2": FieldBlock("S2", d2, 5 + 2*d1),
    }

def _bilinear(left, matrix, right=None):
    if right is None: right = left
    if matrix.shape != (left.dimension, right.dimension):
        raise ValueError("matrix/field dimension mismatch")
    rt2 = sp.sqrt(2); out={}
    for p in range(left.dimension):
        le=((left.real_index(p),1/rt2),(left.imag_index(p),-sp.I/rt2))
        for q in range(right.dimension):
            m=sp.simplify(matrix[p,q])
            if m==0: continue
            re=((right.real_index(q),1/rt2),(right.imag_index(q),sp.I/rt2))
            for i,ci in le:
                for j,cj in re:
                    k=tuple(sorted((i,j)))
                    out[k]=sp.simplify(out.get(k,0)+m*ci*cj)
    return {k:sp.simplify(v) for k,v in out.items() if v!=0}

def _mul_quad(a,b):
    poly={}
    for ka,va in a.items():
        for kb,vb in b.items():
            k=tuple(sorted(ka+kb))
            poly[k]=sp.simplify(poly.get(k,0)+va*vb)
    entries={}
    for k,c in poly.items():
        fac=1
        for n in Counter(k).values(): fac*=factorial(n)
        v=sp.simplify(c*fac)
        if v!=0: entries[k]=v
    return SparseQuarticTensor(entries)

def _sum_tensors(ts: Iterable[SparseQuarticTensor]):
    e={}
    for t in ts:
        for k,v in t.nonzero_items():
            e[k]=sp.simplify(e.get(k,0)+v)
    return SparseQuarticTensor(e)

def _singlet(a,b):
    return _mul_quad(_bilinear(a,sp.eye(a.dimension)), _bilinear(b,sp.eye(b.dimension)))

def _adjoint(a,b):
    if a.dimension==1 or b.dimension==1: return SparseQuarticTensor({})
    ga=su2_complex_generators(a.dimension); gb=su2_complex_generators(b.dimension)
    return _sum_tensors(_mul_quad(_bilinear(a,ta),_bilinear(b,tb)) for ta,tb in zip(ga,gb,strict=True))

def _exchange(a,b):
    if a.dimension != b.dimension: raise ValueError("exchange requires equal dimensions")
    I=sp.eye(a.dimension)
    return _mul_quad(_bilinear(a,I,b),_bilinear(b,I,a))

def _sector_names(seed,prefix):
    """Return every Matchete coupling belonging to one mixed-quartic sector.

    Matchete naming is representation dependent.  Some sectors use only
    prefixInv1, prefixInv2, ... while others keep the singlet contraction under
    the bare name ``prefix`` and label only the additional invariant(s).
    """
    names=set()
    pat=re.compile(
        r"Coupling\[(?P<name>"+re.escape(prefix)+r"(?:Inv\d+)?),"
    )
    for rec in seed.get("ScalarQuarticTerms",[]):
        for m in pat.finditer(rec["TermInputForm"]):
            names.add(m.group("name"))

    def order(name):
        if name == prefix:
            return (0, 0)
        m = re.search(r"Inv(\d+)$", name)
        return (1, int(m.group(1)) if m else 999)

    return sorted(names,key=order)

def _coordinates(target,basis):
    G=sp.Matrix([[tensor_inner_product(x, y) for y in basis] for x in basis])
    rhs=sp.Matrix([tensor_inner_product(x, target) for x in basis])
    if G.det()==0: raise ValueError("Matchete basis linearly dependent")
    coeff=[sp.simplify(x) for x in G.inv()*rhs]
    residual=0
    keys=set(target.entries)
    for t in basis: keys |= set(t.entries)
    for k in keys:
        rec=sp.simplify(sum(c*t[k] for c,t in zip(coeff,basis,strict=True)))
        if sp.simplify(target[k]-rec)!=0: residual+=1
    return coeff,residual

def _rank(ts):
    if not ts: return 0
    return int(sp.Matrix([[tensor_inner_product(a, b) for b in ts] for a in ts]).rank())

def _physical_sector(sector, blocks):
    if sector=="H1": a,b=blocks["H"],blocks["S1"]
    elif sector=="H2": a,b=blocks["H"],blocks["S2"]
    elif sector=="12": a,b=blocks["S1"],blocks["S2"]
    else: raise KeyError(sector)
    labels=["Singlet"]; tensors=[_singlet(a,b)]; note=None
    adj=_adjoint(a,b)
    if adj.entries: labels.append("Adj"); tensors.append(adj)
    if sector=="12" and a.dimension==b.dimension==3:
        labels.append("CrossExchangeCandidate"); tensors.append(_exchange(a,b))
        note="CrossExchangeCandidate=(S1^dagger S2)(S2^dagger S1); identify with RGBeta Cross only after convention check."
    return labels,tensors,note

def build_basis_map(seed_path: Path, rgbeta_path: Path):
    seed=json.loads(seed_path.read_text(encoding="utf-8"))
    rgbeta=json.loads(rgbeta_path.read_text(encoding="utf-8-sig"))
    md=rgbeta["metadata"]; d1=int(md["dS1"]); d2=int(md["dS2"])
    full=load_and_build_eft1_quartic_tensor(seed_path,rgbeta_path)
    blocks=_blocks(d1,d2)
    result={"status":"Success","seed":str(seed_path),"rgbeta":str(rgbeta_path),
            "dimensions":{"dS1":d1,"dS2":d2},
            "definition":"PhysicalTensor[alpha]=sum_i A[i,alpha] MatcheteTensor[i]; lambda_M=A lambda_P; R_P=R_M A.",
            "sectors":{}}
    for sector,prefix in (("H1","lambdaH1"),("H2","lambdaH2"),("12","lambda12")):
        mnames=_sector_names(seed,prefix)
        mbasis=[basis_tensor(full, n, identify_conjugate=True) for n in mnames]
        plabels,pbasis,note=_physical_sector(sector,blocks)
        kept=[(l,t) for l,t in zip(plabels,pbasis,strict=True) if t.entries]
        plabels=[x[0] for x in kept]; pbasis=[x[1] for x in kept]
        if len(mbasis)!=len(pbasis):
            result["status"]="DimensionMismatch"
            result["sectors"][sector]={"matchete_basis":mnames,"physical_basis":plabels,
                "matchete_rank":_rank(mbasis),"physical_rank":_rank(pbasis),
                "error":f"Basis dimensions differ: {len(mbasis)} Matchete vs {len(pbasis)} physical.","note":note}
            continue
        cols=[]; residuals={}
        for lab,t in zip(plabels,pbasis,strict=True):
            c,r=_coordinates(t,mbasis); cols.append(c); residuals[lab]=r
        A=sp.Matrix.hstack(*(sp.Matrix(c) for c in cols))
        det=sp.simplify(A.det())
        ok=all(v==0 for v in residuals.values())
        if not ok: result["status"]="ResidualFailure"
        result["sectors"][sector]={
            "status":"Success" if ok else "ResidualFailure",
            "matchete_basis":mnames,"physical_basis":plabels,
            "matchete_rank":_rank(mbasis),"physical_rank":_rank(pbasis),
            "A_rows_matchete_cols_physical":[[_text(A[i,j]) for j in range(A.cols)] for i in range(A.rows)],
            "det_A":_text(det),"coordinate_residual_nonzero_components":residuals,"note":note}
    return result

def basismap_main():
    p=argparse.ArgumentParser()
    p.add_argument("quartic_seed_json",type=Path)
    p.add_argument("--rgbeta",type=Path,default=None)
    p.add_argument("--output",type=Path,default=None)
    a=p.parse_args()
    rgb=a.rgbeta or a.quartic_seed_json.with_name("eft1_rgbeta_rge.json")
    payload=build_basis_map(a.quartic_seed_json,rgb)
    print(f"dS1={payload['dimensions']['dS1']}, dS2={payload['dimensions']['dS2']}")
    print("Convention: lambda_M = A lambda_P; R_P = R_M A\n")
    for sector,info in payload["sectors"].items():
        print(f"[{sector}]")
        if "error" in info:
            print("  ERROR:",info["error"]); print("  Matchete:",info["matchete_basis"]); print("  physical:",info["physical_basis"]); print(); continue
        print("  Matchete basis:",", ".join(info["matchete_basis"]))
        print("  Physical basis:",", ".join(info["physical_basis"]))
        print("  A (rows=Matchete, cols=physical):")
        for n,row in zip(info["matchete_basis"],info["A_rows_matchete_cols_physical"],strict=True):
            print(f"    {n:<18} "+"  ".join(f"{x:>18}" for x in row))
        print("  det(A) =",info["det_A"])
        print("  residuals =",json.dumps(info["coordinate_residual_nonzero_components"]))
        if info.get("note"): print("  note:",info["note"])
        print()
    if a.output:
        a.output.parent.mkdir(parents=True,exist_ok=True)
        a.output.write_text(json.dumps(payload,indent=2),encoding="utf-8")
        print("JSON summary:",a.output)
    print("OVERALL:",payload["status"])
    return 0 if payload["status"]=="Success" else 1


# ---------------------------------------------------------------------------
# Scalar-quartic invariant-completeness validation
# Former source: ValidateScalarQuarticInvariantCompleteness.py
# ---------------------------------------------------------------------------

"""Validate completeness of mixed scalar-quartic invariant sectors in a Matchete seed.

For two SU(2) multiplets A and B with dimensions d_A and d_B, the quartics

    A^\dagger A B^\dagger B

span one invariant for every common total isospin J appearing in

    A^\dagger \otimes A  and  B^\dagger \otimes B.

For SU(2), j \otimes j contains J = 0,1,...,2j once, hence

    N_invariants(A,B) = min(d_A, d_B).

This gives the expected dimensions:
    H(2)-singlet(1): 1
    H(2)-doublet(2): 2
    H(2)-triplet(3): 2
    doublet-doublet: 2
    triplet-triplet: 3
    singlet-triplet: 1

The script compares this representation-theory count against the distinct
couplings actually exported in ScalarQuarticTerms.  It deliberately does not
assume coupling names such as Inv1/Inv2 beyond extracting the unique coupling
symbol from each quartic term.
"""



FIELD_DIMS_KEY = {
    "H": "H",
    "NewScalar1": "S1",
    "NewScalar2": "S2",
}


def _matching_bracket(text: str, open_index: int) -> int:
    depth = 0
    for i in range(open_index, len(text)):
        if text[i] == "[":
            depth += 1
        elif text[i] == "]":
            depth -= 1
            if depth == 0:
                return i
    raise ValueError("unbalanced brackets")


def _split_top_level(text: str) -> list[str]:
    out = []
    start = 0
    sq = cu = pa = 0
    for i, ch in enumerate(text):
        if ch == "[":
            sq += 1
        elif ch == "]":
            sq -= 1
        elif ch == "{":
            cu += 1
        elif ch == "}":
            cu -= 1
        elif ch == "(":
            pa += 1
        elif ch == ")":
            pa -= 1
        elif ch == "," and sq == 0 and cu == 0 and pa == 0:
            out.append(text[start:i].strip())
            start = i + 1
    out.append(text[start:].strip())
    return out


def _coupling_name(term: str) -> str:
    start = term.find("Coupling[")
    if start < 0:
        raise ValueError("quartic term contains no Coupling[...]")
    open_index = start + len("Coupling")
    close = _matching_bracket(term, open_index)
    return _split_top_level(term[open_index + 1:close])[0]


def _fields(term: str) -> tuple[str, ...]:
    names = []
    pos = 0
    while True:
        start = term.find("Field[", pos)
        if start < 0:
            break
        open_index = start + len("Field")
        close = _matching_bracket(term, open_index)
        args = _split_top_level(term[open_index + 1:close])
        if len(args) >= 2 and args[1] == "Scalar" and args[0] in FIELD_DIMS_KEY:
            names.append(args[0])
        pos = close + 1
    return tuple(names)


def _sector(fields: tuple[str, ...]) -> str | None:
    c = Counter(fields)
    if c == Counter({"H": 2, "NewScalar1": 2}):
        return "H1"
    if c == Counter({"H": 2, "NewScalar2": 2}):
        return "H2"
    if c == Counter({"NewScalar1": 2, "NewScalar2": 2}):
        return "12"
    return None


def expected_invariant_count(d_a: int, d_b: int) -> int:
    return min(int(d_a), int(d_b))


def validate(seed_path: Path, d_s1: int, d_s2: int) -> dict:
    seed = json.loads(seed_path.read_text(encoding="utf-8"))
    dims = {"H": 2, "S1": int(d_s1), "S2": int(d_s2)}
    sector_fields = {
        "H1": ("H", "S1"),
        "H2": ("H", "S2"),
        "12": ("S1", "S2"),
    }

    couplings = {name: set() for name in sector_fields}

    for record in seed.get("ScalarQuarticTerms", []):
        term = str(record.get("TermInputForm", ""))
        sec = _sector(_fields(term))
        if sec is not None:
            couplings[sec].add(_coupling_name(term))

    rows = []
    for sec, (a, b) in sector_fields.items():
        expected = expected_invariant_count(dims[a], dims[b])
        observed_names = sorted(couplings[sec])
        observed = len(observed_names)
        rows.append({
            "sector": sec,
            "fields": [a, b],
            "dimensions": [dims[a], dims[b]],
            "expected_invariant_count": expected,
            "observed_distinct_coupling_count": observed,
            "observed_couplings": observed_names,
            "missing_count": max(0, expected - observed),
            "status": "PASS" if expected == observed else "FAIL",
        })

    return {
        "status": "Success" if all(r["status"] == "PASS" for r in rows) else "Incomplete",
        "seed": str(seed_path),
        "dimensions": {"dH": 2, "dS1": int(d_s1), "dS2": int(d_s2)},
        "representation_theory_rule": "N(A^dagger A B^dagger B)=min(d_A,d_B) for SU(2)",
        "sectors": rows,
    }


def invcomplete_main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("quartic_seed_json", type=Path)
    p.add_argument("dS1", type=int)
    p.add_argument("dS2", type=int)
    p.add_argument("--output", type=Path, default=None)
    a = p.parse_args()

    payload = validate(a.quartic_seed_json, a.dS1, a.dS2)

    print(f"seed: {a.quartic_seed_json}")
    print("rule: N_invariants = min(d_A, d_B)")
    print()
    print(f"{'sector':<7} {'dims':<8} {'expected':>8} {'observed':>8} {'missing':>8}  couplings")
    print("-" * 80)
    for row in payload["sectors"]:
        dims = "x".join(map(str, row["dimensions"]))
        print(
            f"{row['sector']:<7} {dims:<8} "
            f"{row['expected_invariant_count']:>8} "
            f"{row['observed_distinct_coupling_count']:>8} "
            f"{row['missing_count']:>8}  "
            f"{row['observed_couplings']}  {row['status']}"
        )

    if a.output:
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print()
        print(f"JSON summary: {a.output}")

    print(f"OVERALL: {payload['status']}")
    return 0 if payload["status"] == "Success" else 1


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Consolidated group-factor research/validation driver."
    )
    parser.add_argument("mode", choices=[
        "basis-map",
        "invariant-completeness",
    ])
    args, remaining = parser.parse_known_args()
    sys.argv = [sys.argv[0], *remaining]

    if args.mode == "basis-map":
        basismap_main()
    elif args.mode == "invariant-completeness":
        invcomplete_main()


if __name__ == "__main__":
    main()
