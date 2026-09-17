from __future__ import annotations
import argparse, json, re, sys
from collections import Counter
from dataclasses import dataclass
from math import factorial
from pathlib import Path
from typing import Iterable
import sympy as sp

PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from RGE.general.GaugeGenerators import su2_complex_generators
from RGE.running.EFT1QuarticAdapter import SparseQuarticTensor, load_and_build_eft1_quartic_tensor

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

def _weight(key):
    c = Counter(key); w = factorial(4)
    for n in c.values(): w //= factorial(n)
    return w

def _inner(a,b):
    return sp.simplify(sum(_weight(k)*sp.conjugate(a[k])*b[k] for k in set(a.entries)|set(b.entries)))

def _basis_tensor(full, name, identify_conjugate=True):
    s = sp.Symbol(name); cs = sp.conjugate(s); entries={}
    for k,raw in full.nonzero_items():
        e = sp.expand(raw)
        if identify_conjugate: e = sp.expand(e.xreplace({cs:s}))
        c = sp.simplify(e.coeff(s))
        if c != 0: entries[k]=c
    return SparseQuarticTensor(entries)

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
    G=sp.Matrix([[_inner(x,y) for y in basis] for x in basis])
    rhs=sp.Matrix([_inner(x,target) for x in basis])
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
    return int(sp.Matrix([[_inner(a,b) for b in ts] for a in ts]).rank())

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
        mbasis=[_basis_tensor(full,n,True) for n in mnames]
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

def main():
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

if __name__=="__main__":
    raise SystemExit(main())
