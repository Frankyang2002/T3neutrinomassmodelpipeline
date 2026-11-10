"""Read-only SU(2) Ward audit of individual Matchete source operators at MF."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from Numerical.parallel_t3.wilson.basis import build_basis
from Numerical.parallel_t3.wilson.coordinates import embed_formal_llss

def norm(x):return float(np.linalg.norm(x))
def ward(c,f,t):
    return (np.einsum("ki,kjab->ijab",f,c,optimize=True)
           +np.einsum("kj,ikab->ijab",f,c,optimize=True)
           +np.einsum("ca,ijcb->ijab",t,c,optimize=True)
           +np.einsum("cb,ijac->ijab",t,c,optimize=True))
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--capture",type=Path,default=Path("output/parallel_t3/all16_ward_capture01"))
    a=ap.parse_args();root=Path.cwd();basis=build_basis("A_-4")
    li=[i for pair in basis.leptons for i in pair]
    hi=list(range(4,basis.ns))
    gens=[(f[np.ix_(li,li)],t[np.ix_(hi,hi)]) for f,t in basis.generators()["SU2"]]
    for d in sorted(a.capture.glob("fit16_*_ward_diagnostic.json")):
        token=d.name.removesuffix("_ward_diagnostic.json")
        folder=root/"output/parallel_t3/background"/token/"A_m4"
        f=folder/"fermion_boundary_formal.npz"
        m=folder/"fermion_boundary_manifest.json"
        if not f.is_file():raise FileNotFoundError(f)
        with np.load(f,allow_pickle=False) as z:
            raw=z["C_LLSS_formal"].copy()
            by=z["C_LLSS_by_source_formal"].copy()
            indices=z["source_indices"].tolist()
        if by.ndim!=7 or by.shape[1:]!=raw.shape:
            raise ValueError(f"Unexpected source geometry {by.shape} versus {raw.shape}")
        if len(indices)!=len(by):raise ValueError("Source index count mismatch")
        difference=norm(by.sum(axis=0)-raw)/max(norm(raw),1e-100)
        if difference>1e-9:raise ValueError("By-source sum disagrees with saved boundary")
        def metrics(x):
            c=embed_formal_llss(x,basis,np.eye(basis.ns-4),1.)
            c=c[np.ix_(li,li,hi,hi)]
            return c,[norm(ward(c,tf,ts))/max(norm(c),1e-100) for tf,ts in gens]
        total,total_wards=metrics(raw)
        rows=[]
        for idx,term in zip(indices,by):
            c,w=metrics(term)
            rows.append({"source_operator_index":int(idx),"norm":norm(c),
                         "relative_norm":norm(c)/max(norm(total),1e-100),
                         "SU2_ward_relative":w,
                         "overlap_with_total":float(abs(np.vdot(c,total))
                           /max(norm(c)*norm(total),1e-200))})
        rows.sort(key=lambda x:x["norm"],reverse=True)
        meta=json.loads(m.read_text(encoding="utf-8-sig")) if m.is_file() else {}
        report={"trial":token,"status":"READ_ONLY_NOT_A_REPAIR",
                "B_source_file":str(f),"by_source_sum_relative_error":difference,
                "total_SU2_ward_relative":total_wards,
                "source_operators":rows,
                "seed_file":meta.get("source_matchete_seed"),
                "seed_sha256":meta.get("source_matchete_seed_sha256"),
                "source_indices_manifest":meta.get("source_operator_indices"),
                "producer_chain":[
                    "Numerical/parallel_t3/background/fermion.py:evaluate_formal_llss",
                    "Numerical/running/GeneralT3OrderedTensor.py:project_ordered_term",
                    "Numerical/running/GeneralT3DynamicFermionWilson.py:_dense_operator"],
                "warning":"A large source Ward residual does not alone establish which original Matchete operator or CG convention is incorrect."}
        dest=a.capture/f"{token}_per_source_SU2_ward.json"
        dest.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
        print("\n",token,"source sum error:",difference,"total Ward:",total_wards)
        for row in rows:print("source",row["source_operator_index"],
              "norm",row["norm"],"SU2 Ward",row["SU2_ward_relative"])
        print("report:",dest)
if __name__=="__main__":main()
