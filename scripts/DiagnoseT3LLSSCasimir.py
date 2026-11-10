"""Read-only SU(2) Casimir decomposition of retained LLSS boundaries."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from Numerical.parallel_t3.wilson.rough_stream2 import load_canonical

def norm(a):return float(np.linalg.norm(a))
def action(c,f,t):
    return (np.einsum("ki,kjab->ijab",f,c,optimize=True)
           +np.einsum("kj,ikab->ijab",f,c,optimize=True)
           +np.einsum("ca,ijcb->ijab",t,c,optimize=True)
           +np.einsum("cb,ijac->ijab",t,c,optimize=True))
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--capture",type=Path,default=Path("output/parallel_t3/all16_ward_capture01"))
    args=p.parse_args(); root=Path.cwd()
    for d in sorted(args.capture.glob("fit16_*_ward_diagnostic.json")):
        token=d.name.removesuffix("_ward_diagnostic.json")
        case=load_canonical(root,root/"output/parallel_t3/rough_inputs"/token/"a/A_m4",
                            "A_-4",c5_boundary="require",assume_formal_identity=True)
        b=case.basis
        if tuple(b.dims)!=(1,3,2):raise ValueError("A dimensions changed; Casimir polynomial must be updated")
        ls=[i for pair in b.leptons for i in pair]
        hs=list(range(4,b.ns))
        c=case.initial[np.ix_(ls,ls,hs,hs)]
        gens=[(f[np.ix_(ls,ls)],t[np.ix_(hs,hs)]) for f,t in b.generators()["SU2"]]
        def cas(x):
            return sum((action(action(x,f,t),f,t) for f,t in gens),np.zeros_like(x))
        # A: LL = 1/2 tensor 1/2 -> 0+1; S1=0, S2=1.
        # Therefore all LLSS sectors have total J in {0,1,2}.
        # Casimir eigenvalues J(J+1) = {0,2,6}.
        q=cas(c)
        q2=cas(q)
        invariant=c-(2/3)*q+q2/12  # (1-Q/2)(1-Q/6)
        violating=c-invariant
        total=max(norm(c),1e-100)
        ward=[norm(action(invariant,f,t))/max(norm(invariant),1e-100) for f,t in gens]
        orth=abs(np.vdot(invariant,violating))/max(norm(c)**2,1e-200)
        sector={}
        # Local heavy scalar indices: S1 (2 real), S2 (6 real).
        for name,aa,bb in (("S1_S1",range(0,2),range(0,2)),
                           ("S1_S2",range(0,2),range(2,8)),
                           ("S2_S1",range(2,8),range(0,2)),
                           ("S2_S2",range(2,8),range(2,8))):
            ix=np.ix_(range(6),range(6),list(aa),list(bb))
            x=c[ix]; v=violating[ix]; inv=invariant[ix]
            sector[name]={"norm":norm(x),"invariant_norm":norm(inv),
                          "violating_norm":norm(v),
                          "violating_squared_fraction":norm(v)**2/max(norm(x)**2,1e-200)}
        report={"trial":token,"status":"READ_ONLY_NOT_A_REPAIR",
                "SU2_invariant_squared_fraction":norm(invariant)**2/total**2,
                "SU2_violating_squared_fraction":norm(violating)**2/total**2,
                "SU2_projected_Ward_relative":ward,
                "invariant_violating_orthogonality":orth,
                "reconstruction_relative":norm(invariant+violating-c)/total,
                "sectors":sector,
                "interpretation":"Exact Casimir projector for A with total spins J=0,1,2. Projection is diagnostic only; never substitute it into matching without a derived UV contraction."}
        dest=args.capture/f"{token}_LLSS_casimir_sectors.json"
        dest.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
        print("\n",token,"SU2 invariant fraction:",report["SU2_invariant_squared_fraction"],
              "violating fraction:",report["SU2_violating_squared_fraction"])
        print("projected Ward:",ward,"orthogonality:",orth)
        for name,x in sector.items():print(name,x)
        print("report:",dest)
if __name__=="__main__":main()
