"""Read-only LLSS weak-index and source/canonical audit for retained A(-4) trials."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from Numerical.parallel_t3.wilson.basis import build_basis
from Numerical.parallel_t3.wilson.rough_stream2 import load_canonical

def ward(c, f, s):
    return (np.einsum("ki,kjab->ijab", f, c, optimize=True)
            + np.einsum("kj,ikab->ijab", f, c, optimize=True)
            + np.einsum("ca,ijcb->ijab", s, c, optimize=True)
            + np.einsum("cb,ijac->ijab", s, c, optimize=True))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--capture",type=Path,default=Path("output/parallel_t3/all16_ward_capture01"))
    args=ap.parse_args()
    root=Path.cwd()
    cap=args.capture
    if not cap.is_dir(): raise FileNotFoundError(cap)
    for diagnostic in sorted(cap.glob("fit16_*_ward_diagnostic.json")):
        token=diagnostic.name.removesuffix("_ward_diagnostic.json")
        if not token.startswith("fit16_"): continue
        base=root/"output/parallel_t3/rough_inputs"/token/"a/A_m4"
        case=load_canonical(root,base,"A_-4",c5_boundary="require",assume_formal_identity=True)
        b=case.basis
        c=case.initial
        f=[i for pair in b.leptons for i in pair]
        heavy=list(range(4,b.ns))
        llss=c[np.ix_(f,f,heavy,heavy)]
        generators=b.generators()["SU2"]
        rows=[]
        for k,(tf,ts) in enumerate(generators):
            delta=ward(c,tf,ts)[np.ix_(f,f,heavy,heavy)]
            by_weak={}
            for i in range(2):
                for j in range(2):
                    sub=delta[np.ix_(list(range(i,6,2)),list(range(j,6,2)),list(range(len(heavy))),list(range(len(heavy))))]
                    by_weak[f"L{i+1}L{j+1}"]=float(np.linalg.norm(sub))
            flat=np.argsort(abs(delta).ravel())[::-1][:8]
            largest=[]
            for ix in flat:
                a1,a2,s1,s2=np.unravel_index(int(ix),delta.shape)
                largest.append({"lepton_labels":[b.fermion_labels[f[a1]],b.fermion_labels[f[a2]]],
                                "scalar_labels":[b.scalar_labels[heavy[s1]],b.scalar_labels[heavy[s2]]],
                                "ward_abs":float(abs(delta[a1,a2,s1,s2]))})
            rows.append({"generator":f"SU2_T{k+1}","ward_norm":float(np.linalg.norm(delta)),
                         "relative_to_LLSS":float(np.linalg.norm(delta)/max(np.linalg.norm(llss),1e-100)),
                         "weak_pair_ward_norms":by_weak,"largest_violations":largest})
        # Compare retained Wilson boundary with the independently loaded canonical source.
        wpath=root/"output/parallel_t3/rough_wilson"/token/"A_m4/wilson.json"
        meta=json.loads(wpath.read_text(encoding="utf-8"))
        with np.load(wpath.parent/meta["tensor_npz"],allow_pickle=False) as z:
            saved=z["C_psi2phi2_tree_mf"]
        discrepancy=float(np.linalg.norm(c-saved)/max(np.linalg.norm(saved),1e-100))
        report={"status":"READ_ONLY_NOT_A_REPAIR","trial":token,
                "boundary_source":case.source["boundary"],
                "LLSS_norm":float(np.linalg.norm(llss)),
                "canonical_vs_retained_boundary_relative":discrepancy,
                "SU2_LLSS":rows,
                "note":"An assumed identity heavy map or unsigned spinor factor is not certification. Inspect source mapping before modifying coefficients."}
        target=cap/f"{token}_LLSS_source_ward_audit.json"
        target.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
        print("\n",token,"source:",report["boundary_source"])
        print("LLSS norm:",report["LLSS_norm"],"source mismatch:",discrepancy)
        for item in rows:
            print(item["generator"],"relative:",item["relative_to_LLSS"])
            print("  weak pairs:",item["weak_pair_ward_norms"])
            print("  largest:",item["largest_violations"][0] if item["largest_violations"] else None)
        print("report:",target)

if __name__=="__main__":
    main()
