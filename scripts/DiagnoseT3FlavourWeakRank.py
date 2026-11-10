"""Corrected read-only rank-one and cross-trial LLSS factorization audit."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from Numerical.parallel_t3.wilson.rough_stream2 import load_canonical

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--capture",type=Path,default=Path("output/parallel_t3/all16_ward_capture01"))
    args=p.parse_args()
    root=Path.cwd()
    weak={}
    for diag in sorted(args.capture.glob("fit16_*_ward_diagnostic.json")):
        token=diag.name.removesuffix("_ward_diagnostic.json")
        folder=root/"output/parallel_t3/rough_inputs"/token/"a/A_m4"
        case=load_canonical(root,folder,"A_-4",c5_boundary="require",assume_formal_identity=True)
        b=case.basis
        ls=[i for pair in b.leptons for i in pair]
        hs=list(range(4,b.ns))
        full=case.initial[np.ix_(ls,ls,hs,hs)]
        matrix=full.reshape(3,2,3,2,len(hs),len(hs)).transpose(0,2,1,3,4,5).reshape(9,-1)
        u,s,vh=np.linalg.svd(matrix,full_matrices=False)
        flavour=(u[:,0]*s[0])
        weak_row=vh[0].copy()  # V dagger, NOT its complex conjugate
        recon=np.outer(flavour,weak_row)
        relative=float(np.linalg.norm(matrix-recon)/max(np.linalg.norm(matrix),1e-100))
        # Independent check: project all rows onto the normalized weak row.
        flavour_direct=matrix@weak_row.conj()
        independent=float(np.linalg.norm(matrix-np.outer(flavour_direct,weak_row))
                          /max(np.linalg.norm(matrix),1e-100))
        weak[token]=weak_row
        report={"trial":token,"status":"READ_ONLY_NOT_A_REPAIR",
                "relative_singular_values":(s/s[0]).tolist() if s[0] else [],
                "rank_one_relative_error":relative,
                "independent_row_projection_relative_error":independent,
                "flavour_matrix_real":flavour.reshape(3,3).real.tolist(),
                "flavour_matrix_imag":flavour.reshape(3,3).imag.tolist(),
                "weak_tensor_real":weak_row.reshape(2,2,len(hs),len(hs)).real.tolist(),
                "weak_tensor_imag":weak_row.reshape(2,2,len(hs),len(hs)).imag.tolist(),
                "note":"The SVD right singular row is vh[0] without conjugation. A universal tensor does not certify the UV matching contraction."}
        dest=args.capture/f"{token}_flavour_weak_factorization_corrected.json"
        dest.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
        print(token,"sigma2/sigma1:",report["relative_singular_values"][1],
              "rank-one error:",relative,"independent error:",independent)
        print("report:",dest)
    keys=list(weak)
    for i in range(len(keys)):
        for j in range(i+1,len(keys)):
            overlap=float(abs(np.vdot(weak[keys[i]],weak[keys[j]])))
            print("weak tensor overlap",keys[i],keys[j],":",overlap,
                  "phase-independent mismatch:",float(np.sqrt(max(0.,2-2*min(1.,overlap)))))
if __name__=="__main__":main()
