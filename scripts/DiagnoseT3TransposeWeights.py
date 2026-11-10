"""Read-only LLSS T3 weight decomposition with correctly transformed coefficient indices."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from Numerical.parallel_t3.wilson.rough_stream2 import load_canonical

def norm(a): return float(np.linalg.norm(a))
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--capture",type=Path,default=Path("output/parallel_t3/all16_ward_capture01"))
    args=ap.parse_args()
    root=Path.cwd()
    for old in sorted(args.capture.glob("fit16_*_ward_diagnostic.json")):
        token=old.name.removesuffix("_ward_diagnostic.json")
        folder=root/"output/parallel_t3/rough_inputs"/token/"a/A_m4"
        case=load_canonical(root,folder,"A_-4",c5_boundary="require",assume_formal_identity=True)
        b=case.basis
        lidx=[i for p in b.leptons for i in p]
        hidx=list(range(4,b.ns))
        c=case.initial[np.ix_(lidx,lidx,hidx,hidx)]
        ff,ss=b.generators()["SU2"][2]
        f=ff[np.ix_(lidx,lidx)]
        t=ss[np.ix_(hidx,hidx)]
        # Ward action on coefficient indices is T^T, not T.
        # Eigenvector matrices U diagonalise the transpose action.
        wf,uf=np.linalg.eigh(f.T)
        ws,us=np.linalg.eigh(t.T)
        # Transform all four tensor indices with U^dagger.
        z=np.einsum("Ii,Jj,Aa,Bb,ijab->IJAB",uf.conj().T,uf.conj().T,
                    us.conj().T,us.conj().T,c,optimize=True)
        weights=(wf[:,None,None,None]+wf[None,:,None,None]
                 +ws[None,None,:,None]+ws[None,None,None,:])
        forbidden=np.abs(weights)>1e-8
        frac=float(np.sum(np.abs(z[forbidden])**2)/max(np.sum(np.abs(z)**2),1e-200))
        # Exact T3 Ward in eigenbasis must equal weights * transformed tensor.
        direct=(np.einsum("ki,kjab->ijab",f,c,optimize=True)
                +np.einsum("kj,ikab->ijab",f,c,optimize=True)
                +np.einsum("ca,ijcb->ijab",t,c,optimize=True)
                +np.einsum("cb,ijac->ijab",t,c,optimize=True))
        transformed_ward=np.einsum("Ii,Jj,Aa,Bb,ijab->IJAB",
                                  uf.conj().T,uf.conj().T,us.conj().T,us.conj().T,
                                  direct,optimize=True)
        check=norm(transformed_ward-weights*z)/max(norm(direct),1e-100)
        indices=np.argwhere(forbidden)
        top=sorted(indices,key=lambda idx:abs(z[tuple(idx)]),reverse=True)[:12]
        report={
            "trial":token,"status":"READ_ONLY_NOT_A_REPAIR",
            "LLSS_norm":norm(c),
            "T3_ward_norm_relative":norm(direct)/max(norm(c),1e-100),
            "T3_weight_basis_covariance_error":check,
            "T3_forbidden_squared_norm_fraction":frac,
            "fermion_T3_eigenvalues":wf.tolist(),
            "scalar_T3_eigenvalues":ws.tolist(),
            "largest_forbidden":[{"weight":float(weights[tuple(idx)]),
                "coefficient_abs":float(abs(z[tuple(idx)])),
                "eigenbasis_indices":idx.tolist()} for idx in top],
            "source":case.source["boundary"],
            "source_folder":str(folder),
            "note":"T3 charge basis uses eigenvectors of the transposed generators acting on coefficient indices. Does not identify the original matching code or certify its representation."
        }
        dest=args.capture/f"{token}_T3_transpose_weight_audit.json"
        dest.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
        print("\n",token,"T3 forbidden fraction:",frac,
              "Ward:",report["T3_ward_norm_relative"],
              "basis consistency:",check)
        print("source:",report["source"])
        print("largest forbidden:",report["largest_forbidden"][:3])
        print("report:",dest)

if __name__=="__main__":main()
