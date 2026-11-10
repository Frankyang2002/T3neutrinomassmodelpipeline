"""Read-only SU2 generator-algebra and LLSS weight-selection audit."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from Numerical.parallel_t3.wilson.rough_stream2 import load_canonical

def ward(c,t,s):
    return (np.einsum("ki,kjab->ijab",t,c,optimize=True)
            +np.einsum("kj,ikab->ijab",t,c,optimize=True)
            +np.einsum("ca,ijcb->ijab",s,c,optimize=True)
            +np.einsum("cb,ijac->ijab",s,c,optimize=True))

def norm(x): return float(np.linalg.norm(x))
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--capture",type=Path,default=Path("output/parallel_t3/all16_ward_capture01"))
    args=p.parse_args()
    root=Path.cwd()
    for d in sorted(args.capture.glob("fit16_*_ward_diagnostic.json")):
        token=d.name.removesuffix("_ward_diagnostic.json")
        case=load_canonical(root,root/"output/parallel_t3/rough_inputs"/token/"a/A_m4",
                            "A_-4",c5_boundary="require",assume_formal_identity=True)
        b=case.basis; inds=[i for pair in b.leptons for i in pair]; heavy=list(range(4,b.ns))
        c=case.initial[np.ix_(inds,inds,heavy,heavy)]
        gen=b.generators()["SU2"]
        fs=[f[np.ix_(inds,inds)] for f,_ in gen]
        ss=[s[np.ix_(heavy,heavy)] for _,s in gen]
        denom=max(norm(c),1e-100)
        # Validate generator Lie algebra, not just Ward output.
        def algebra(g):
            signs={}
            for sign in (1,-1):
                errs=[]
                for a,b1,k in ((0,1,2),(1,2,0),(2,0,1)):
                    comm=g[a]@g[b1]-g[b1]@g[a]
                    errs.append(norm(comm-sign*1j*g[k])/max(norm(g[k]),1e-100))
                signs[str(sign)]=errs
            return signs
        # Test whether a consistent conjugate/dual action resolves the tensor.
        # These are diagnostics, NOT sanctioned replacements for the representation.
        variants={}
        for name, fset,sset in (
            ("declared",fs,ss),
            ("scalar_dual_minus_transpose",fs,[-s.T for s in ss]),
            ("fermion_dual_minus_transpose",[-f.T for f in fs],ss),
            ("both_dual_minus_transpose",[-f.T for f in fs],[-s.T for s in ss]),
        ):
            variants[name]=[norm(ward(c,f,s))/denom for f,s in zip(fset,sset)]
        # T3 eigenspace selection: diagonalise the Hermitian scalar generator.
        eig,vec=np.linalg.eigh(ss[2])
        transformed=np.einsum("ijab,aA,bB->ijAB",c,vec,vec,optimize=True)
        # Since C transforms with two scalar LOWER indices, use transpose map:
        # sum_a V[a,A] C_ab V[b,B] (no conjugation).
        forbidden=[]
        total=0.0
        for i in range(6):
            for j in range(6):
                for a in range(len(eig)):
                    for bb in range(len(eig)):
                        weight=float(fs[2][i,i].real+fs[2][j,j].real+eig[a]+eig[bb])
                        v=abs(transformed[i,j,a,bb])**2
                        total+=v
                        if abs(weight)>1e-7 and v>0:
                            forbidden.append((float(v),weight,i,j,a,bb))
        forbidden.sort(reverse=True)
        report={
            "trial":token,
            "status":"READ_ONLY_NOT_A_REPAIR",
            "LLSS_norm":denom,
            "generator_algebra_fermion":algebra(fs),
            "generator_algebra_scalar":algebra(ss),
            "SU2_ward_variants":variants,
            "T3_weight_forbidden_fraction":float(sum(x[0] for x in forbidden)/max(total,1e-200)),
            "largest_T3_forbidden":[{"abs_coefficient":float(np.sqrt(v)),"total_T3_weight":w,
                 "leptons":[b.fermion_labels[inds[i]],b.fermion_labels[inds[j]]],
                 "scalar_T3_eigenvalues":[float(eig[a]),float(eig[bb])]}
                 for v,w,i,j,a,bb in forbidden[:12]],
            "warning":"Dual-action variants are counterfactual convention probes; lower Ward residual does not establish a valid basis change."
        }
        out=args.capture/f"{token}_SU2_weight_convention_audit.json"
        out.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
        print("\n",token)
        print("fermion algebra:",report["generator_algebra_fermion"])
        print("scalar algebra:",report["generator_algebra_scalar"])
        print("Ward convention probes:",variants)
        print("T3 forbidden fraction:",report["T3_weight_forbidden_fraction"])
        print("largest forbidden:",report["largest_T3_forbidden"][:2])
        print("report:",out)
if __name__=="__main__":main()
