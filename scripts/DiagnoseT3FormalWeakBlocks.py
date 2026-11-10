"""Read-only audit of original formal LLSS weak-index blocks."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from Numerical.parallel_t3.wilson.rough_stream2 import load_canonical

def n(x):return float(np.linalg.norm(x))
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--capture",type=Path,default=Path("output/parallel_t3/all16_ward_capture01"))
    a=p.parse_args()
    root=Path.cwd()
    for d in sorted(a.capture.glob("fit16_*_ward_diagnostic.json")):
        token=d.name.removesuffix("_ward_diagnostic.json")
        folder=root/"output/parallel_t3/rough_inputs"/token/"a/A_m4"
        case=load_canonical(root,folder,"A_-4",c5_boundary="require",assume_formal_identity=True)
        b=case.basis
        # Source tensor is in the retained canonical input; identify exact NPZ and key.
        source=None; source_path=None
        for name in ("canonical_input.npz","rough_input.npz","boundary.npz",
                     "fermion_boundary_canonical.npz","fermion_boundary_formal.npz"):
            f=folder/name
            if not f.is_file():continue
            with np.load(f,allow_pickle=False) as z:
                if "C_LLSS_formal_muF" in z:
                    source=z["C_LLSS_formal_muF"].copy();source_path=f;break
        if source is None:raise FileNotFoundError(f"Formal tensor missing in {folder}")
        if source.shape[:4]!=(3,3,2,2):
            raise ValueError(f"Unexpected source shape {source.shape}")
        blocks={}
        for i in range(2):
            for j in range(2):
                x=source[:,:,i,j,:,:]
                blocks[f"L{i+1}L{j+1}"]={"norm":n(x),
                    "squared_norm_fraction":n(x)**2/max(n(source)**2,1e-200)}
        # Compare weak patterns for each flavour pair. No coefficient edits.
        pairs=[]
        for p1 in range(3):
            for p2 in range(p1,3):
                x=source[p1,p2]
                pairs.append({"flavours":[p1+1,p2+1],
                    "weak_norms":[[n(x[i,j]) for j in range(2)] for i in range(2)],
                    "weak_block_gram_real":[[
                        float(np.vdot(x[i,j],x[k,l]).real)
                        for k,l in ((0,0),(0,1),(1,0),(1,1))]
                        for i,j in ((0,0),(0,1),(1,0),(1,1))]})
        # Preserve source metadata for tracing the producer.
        manifests={}
        for f in folder.glob("*.json"):
            try:
                v=json.loads(f.read_text(encoding="utf-8-sig"))
                manifests[f.name]={k:v.get(k) for k in
                    ("schema","format","source","producer","generator","model",
                     "gamma_cc_to_weyl_signed_factor","C_full_is_canonical")
                    if k in v}
            except (ValueError,UnicodeError):pass
        report={"trial":token,"status":"READ_ONLY_NOT_A_REPAIR",
                "formal_tensor_file":str(source_path),
                "formal_tensor_shape":list(source.shape),
                "formal_tensor_norm":n(source),
                "weak_blocks":blocks,"flavour_pair_weak_patterns":pairs,
                "manifests":manifests,
                "note":"Block norms do not themselves prove which contraction is wrong. Use the source producer and its SU2 Clebsch-Gordan conventions to establish the correction."}
        dest=a.capture/f"{token}_formal_weak_blocks.json"
        dest.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
        print("\n",token,"formal source:",source_path)
        print("weak blocks:",blocks)
        print("manifest provenance:",manifests)
        print("report:",dest)
if __name__=="__main__":main()
