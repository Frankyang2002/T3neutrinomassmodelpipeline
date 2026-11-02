"""Validate and merge Chat1 certified canonical tensor + EFT1 couplings into one file."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0,str(ROOT))
from RGE.running.intermediate.CertifiedFermionBoundaryAdapter import adapt_chat1_boundary


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--chat1", type=Path, required=True)
    p.add_argument("--eft1-couplings", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a=p.parse_args(argv)
    with a.chat1.open(encoding="utf-8-sig") as f:
        chat1=json.load(f)
    with a.eft1_couplings.open(encoding="utf-8-sig") as f:
        eft1=json.load(f)
    result=adapt_chat1_boundary(chat1,eft1)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open("x", encoding="utf-8") as f:
        json.dump(result,f,indent=2,allow_nan=False)
        f.write("\n")
    print("Created certified boundary",a.output)
    return 0

if __name__=="__main__":
    raise SystemExit(main())
