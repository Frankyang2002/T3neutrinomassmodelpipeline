from __future__ import annotations
import argparse, json, shutil, subprocess, sys
from dataclasses import dataclass
from pathlib import Path

ROOT=Path(__file__).resolve().parent
W=ROOT/"wolfram"
OUT=W/"output"
RUN=W/"RunModel.wl"
EFT_ORDER=5
LOOP_ORDER=1

# Curated physically motivated T3 benchmarks; no arbitrary 52-point grid.
INTERESTING=[("A",0),("B",-1),("C",-1),("D",-2),("E",0)]
# One representative of every known T3 SU(2) class.  D(-2) is the
# scalar-line-swapped/conjugate benchmark corresponding to A(0).
SMOKE=[("B",-1),("C",-1),("A",0),("D",-2),("E",0)]
EXTENDED=[("A",0),("A",-2),("B",-1),("C",-1),("D",-2),("E",0),("E",-2)]

@dataclass
class Record:
    cls:str
    alpha:int
    rc:int
    summary:dict
    path:Path

def label(cls:str,alpha:int)->str:
    a=f"m{abs(alpha)}" if alpha<0 else f"p{alpha}"
    return f"T3_{cls}_alpha_{a}"

def run_one(cls:str,alpha:int,verbose:bool=False)->Record:
    path=OUT/label(cls,alpha)
    shutil.rmtree(path,ignore_errors=True); path.mkdir(parents=True,exist_ok=True)
    alpha_arg = f"m{abs(alpha)}" if alpha < 0 else f"p{alpha}"
    cmd=["wolframscript","-file",str(RUN),str(path),str(EFT_ORDER),str(LOOP_ORDER),cls,alpha_arg]
    print(f"Running T3-{cls}, alpha={alpha} ...", flush=True)
    p=subprocess.run(cmd,cwd=W,capture_output=True,text=True,check=False)
    (path/"wolfram_stdout.log").write_text(p.stdout,encoding="utf-8")
    (path/"wolfram_stderr.log").write_text(p.stderr,encoding="utf-8")
    if verbose:
        if p.stdout: print(p.stdout)
        if p.stderr: print(p.stderr,file=sys.stderr)
    elif p.returncode != 0:
        # On failure, surface enough context to diagnose without dumping every
        # successful Matchete stage in normal runs.
        if p.stdout:
            print("\n".join(p.stdout.splitlines()[-40:]))
        if p.stderr:
            print("\n".join(p.stderr.splitlines()[-20:]),file=sys.stderr)
    sp=path/"comparison_summary.json"
    summary=json.loads(sp.read_text(encoding="utf-8")) if sp.exists() else {"BuildStatus":"ProcessFailed","MatchingStatus":"NotRun"}
    return Record(cls,alpha,p.returncode,summary,path)

def main()->int:
    ap=argparse.ArgumentParser()
    mode=ap.add_mutually_exclusive_group()
    mode.add_argument("--smoke",action="store_true",help="5 representative T3 models")
    mode.add_argument("--extended",action="store_true",help="7 viable/equivalent benchmark points")
    ap.add_argument("--verbose",action="store_true",help="show full Wolfram/Matchete output")
    ns=ap.parse_args()
    points=SMOKE if ns.smoke else EXTENDED if ns.extended else INTERESTING
    OUT.mkdir(parents=True,exist_ok=True)
    print(f"T3 scan mode: {'smoke' if ns.smoke else 'extended' if ns.extended else 'interesting'}; {len(points)} model(s).")
    records=[run_one(c,a,verbose=ns.verbose) for c,a in points]
    print("\n"+"="*72+"\nT3 MODEL SUMMARY\n"+"="*72)
    ok=0
    for r in records:
        s=r.summary
        success=s.get("BuildStatus")=="Success" and s.get("MatchingStatus")=="Success"
        ok+=success
        print(f"T3-{r.cls} alpha={r.alpha}: build={s.get('BuildStatus')}, match={s.get('MatchingStatus')}, T3={s.get('T3IngredientsPresent')}, Weinberg={s.get('WeinbergOperatorPresent')}")
        if s.get("WeinbergOperatorPresent"):
            extraction=s.get("WeinbergExtractionStatus", "Unknown")
            hcount=s.get("WeinbergHolomorphicTermCount", 0)
            hccount=s.get("WeinbergConjugateTermCount", 0)
            if extraction == "Success":
                print(f"  C5 extraction: Success ({hcount} holomorphic + {hccount} HC terms)")
                c5_file=s.get("WeinbergCoefficientFile")
                if c5_file:
                    print(f"  C5: {r.path/c5_file}")
            else:
                print(f"  Weinberg terms: {s.get('WeinbergTermCount', 0)}; C5: {extraction}")
            if ns.verbose:
                raw_file=s.get("WeinbergRawFile")
                if raw_file:
                    print(f"  raw: {r.path/raw_file}")
        if ns.verbose:
            print("  accepted: "+", ".join(s.get("AcceptedInteractions",[])))
    aggregate=OUT/"t3_model_comparison.json"
    aggregate.write_text(json.dumps([r.summary for r in records],indent=2),encoding="utf-8")
    print(f"\n{ok}/{len(records)} completed build+matching.\nAggregate: {aggregate}")
    return 0 if ok==len(records) else 1

if __name__=="__main__": raise SystemExit(main())
