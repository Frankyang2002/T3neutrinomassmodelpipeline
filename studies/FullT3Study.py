"""Complete 16-model numerical T3 study."""
from __future__ import annotations
from copy import deepcopy
import argparse, json, shutil, subprocess, sys, time
from pathlib import Path
from typing import Any
from Numerical.InteractiveModelComparison import write_dashboard

PROJECT_ROOT=Path(__file__).resolve().parents[1]
PIPELINE_SCRIPT=PROJECT_ROOT/"pipeline.py"
OUTPUT_ROOT=PROJECT_ROOT/"output"/"full"
REPORT_ROOT=PROJECT_ROOT/"Reports"/"output"/"full"
DEFAULT_TEMPLATE=PROJECT_ROOT/"configs"/"t3_numerical_default_B_alpha_m1.json"
ACTIVE_CONFIG_DIR=PROJECT_ROOT/"configs"/"generated_models"
CONFIG_SET_ROOT=PROJECT_ROOT/"configs"/"generated_model_sets"

MODELS=(
 ("A",1,3,2,-4),("A",1,3,2,-2),("A",1,3,2,0),
 ("B",2,2,1,-3),("B",2,2,1,-1),("B",2,2,1,1),
 ("C",2,2,3,-3),("C",2,2,3,-1),("C",2,2,3,1),
 ("D",3,1,2,-2),("D",3,1,2,0),("D",3,1,2,2),
 ("E",3,3,2,-4),("E",3,3,2,-2),("E",3,3,2,0),("E",3,3,2,2),
)
COMPARISON_SCENARIOS={
 "smallY_smallL":(0.005,0.1),
 "smallY_largeL":(0.005,1.0),
 "largeY_smallL":(0.5,0.1),
 "largeY_largeL":(0.5,1.0),
}
OPTIONAL_REAL_QUARTICS=("lambdaH1Adj","lambdaH2Adj","lambdaS1Adj","lambdaS2Adj","lambda12Adj","lambda12Cross")
SPECIAL_DOUBLET_COMPLEX_QUARTICS=("lambdaHHdagS2S2","lambdaHHdagS1barS1bar","lambdaS1bar2S2bar2","lambdaS1barS2S2bar2","lambdaS1S1bar2S2bar","lambdaHHdagS1barS2barCross")
COMPARISON_Y1_TEXTURE=((1,.2,-.1),(.2,.8,.3),(-.1,.3,.6))
COMPARISON_Y2_TEXTURE=((.7,-.3,.2),(.4,1,-.2),(.1,.3,.8))

def _model_key(ds1:int,ds2:int,df:int,alpha:int)->str:
 a=f"p{alpha}" if alpha>=0 else f"m{abs(alpha)}"; return f"T3_dS1_{ds1}_dS2_{ds2}_dF_{df}_alpha_{a}"
def _load_json(path:Path)->dict[str,Any]:
 v=json.loads(path.read_text(encoding="utf-8"))
 if not isinstance(v,dict): raise ValueError(f"{path} must contain a JSON object")
 return v
def _write_json(path:Path,value:dict[str,Any])->None:
 path.parent.mkdir(parents=True,exist_ok=True); tmp=path.with_name(path.name+".tmp")
 tmp.write_text(json.dumps(value,indent=2)+"\n",encoding="utf-8"); tmp.replace(path)
def _common_args(args:argparse.Namespace)->list[str]:
 r=[]
 if getattr(args,"debug_reports",False):r.append("--debug-reports")
 if getattr(args,"force",False):r.append("--force")
 if getattr(args,"threshold",None):
  for g in args.threshold:r+=["--threshold",*g]
 if getattr(args,"threshold_scale",None):
  for x in args.threshold_scale:r+=["--threshold-scale",str(x)]
 return r
def _run(command:list[str])->int:
 print("\n> "+" ".join(command),flush=True); return int(subprocess.run(command,cwd=PROJECT_ROOT).returncode)
def _snapshot_active(dst:Path)->None:
 if dst.exists():shutil.rmtree(dst)
 shutil.copytree(ACTIVE_CONFIG_DIR,dst) if ACTIVE_CONFIG_DIR.exists() else dst.mkdir(parents=True,exist_ok=True)
def _restore_active(src:Path)->None:
 if ACTIVE_CONFIG_DIR.exists():shutil.rmtree(ACTIVE_CONFIG_DIR)
 shutil.copytree(src,ACTIVE_CONFIG_DIR) if src.exists() else ACTIVE_CONFIG_DIR.mkdir(parents=True,exist_ok=True)
def _normalise_representation_quartics(o:dict[str,Any],ds1:int,ds2:int,alpha:int)->None:
 req={"lambdaH1Adj":ds1>1,"lambdaH2Adj":ds2>1,"lambdaS1Adj":ds1==3,"lambdaS2Adj":ds2==3,"lambda12Adj":ds1>1 and ds2>1,"lambda12Cross":ds1==3 and ds2==3}
 for n in OPTIONAL_REAL_QUARTICS:
  o.setdefault(n,0.0) if req[n] else o.pop(n,None)
 special=ds1==2 and ds2==2 and alpha==-1
 for n in SPECIAL_DOUBLET_COMPLEX_QUARTICS:
  o.setdefault(n,{"real":0.0,"imag":0.0}) if special else o.pop(n,None)
def _retarget_config(template:dict[str,Any],ds1:int,ds2:int,df:int,alpha:int)->dict[str,Any]:
 c=deepcopy(template); c["representation"]={"d_s1":ds1,"d_s2":ds2,"d_f":df,"alpha":alpha,"shared_scalar":False}
 b=c.get("base_state")
 if not isinstance(b,dict):raise ValueError("Numerical config requires base_state object.")
 o=b.get("ordinary")
 if not isinstance(o,dict):raise ValueError("Numerical config requires base_state.ordinary object.")
 _normalise_representation_quartics(o,ds1,ds2,alpha); return c
def _normalise_saved_model_config(path:Path,ds1:int,ds2:int,df:int,alpha:int)->None:
 if path.is_file():_write_json(path,_retarget_config(_load_json(path),ds1,ds2,df,alpha))
def _scaled_real_texture(t,scale):return [[scale*float(v) for v in row] for row in t]
def _fixed_comparison_config(template,ds1,ds2,df,alpha,yukawa,scalar):
 c=_retarget_config(template,ds1,ds2,df,alpha);o=c["base_state"]["ordinary"]
 o["y1"]["real"]=_scaled_real_texture(COMPARISON_Y1_TEXTURE,yukawa);o["y2"]["real"]=_scaled_real_texture(COMPARISON_Y2_TEXTURE,yukawa)
 o["y1"]["imag"]=[[0.0]*3 for _ in range(3)];o["y2"]["imag"]=[[0.0]*3 for _ in range(3)];o["lambdaT3"]={"real":scalar,"imag":0.0}
 search=c.setdefault("benchmark_search",{});search["enabled"]=False;search["use_current_model_representation"]=False;c.setdefault("sensitivity",{})["enabled"]=False
 c["comparison_point"]={"yukawa_scale":yukawa,"lambdaT3_real":scalar,"representation_specific_extra_quartics":"required optional quartics fixed to 0","texture":"fixed real non-diagonal y1/y2 textures scaled by common Y; Im(y1)=Im(y2)=0; Im(lambdaT3)=0","y1_texture":[list(r) for r in COMPARISON_Y1_TEXTURE],"y2_texture":[list(r) for r in COMPARISON_Y2_TEXTURE]}
 return c
def _dashboard(raw:Path,report:Path,title:str)->str|None:
 try:
  p=write_dashboard(raw,report/"interactive_comparison.html",title);print(f"Interactive comparison: {p}");return str(p)
 except FileNotFoundError as e:print(f"Interactive comparison skipped: {e}");return None

def run_full_study(args:argparse.Namespace)->int:
 tp=Path(args.numerical) if getattr(args,"numerical",None) else DEFAULT_TEMPLATE
 if not tp.is_absolute():tp=PROJECT_ROOT/tp
 template=_load_json(tp);OUTPUT_ROOT.mkdir(parents=True,exist_ok=True);REPORT_ROOT.mkdir(parents=True,exist_ok=True)
 optimal_store=CONFIG_SET_ROOT/"optimal";comparison_store=CONFIG_SET_ROOT/"comparison";common=_common_args(args);started=time.time();results=[];status=0
 print("="*72);print("FULL 16-MODEL T3 NUMERICAL STUDY");print("="*72)
 print("Modes: per-model optimal benchmark + four common-parameter comparisons")
 print("Comparison scales: Y={0.005,0.5}, lambdaT3={0.1,1.0}")
 print("Mass thresholds: MF=100 TeV, MS=1 TeV; UV=1e7 GeV; low=100 GeV")
 optimal_started=time.time();optimal_rc=0;inp=optimal_store/"input"
 if inp.exists():shutil.rmtree(inp)
 inp.mkdir(parents=True,exist_ok=True)
 for _,ds1,ds2,df,alpha in MODELS:
  key=_model_key(ds1,ds2,df,alpha);clean=_retarget_config(template,ds1,ds2,df,alpha);search=clean.setdefault("benchmark_search",{});search["enabled"]=True;search["use_current_model_representation"]=False
  cp=inp/f"{key}.json";_write_json(cp,clean);saved=ACTIVE_CONFIG_DIR/f"{key}.json"
  if saved.is_file() and not getattr(args,"reset_numerical_configs",False):_normalise_saved_model_config(saved,ds1,ds2,df,alpha)
  cmd=[sys.executable,str(PIPELINE_SCRIPT),"--dims",str(ds1),str(ds2),str(df),"--alpha",str(alpha),"--study",f"full/optimal/{key}","--numerical",str(cp),*common]
  if getattr(args,"reset_numerical_configs",False):cmd.append("--reset-numerical-configs")
  rc=_run(cmd);optimal_rc=max(optimal_rc,rc);status=max(status,rc)
 _snapshot_active(optimal_store);dash=_dashboard(OUTPUT_ROOT/"optimal",REPORT_ROOT/"optimal","T3 optimal-benchmark comparison (16 models)")
 results.append({"mode":"optimal","return_code":optimal_rc,"runtime_seconds":time.time()-optimal_started,"dashboard":dash})
 for scenario,(yukawa,scalar) in COMPARISON_SCENARIOS.items():
  t=time.time();scenario_rc=0;inp=comparison_store/scenario/"input";resolved=comparison_store/scenario/"resolved"
  for d in (inp,resolved):
   if d.exists():shutil.rmtree(d)
   d.mkdir(parents=True,exist_ok=True)
  for _,ds1,ds2,df,alpha in MODELS:
   key=_model_key(ds1,ds2,df,alpha);cp=inp/f"{key}.json";_write_json(cp,_fixed_comparison_config(template,ds1,ds2,df,alpha,yukawa,scalar))
   cmd=[sys.executable,str(PIPELINE_SCRIPT),"--dims",str(ds1),str(ds2),str(df),"--alpha",str(alpha),"--study",f"full/comparison/{scenario}/{key}","--numerical",str(cp),"--reset-numerical-configs",*common]
   rc=_run(cmd);scenario_rc=max(scenario_rc,rc);status=max(status,rc);g=ACTIVE_CONFIG_DIR/f"{key}.json"
   if g.is_file():shutil.copy2(g,resolved/g.name)
  dash=_dashboard(OUTPUT_ROOT/"comparison"/scenario,REPORT_ROOT/"comparison"/scenario,f"T3 common-parameter comparison: {scenario} (Y={yukawa:g}, λT3={scalar:g})")
  results.append({"mode":"comparison","scenario":scenario,"yukawa":yukawa,"lambdaT3":scalar,"return_code":scenario_rc,"runtime_seconds":time.time()-t,"dashboard":dash})
 _restore_active(optimal_store)
 summary={"status":"Success" if status==0 else "Failed","ordinary_model_count":len(MODELS),"shared_scalar_models_included":False,
 "comparison_definition":{"yukawa_texture":"fixed real non-diagonal y1/y2 textures scaled by common Y","scalar_parameter":"Re(lambdaT3), Im(lambdaT3)=0","representation_specific_extra_quartics":"required optional quartics fixed to 0","scenarios":{k:{"yukawa":v[0],"lambdaT3":v[1]} for k,v in COMPARISON_SCENARIOS.items()}},
 "results":results,"runtime_seconds":time.time()-started}
 _write_json(OUTPUT_ROOT/"full_numerical_summary.json",summary);print("\n"+"="*72);print("FULL NUMERICAL SUMMARY");print("="*72);print(f"Status: {summary['status']}");print(f"Models: {len(MODELS)} ordinary neutral-compatible");print(f"Summary: {OUTPUT_ROOT/'full_numerical_summary.json'}");print(f"Interactive reports: {REPORT_ROOT}");return status
