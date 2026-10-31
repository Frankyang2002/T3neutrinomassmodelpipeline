"""Run independent complex Majorana action witnesses; no production writes."""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Numerical.experimental.T3BCMajoranaActionBatch import run

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--trials',type=int,default=64)
    p.add_argument('--native',type=Path,default=ROOT/'output/studies/t3_bc_majorana_chiral_rotation_v1/native_matchete.json')
    p.add_argument('--output',type=Path,default=ROOT/'output/studies/t3_bc_majorana_action_batch_v1/summary.json')
    a=p.parse_args()
    result=run(trials=a.trials,native=a.native if a.native.exists() else None)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(result['status'])
    print('max covariance residual:',max(result['max_covariance_residuals'].values()))
    print('negative control:',result['negative_control_min_residuals'])
    print('native:',result['native']['status'])
    print('physical C5 certified:',result['physical_C5_certified'])
    print('output:',a.output)
if __name__=='__main__':main()
