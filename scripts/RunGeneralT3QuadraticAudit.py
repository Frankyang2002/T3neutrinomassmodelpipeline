"""Report the eight optional scalar quadratic mass structures (basis prototype)."""
from __future__ import annotations

import argparse,json
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from Numerical.experimental.GeneralT3ScalarQuadraticMasses import (load_manifest, quadratic_mass_hessian)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ms1-gev',type=float,default=1350)
    p.add_argument('--ms2-gev',type=float,default=1780)
    p.add_argument('--extra-mu2-real-gev2',type=float,default=0)
    p.add_argument('--extra-mu2-imag-gev2',type=float,default=0)
    a=p.parse_args()
    man=load_manifest(ROOT/'docs/general_scalar_trial_operator_manifest.json')
    result={}
    for name in sorted(man['models']):
        result[name]=quadratic_mass_hessian(man,name,a.ms1_gev**2,a.ms2_gev**2,
                                           complex(a.extra_mu2_real_gev2,a.extra_mu2_imag_gev2)) if (
                                           man['models'][name]['quadratic_complex_couplings_not_enabled']
                                           or (a.extra_mu2_real_gev2==0 and a.extra_mu2_imag_gev2==0)) else None
    result={k:v for k,v in result.items() if v is not None}
    path=ROOT/'output/general_t3_stage20/quadratic_mass_hessian_prototype.json'
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps({'status':'PrototypeNotThresholdMatched','models':result},indent=2)+'\n',encoding='utf-8')
    print('Quadratic prototypes:',len(result))
    print('Models with extra quadratic:',sum(bool(v['additional_quadratic_structures']) for v in result.values()))
    print('Report:',path)
if __name__=='__main__':main()
