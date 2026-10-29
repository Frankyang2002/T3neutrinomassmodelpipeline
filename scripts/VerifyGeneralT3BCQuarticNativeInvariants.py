"""Fail-closed comparison of independent Stage24 Wolfram and Python B/C audits.

python -m scripts.VerifyGeneralT3BCQuarticNativeInvariants
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from itertools import product
import sympy as sp
from scripts.RGBetaTensorAlgebra import exact
from scripts.BCQuarticRemainderCovariance import PROPOSED, OLD, MODELS, components, run as run_python

DEFAULT=Path('output/general_t3_stage22')


def verify(py, wl):
    if py.get('totals',{}).get('old_non_covariant_slots')!=8 or py.get('physical_coupling_mapping_certified') is not False:
        raise ValueError('Invalid Stage24 Python evidence')
    if wl.get('format')!='stage24-live-rgbeta-bc-other-quartic-covariance-v1' or wl.get('physical_coupling_mapping_certified') is not False:
        raise ValueError('Invalid Wolfram Stage24 export')
    if set(wl.get('models',{})) != set(MODELS):raise ValueError('B/C live model count mismatch')
    n=0
    for model in MODELS:
        w=wl['models'][model]
        if set(w)!=set(PROPOSED):raise ValueError('Unexpected B/C live quartic name inventory')
        for sym,item in w.items():
            expected_new=[components(PROPOSED[sym])[i] for i in product(range(2),repeat=4)]
            expected_old=[components(OLD[sym])[i] for i in product(range(2),repeat=4)]
            for key,values,flag in [('native',expected_old,False),('proposed',expected_new,True)]:
                audit=item[key]
                if audit.get('su2_covariant') is not flag or len(audit.get('nonzero_residual_components_by_generator',[]))!=3:
                    raise ValueError(model+'/'+sym+': unexpected live generator result')
                residuals=audit['nonzero_residual_components_by_generator']
                if (residuals==[0,0,0]) is not flag:
                    raise ValueError(model+'/'+sym+': live residual counts incoherent')
                actual=audit.get('component_tensor_flat_inputform')
                if not isinstance(actual,list) or len(actual)!=16 or any(sp.simplify(exact(v)-e)!=0 for v,e in zip(actual,values)):
                    raise ValueError(model+'/'+sym+': Wolfram/Python component tensor mismatch: '+key)
            n+=1
    return {'models_verified':2,'native_failed_expected':n,'proposed_su2_pass':n,
            'exact_python_wolfram_components_agree':True,
            'physical_coupling_mapping_certified':False}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir',type=Path,default=DEFAULT)
    p.add_argument('--wolfram-export',type=Path,default=DEFAULT/'stage24_bc_wolfram_native_covariance.json')
    args=p.parse_args()
    try:
        if not args.wolfram_export.is_file():raise ValueError('Missing live Wolfram Stage24 export')
        py=run_python(args.input_dir)
        wl=json.loads(args.wolfram_export.read_text(encoding='utf-8-sig'))
        result=verify(py,wl)
    except (ValueError,KeyError,TypeError,IndexError,OSError) as exc:
        p.exit(2,'Stage24 live check FAILED: '+str(exc)+'\n')
    print('Stage24 live B/C proposed invariant check:',result)
    print('Physical coupling mapping: UNCERTIFIED; production unchanged')

if __name__=='__main__':main()
