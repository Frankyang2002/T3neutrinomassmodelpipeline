"""Fail-closed Stage24 local post-fix live RGBeta SU(2) comparison.

This verifies 8 registered native quartic contractions against independently
specified component tensors and the original Matchete quartic tensor basis.
It is NOT a physical coupling-normalisation certificate.

python -m scripts.VerifyGeneralT3Stage24BCPostFix
"""
from __future__ import annotations
import argparse
import json
from itertools import product
from pathlib import Path
import sympy as sp

from scripts.Stage24BCRemainderCovariance import (
    PROPOSED, FIELDS, MODELS, components, group_blocks,
)
from scripts.Stage22TensorAlgebra import solve_basis, exact

ROOT = Path('output/general_t3_stage22')
FORMAT = 'stage24-live-rgbeta-bc-postfix-covariance-v1'


def verify(wl: dict, matchete: dict) -> dict:
    if wl.get('format') != FORMAT or wl.get('physical_coupling_mapping_certified') is not False:
        raise ValueError('Wrong Stage24 post-fix Wolfram source phase or certification status')
    if wl.get('native_non_covariant_slots') != 0 or wl.get('native_covariant_slots') != 8:
        raise ValueError('Native Stage24 covariance inventory incomplete')
    if set(wl.get('models', {})) != set(MODELS):
        raise ValueError('Expected B_-1 and C_-1 post-fix native models')
    if not isinstance(matchete, dict) or any(m not in matchete.get('models',{}) for m in MODELS):
        raise ValueError('Matchete models missing')
    n = 0
    for model in MODELS:
        records = wl['models'][model]
        if set(records) != set(PROPOSED):
            raise ValueError(model+': incomplete native quartic symbol inventory')
        mt = matchete['models'][model]['tensors']
        bycount = {}
        for row in mt:
            k = ''.join(str(v) for v in row['counts'])
            bycount.setdefault(k, []).append(row)
        for symbol, entry in records.items():
            ref = [components(PROPOSED[symbol])[idx]
                   for idx in product(range(2), repeat=4)]
            for tag in ('native', 'proposed'):
                audit = entry.get(tag, {})
                if audit.get('su2_covariant') is not True or audit.get('nonzero_residual_components_by_generator') != [0,0,0]:
                    raise ValueError(model+'/'+symbol+': noncovariant registered '+tag)
                actual = audit.get('component_tensor_flat_inputform')
                if not isinstance(actual,list) or len(actual) != len(ref):
                    raise ValueError(model+'/'+symbol+': incomplete '+tag+' components')
                if any(sp.simplify(exact(v) - wanted) != 0 for v,wanted in zip(actual,ref)):
                    raise ValueError(model+'/'+symbol+': incorrect native '+tag+' tensor')
            key = ''.join(str(FIELDS[symbol].count(f)) for f in ('H','Hdag','S1','S1dag','S2','S2dag'))
            mats = bycount.get(key,[])
            if len(mats)!=1 or mats[0]['species'] != FIELDS[symbol]:
                raise ValueError(model+'/'+symbol+': Matchete reference field order or rank changed')
            candidate = {
                'counts': mats[0]['counts'], 'dimensions': mats[0]['dimensions'],
                'species': mats[0]['species'], 'symbol': symbol,
                'components': entry['native']['component_tensor_flat_inputform']
            }
            geom = solve_basis(mats,[candidate])
            if geom.get('r_in_m') != [['1']] or geom.get('maximum_exact_component_residual') != '0':
                raise ValueError(model+'/'+symbol+': exact native Matchete geometry changed')
            n += 1
    if n!=8:
        raise ValueError('Wrong post-fix slot count')
    return {'models_verified':2,'native_su2_pass':8,
            'native_exact_matchete_tensor_matches':8,
            'global_field_component_basis_certified':False,
            'physical_coupling_mapping_certified':False}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir',type=Path,default=ROOT)
    p.add_argument('--wolfram-export',type=Path,
                   default=ROOT/'stage24_bc_postfix_wolfram_native_covariance.json')
    a=p.parse_args()
    try:
        wl=json.loads(a.wolfram_export.read_text(encoding='utf-8-sig'))
        m=json.loads((a.input_dir/'matchete_tensors.json').read_text(encoding='utf-8-sig'))
        status=verify(wl,m)
    except (ValueError,KeyError,TypeError,IndexError,OSError) as exc:
        p.exit(2,'Stage24 post-fix verification FAILED: '+str(exc)+'\n')
    print('Stage24 post-fix live registered RGBeta verification:',status)
    print('Physical coupling normalisation: UNCERTIFIED')

if __name__=='__main__':
    main()
