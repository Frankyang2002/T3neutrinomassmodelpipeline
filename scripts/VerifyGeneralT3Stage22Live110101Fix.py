"""Verify the new live Wolfram registered 110101 quartic invariants in B_-1/C_-1.

Run the pre-existing ExportGeneralT3Stage22SU2NativeCovariance.wl after installing
T3RGBetaModel.wl. This checker accepts ONLY post-fix live source evidence.
Geometry is not a physical coupling-normalisation certificate.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from scripts.Stage22RGBeta110101Covariance import covariance, ORIGINAL_FIRST, ORIGINAL_REUSED, EXPECTED, PROPOSED, components, INDEX_GRID


def check(data: dict) -> dict:
    if data.get('format') != 'stage22-live-rgbeta-su2-covariance-v1':
        raise ValueError('Unrecognised live Wolfram covariance evidence format')
    if data.get('physical_coupling_mapping_certified') is not False or data.get('potential_prefactors_verified') is not False:
        raise ValueError('Unexpected physical certification in live diagnostic')
    if set(data.get('models',{})) != {'B_-1','C_-1'}:
        raise ValueError('Expected exactly B_-1/C_-1')
    for label in ('B_-1','C_-1'):
        m=data['models'][label]
        if m.get('native_registered_symbols') != [ORIGINAL_FIRST,ORIGINAL_REUSED]:
            raise ValueError(f'{label}: unexpected registered coupling inventory')
        expr=m.get('native_expression_inputform')
        if not isinstance(expr,list) or len(expr)!=2:
            raise ValueError(f'{label}: missing two native invariant expressions')
        e=[''.join(s.split()) for s in expr]
        if e != [EXPECTED[ORIGINAL_FIRST],PROPOSED]:
            raise ValueError(f'{label}: corrected native expressions not registered')
        audits=m.get('native_covariance_audits')
        if not isinstance(audits,list) or len(audits)!=2:
            raise ValueError(f'{label}: missing native covariance tests')
        for i,expression in enumerate(expr):
            audit=audits[i]
            if audit.get('su2_covariant') is not True or audit.get('nonzero_residual_components_by_generator') != [0,0,0]:
                raise ValueError(f'{label}: native invariant {i} is not SU(2)-covariant')
            expected=[str(components(expression)[idx]) for idx in INDEX_GRID]
            if audit.get('component_tensor_flat_inputform') != expected:
                raise ValueError(f'{label}: live component tensor differs from the exact Python tensor')
    return {'models_verified':2,'registered_symbols_verified':4,
            'all_native_generator_residuals_zero':True,
            'physical_coupling_mapping_certified':False}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=Path('output/general_t3_stage22/stage22_110101_wolfram_native_covariance.json'))
    args=parser.parse_args()
    result=check(json.loads(args.input.read_text(encoding='utf-8-sig')))
    print('Stage22 live B/C corrected SU(2) invariant check:',result)
    print('Physical coupling normalisation: UNCERTIFIED')


if __name__=='__main__':main()
